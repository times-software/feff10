! Fixed problem with number density in subroutine opcons and problem with
! format 36 that gave too little precision to the energies in opcons.dat.
! MPP 12/17/03

! This is the 4th and last piece of the new modular version of the code.
! MPP 8/11/03

! Kevin Jorissen 2014
! turning this into a "FEFF9 module" of sorts, for use with JFEFF
! Changing name from "step4" to "fullspectrum"
! Changing "core.inp" to "fullspectrum.inp"
! step1 and step2 can be run by JFEFF (in Java) (step3 hasn't existed in ages)

      program fullspectrum
      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu, init_dimensions
      use fullspectrum_inp
      use par
      use errorfile
      implicit none
      include 'HEADERS/params.h'


      character*2 l1*75, l2*75, l3*75,flnm2*20,cmd*80,cmnt 
      character*2 edname
      character*2, allocatable :: alledg(:,:)
      integer  ie, istrln, ihole
      character*3, allocatable :: cmpnm(:)
      integer i, j, k,iepts, ncomps, cntrl(1:6),                        &
     &        ios,ierr,ndos,l, nhi, nlo, z
      integer, allocatable :: numedg(:),nz(:),errno (:,:)
      real*8, allocatable :: numden(:)
      real*8 emin, emax, emaxi, myemin, myemax, myestp, am, rho,          &
     &     tau,ndrude, omega(fullpts),eps2(fullpts),       &
     &     eps2bk(fullpts), edge(fullpts), edgebk(fullpts),             &
     &     edos(maxpts), dos(0:lmax,maxpts),curdos(maxpts),efermi,      &
     &     numsum, neff, norm, omgcnv(numcnv),stpcnv,edgcnv(numcnv),    &
     &     x,y, eps1(fullpts), temp,repart(fullpts),impart(fullpts),    &
     &     sigma(fullpts), valeps2(fullpts)
      double precision tempdp
      character*512 slog, str, str2, str3, root,flnm, path,file1,file2
!     these logical vars are false if reading some file needed to 
!     compute for an edge fails
      logical rdokc(mxedgs),rdokv(mxedgs),drude,flag, val, eels
      logical, allocatable :: finest (:,:),conv(:,:), ldos(:)
      complex eps(fullpts),epsbk(fullpts) !for a single edge
      complex efree(fullpts) !for free electrons (drude model)
      complex epstot(fullpts),epstot0(fullpts) !summed over all edges
      !these only used for kk test; can be removed when satisfied w/ kk
      real*8 gam,omega0,width
      complex*16 epsi(fullpts)
      !Kevin Jorissen to bypass Hamaker test:
      logical,parameter :: dohamaker=.false.

      call par_begin
      if (worker) goto 400
      call OpenErrorfileAtLaunch('ff2x')	

      call init_dimensions

      allocate(numden(nphx),numedg(nphx),nz(nphx),errno (2*mxedgs,nphx),cmpnm(nphx),alledg(2*mxedgs,nphx), &
          finest (2*mxedgs,nphx),conv(2*mxedgs,nphx), ldos(nphx) )


      call fullspectrum_read
      if (mFullSpectrum .le. 0)  goto 400

!     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='logfullspectrum.dat', status='unknown', iostat=ios)
      call chopen (ios, 'logfullspectrum.dat', 'fullspectrum')
!     open file to store partial f-sumrule for each edge
      open (unit=15, file='osc_str.dat', status='unknown', iostat=ios)
      call chopen (ios, 'osc_str.dat', 'step4')
      write (unit=15, fmt="(a)") '# component  edge  n_eff'
      write (unit=15, fmt="(a)") ' '
      open (unit=16, file='xmu.dat',iostat=ios)
      call chopen (ios, 'xmu.dat', 'step4')
      write (unit=16, fmt="(a)") '# component  edge  n_eff'
      write (unit=16, fmt="(a)") ' '

!     initialize arrays
      eps2(:)=0.0
      eps2bk(:)=0.0
      edge(:)=0.0
      edgebk(:)=0.0
      sigma(:)=0.0

!     read input file
      call rdop (cntrl,emin,emax,iepts,ncomps,cmpnm,nz,numden,numedg,alledg,conv,ldos,finest,drude,tau,ndrude,val,eels)

!     construct energy grid omega for calculation 
      call egrid (emin, emax, nz, ncomps, omega, iepts)

!     separate linear grid for convolution
      emin=omega(1)
      stpcnv=(ecnv-emin)/(numcnv-1)
      do i=1,numcnv
        omgcnv(i)=emin+(i-1)*stpcnv
      enddo 

!*******************************  loop over components collecting the contribution to epsilon from each distinguishable site.
      do j=1,ncomps

!       store root path for this component. root will be something like 'edges/Cu2/'
        str='edges/'
        str2=cmpnm(j)
        call concat(str,str2,str3,k)
        str='/'
        call concat(str3,str,root,k)

        !if we are doing any convolution, read LDOS around the current atomic site
        if (ldos(j)) then
          str='ldos/ldos00.dat'
          call concat(root,str,path,k)
          call rdldos (path,flag,ndos,efermi,edos,dos)
          if (.not. flag) then !reading of ldos failed, so don't convolve
            write (slog,fmt="('    ',a4, 'edges will not be convolved.')") cmpnm(j)
            conv(1:2*mxedgs,j)=.false.
          endif
        endif


        !if no number density was specified in the input file, estimate it from pot.bin 
        if (numden(j).le.0.0) then
          !store pathname to first edge directory for component j
          i=1
          str=alledg(i,j)
          call concat(root,str,path,k)
          !need to pass atomic number so rddens can count the number of atoms of that species
          k=nz(j)
          !rddens estimates number density
          call rddens (path,k,tempdp)
          numden(j)=dble(tempdp)
          write (slog, fmt="('number density of ',a3)") cmpnm(j)
          call wlog(slog)
          write(slog,fmt= "('estimated to be ', e20.10,' inverse cubic angstroms')") numden(j)/bohr**3
          call wlog(slog)
        end if



!********************** loop over initial states around the jth type of site.
        do i=1,numedg(j)
          norm=1.0
!         store full path for this edge
          str=alledg(i,j)
          call concat(root,str,path,k)
!         read abs. data from xmu.dat files into eps and epsbk.
!         eps, epsbk contain the complex scattering factor f after this call.
          call addedg(path,omega, iepts, eps, epsbk, neff, ierr)
          errno(i,j)=ierr
          do k=1,2000
             write(unit=77,fmt="(20e20.10)") omega(k),eps(k),epsbk(k)
          enddo

          !up to this point our eps arrays have actually held f
          !now we convert f --> epsilon
          !store real and im parts in real arrays b/c of problems using complex arithmetic when re part holds NaN
          do k=1,iepts
            repart(k)=-1*4*pi*numden(j)*dble(eps(k))/omega(k)**2
            impart(k)=-1*4*pi*numden(j)*aimag(eps(k))/omega(k)**2
            sigma(k)=sigma(k)-aimag(eps(k))/alpinv/bohr**2/omega(k)
          enddo

          !back to cmplx variables
          do k=1,iepts
            eps(k)=repart(k)+coni*impart(k)
          enddo
          !repeat the charade for the atomic background
          !store real and im parts in real arrays b/c of problems using complex arithmetic when re part holds NaN
          do k=1,iepts
            repart(k)=-1*4*pi*numden(j)*dble(epsbk(k))/omega(k)**2
            impart(k)=-1*4*pi*numden(j)*aimag(epsbk(k))/omega(k)**2
          enddo

          !back to cmplx variables
          do k=1,iepts
            epsbk(k)=repart(k)+coni*impart(k)
          enddo

          !load im part into separate arrays to use old subroutines from KK version of code.
          do k=1,iepts
            edge(k)=aimag(eps(k))
            edgebk(k)=aimag(epsbk(k))
          enddo

!         see how many effective electrons in this edge before cnvlt
          numsum=numden(j)
          call qsum(numsum,edge,omega,iepts,neff)
          edname=alledg(i,j)
          call setedg(edname,ihole)
          !report n_eff for this edge in osc_str.dat
          write (unit=15, fmt="(a11,a6,i4,f8.3)") cmpnm(j),edname,ihole,neff
          write (unit=16, fmt="('# ',a11,a6,f8.3)") cmpnm(j),alledg(i,j),neff

          !add contribution from current initial state to total dielectric function. 
          do k=1,iepts
            epstot(k)=epstot(k)+eps(k)
            epstot0(k)=epstot0(k)+epsbk(k)
          enddo

        enddo
!********************** end of loop over initial states around the jth type of site.

      enddo
!**************************************** end of big loop over components


! ****************** report on succes
      !now that all edges are read, report via the log which ones were problematic.
      !See if any edges failed; if so report them.
      j=0
      do i=1,ncomps
        do ie=1,numedg(i)
           if (errno(ie,i).ne.0) j=j+1
        enddo
      enddo
      if(j.ne.0) call wlog('We attempted to include the following edges, but failed.')
      do i=1,ncomps
        str=cmpnm(i)
        j=istrln(str)
        do ie=1,numedg(i)
           if (errno(ie,i).ne.0) then
             slog='    '//str(1:j)//': '//alledg(ie,i)
             call wlog(slog)
           endif
        enddo
      enddo
      !Even if all edges were read succesfully, report the number of edges read succesfully.
      k=0
      j=0
      do i=1,ncomps
        j=j+numedg(i) !j holds total number of edges.
        do ie=1,numedg(i)
           if (errno(ie,i).eq.0) k=k+1 !k holds number of good edges
        enddo
      enddo
      write(slog,fmt="('Succeeded in computing ',i2,' out of ',i2,' edges.')") k,j
      call wlog(slog)
      if (k.eq.0) then
        slog='No edges were computed succesfully, so I am stopping without writing output.'
        call wlog (slog)
        stop
      else
        write(slog,fmt="('k: ',i2)")k
        call wlog(slog)
      end if
!********************** end of report on success


!********************* make and write output

!**************** xmu.val
      if (val) then
        !loop over components collecting the contribution to epsilon from each distinguisable site.

        do j=1,ncomps
!         store root path for this component. root will be something like 'edges/Cu2/'
          str='edges/'
          str2=cmpnm(j)
          call concat(str,str2,str3,k)
          str='/valence/xmu.val'
          call concat(str3,str,path,k)
          write (slog,fmt="('path: ',a400)") path
          call wlog(slog)

          !read valence code results and add them here
          call wlog ('adding in valence response from '//path)
          call rdval(numden(j),iepts,omega,path,edge)
          do ie=1,iepts
            epstot(ie)=epstot(ie)+coni*edge(ie)
          enddo
        enddo
      end if


!**************** xmu.dat
      !write header to our 'fake' xmu.dat file
      write(unit=16,fmt="( a2, 1x, i4, '/', i4, ' paths used')") '# ',0,0
      write(unit=16,fmt="(a2, ' xsedge+ 50, used to normalize mu ', 1pe20.4)") '# ',1.0
      write(unit=16,fmt="(a2)") '# '
      write(unit=16,fmt="(a2)") '# '

!**************** eps.dat
      !epstot, epstot0 now hold the entire contribution from the bound charges to the dielectric function. 
      !write them to a file.
      open (file='eps.dat', unit=12)
      do i=2,iepts-1
        write (unit=12, fmt="(6e20.10)") omega(i),epstot(i),epstot0(i),sigma(i)
        write (unit=16, fmt="(6e20.10)") omega(i)*hart,omega(i)*hart,sqrt(2*omega(i))/bohr,sigma(i),sigma(i)
      enddo
      close (unit=12)

!**************** opcons.dat and opconsKK.dat
      if (cntrl(6).eq.1) then

!**************** drude.dat
        if (drude) then !compute drude term
          ndrude=ndrude*numden(1)
          call drdtrm(omega,efree,iepts,tau,ndrude)
          call wlog ('drude term added; written in drude.dat')
        else
          efree(1:iepts)=0.0
        endif 
        !efree either holds the drude term or 0s

!**************** opcons.dat
        !load eps with contribution from free and bound charges to
        !compute the optical constants. This calculation uses FEFF
        !results for both real and im parts of epsilon, augmented by a
        !drude term if the user specified one in the input file.
        eps(1:iepts)=epstot(1:iepts)+efree(1:iepts)
        !store pathname to first edge directory for component j
        i=1
        j=ncomps
        str=alledg(i,j)
        call concat(root,str,path,k)
        str=path !save pathname for call to write opcons0.dat
        !write out optical constants to file
        flnm='opcons.dat'
        call opcons(omega, eps, iepts,flnm, path)

!**************** opconsKK.dat
        !run KK transforms to test FEFF calculation of f'.
        !first load eps2 with our final epsilon_2 ...
        eps2(1:iepts)=aimag(epstot(1:iepts))
        !then call the KK transform, the transformed data is passed back in eps1
        slog='Kramers-Kronig transforms...' 
        call wlog(slog)
        call kk(omega,eps2,iepts,eps1)
        !since we have lost confidence in FPRIME/DANES we want to have
        !OC's computed from abs. cross section via KK tranform. Do it here.
        do k=1,iepts
          eps(k)=eps1(k)+coni*eps2(k)+efree(k)
        enddo
        str='edges/'
        str2=cmpnm(1)
        call wlog('reading header from pot.bin for this comp: '//str2)
        call concat(str,str2,str3,k)
        str='/'
        call concat(str3,str,str2,k)
        str=alledg(1,1)
        call concat(str2,str,path,k)
        flnm='opconsKK.dat'
        call opcons(omega, eps, iepts,flnm, path)

!**************** opcons0.dat
        do k=1,iepts
          eps2(k)=aimag(epstot0(k))
        enddo
        slog='Kramers-Kronig transforms (atomic background)...' 
        call wlog(slog)
        call kk(omega,eps2,iepts,eps1)
        do k=1,iepts
          eps(k)=eps1(k)+coni*eps2(k) +efree(k)
        enddo
        !write out background to file
        call wlog('path: '//path)
        flnm='opcons0.dat'
        call opcons(omega, eps, iepts,flnm, path)

!**************** sumrules.dat
        slog='Calculating sumrules...' 
        call wlog(slog)
!       find the minimum number density for all components
!       we take this to be the number density of chemical formulas
        numsum=numden(1)
        do i=2,ncomps
          if (numsum.gt.numden(i)) numsum=numden(i)
        enddo 
        file1='opconsKK.dat'
        file2='sumrules.dat'
        call sumrules(numsum,file1,file2)

!**************** hamaker.dat
        if (dohamaker) then  !switch added Kevin Jorissen b/c hamaker crashes sometimes
!       !debugging test of hamaker const. code; computes the hamaker
!       !transform of a step-function epsilon; hamakertest.gnuplot plots
!       !the results against the analytic answer.
!       do k=1,iepts
!         if(omega(k).ge.1000.0.and.omega(k).le.1200.0) then
!           eps(k)=1.0+coni
!         else
!           eps(k)=0.0
!         end if
!       end do

        call hamaker(omega,eps,iepts,epsi)
        open (unit=89,file="hamaker.dat")
        do k=1,iepts
          write(unit=89,fmt="(20e20.10)") omega(k),epsi(k)
        enddo
        close(unit=89)
        else
        call wlog('Hamaker calculation skipped - change dohamaker in step4.f to switch it back on.')
        endif !dohamaker

      else
        call wlog('Skipping computation of optical constants')
      end if
        
      close (unit=11) !log file logfullspectrum.dat

  400 call par_barrier
      call par_end
	  if(master)call WipeErrorfileAtFinish

      stop
      end

