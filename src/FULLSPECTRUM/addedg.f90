      subroutine addedg (path,  omega, npts,                            &
     &                   sctfc ,sctfc0, neff, ierr)
      !Returns the contribution to the scattering factor f and its atomic
      !background from the edge whose xmu.dat files are found at the path
      !pathname stored in path. Looks for xmu.dat files at
      !path//'fms_im/xmu.dat', path//'fms_re/xmu.dat',
      !path//'path_im/xmu.dat', path//'fprime/xmu.dat', and  combines
      !them into a single spectrum which is returned in sctfc (full
      !signal including any fine structure stored in the xmu.dat files)
      !and sctfc0 (atomic background).  A smooth transition is made from
      !FMS to path expansion calculations between photoelectron momenta
      !khi, klo. The effective number of electrons (estimated by the
      !eps_2 sumrule) is passed back in neff.

      !inputs:
      !   path -- path to dir containing calculations for a given edge
      !   numden -- number density of this type of absorber for
      !             conversion from atomic to macroscopic quantities
      !   omega -- energy grid for epsilon
      !   npts -- number of points in omega
      !   khi -- highest photo-electron momentum to compute via FMS
      !          (inverse bohrs)
      !   klo -- lowest photo-electron momentum to compute via path
      !           expansion. 
      !   nz  -- atomic number of abs. for calculation of f_0 term in
      !          scattering amplitude.
      !outputs:
      !   sctfc -- array holding complex (scalar) scattering factor
      !   sctfc0 -- atomic background of sctfc
      !   ierr -- error flag: returned 0 if everything is OK, positive
      !           otherwize

!     implicit real (a-h,o-z)
      use constants
      implicit none

      include 'HEADERS/params.h'
      real*8 trans,emin
      real*8 eps1,eps1bk,eps2,eps2bk !these hold real or im part for one freq.
      complex sctfc(fullpts),sctfc0(fullpts) !these hold cmplx sctfc for all omega
      !nebk holds the atomic background from the xanes/danes/exafs calculations.
      complex nebk(fullpts)
      integer npts,i,k,indexi,indexj,j,nz,                              &
     &        calc1_re, calc2_re, ierr,calc1_im, calc2_im,              &
     &        trstart_im(nexts+2),trend_im(nexts+2),                    &
     &        trstart_re(nexts+2),trend_re(nexts+2), istrln
      external istrln
      real*8 ehi(2), omega(fullpts),y,y2,y0,y02,frac,numden, elo(2)
      character*512 path, slog,str, fname
      logical conv,rdok,flag
      real*8 col1(maxpts), col2(maxpts),col3(maxpts),                     &
     &     col4(maxpts), col5(maxpts),                                  &
     &     mu(maxpts,nexts+2), & !array to hold results read from xanes/exafs xmu.dat&
     &     fprime(maxpts,nexts+2), &  !array to hold results read from danes xmu.dat&
     &     mu0(maxpts,nexts+2), & !array to hold atomic background of xanes/exafs&
     &     fprm0(maxpts,nexts+2), & !array to hold atomic background of fprime&
     &     em_re(maxpts,nexts+2), & !array to hold energy grids for data in fprime, fprm0&
     &     em_im(maxpts,nexts+2)  !array to hold energy grids for data in mu, mu0
      real*8 exten(maxpts,nexts),                                         &
     &     scrtch(maxpts),scrtc2(maxpts), theta, f0, neff,              &
     &     slope1,slope2,icept1,icept2,y1,yt1,yt0,x1,x0
      !ilo, ihi index the interval containing the fine structure
      !calculation. indx indexes which method was used to calculate each
      !point in the atomic background of the im. part.
      integer ilo(2),  ihi(2), mult, ovrlap, indx(fullpts)
      parameter (mult = 10)
      integer icalc(fullpts)
      real*8 fp0


!     initialize arrays to 0
      do i=1,maxpts
        do j=1,nexts+2
          em_re(i,j)=0.0
          em_im(i,j)=0.0
          mu(i,j)=0.0
          mu0(i,j)=0.0
          fprime(i,j)=0.0
          fprm0(i,j)=0.0
        enddo
      enddo
      do i=1,fullpts
        sctfc(i)=0.0
        sctfc0(i)=0.0
        nebk(i)=0.0
        indx(i)=0
      enddo

      !read atomic background from fprime calculations
      call rdbkg(path,omega,npts,sctfc0,neff,fp0,ierr)
      !sctfc0 now contains f from fprime calculation on the output grid

      !write out intermediate results as a check
      open (unit=33,file='background.dat')
      do i=1,npts
        write(unit=33,fmt="(5e20.10)") hart*omega(i),sctfc0(i)
      enddo
      close (unit=33)

      call rdst(path,omega,npts, sctfc, nebk, ehi, elo, ierr )
      !returns fine structure on output grid. ehi(1:2), elo(1:2) define the
      !intervals over which fine structure was read for the real and im
      !parts of the scattering factor.

      !write out intermediate results as a check
      open (unit=33,file='fine_st.dat')
      do i=1,npts
        write(unit=33,fmt="(5e20.10)") hart*omega(i),sctfc(i),nebk(i)
      enddo
      close (unit=33)

      !find the begining and end of the fine structure calculations.
      !We will use the endpoints to smooth the transitions to/from
      !FPRIME background calculations to fine structure calculations
      !with XANES, etc.
      do j=1,2
        ilo(j)=0
        ihi(j)=0
      enddo
      do i=1,npts

        !needed to match sign convention used by feff for the
        !scattering factor. see Phys. Rev. B 62, 2437-2445 for formula
        !and comments on the sign.
        sctfc0(i)=conjg(sctfc0(i))
        sctfc(i)=conjg(sctfc(i))
        nebk(i)=conjg(nebk(i))


        !update the indices; j indexes the real and im parts of the
        !spectrum; ihi(j) indexes the highest freq. for which we have a
        !fines structure calculation; ilo(j) gives the start of the
        !fine structure calculation.
        do j=1,2
          if(omega(i).le.elo(j)) then
            ilo(j)=i
          end if
          if(omega(i).le.ehi(j)) then
            ihi(j)=i
          end if
        enddo

      enddo

      !decide how many points should be in the smooth transitions
      !between the various calculations
      ovrlap=5 !use 5 points if loop below does not set ovrlap
      do i=ilo(2),ihi(2)
        if (omega(i).le.elo(2)+trsize) then
          ovrlap=i-ilo(2)+1
        end if
      enddo
!     ovrlap=5
     
  
      if (omega(1).le.elo(1)) then 
        !use background for total signal in interval before edge
        do i=1,ilo(1)
          sctfc(i)=dble(sctfc0(i))+coni*aimag(sctfc(i))
        enddo 
      endif
      if (omega(1).le.elo(2)) then 
        do i=1,ilo(2)
!         sctfc(i)=coni*aimag(sctfc0(i))+real(sctfc(i))
          !as a hack, enforce a sharp step function at the absorption
          !edge since FPRIME calculation lacks broadening, which can
          !screw up the spectrum in the pre-edge region.
          sctfc(i)=dble(sctfc(i)) 
          sctfc0(i)=dble(sctfc0(i)) 
          icalc(i)=1
        enddo 
      endif

      !for the background, smoothly transition from the FPRIME
      !background to the background from XANES/DANES/EXAFS. We use a
      !smaller transition interval than used for the fine structure (use
      !j = ovrlap/5 points.
      j=ovrlap/5 
      if (omega(1).le.elo(1)) then 
        do i=ilo(1),ilo(1)+j
          !set mixing angle for transition
          frac=dble((i-ilo(1)))/j
          theta=frac*pi/2

          !update background with near edge calculation
          sctfc0(i)=coni*aimag(sctfc0(i))+                              &
     &             dble(sctfc0(i))*cos(theta)**2 +                      &
     &             dble(nebk(i)  )*sin(theta)**2   
        enddo
        !use the good background for the remainder of the transition
        !interval 
        do i=ilo(1)+j, ilo(1)+ovrlap
          sctfc0(i)=coni*aimag(sctfc0(i))+ dble(nebk(i))
        enddo 
      endif
      if (omega(1).le.elo(2)) then 
        do i=ilo(2),ilo(2)+j
          !set mixing angle for transition
          frac=dble((i-ilo(2)))/j
          theta=frac*pi/2

          !update background with near edge calculation
          sctfc0(i)=dble(sctfc0(i))+                                    &
     &             (aimag(sctfc0(i))*cos(theta)**2 +                    &
     &              aimag(nebk(i)  )*sin(theta)**2)*coni
          indx(i)=2
          icalc(i)=2
!         write(slog,fmt="(e20.10,i3)") omega(i)*hart, icalc(i)
!         call wlog(slog)
        enddo
        !use the good background for the remainder of the transition
        !interval 
        do i=ilo(2)+j, ilo(2)+ovrlap
          sctfc0(i)=coni*aimag(nebk(i))+ dble(sctfc0(i))
          indx(i)=3
        enddo 
      endif
    

      if (omega(1).le.elo(1)) then 
        !smoothly transition into fine structure calculation.
        do i=ilo(1),ilo(1)+ovrlap

          !set mixing angle for transition
          frac=dble((i-ilo(1)))/ovrlap
          theta=frac*pi/2

          !update full signal with near edge calculation
          sctfc(i)=coni*aimag(sctfc(i))+                                &
     &             dble(sctfc0(i))*cos(theta)**2 +                      &
     &             dble(sctfc (i))*sin(theta)**2   
        enddo
      end if
      if (omega(1).le.elo(2)) then 
        do i=ilo(2),ilo(2)+ovrlap

          !set mixing angle for transition
          frac=dble(i-ilo(2))/ovrlap
          theta=frac*pi/2

          !update background and full signal with near edge calculation
          sctfc(i)=dble(sctfc(i))+coni*(                                &
     &             aimag(sctfc0(i))*cos(theta)**2 +                     &
     &             aimag(sctfc (i))*sin(theta)**2 )  
          icalc(i)=3
        enddo
      endif

      if (ihi(1)-ovrlap.gt.0) then
        !smoothly transition out of fine structure calculation.
        do i=ihi(1)-ovrlap, ihi(1)

          !set mixing angle for this transition
          frac=dble(ihi(1)-i)/ovrlap
          theta=frac*pi/2

          !update background and full signal
          sctfc0(i)=coni*aimag(sctfc0(i))+                              &
     &             dble(sctfc0(i))*cos(theta)**2 +                      &
     &             dble(nebk(i)  )*sin(theta)**2   
          sctfc(i)=coni*aimag(sctfc(i))+                                &
     &             dble(sctfc0(i))*cos(theta)**2 +                      &
     &             dble(sctfc (i))*sin(theta)**2   
        enddo
      end if
      !For the imaginary part, smoothly transition out of fine structure
      !calculation over an interval (at most)  mult times as large as
      !the the other seams in the calculation to hide any discrepency
      !between the tails calculated with FPRIME and EXAFS/DANES.
      !j indexes the begining of this transition
      j=max(ilo(2)+ovrlap,ihi(2)-ovrlap*mult)
      do i=j,ihi(2)

        !set mixing angle
        frac=dble(ihi(2)-i)/(ihi(2)-j)
        theta=frac*pi/2

        !update background and full signal
        sctfc(i)=dble(sctfc(i))+coni*(                                  &
     &           aimag(sctfc0(i))*cos(theta)**2 +                       &
     &           aimag(sctfc (i))*sin(theta)**2 )  
        sctfc0(i)=dble(sctfc0(i))+coni*(                                &
     &           aimag(sctfc0(i))*cos(theta)**2 +                       &
     &           aimag(nebk(i)  )*sin(theta)**2 )  
        indx(i)=5
        icalc(i)=4
      enddo

      !In between the seams we want to use the background calculated by
      !XANES/EXAFS/DANES cards, so load it here.
      if (ilo(1)+ovrlap-1.gt.0) then
        do i=ilo(1)+ovrlap-1, ihi(1)-ovrlap+1
          sctfc0(i)=coni*aimag(sctfc0(i))+ dble(nebk(i))
        enddo
      end if
      if (ilo(2)+ovrlap-1.gt.0) then
        !j indexes the begining of the transition at the end of the
        !calculations for this edge
        do i=ilo(2)+ovrlap, j
          sctfc0(i)=coni*aimag(nebk(i))+ dble(sctfc0(i))
          indx(i)=6
          icalc(i)=5
        enddo
      end if

      !use background for total signal in interval beyound fine
      !structure calculations.
      do i=max(ihi(1),1),npts
        sctfc(i)=dble(sctfc0(i))+coni*aimag(sctfc(i))
      enddo 
      do i=max(ihi(2),1),npts
        sctfc(i)=coni*aimag(sctfc0(i))+dble(sctfc(i))
        indx(i)=7
        icalc(i)=6
      enddo 

      !shift the scaterring factor so that f(0)=0 as required for a
      !insulator. Conducting states will be handeled differently.
      do i=1,npts
        sctfc(i)=sctfc(i)-fp0
        sctfc0(i)=sctfc0(i)-fp0
      enddo



      !diagnostic file
      j=istrln(path)
      call wlog(path(j-1:j+1))
      if (path(j-1:j+1).eq.'/K') then
        slog='tmp_K.dat'
      else
        write (slog,fmt="('tmp_',a2,'.dat')") path(j-1:j+1)
      end if
      call wlog ('writing file: '//slog)
      open (unit=12, file=slog)
      write(12,fmt="('# ',a)") path(1:80)
      do i=1,npts
        write(unit=12,fmt="(5e20.10,' ',i2)")                           &
     &     omega(i),sctfc(i)/alpinv/bohr**2/omega(i),sctfc0(i),indx(i)
      enddo
      close (unit=12)
      

      end
