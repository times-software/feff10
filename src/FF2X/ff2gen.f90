!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ff2gen.f90,v $:
! $Revision: 1.18 $
! $Author: jorissen $
! $Date: 2012/05/15 21:29:59 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,&
     &                   emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv,&
     &                   title, ntitle)

      use dimsmod, only: npx=>npx_ff2x, nheadx, nex, nphx=>nphu, legtot
	  use constants
      implicit double precision (a-h, o-z)

!     header from xsect.dat
      dimension ltitle(nheadx)
      character*80 title(nheadx)
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex)
	  character*1 dummy !KJ 7-09 to get rid of the '#' in header lines
      integer istrln
      external istrln

      open (unit=8, file='xsect.dat', status='old', iostat=ios)  !KJ  7-09 changed xsect.bin to xsect.dat
!     read xsect.dat
      ntitle = nheadx
      call rdhead (8, ntitle, title, ltitle)
!     read method for xsec calculation
      read(8,*)  dummy,s02p, erelax, wp, edgep, emu  !KJ added '#' 7-09
      if (mbconv .gt.0 .or. s02.le.0.1) s02=s02p
!     read gamach (in eV) for use in atan at absorption edge
!     and convert to code units
      read(8,*)  dummy,gamach, ne1, ik0  !KJ added '#' 7-09
      gamach = gamach / hart
!     skip label and read after it
      read(8,*)
      i = 1
  300    read(8,*,end=310)  ereal, eimag, xsnorm(i), dum1, dum2
         xsec(i) = dum1 + coni*dum2
!        xsect.bin is in eV and invA, convert to code units here
         emxs(i) = (ereal + coni*eimag) / hart
         xkxs(i) = getxk (dble(emxs(i)) - edgep)
         omega(i) = dble(emxs(i)) - edgep + emu
         nxsec = i
         i = i + 1
         if (i.le.nex) goto 300
  310 continue
      close(unit=8)
 
      return
      end

      subroutine wrhead (iunit, nhead, head, dwcorr, s02,               &
     &  tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

      use dimsmod, only: nheadx
	  use constants
	  use global_inp, qweights=>qw,qvecs=>qs
	  implicit double precision (a-h, o-z)

      include '../HEADERS/vers.h'
      parameter (eps4 = 1.0e-4)
      character*80  head(nheadx)
      logical dwcorr
      character*2 coment
      parameter (coment='# ')
	  integer iq

!KJ 3-3012  I've "dumbed down" this routine so it always writes the same info to the header,
!           even when a card wasn't used.
!           In particular, Vr and Vi (EXCHANGE), vr and vi (CORRECTIONS), and ThetaE and ThetaD (DEBYE) are now always written.

!     write miscellaneous staff into headers
!     add feff version to the first line
      ll = istrln(head(1))
      if (ll .lt. 0)  then
        head(1)= 'Untitled'
        ll = istrln(head(1))
      endif
      write(iunit,310) coment, head(1)(1:), vfeff
  310 format (a2, a55, t66, a12)

!     the rest of the title
      do 330  ihead = 2, nhead
         ll = istrln(head(ihead))
         if (ll .gt. 0)  then
            write(iunit,320) coment, head(ihead)(1:ll)
         endif
  320    format (a2, a)
  330 continue
!      if (dwcorr)  then
         write(iunit,340)  coment, s02, tk, thetad, sig2g
  340    format (a2,' S02=', f5.3, '  Temp=', f7.2,'  Debye_temp=',f7.2,&
     &        '  Global_sig2=', f8.5)
!      else
!         write(iunit,341)  coment, s02, sig2g
!  341    format (a2, ' S02=', f5.3,                                     &
!     &   '                                        Global_sig2=', f8.5)
!      endif
      if (alphat .gt. zero)  then
         write(iunit,321)  coment, alphat
  321    format (a2, ' 1st and 3rd cumulants, alphat = ', 1pe20.4)
      endif

!      if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
         write(iunit,342) coment, vrcorr*hart, vicorr*hart
!      endif
  342 format (a2, ' Energy zero shift, vr, vi ', 1p, 2e14.5)

      if (critcw .gt. 0)  write(iunit,350) coment, critcw
  350 format (a2, ' Curved wave amplitude ratio filter ', f7.3, '%')
      write(iunit,360) coment
  360 format (a2, '    file         sig2 tot  cw amp ratio   deg',      &
     &        '  nlegs   reff  inp sig2')
	 
	 if (do_nrixs .eq. 0) return
	 !KJ Not sure if the below (using elpty) makes sense in feff91; copied from withnqwithmdff and I don't care enough to test it ...

       if (elpty.eq.0) then 
   370    format (a2, 'The momentum transfer was ', 3f7.3)
          write(iunit,370) coment,(xivec(i),i=1,3)
       else if (qaverage .and. elpty.lt.0.0d0) then
   371    format (a2, 'This was a spherically averaged calculation')
          write(iunit,371) coment
   372    format (a2, 'The magnitude of the momentum transfer was', f7.3)
          write(iunit,372) coment,xivec(3)
       else if ((.not.qaverage) .and. elpty.lt.0.0d0) then
   373    format (a2, 'This was a spherically averaged calculation')
          write(iunit,373) coment
   374    format (a2, 'The momentum transfers and weights were')
          write(iunit,374) coment
   375    format (a2,5f7.3) !KJ
          do iq=1,nq
             write(iunit,375) coment,qvecs(iq,3),qweights(iq)  !KJ 12-2011 had to switch "iq" and "3"
          end do
       else
   376    format (a2, 'The momentum transfers and weights were')
          write(iunit,374) coment
   377    format (a2,4f7.3)
          do iq=1,nq
             write(iunit,375) coment,(qvecs(iq,i),i=1,3),qweights(iq)  !KJ 12-2011 had to switch "iq" and "i"
          end do
       end if 	 
	 
!     stop writing misc. staff to files

      return
      end


      subroutine dwadd (ntotal,nptot,idwopt,ip,index,crit,critcw,sig2g, &
     &  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad,&
     &  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, nkx, xk, xk0,&
     &  xkp, cchi, iabs, nabs, ispec, ipr4, nhead,                      &
     &  head, vrcorr, vicorr, nused)

! Used to do the DW factors from the dynamical matrix
      use m_DMDW

      use dimsmod, only: npx=>npx_ff2x, legtot, nphx=>nphu, nex, nheadx
	  use constants
      implicit double precision (a-h, o-z)

      parameter (eps4 = 1.0e-4)
      character*80  head(nheadx)
      !parameter (npx=1200) !now in dimsmod
!     indices of paths to do, read from list.dat
      dimension ip(npx)

      real sig2u(npx)

      parameter (nfinex = nex*100) ! 
	  !KJ 11-2011 The following 4 arrays were defined as (nfinex).  However the compiler doesn't like this.  Furthermore, they are only used
	  ! over (1:nkx) where nkx is a calling argument; it is (afaik) given as "ne" which one expects to be smaller than "nex" in dimsmod.
	  ! Hence I am changing the definition of these arrays to match that in the calling routine.
	  !KJ 05-2012 Hmm, actually not a good solution.  Problem: nkx=ne in ff2xmu, but nkx=nfinex in ff2chi.  What to do?  --> try setting them to size nkx.
!      complex*16 cchi(nex), ccpath(nex), ccc, ckp
!      dimension xkp(nex), xk0(nex)
!      complex*16 cchi(nfinex), ccpath(nfinex), ccc, ckp
!      dimension xkp(nfinex), xk0(nfinex)
      complex*16 cchi(nkx), ccpath(nkx), ccc, ckp
      dimension xkp(nkx), xk0(nkx)

!     to keep Im part of cchi 11.18.97 ala
      complex*16 dw, dw1, dw3

      logical dwcorr
      dimension rattmp(3,0:legtot)
      dimension iztmp(0:legtot)
      character*512 slog
      character*12 fname
      real rnrmav
      dimension iz(0:nphx)
!     central atom phase shift at l0
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx)
      real achi(nex,npx), phchi(nex,npx)
      dimension sig2x(0:nphx, 0:nphx)
      character*2 coment
      parameter (coment='# ')

! Variables used by the DMDW module
      type(Lanczos_Info)                 :: Lanc_In
      type(dym_Info)                     :: dym_In
! Modified by FDV
! These are not used anymore, substituted by compact output
!     real*8, dimension(:), allocatable  :: DMDW_sig
!     real*8, dimension(:), allocatable  :: DMDW_vfe
!     real*8, dimension(:), allocatable  :: DMDW_mef
      type(DW_Out_Info)                  :: DW_Out
      integer, dimension(:), allocatable :: DMDW_Path
      type(Error_Info)                   :: Err
      integer                            :: iCount
      integer                            :: DMDW_IO_Flag_Save
      real*8, Parameter  :: DMDW_Too_Close = 0.01

! For debug purposes
      type(Paths_Info)   :: Paths_In

!     Keep stats on paths used to make chi
      nused = 0
      xkref = dble(ck(1)**2) - xk(1)*abs(xk(1)) 
! If we are using Dynamical Matrix DW factors, read the information only once
      if ( idwopt == 5 ) then
  
! Open the dmdw.inp file
        call DMDW_Open_I(Err)
  
        if ( Err%Flag ) then
          write(6,fmt='(a)') trim(Err%Message)
          stop
        end if

! Read Lanczos information
        call Read_Lanczos_Info(Lanc_In)

! Read the path info in dmdw.inp, for debug purposes only
!       call Read_Paths_Info(Paths_In)

! Read the Dynamical Matrix
        call Read_dym_Info(Lanc_In%dym_file,dym_In)

! Create the full dynamical matrix
        call Make_DM(dym_In)

! Calculate the transformation to internal coordinates (eliminates
! rotations and translations)
! NOTE: The subroutine is still missing parts of the code and shouldn't be used.
!       If used, the TrfD is deallocated before return and any attempt to use it
!       will result in a segfault.
! NEW NOTE: TrfD works partially now, at least the first 6 vectors (2
!           trans, 3 rot). I use it to project out these modes to make
!           the DW factors more stable.
        call Make_TrfD(dym_In)

! Debug
!       call Print_Header(Lanc_In,Paths_In,dym_In)

! Close the dmdw.inp file
        call DMDW_Close_I

! Allocate the auxiliary variable
!       allocate(DMDW_sig(Lanc_In%nT),&
!                DMDW_vfe(Lanc_In%nT),&
!                DMDW_mef(Lanc_In%nT))

! Debug
!      stop

!             do iileg=1,nleg(ipath)
!               write(6,fmt='(3f14.8)') rattmp(:,iileg)
!             end do
!             write(6,fmt='(a)') 'xxxx'
     end if

!     open the files for sigrm and sigem
      if (idwopt.eq.1) then
         iem = 111
         open(unit=iem,file='s2_em.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_em.dat', 'sigem')
      elseif (idwopt.eq.2) then
         irm1 =111
         open(unit=irm1,file='s2_rm2.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm2.dat', 'sigrm')
         irm2 = 112
         open(unit=irm2,file='s2_rm1.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm1.dat', 'sigrm')
      endif
      if (alphat .gt. 0) then
        icum = 113
        open(unit=icum, file='cum.dat', status='unknown', iostat=ios)
        call chopen (ios, 'cum.dat', 'sig3')
        Write(icum, 363)
  363  format('# first and third icumulant for single scattering paths')
        write(icum,364) thetae, alphat
        write(icum,365)
  364   format ('# Einstein-Temp. =',f9.2 ,'   ', 'alpha=',f9.5)
  365   format ('#       file   sig1    sig2    sig3 ')
      endif

      if (idwopt.ge.1) then
!        initialize statistics for max DW for sigrm
         sig2mx=0
         do 400 iph1=0,nphx
         do 400 iph2=0,nphx
  400    sig2x(iph1, iph2) = 0
      endif

      if(ntotal.gt.0 ) then
            if (idwopt.eq.0) then 
              call wlog('Applying Debye-Waller factors using a Correlated Debye model.')
            elseif (idwopt.eq.1) then
              call wlog('Applying Debye-Waller factors using the Equation-of-Motion method.')
	        elseif (idwopt.eq.2) then
              call wlog('Applying Debye-Waller factors using the Recursion method.')
            elseif (idwopt.eq.3) then  !KJ 7/06 added this section
              call wlog('Applying Debye-Waller factors using the Classical Debye model.')
            elseif (idwopt.eq.4) then  !KJ 7/06 added this section
              call wlog('Applying Debye-Waller factors using the sig.dat file.')
            elseif (idwopt.eq.5) then  ! FDV
              call wlog('Applying Debye-Waller factors using the ab-initio Dynamical Matrix model.')
            endif
      endif

!     cycle over all paths in the list
      do 560  ilist = 1, ntotal
!        find index of path
         do 410  j = 1, nptot
            if (ip(ilist) .eq. index(j))  then
               ipath = j
               goto 430
            endif
  410    continue
         write(slog,420)  ilist, ip(ilist)
  420    format (' did not find path i, ip(i) ', 2i10)
         call wlog(slog)
  430    continue
!        Path ipath is the path from feff.bin that corresponds to
!        the path ilist in list.dat.  The index of the path is
!        ip(ilist) and index(ipath).

!        Use this path if it passes critcw filter
         if (crit(ipath) .lt. critcw)  goto 550

!        do debye-waller factors, get sig2d from correlated debye 
!        model if required
!        A note about units:  sig2g, sig2u() and sig2d are all in
!        Angs**2.  Convert to code units after we've summed them.
         sig2 = sig2g + sig2u(ilist)
         if (dwcorr .and. idwopt.ge.0)  then
!           note that stuff from feff.bin is single precision and
!           mostly in multidim. arrays.  sigms is double precision
!           and its arrays are dimensioned for a single path, so
!           use tmp variables to call it.  tk, thetad and sig2d are 
!           all dp, and therefore OK.  Also note that sigms takes
!           inputs in angstroms, except for rs which is in bohr.
            rs = rnrmav
            do 460  ileg = 1, nleg(ipath)
               iztmp(ileg) = iz(ipot(ileg,ipath))
               do 450  j = 1, 3
                  rattmp(j,ileg) = rat(j,ileg,ipath) * bohr
  450          continue
  460       continue

            iztmp(0) = iztmp(nleg(ipath))
            do 470  j = 1,3
               rattmp(j,0) = rattmp(j,nleg(ipath))
  470       continue
            if (idwopt.eq.0) then 
!             use CD model
              call sigms (tk, thetad, rs, legtot, nleg(ipath),rattmp, iztmp, sig2d)
              !call wlog('Applying Debye-Waller factors using a Correlated Debye model.')
! Debug (FDV)
!             print *, 'sig2d: ', sig2d
            elseif (idwopt.eq.1) then 
!             use EM method
              call sigem(sig2mx,sig2x,iem,tk,ipath,nleg(ipath),rattmp,sig2d)
              !call wlog('Applying Debye-Waller factors using an Einstein model.')
	        elseif (idwopt.eq.2) then
			! use RM method
              call sigrm (sig2mx,sig2x,irm1,irm2,tk,ipath,nleg(ipath),rattmp,sig2d)
              !call wlog('Applying Debye-Waller factors using the Recursion Method.')
            elseif (idwopt.eq.3) then  !KJ 7/06 added this section
!             use CL model
              call sigcl (tk, thetad, rs, legtot, nleg(ipath),rattmp, iztmp, sig2d)   
              !call wlog('Applying Debye-Waller factors using the CL model.')
            elseif (idwopt.eq.5) then  ! FDV
              !call wlog('Applying Debye-Waller factors using the ab-initio Dynamical Matrix model.')
!             Use DMDW approach
!             write(6,fmt='(a,3f14.8)') 'tk, thetad, rs', tk, thetad, rs
!             write(6,fmt='(a,2i6)') 'ipath, legtot', ipath, legtot
!             write(6,fmt='(a,2i6)') 'ipath, nleg(ipath)', ipath, nleg(ipath)
! Find the indices of the path atoms in the dynamical matrix structure
              allocate(DMDW_Path(0:nleg(ipath)-1))
! Debug
!             write(6,fmt='(3F20.10)'), transpose(dym_In%xyz(1:5,:))
!             print *, 'rattmp'
!             print *, rattmp(:,0:nleg(ipath)-1)/bohr
!             stop
              do iileg=0,nleg(ipath)-1
                iCount = 0
                do iAt=1,dym_In%nAt
!                 print *, dym_In%xyz(iAt,:)
                  D2_iAt_iileg = sqrt(sum((dym_In%xyz(iAt,:) - rattmp(:,iileg)/bohr)**2))
                  if ( D2_iAt_iileg < DMDW_Too_Close ) then
!                   write(6,fmt='(2i4)') iileg, iAt
! Modified by FDV
! The format inside Calc_DW requires the paths to be actual atomic indices
! starting from 1, not 0. So we use that here.
!                   DMDW_Path(iileg) = iAt - 1
                    DMDW_Path(iileg) = iAt
                    iCount = iCount + 1
                  end if
                end do
                if ( iCount > 1 ) then
                  write(6,fmt='(a)') &
                        'ERROR: DMDW found more than one possible atom'
                  write(6,fmt='(a,i5)') 'Path:', ipath
                  write(6,fmt='(a,i5)') 'Leg: ', iileg
                  stop
                end if
                if ( iCount == 0 ) then
                  write(6,fmt='(a)') &
                        'ERROR: DMDW could not match the path structure with the one in the dym file.', &
                        '       Please check that the input and dym files are a matching pair.'
                  write(6,fmt='(a,i5)') 'Path:', ipath
                  write(6,fmt='(a,i5)') 'Leg: ', iileg
                  stop
                end if
!               write(6,fmt='(3f14.8)') rattmp(:,iileg)
              end do
! Debug (FDV)
!             print *, DMDW_Path
!             write(6,fmt='(a)') 'xxxx'
! We want Calc_DW to run silently, so we set the IO flag to 0
              DMDW_IO_Flag_Save = Lanc_In%IOFlag
              Lanc_In%IOFlag = 0
! Changed by FDV
! Updating to new interface
!             call Calc_DW(Lanc_In, dym_In, DMDW_Path, nleg(ipath), &
!                          DMDW_sig, DMDW_vfe, 1, DMDW_mef)
! Debug (FDV)
!             print *, 'Before Calc_DW'
              call Calc_DW(Lanc_In, dym_In, DMDW_Path, nleg(ipath), &
                           0, 0, DW_Out)
! Debug (FDV)
!             print *, 'After Calc_DW'
!             sig2d = DMDW_sig(1)
              sig2d = DW_Out%s2(1)
! Now we reset the flag to its previous value
              Lanc_In%IOFlag = DMDW_IO_Flag_Save
              deallocate(DMDW_Path)
            else 
			  call wlog('Unknown Debye Waller method requested in dwadd.  Exiting.')
			  stop 'ERROR'
            endif
            sig2 = sig2 + sig2d
         endif
         sig2 = sig2 / (bohr**2)

!        Do first and third cumulants
         sig1 = 0
         sig3 = 0
         if (alphat .gt. zero  .and. nleg(ipath) .eq. 2)  then
           if (thetae.le.0.d0) then
!            call sig3  to get sig1 and sig3 for single scattering paths
!           use reff(ipath) for r, note that reff is single precision
             iz1 = iztmp(nleg(ipath))
             iz2 = iztmp(1)
             call sigte3(iz1, iz2, sig2, alphat, thetad, reff(ipath), sig1, sig3)
           else
!            this gets sig1 and sig3 for single scattering paths
!            using Morse potential
             call sigm3(sig1, sig2, sig3, tk, alphat, thetae)
           endif
           write(icum,475) index(ipath),  sig1 * bohr, sig2*(bohr**2), sig3*(bohr**3)
  475      format( i10,f9.5,f9.5,' ',f9.7)
         endif

!        put the debye-waller factor and other cumulants into 
!        achi and phchi
         if (mbconv .gt. 0) s02 = 1.0
         do 480  i = 1, ne1
            dw = exp(-2 * sig2 * ck(i)**2)
            dw1 = exp (2 * coni * ck(i) * sig1)
            dw3 = exp ((-4 * coni * ck(i)**3 * sig3) / 3)
            dw = dw * dw1 * dw3
            phdw = 0.0
            if (abs(dw).gt.0) phdw = atan2 (dimag(dw), dble(dw))
            achi(i,ipath) = achi(i,ipath) * abs(dw) * s02 * deg(ipath)
            phchi(i,ipath) = phchi(i,ipath) + phdw
  480    continue
!        make sure no 2pi jumps in phase
         do 490  i = 2, ne1
!           phchi is single precision, so use tmp variables
            curr = phchi (i, ipath)
            old = phchi (i-1, ipath)
            call pijump (curr, old)
            phchi (i, ipath) = curr
  490    continue
         do 500  ik = 1, nkx
            call terp1 (xk, achi(1,ipath),  ne1, xk0(ik), achi0)
            call terp1 (xk, phchi(1,ipath), ne1, xk0(ik), phchi0)
            ccpath(ik) =                                                &
     &        achi0 * exp (coni * (2 * xk0(ik) * reff(ipath) + phchi0))
!           note that this already includes s02, deg, sig2, etc.
!           sum total complex chi
            cchi(ik) = cchi(ik) + ccpath(ik)
  500    continue
         nused = nused + 1
         if (iabs.eq.nabs) then
!           Put path into chi.dat, xmu.dat as required
            if (abs(sig2u(ilist)) .gt. 0.000001)  then
              write(3,515)  coment, index(ipath), sig2*(bohr**2),       &
     &          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr, &
     &          sig2u(ilist)
              write(8,515)  coment, index(ipath), sig2*(bohr**2),       &
     &          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr, &
     &            sig2u(ilist)
            else
              write(3,515) coment, index(ipath), sig2*(bohr**2),        &
     &          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
              write(8,515) coment, index(ipath), sig2*(bohr**2),        &
     &          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
            endif
  515       format(a2, 1x, i10, 5x, f9.5, 2f10.2, i6, f9.4, f9.5)
         endif

!        write out a chinnnn.dat for this path, if necessary.
         if (ipr4 .ge. 2 .and. iabs.eq.nabs .and. ispec.eq.0)  then
!           make filename chipnnnn.dat
            write(fname,520)  index(ipath)
  520       format('chip', i4.4, '.dat')
            open (unit=9, file=fname, status='unknown',iostat=ios)
            call chopen (ios, fname, 'ff2chi')
            do 530  ihead = 1, nhead
               lhead = istrln(head(ihead))
               if (lhead .gt. 0)  then
                  write(9,320) head(ihead)(1:lhead)
  320             format (a)
               endif
  530       continue
            if (dwcorr)  then
               write(9,340)  s02, tk, thetad, sig2g
  340          format (' S02', f7.3, '  Temp', f8.2,'  Debye temp',f8.2,&
     &        '  Global sig2', f9.5)
            else
               write(9,341)  s02, sig2g
  341          format (' S02', f7.3,                                    &
     &      '                                        Global sig2', f9.5)
            endif
            if (alphat .gt. zero)  then
               write(9,321)  alphat
  321          format (' 1st and 3rd cumulants, alphat = ', 1pe20.4)
            endif

            if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
               write(9,342)  vrcorr, vicorr
  342          format (' Energy zero shift, vr, vi ', 1p, 2e14.5)
            endif
            write(9,*) 'Debye-waller factor ', sig2, sig3

            write(9,610)
  610       format (1x, 71('-'))
            write(9,535)
  535       format ('       k         chi           mag          ',     &
     &              'phase        phase-2kr  @#')
            do 540  i = 1, nkx
               ckp = sqrt (xkp(i)*abs(xkp(i)) + xkref)
!              it would be better to use interpolation for ckp
!              fix later if complaints about chipnnn.dat files, ala
               xlam0 =  - dimag(ckp)
               ccc = ccpath(i) * exp(2 * reff(ipath) * xlam0)
               phase = 0
               if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
               if (i .gt. 1)  call pijump (phase, phase0)
               phase0 = phase
               write(9,630)  xkp(i)/bohr, dimag(ccc), abs(ccc), phase,  &
     &                       phase-2*xk0(i)*reff(ipath)
  630          format (1x, f10.4, 3x, 4(1pe13.6,1x))
  540       continue
            close (unit=9)
         endif

  550    continue
  560 continue

!  if(allocated(DMDW_sig)) then ; deallocate(DMDW_sig) ; endif
!     close files opened for sigem and sigrem
      if (idwopt.eq.1) then
        close (unit=iem)
      elseif (idwopt.eq.2) then
        close (unit=irm1)
        close (unit=irm2)
      endif
      if (alphat .gt. 0) then
        close (unit=icum)
      endif

      return
      end
