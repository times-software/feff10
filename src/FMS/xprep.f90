!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xprep.f90,v $:
! $Revision: 1.15 $
! $Author: jorissen $
! $Date: 2012/04/20 22:33:36 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xprep(iph0, idwopt, nat, inclus, npot,      &
     &            iphat, rmax, rat,                    &
     &            izx, rnrmav, temper, thetad, sig2,   &
     &            minv, rdirec )
  use IOMod
  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx, legtot, natxx, nclusx
  use par
  use constants

  use rotx
  use lnlm
  use xstruc
  use t3j
  
! Used to do the DW factors from the dynamical matrix
  use m_DMDW

  implicit none
  !--------------------------------------------------------------------
  ! This subroutine prepares geometrical information for a subsequent
  ! call to the fms full multiple scattering package.  Information is
  ! passed to fms via common blocks (which live in header files during
  ! development).
  !   !!!!!!!!!!!!!!!!!!KJ I KILLED THE COMMON BLOCKS AND .H FILES - GONE!!!!!
  ! These header files are required:  dim.h xparam.h  xstruc.h
  ! dim.h and xparam.h must be included in the calling routine
  ! This routine calls wlog, so must be compiled with it
  !
  ! This is the main file of xpreppack for use with the full multiple
  ! scattering package (fmspack).  The calling protocol for xpreppack
  ! and fmspack is;
  !
  !      include '../HEADERS/dim.h'
  !      include 'xparam.h'
  !      ...
  !      call xprep(iph0, idwopt, nat, inclus, npot, iphat, rmax, rat,
  ! $        izx, rnrmav, temper, thetad)
  !      energy loop {
  !         ...
  !         call fms(nsp, inclus, npot, ck, lipotx, xphase, ik, iverb, gg)
  !         ... }
  !
  ! xpreppack contains the following routines:
  !   xprep:  main routine of xpreppack
  !   getang: determine angles between the z axis and all pairs of atoms
  !   rotxan: get all rotation matrix elements for the cluster
  !   rotint: initialize arrays used in the construction of rotation
  !           matrices
  !   sortat: organize atoms and potentials lists for computational and
  !             organizational efficiency
  !   atheap: heap sort extended cluster by distance from central atom
  !   xanlm:  get all legendre normalization factor
  !   factst: part of legendre factor computation
  !
  ! xpreppack currently supports use of the Debye-Waller factors for
  !   estimating the effect of thermal disorder on the xanes spectrum.
  !   It does this by filling a matrix with the pairwise mean square
  !   displacement between atoms.  Other forms of this calculation may
  !   be included in the future.  Note that it is strictly impossible
  !   to correctly model disorder in the MS scattering contribution to
  !   the spectrum when using the FMS technique. 
  !--------------------------------------------------------------------
  !  input:
  !     iph0:   potential index for DOS calculations (added by ala to
  !             handle other-than-the-central atom) (iph0=0 for the
  !             absorbing atom)
  !     nat:    number of atoms in extended cluster
  !     npot:   number of unique potentials in cluster
  !     iphat:  (natxx) potential index for each atom in extended
  !             cluster
  !     rmax:   radial size of cluster
  !     rat:    (3, natxx) coordiantes of each atom in extended cluster
  !             as read from geometry file
  !  input for correlated debye model:
  !     izx:    (natxx) Z number of each atom in the cluster
  !     rnrmav: average norman radius in cluster (from pahse.bin)
  !     temper: sample temperature
  !     thetad: Debye temperature
  !
  !  output:
  !     inclus: number of atoms in cluster (inclus <= nat)
  !
  !  output (all via commmon blocks in xstruc.h):
  !     xphi:  matrix of angles between the z axis and pairs of atoms
  !     xrat:  xyz coordinates of the atoms in the cluster, the first
  !            npot+1 entries are examples of each unique potential
  !     iphx:  potential indeces of each atom in the cluster, ordered
  !            like xrat
  !     drix:  huge matrix containing all rotation matrix elements
  !            needed for computation of free elctron propagators
  !     xnlm:  matrix of legendre polynomial normalization factors
  !--------------------------------------------------------------------

  integer, intent(in) :: iphat(natxx),idwopt,nat,npot, iph0
  real,    intent(in) :: rat(3,natxx), rmax, sig2, rdirec
  integer, intent(in) :: minv

  ! Input for correlated debye model:
  real,    intent(in) :: temper, thetad, rnrmav
  integer, intent(in) :: izx(0:nphasx)

  integer, intent(out) :: inclus

  ! Added to satisfy implicit none
  integer :: iat,iat1,iat2,i,j,k,iph1,iph2 !Loop indecies
  integer :: ipair,l1,mm,isp1! Loop indecies
  integer :: j1,j2,j3p,j3m,m1,m2,icen
  integer :: irm1,irm2,ipath0,lplus1,mplus1,npair
  real    :: rr,xr12,yr12,zr12,rr12,rmaxem,rmax2,xbeta,Tmp
  integer :: iem, ios, irm

  integer :: iphat2(natxx), izpair(0:2)
  real :: rat2(3,natxx)
  real*8 :: ra(natxx)
  character*78 :: line

  !     sigms is written in double precision.  these are the variables
  !     that it uses
  real*8 :: dtemp, dthet, drs, dsigsq, pair(3,0:2)
  real*8 :: sig2mx, sig2x(0:nphx,0:nphx)
  !     iwarn - needed to write warning just one time
  integer iwarn
  save iwarn
  data  iwarn /0/
  logical,parameter :: WarnWhenClusterReducedInRadius = .false.  !Turn off these warnings -- not very helpful to the user anyway

  double precision cwig3j
  external cwig3j

  !     Josh - logical flag for reading in dw factors
  logical readdw
  integer iTmp1, iTmp2
  !     Josh END

! Variables used by the DMDW module
  type(Lanczos_Info) :: Lanc_In
  type(dym_Info)     :: dym_In
! Modified by FDV
! These are not used anymore, substituted by compact output
! real*8, dimension(:), allocatable :: DMDW_sig
! real*8, dimension(:), allocatable :: DMDW_vfe
! real*8, dimension(:), allocatable :: DMDW_mef
  type(DW_Out_Info)                 :: DW_Out
  integer            :: DMDW_IO_Flag_Save
  type(Error_Info)   :: Err

! For debug purposes
! type(Paths_Info)   :: Paths_In

  real*8 ratdw(3,0:legtot) !KJ 11-2011 to fix sloppy procedure calls

  if(.not.WarnWhenClusterReducedInRadius) iwarn=1
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
!   call Read_Paths_Info(Paths_In)

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
!   call Print_Header(Lanc_In,Paths_In,dym_In)

! Close the dmdw.inp file
    call DMDW_Close_I

! Allocate the auxiliary variable
!   allocate(DMDW_sig(Lanc_In%nT), &
!            DMDW_vfe(Lanc_In%nT), &
!            DMDW_mef(Lanc_In%nT))

! Debug
!   stop

  end if

  !  initialize geometrical arrays
  do i=1,nclusx
     do j=1,nclusx
        xphi(j,i) = zero
     enddo

     do j=1,3
        xrat(j,i) = zero
     enddo

     iphx(i) = 0
  enddo
  inclus = 0



  ! --- find the central atom, ipot=iph0 (iph0=0 for the absorbing atom)
  icen = 0
  do i=1,nat
     iphat2(i) = iphat(i)
     if (iphat(i).eq.iph0) then
        if (icen.eq.0) then
           icen = i
        elseif (iph0.eq.0) then
           call wlog('* * * ERROR!  More than one atom in the extended cluster has ipot=0')
           call wlog('      You may only have one central atom.')
           call wlog('      Stopping in xprep.')
           call par_stop('XPREP-1')
        endif
     endif
  enddo
  ! --- make sure central atom is at (0,0,0)
  do i=1,nat
     rat2(1,i) = rat(1,i)-rat(1,icen)
     rat2(2,i) = rat(2,i)-rat(2,icen)
     rat2(3,i) = rat(3,i)-rat(3,icen)
  enddo
  ! --- sort the atoms from extended cluster by distance from central
  !     atom.
!KJ  call atheap(nat, rat2, iphat2, ra) !KJ 7-09 disabled b/c inaccurate and annoying and redundant!!

  ! --- define cluster from extended cluster by as those closer than
  !     rmax to central atom
  inclus=0
  rmax2 = rmax**2
  NAT_LOOP: do i=1,nat
     rr = (rat2(1,i)**2 + rat2(2,i)**2 + rat2(3,i)**2)
     if (rr .gt. rmax2) then
        inclus = i-1
        exit NAT_LOOP
     endif
  enddo NAT_LOOP
  if (inclus.eq.0) inclus=nat

  ! --- sanity check size of cluster
  if (inclus.gt.nclusx) then
     if (iwarn.eq.0) then
        call wlog('* * * WARNING preparing cluster for FMS calculation.')
        write(line,400) inclus
400     format('      You specified a cluster of ', i3, ' atoms for the FMS calculation.')
        call wlog(line)
        write(line,410)nclusx
        call wlog(line)
410     format('      This exceeds the hard wired limit of ', i3, ' atoms.')
        write(line,420)nclusx
        call wlog(line)
420     format('      The cluster size was reset to ', i3, ' and the calculation will continue.')
        iwarn = 1
     endif
     inclus = nclusx
  endif

  ! --- make the first few entries in xrat represent each of the
  !     unique potentials, sorting around the chosen center
  !     (iph0=0 for the absorbing atom)
  !     call sortat(iph0, inclus, npot, iphat2, iphx, rat2, xrat)
  do iat = 1, inclus
     iphx(iat) = iphat2(iat)
     xrat(1,iat) = real (rat2(1,iat))
     xrat(2,iat) = real (rat2(2,iat))
     xrat(3,iat) = real (rat2(3,iat))
  enddo


  ! --- Calculate and store rotation matrix elements and phi angles
  !     the k loop calculates the forward then the backward rotation
  !     for an atom pair (ij). k = 0-->forward, 1-->backward
  call rotint
  lplus1 = lx+1
  mplus1 = lx+1
  INCLUSI_LOOP: do i=1,inclus
     INCLUSJ_LOOP: do j=1,inclus
        rr = (xrat(1,i)-xrat(1,j))**2 + (xrat(2,i)-xrat(2,j))**2 + (xrat(3,i)-xrat(3,j))**2
        !         if (rr.gt. rdirec**2) goto 140

        call getang(nclusx, xrat, i, j, xbeta, xphi(i,j))
        if (i.eq.j) cycle INCLUSJ_LOOP
        do k=0,1
           if (k.eq.1) xbeta = (-1) * xbeta
           call rotxan(lplus1, mplus1, xbeta, i, j, k)
        enddo
     enddo INCLUSJ_LOOP
  enddo INCLUSI_LOOP


  ! --- calculate spherical harmonic normalization factors
  call xanlm(lplus1,mplus1)

  ! --- calculate array of correlated debye waller factors
  !     initialize
  do iat2=1,nclusx
     do iat1=1,nclusx
        sigsqr(iat1,iat2) = zero
     enddo
  enddo

  iem = 111
  irm1 =113   
  irm2 = 112

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

  !     open files for sigrm and sigem
  if (idwopt.eq.1.and.master) then
     open(unit=iem,file='s2_em.dat',status='unknown', iostat=ios)
     call chopen (ios, 's2_em.dat', 'sigem')
  endif
  if (idwopt.ge.1.and.idwopt.lt.4) then
     if(master) then
        open(unit=irm1,file='s2_rm2.dat',status='unknown', iostat=ios)
        call chopen (ios, 's2_rm2.dat', 'sigrm')
        open(unit=irm2,file='s2_rm1.dat',status='unknown', iostat=ios)
        call chopen (ios, 's2_rm1.dat', 'sigrm')
     endif

     !        initialize statistics for max DW factors and set rmaxem
     sig2mx=0
     do iph1=0,nphx
        do iph2=0,nphx
           sig2x(iph1,iph2) = 0
        enddo
     enddo
     rmaxem = 20.0
     do i = 2,inclus
        rr = 0
        do j=1,3
           rr = rr + (xrat(j,i)-xrat(j,1))**2
        enddo
        rr = sqrt(rr)
        if (rr.lt.rmaxem) rmaxem = rr
     enddo
     rmaxem = max(2.2*rmaxem, 5.0/bohr)
  endif

  npair = 0
  IAT1_LOOP: do iat1=1,inclus-1
     IAT2_LOOP: do iat2=iat1+1, inclus
        rr = (xrat(1,iat1)-xrat(1,iat2))**2 +                         &
             &    (xrat(2,iat1)-xrat(2,iat2))**2 +(xrat(3,iat1)-xrat(3,iat2))**2
        !         if (rr.gt. rdirec**2) goto 240

        if (idwopt.ge.0) then
           do ipair=1,3
              pair(ipair,0) = dble(xrat(ipair, iat1)*bohr)
              pair(ipair,1) = dble(xrat(ipair, iat2)*bohr)
              pair(ipair,2) = dble(pair(ipair,0))
           enddo

           izpair(0) = izx(iphx(iat1))
           izpair(1) = izx(iphx(iat2))
           izpair(2) = izpair(0)
           dtemp = dble(temper)
           dthet = dble(thetad)
           drs   = dble(rnrmav)
           ipath0=0
		   ratdw=0.d0 !KJ
		   ratdw(1:3,0:2)=pair(1:3,0:2) !KJ  ratdw is simply a copy of pair with the correct dimensions to match the dummy argument in the sigrm and sigem routines.
           if (idwopt.eq.0) then 
              ! Use CD model
              call sigms(dtemp,dthet,drs,2,2,pair,izpair,dsigsq)
           elseif (idwopt.eq.1) then 
              xr12 = (xrat(1,iat1) - xrat(1,iat2))**2
              yr12 = (xrat(2,iat1) - xrat(2,iat2))**2
              zr12 = (xrat(3,iat1) - xrat(3,iat2))**2
              rr12 = sqrt(xr12 +yr12 +zr12)
              if (rr12.le.rmaxem) then
                 ! Use EM method
                 npair = npair + 1
                 if (mod(npair,100).eq.0) then
                    write (line, 337) npair
337                 format('    Doing DW factors via EM method for the pair number ', i5)
                    call wlog(line)
                 endif
                 call sigem(sig2mx,sig2x,iem,dtemp,ipath0,2,ratdw,dsigsq) !KJ ratdw replaces pair
              else
                 ! Use RM method
                 call sigrm(sig2mx,sig2x,irm1,irm2,dtemp,ipath0,2,ratdw,dsigsq) !KJ ratdw replaces pair
              endif
           elseif (idwopt.eq.3) then  !KJ 7/06 added this section
              !               use classical model
              call sigcl(dtemp,dthet,drs,2,2,pair,izpair,dsigsq)
           elseif (idwopt.eq.4) then
              ! User defined dw-factors defined in sig2.dat
              CALL ReadData('sig2.dat', Int1 = iTmp1, Int2 = iTmp2, & 
                   &        Real3 = sigsqr(iat1,iat2), Real4 = Tmp)
              dsigsq = DBLE(sigsqr(iat1,iat2))
           elseif (idwopt.eq.5) then
! We want Calc_DW to run silently, so we set the IO flag to 0
              DMDW_IO_Flag_Save = Lanc_In%IOFlag
              Lanc_In%IOFlag = 0
! Changed by FDV
! Updating to new interface
!             call Calc_DW(Lanc_In, dym_In, &
!                  (/ iat1, iat2 /)-1, 2, &
!                  DMDW_sig, DMDW_vfe, 1, DMDW_mef)
              call Calc_DW(Lanc_In, dym_In, &
                   (/ iat1, iat2 /), 2, &
                   0, 0, DW_Out)
!             dsigsq = DMDW_sig(1)
              dsigsq = DW_Out%s2(1)
! Now we reset the flag to its previous value
              Lanc_In%IOFlag = DMDW_IO_Flag_Save
!             print *, 'DW: ', iat1, iat2, dsigsq
           else
              !               use RM
              call sigrm(sig2mx,sig2x,irm1,irm2,dtemp,ipath0,2,ratdw,dsigsq) !KJ ratdw replaces pair
           endif
           sigsqr(iat1,iat2) = real(dsigsq)
           !            Josh - Temporary. write to sigsqr.dat Only works on one
           !            processor.
           !CALL WriteData('sig2FEFF.dat', Int1 = iat1, Int2 = iat2, & 
           !     & Real3 = sigsqr(iat1,iat2), &
           !     & Real4 = REAL(bohr*SQRT(rat2(1,iat2)**2 + &
           !     & rat2(2,iat2)**2 + rat2(3,iat2)**2)))
        endif
        sigsqr(iat1,iat2) = sigsqr(iat1,iat2) + sig2
        sigsqr(iat2,iat1) = sigsqr(iat1,iat2)
     enddo IAT2_LOOP
  enddo IAT1_LOOP

  !     close output for sigem sigrm
  if (master) then
     if (idwopt.eq.1) close (unit=iem)
     if (idwopt.ge.1) close (unit=irm1)
     if (idwopt.ge.1) close (unit=irm2)
  endif

  !     Calculate Clebsch-Gordon coefficients <LS|J>
  do l1 = 0, lx
     do mm = -l1, l1
        do isp1 = 1, 2
           j1 = 2 * l1
           j2 = 1
           j3p = j1 + 1
           j3m = j1 - 1
           m1 = 2*mm
           m2 = 2*isp1 - 3
           !  j = l+1/2
           t3jp( l1, mm, isp1) = sqrt( j3p + 1.0e0 ) *                     &
                &                 real( cwig3j( j1, j2, j3p, m1, m2, 2) )
           if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0)                            &
                &          t3jp( l1, mm, isp1) = - t3jp( l1, mm, isp1)

           !  j = l-1/2
           t3jm( l1, mm, isp1) = sqrt( j3m + 1.0e0 ) *                     &
                &                 real( cwig3j( j1, j2, j3m, m1, m2, 2) )
           if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0)                            &
                &          t3jm( l1, mm, isp1) = - t3jm( l1, mm, isp1)
        enddo
     enddo
  enddo

! If we are using Dynamical Matrix DW factors, do som clean up
  if ( idwopt == 5 ) then

! Allocate the auxiliary variable
!   deallocate(DMDW_sig)

  end if

  return
  !  end subroutine xprep
end subroutine xprep
