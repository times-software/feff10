!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: yprep.f90,v $:
! $Revision: 1.6 $
! $Author: jorissen $
! $Date: 2011/11/18 02:13:31 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine yprep(iph0, nat, inclus, npot, iphat, rmax, rat)
  !    yprep is the same as xprep for negative idwopt
  !    simlifies calls in SCF and LDOS where DW factors should not enter
  use DimsMod, only: nphx=>nphu, istatx, nspx=>nspu, lx, nphasx, legtot, natxx, nclusx
  use constants

  use rotx
  use lnlm
  use xstruc

  implicit none

  integer, intent(in) :: iph0,nat
  integer, intent(in) :: iphat(natxx)
  real,    intent(in) :: rat(3,natxx),rmax

  integer, intent(out) :: inclus
  
  ! These are not used in this routine; are the intents correct???
  integer, intent(in) :: npot

  integer :: iphat2(natxx), izpair(0:2)
  real    :: rat2(3,natxx)
  real*8  :: ra(natxx)

  character*78 line

  ! Sigms is written in double precision.  these are the variables
  ! that it uses
  real*8 :: dtemp, dthet, drs, dsigsq, pair(3,0:2)
  real*8 :: sig2mx, sig2x(0:nphx,0:nphx)
  ! iwarn - needed to wrtite waqrning just one time
  integer iwarn
  save iwarn
  data  iwarn /0/
  logical,parameter :: WarnWhenClusterReducedInRadius = .false.  !Turn off these warnings -- not very helpful to the user anyway

  ! Added to satisfy implicit none:
  integer :: i,j,k,iat,iat1,iat2  ! Loop indecies
  integer :: lplus1,mplus1,icen
  
  real :: rr,rmax2,xbeta

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
  if(.not.WarnWhenClusterReducedInRadius) iwarn=1

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
           call par_stop('YPREP-1')
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
  call atheap(nat, rat2, iphat2, ra) !KJ 7-09 disabled b/c inaccurate and annoying and redundant

  ! --- define cluster from extended cluster by as those closer than
  !     rmax to central atom
  inclus=0
  rmax2 = rmax**2
  NAT_LOOP: do i=1,nat
     rr = (rat2(1,i)**2 + rat2(2,i)**2 + rat2(3,i)**2)
     if (rr.gt.rmax2) then
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
400     format('      You specified a cluster of ', i3,' atoms for the FMS calculation.')
        call wlog(line)
        write(line,410)nclusx
        call wlog(line)
410     format('      This exceeds the hard wired limit of ', i3,' atoms.')
        write(line,420)nclusx
        call wlog(line)
420     format('      The cluster size was reset to ', i3,' and the calculation will continue.')
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
  INCLUSI_LOOP: do  i=1,inclus
     INCLUSJ_LOOP: do j=1,inclus
        rr = (xrat(1,i)-xrat(1,j))**2 + (xrat(2,i)-xrat(2,j))**2      &
             &       + (xrat(3,i)-xrat(3,j))**2
        !         if (rr.gt.rdirec**2) goto 140

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

  do iat2=1,nclusx
     do iat1=1,nclusx
        sigsqr(iat1,iat2) = zero
     enddo
  enddo


  return
end subroutine yprep
