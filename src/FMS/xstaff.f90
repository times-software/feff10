!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xstaff.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine atheap(nat, rat, iphat, ra)

  !--------------------------------------------------------------
  !  copyright 1993 university of washington         bruce ravel
  !  modified by alexei ankudinov in march 1999
  !--------------------------------------------------------------
  
  implicit none
  !-------------------------------------------------------------------
  !  heapsort adapted from numerical recipes.  sort atoms by distance.
  !  all the pesky little do loops are for transferring rows
  !  of temp into toss.
  !-------------------------------------------------------------------
  !  alexei ankudinov: needed to avoid unnecessary permutations when atoms
  !  are at the same distance from the central atom, in order to comply 
  !  feff document: the sample atom should be the nearest to absorber or
  !  first in the list amongequidistant
  !  Add small contribution 10**-8 * number to the sorting variable ra
  !  in order to achieve this.
  !-------------------------------------------------------------------
  !  natx:   dimension parameter from calling program
  !-------------------------------------------------------------------
  integer, intent(in) :: nat
  real,    intent(inout) :: rat(3,nat)
  integer, intent(inout) :: iphat(nat)
  real*8,    intent(out) :: ra(nat)

  real :: toss(3)
  real*8 ::dum

  integer :: index,i,j,l,ir,itoss

  if (nat.lt.2) return

  l=0
  do i=1,nat
     ra(i) = dble(rat(1,i))**2 + dble(rat(2,i))**2 + dble(rat(3,i))**2 + dble(i)*1.d-6
     ! Small addition at to prefer the old ordering
!     write(*,'(a,i5,1x,f20.10,1x,L2,4(1x,f20.10))') 'ATHEAP - i,ra,rat',i,ra(i),ra(i).lt.ra(max(i-1,1)),rat(1:3,i)
     if (l.eq.0 .and.i.gt.1) then
        if (ra(i).lt.ra(i-1)) l=1
     endif
  enddo

!  write(*,*) 'ATHEAP -l=',l 
  ! Check if array is already in order
  if (l.eq.0) return

  l  = nat/2+1
  ir = nat

110 continue
  if (l.gt.1) then
     l = l-1
     do index=1,3
        toss(index)=rat(index,l)
     enddo
     itoss = iphat(l)
     dum = ra(l)
  else
     do index=1,3
        toss(index) = rat(index,ir)
     enddo
     itoss = iphat(ir)
     dum = ra(ir)
     do index=1,3
        rat(index,ir) = rat(index,1)
     enddo

     iphat(ir) = iphat(1)
     ra(ir) = ra(1)
     ir=ir-1
     if (ir.eq.1) then
        do index=1,3
           rat(index,1)=toss(index)
        enddo
        iphat(1) = itoss
        ra(1) = dum
        !              sort is finished
        goto 300
     endif
  endif
  i=l
  j=l+l

160 if (j.le.ir) then
     if (j.lt.ir) then
        if ( ra(j) .lt. ra(j+1) ) then
           j  = j + 1
        endif
     endif

     if ( dum .lt. ra(j) ) then
        do index=1,3
           rat(index,i) = rat(index,j)
        enddo
        iphat(i) = iphat(j)
        ra(i) = ra(j)
        i=j
        j=j+j
     else
        j=ir+1
     endif
     goto 160
  endif

  do index=1,3
     rat(index,i) = toss(index)
  enddo
  iphat(i) = itoss
  ra(i) = dum

  goto 110
300 continue

  return
end subroutine atheap

subroutine getang(nclusx, rat, i, j, theta, phi)

  !------------------------------------------------------------------
  !  determine theta and phi polar angles of the vector between two
  !  atom positions
  !
  !  inputs
  !    rat:   (3,nclusx) x,y,z of all atoms in cluster
  !    i, j:  indices of atoms at ends of vector Ri-Rj
  !
  !  outputs
  !    theta: polar angle theta of vector Ri-Rj
  !    phi:   polar angle phi of vector Ri-Rj
  !------------------------------------------------------------------
  use constants, only: zero


  implicit none
  integer, intent(in) :: nclusx,i,j
  real, intent(in) :: rat(3,nclusx)
  real, intent(out) :: theta, phi

  real :: x,y,z,r

  real,parameter:: tiny=1.e-7, pi=3.141592654

  x = rat(1,i) - rat(1,j)
  y = rat(2,i) - rat(2,j)
  z = rat(3,i) - rat(3,j)
  r = sqrt(x**2 + y**2 + z**2)

  !  this fails to calculate phi correctly for, as an example,
  !  x=0.5e-7 and y=2e-7.  However, those numbers are below the
  !  precision of the numbers stored in potph.bin.

  phi = zero
  theta  = zero
  if (i.ne.j) then
     !           phi = atan2(y,x)
     !        all of these conditionals will do the work for a machine that
     !        cannot correctly handle a zero value for the second argument
     !        of atan2
     if (abs(x).lt.tiny) then
        if (abs(y).lt.tiny) then
           phi = zero
        elseif (y.gt.tiny) then
           phi = pi/2
        else
           phi = -pi/2
        endif
     else
        phi = atan2(y,x)
     endif
     if (r.gt.tiny) then
        if (z.le.-r) then
           theta = pi
        elseif ( z.lt.r) then
           theta = acos(z/r)
        endif
     endif
  endif

  return
end subroutine getang


!====================================================================
subroutine rotxan (lxp1, mxp1, betax, i, j, k)
  use DimsMod, only: lx
  use constants, only: zero

  use rotx
  use lnlm
  use xstruc

  implicit none
  !     input:  lxp1, mxp1: lmax+1 & mmax+1, largest L states in matrix
  !             betax is the rotation angle
  !             i and j are the indeces of the atoms, thus denote
  !                 which pair of atoms this is the rotation matrix for
  !             k=0 for forward rotation, k=1 for backward rotation
  !     output: drix(L,k,j,i) in common /rotx/
  !+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
  !     adapted by BR from rot3i, version for genfmt by SIZ
  !        new data structure for rotation matrices to accomodate
  !        xanes calculation
  !+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
  !     subroutine rot3 calculates rotation matrices for l = 0,lxp1-1

  !     subroutine rot3 calculates the beta dependence of rotation
  !     matrix elements using recursion of an iterated version of
  !     formula (4.4.1) in edmonds.
  !
  !     first written:(september 17,1986) by j. mustre
  !     version 2  (17 sep 86)
  !     version 3  (22 feb 87) modified by j. rehr
  !     version for genfmt, modified by s. zabinsky, Sept 1991
  !     Initialized dri0.  Some elements may be used before being
  !        initialized elsewhere -- rot3i needs to be carefully
  !        checked.  S. Zabinsky, April 1993
  !
  !******************** warning****************************************
  !     lxx must be at least lxp1 or overwriting will occur
  !     nmax must be at least nm or overwriting will occur
  !--------------------------------------------------------------------
  !     notation dri0(l,m,n) =  drot_i(l'm'n')
  !     l = l'+1, n' = n-l, m' = m-l, primes denoting subscripts
  !     thus dri0(1,1,1) corresponds to the rotation matrix with
  !     l' = 0, and n' and m' = 0; dri0(3,5,5) : l' = 2,n' = 2,m' = 2.
  !--------------------------------------------------------------------

  integer, intent(in) :: lxp1,mxp1,i,j,k
  real, intent(in)    :: betax


  !      needed for commented out diagnostic file
  !      logical open

  ! Added to satisfy implicit none
  integer :: isav,il,in,im,l,lm,ln,m,m1,m2,mmx,n,ndm,nm
  real :: t,t1,t3,f1,f2,f3,dlnm
  real :: xc,xs,s

  complex dum
  integer, parameter :: lxx=24
  real, parameter :: pi=3.1415926535897932384626433e0, one=1
  complex, parameter :: coni = (0,1)

  !     dri0 is larger than needed for genfmt, but necessary for
  !     this calculation algorithm.  Copy result into smaller
  !     dri arrays (in common) at end of this routine.
  real :: dri0(lxx+1, 2*lxx+1, 2*lxx+1)

  !#mn{
  !  check whether a rotation matrix for this {beta(ileg),lxp1,mxp1} has
  !  been calculated and saved.  If so, just use the saved value
  ISAV_LOOP: do isav = 1, jsav
     if (betsav(isav).eq.jbmagk) exit ISAV_LOOP
     if (   (lxp1.eq.ldsav(isav)).and.(mxp1.eq.mdsav(isav)).and.      &
          & (abs(betax-betsav(isav)).le.roteps) ) then
        !c    print*, 'using drisav for ', isav, betax, lxp1, mxp1
        do il = 0, lx
           do m1 = -il, il
              do m2 = -il, il
                 drix(m2,m1,il,k,j,i)=cmplx(drisav(m2,m1,il,isav),zero)
              enddo
           enddo
        enddo

        go to 770
     end if
  enddo ISAV_LOOP

  !  initialize dri0
  do in = 1, 2*lxx+1
     do im = 1, 2*lxx+1
        do il = 1, lxx+1
           dri0(il,im,in) = zero
        enddo
     enddo
  enddo


  nm  = mxp1
  ndm = lxp1+nm-1
  xc  = cos(betax/2)
  xs  = sin(betax/2)
  s   = sin(betax)
  dri0(1,1,1) =  1
  dri0(2,1,1) =  xc**2
  dri0(2,1,2) =  s/sqrt(2*one)
  dri0(2,1,3) =  xs**2
  dri0(2,2,1) = -dri0(2,1,2)
  dri0(2,2,2) =  cos(betax)
  dri0(2,2,3) =  dri0(2,1,2)
  dri0(2,3,1) =  dri0(2,1,3)
  dri0(2,3,2) = -dri0(2,2,3)
  dri0(2,3,3) =  dri0(2,1,1)
  do l = 3, lxp1
     ln = 2*l - 1
     lm = 2*l - 3
     if (ln .gt. ndm)  ln = ndm
     if (lm .gt. ndm)  lm = ndm
     do n = 1, ln
        do m = 1, lm
           t1   = (2*l-1-n) * (2*l-2-n)
           t    = (2*l-1-m) * (2*l-2-m)
           f1   = sqrt(t1/t)
           f2   = sqrt( (2*l-1-n) * (n-1) / t )
           t3   = (n-2) * (n-1)
           f3   = sqrt(t3/t)
           dlnm = f1 * xc**2 * dri0(l-1,n,m)
           if (n-1 .gt. 0) dlnm = dlnm - f2*s*dri0(l-1,n-1,m)
           if (n-2 .gt. 0) dlnm = dlnm + f3*xs**2*dri0(l-1,n-2,m)
           dri0(l,n,m) = dlnm
           if (n .gt. (2*l-3))                                         &
                &                  dri0(l,m,n) = (-1)**(n-m) * dri0(l,n,m)
        enddo

        if (n .gt. (2*l-3)) then
           dri0(l,2*l-2,2*l-2) =  dri0(l,2,2)
           dri0(l,2*l-1,2*l-2) = -dri0(l,1,2)
           dri0(l,2*l-2,2*l-1) = -dri0(l,2,1)
           dri0(l,2*l-1,2*l-1) =  dri0(l,1,1)
        endif
     enddo
  enddo



  !     initialize drix
  do il = 0, lx
     do m1 = -lx, lx
        do m2 = -lx, lx
           drix(m2,m1,il,k,j,i) = cmplx(zero,zero)
           drix(m2,m1,il,k,i,i) = cmplx(zero,zero)
        enddo
     enddo
  enddo

  !     Copy result into drix(...,k,j,i) in /rotx/
  do il = 1, lxp1
     mmx = min (il-1, mxp1-1)
     do m1 = -mmx, mmx
        do m2 = -mmx, mmx
           drix(m2, m1, il-1, k, j, i)=cmplx(dri0(il,m1+il,m2+il),zero)
        enddo
     enddo
  enddo

  !#mn{
  !      save dri if there's room
  if (jsav.lt.jsavx) then
     jsav = jsav + 1
     !c   print*, 'saving dri to ',  jsav, betax, lxp1, mxp1
     betsav(jsav) = betax
     ldsav(jsav)  = lxp1
     mdsav(jsav)  = mxp1
     do il = 0, lx
        do m1 = -il, il
           do m2 = -il, il
              drisav(m2,m1,il,jsav) = real(drix(m2,m1,il,k,j,i))
           enddo
        enddo
     enddo

  else
     !c print*, 'not saving dri to ',  betax, lxp1, mxp1
  end if

770 continue
  !#mn}

  !-----test sum rule on d
  !       if (idbg(1).eq.1) then
  !           inquire(file='rotmat.dat', opened=open)
  !           if (.not.open) then
  !               iun = nxtunt(25)
  !               open (iun,file='rotmat.dat',status='unknown')
  !           endif
  !           write(iun,*)'  '
  !           write(iun,*)'atom #s : ',i,j
  !           write(iun,*)  ' il, im, sum, beta'
  !           write(iun,*) ' (drix(il,im,in,k,j,i),in = -il,il)'
  !           do 880 il = 0,lxp1-1
  !             do 870 im = -il,il
  !               sum = 0
  !               do 850 in = -il,il
  !                 term = drix(in,im,il,k,j,i)
  !                 sum = sum+term**2
  !  850           continue
  !               write(iun,860) il,im,sum,betax
  !               write(iun,862) (drix(in,im,il,k,j,i),in = -il,il)
  !  860          format(2i3,1x,f16.12,1x,f8.4)
  !  862          format(5f14.6)
  !  870         continue
  !  880       continue
  ! c          close(iun)
  !       endif
  !-----end test------------------------

  do il = 0, lx
     do m1 = -il, il
        dum = coni * m1 * (xphi(i,j)-pi)
        if (k.eq.1) dum = -dum
        dum = exp( dum )
        do m2 = -il, il
           if (k.eq.1) then
              drix(m2,m1,il,k,j,i) = drix(m2,m1,il,k,j,i) * dum
           else
              drix(m1,m2,il,k,j,i) = drix(m1,m2,il,k,j,i) * dum
           endif
        enddo
     enddo
  enddo

  return
end subroutine rotxan

!====================================================================
subroutine rotint
  use DimsMod, only: lx

  use rotx
  use lnlm
  use xstruc

  implicit none
  integer js,il,m1,m2

  ! Initialize /rotsav/
  jsav = 0
  do js = 1, jsavx
     betsav(js) = jbmagk
     ldsav(js)  = 0
     mdsav(js)  = 0
     do il  = 0, lx
        do m1 = -lx, lx
           do m2 = -lx, lx
              drisav(m2,m1,il,js) = 0
           enddo
        enddo
     enddo
  enddo

  return
end subroutine rotint

subroutine sortat(iph0, nat, npot, iphat, iphx, rat, xrat)

  use DimsMod, only: natxx, nclusx, nphasx
  implicit none
  !--------------------------------------------------------------------
  !  this subroutine sorts the atoms in xrat such that the first npot
  !  entries are each a representative atom of a unique potential.  This
  !  will mean that the upper left corner of the full MS matrix will
  !  contain all of the information needed to compute the fine structure
  !  and all of the electron densities.
  !  NOTA BENE:  the atoms *must* have already been sorted by radial
  !    distance!
  !--------------------------------------------------------------------
  !  input:
  !    iph0:    potential index for central atom in LDOS (added by ala)
  !                (iph0=0 for absorbing atom as the central atom)
  !     nat:    number of atoms in cluster
  !    npot:    number of unique potentials in cluster
  !    iphat:   (nclusx) potential index of each atom in cluster as read
  !             from geometry file
  !    rat:     (3, nclusx) coordinates of each atom in cluster as read
  !             from geometry file
  !  output:
  !    iphx:    (nclusx) potential index of each atom in cluster sorted
  !             so that the first npot+1 entries are examples of each
  !             ipot
  !    xrat:    (3, nclusx) coordinates of each atom in cluster sorted
  !             so that the first npot+1 entries are examples of each
  !             ipot
  !--------------------------------------------------------------------

  integer, intent(in) :: iph0,nat,npot,iphat(natxx)
  real,    intent(in) :: rat(3,natxx)

  real,    intent(out) :: xrat(3,nclusx)
  integer, intent(out) :: iphx(nclusx)

  integer :: ip, ilast, ipoint(0:nphasx)

  ! Added to satisfy implicit none:
  integer :: i,ic,ix,iat,ipp
  integer :: iph,nmin
  real :: xx,yy,zz

  do i=0,nphasx
     ipoint(i) = 0
  enddo

  do ic=1,nat
     iphx(ic) = iphat(ic)
     do ix=1,3
        xrat(ix,ic) = rat(ix,ic)
     enddo
  enddo

  !     (iph0=0 for absorbing atom as the central atom)
  if (iphx(1).ne.iph0) then
     call wlog('* * * ERROR in sortat * * *')
     call wlog('            The first atom in xrat is not '//      &
          &                'the central atom.')
     call wlog('            Complain to Bruce immediately!')
     call par_stop('SORTAT-1')
  endif

  !       if (idbg(4).eq.1) print*,'SORTAT: nat,npot: ',nat,npot
  !       if (idbg(4).eq.1) print*,'SORTAT: xcen,ycen,zcen: ',
  !      $            xcen,ycen,zcen

  ! --- find the example of each unique potential that is closest to the
  !     central atom.  This will presumably be well within the cluster
  !     that was used to compute the overlapped potentials
  ipoint(iph0) = 1
  do ip=0,npot
     if (ip .ne. iph0) then
        do iat=2,nat
           if (iphx(iat).eq.ip .and. ipoint(ip).eq.0) then
              ipoint(ip) = iat
              !                print*,'>>>>> ip, ipoint(ip)', ip, ipoint(ip)
           endif
        enddo

     endif
  enddo

  ! --- now swap the first few atoms with the atoms found above
  IP_LOOP: do ip=0,npot

     ! Some potentials might not be in the xanes cluster
     if (ipoint(ip).eq.0) cycle IP_LOOP
     ! Don't swap two potentials if examples live in the first npot entries
     if (ipoint(ip).le.ip+1) cycle IP_LOOP

     xx  = xrat(1,1+ip)
     yy  = xrat(2,1+ip)
     zz  = xrat(3,1+ip)
     iph = iphx(1+ip)

     xrat(1,1+ip) = xrat(1,ipoint(ip))
     xrat(2,1+ip) = xrat(2,ipoint(ip))
     xrat(3,1+ip) = xrat(3,ipoint(ip))
     iphx(1+ip)  = iphx(ipoint(ip))

     xrat(1,ipoint(ip)) = xx
     xrat(2,ipoint(ip)) = yy
     xrat(3,ipoint(ip)) = zz
     iphx(ipoint(ip))  = iph

     ! added by ala
     ! check that substituted atom was not some ip example
     ! ???BR Jan 16 1998???
     do ipp = ip+1, npot
        if (ipoint(ipp).eq.ip+1) ipoint(ipp) = ipoint(ip)
     enddo
     !       set the correct pointer to ip example
     ipoint(ip) = ip+1

  enddo IP_LOOP

  !     added by ala
  !     Notice that fms will take the last atom of given type ip
  !     from first npot atoms in the list as an example for ip.
  !     Make more permutaions if necesary.
  ilast = -1
  nmin = min (npot+1, nat)
  do ip = 0, npot
     if (ipoint(ip).ne.0) then
        do iat = 1,nmin
           if (iphx(iat).eq.ip) ilast = iat
        enddo
        if (ilast.ne.ipoint(ip)) then
           xx  = xrat(1,ilast)
           yy  = xrat(2,ilast)
           zz  = xrat(3,ilast)

           xrat(1,ilast)= xrat(1,ipoint(ip))
           xrat(2,ilast)= xrat(2,ipoint(ip))
           xrat(3,ilast)= xrat(3,ipoint(ip))

           xrat(1,ipoint(ip)) = xx
           xrat(2,ipoint(ip)) = yy
           xrat(3,ipoint(ip)) = zz
           !           now ipoint(ip) = ilast, but don't need ipoint anymore
        endif
     endif
  enddo

  !       if (idbg(4).eq.1) then
  !           do 220 i=1,npot+1
  !             print *,i,xrat(1,i),xrat(2,i),xrat(3,i),iphx(i)
  !  220      continue
  !       endif
  return
end subroutine sortat

subroutine xanlm(lmaxp1,mmaxp1)

  use rotx
  use lnlm
  use xstruc
  use afctr

  !------------------------------------------------------------------
  !  calculate and store all of the legendre polynomial normalization
  !  factors needed in the problem
  !     xnlm= sqrt ((2l+1)(l-m)!/(l+m)!)
  !  see, for instance, Arfken section 12.6.  Note that this lacks the
  !  factor of sqrt(4*pi)
  !
  !  inputs:
  !     lmaxp1, nmaxp1:  maximun l and m considered in the problem +1
  !                      i.e. lmaxp1 = l_max+1
  !
  !  outputs:
  !     all normalization factors passed in common /xnlm/
  !------------------------------------------------------------------

  implicit none

  integer, intent(in) :: lmaxp1,mmaxp1

  ! Added to satisfy implicit none
  integer :: im,il,l,m,mmxp1
  real    :: cnlm
  
  call xfctst
  do il=1,lmaxp1
     mmxp1 = min(mmaxp1,il)
     do im=1,mmxp1
        l    = il-1
        m    = im-1
        cnlm = (2*l+1) * flg(l-m) / flg(l+m)
        cnlm = sqrt(cnlm) * afac**m
        xnlm(m,l) = cnlm
     enddo
  enddo

  return
end subroutine xanlm

subroutine xfctst
  !  same as feff's factst, but with a different name
  use afctr

  implicit none
  !     program for s3j and s6j symbols obtained from
  !     steve younger of n.b.s.   modified by j.b.mann
  !--------------------------------------------------------------------
  !     a set to 1/64 to prevent overflow on vax
  !     range on  flg set to 0:210, rather than flg(210)
  !--------------------------------------------------------------------
  !BR   This allows calculation of a large factorial (~100) without
  !BR   overflow problems -- factor in a power of a small number then
  !BR   factor it out
  !--------------------------------------------------------------------

  integer :: i

  afac=0.03125
  !afac=0.015625
  flzero = 1.0
  flg(0) = 1.0
  flg(1) = afac
  do i=2,50
     flg(i) = flg(i-1) * i * afac
  enddo

  return
end subroutine xfctst
