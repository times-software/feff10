!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: s02at.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Josh Kas - Changed array dimensions from 30 to 41 (and others) for high Z elements
! according to Pavlo Baranov's changes.
      subroutine s02at(ihole, norb, nk, xnel, ovpint, dval)
      implicit double precision (a-h,o-z)
      double precision  m1(8,8), m2(8,8)
      dimension nk(41), xnel(41), iorb(41), ovpint(41,41)
      external determ

      dval = 1.0
!     loop over possible kappa for existing atoms
      do 100 kap = -5,4

!        initialize matrices and other stuff
         do 10 i = 1,8
         do 10 j = 1,8
            m1(i,j) = 0
  10        m2(i,j) = 0
         do 20 i= 1,8
            iorb(i) = 0
            m1(i,i) = 1.0
  20        m2(i,i) = 1.0
!        morb - number of orbitals with quantum number kappa
         morb = 0
         nhole = 0

!        construct the largest possible matrix for given value of kappa.
         do 40 i = 1, norb
            if (nk(i) .eq. kap) then
               morb = morb + 1
               iorb(morb) = i
               do 50 j = 1, morb
!                 print overlap integrals
!                 print*, kap,' ', iorb(j),' ', iorb(morb), '
!    1                            ovp= ',ovpint(iorb(j), iorb(morb))
   50             m1(j,morb) = ovpint(iorb(j), iorb(morb))
               do 60 j = 1, morb - 1
   60             m1(morb,j) = m1(j,morb)
               
               if (ihole .eq. i) nhole = morb
            endif
   40    continue
         if (morb .eq. 0) goto 100
         dum1 = determ(m1, morb, 8)
         dum1 = dum1**2

         dum3 = determ(m1, morb-1, 8)
         dum3 = dum3**2
         xn = xnel(iorb(morb))
         nmax = 2*abs(kap)
         xnh = nmax - xn
         if (nhole .eq. 0) then 
            dval = dval * dum1**xn * dum3**xnh
         elseif (nhole .eq. morb) then
            dval = dval * dum1**(xn-1) * dum3**(xnh+1)
         else
            call elimin(m1,nhole,m2)
            dum2 = determ(m2,morb,8)
            dum2 = dum2**2
            dum4 = determ(m2,morb-1,8)
            dum4 = dum4**2
            dum5 = (dum4*dum1*xnh + dum2*dum3*xn)/nmax
            dval = dval * dum5 * dum1**(xn-1) * dum3**(xnh-1)
         endif

100   continue

      return
      end

      subroutine elimin(d1,n,d2)
      implicit double precision (a-h,o-z)
      dimension d1(8,8), d2(8,8)

      do 10 i = 1,8
      do 10 j = 1,8
         if (i .ne. n) then
            if (j .ne. n) then
               d2(i,j)=d1(i,j)
            else
               d2(i,j) = 0
            endif
         else
            if (j .ne. n) then
               d2(i,j) = 0
            else
               d2(i,j) = 1.0
            endif
         endif
   10 continue
      return
      end
