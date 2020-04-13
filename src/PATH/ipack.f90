!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ipack.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ipack (iout, n, ipat)

!     Input:  n          number of things to pack, nmax=8
!             ipat(1:n)  integers to pack
!     Output: iout(3)    packed version of n and ipat(1:n)
!
!     Packs n and ipat(1:n) into 3 integers, iout(1:3).  Algorithm
!     packs three integers (each between 0 and 1289 inclusive) into a
!     single integer.  Single integer must be INT*4 or larger, we assume
!     that one bit is wasted as a sign bit so largest positive int
!     is 2,147,483,647 = (2**31 - 1).
!     This version is specifically for the path finder and
!     degeneracy checker.

      dimension iout(3), ipat(n)
      dimension itmp(8)
      parameter (ifac1 = 1290, ifac2 = 1290**2)

      if (n .gt. 8)  call par_stop('ipack n too big')

      do 10  i = 1, n
         itmp(i) = ipat(i)
   10 continue
      do 20  i = n+1, 8
         itmp(i) = 0
   20 continue

      iout(1) = n       + itmp(1)*ifac1 + itmp(2)*ifac2
      iout(2) = itmp(3) + itmp(4)*ifac1 + itmp(5)*ifac2
      iout(3) = itmp(6) + itmp(7)*ifac1 + itmp(8)*ifac2

      return
      end
      subroutine upack (iout, n, ipat)

!     retrieve n and ipat from iout
!     Input:  iout(3)  packed integers
!             n        max number to get, must be .le. 8
!     Output: n        number unpacked
!             ipat(1:n) unpacked integers

      dimension iout(3), ipat(n)
      dimension itmp(8)
      parameter (ifac1 = 1290, ifac2 = 1290**2)

      nmax = n
      if (nmax .gt. 8)  call par_stop('nmax .gt. 8 in upack')

      n = mod (iout(1), ifac1)
      if (n .gt. nmax)  call par_stop('nmax in upack too small')

      itmp(1) = mod (iout(1), ifac2) / ifac1
      itmp(2) = iout(1) / ifac2
      itmp(3) = mod (iout(2), ifac1)
      itmp(4) = mod (iout(2), ifac2) / ifac1
      itmp(5) = iout(2) / ifac2
      itmp(6) = mod (iout(3), ifac1)
      itmp(7) = mod (iout(3), ifac2) / ifac1
      itmp(8) = iout(3) / ifac2

      do 10  i = 1, n
         ipat(i) = itmp(i)
   10 continue

      return
      end
