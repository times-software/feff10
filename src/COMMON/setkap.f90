!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: setkap.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setkap(ihole, kinit, linit)
      implicit double precision (a-h, o-z)

!     Set initial state ang mom and quantum number kappa
!     ihole  initial state from ihole    
!     1      K    1s      L=0 -> linit=0 
!     2      LI   2s      L=0 -> linit=0
!     3      LII  2p 1/2  L=1 -> linit=1
!     4      LIII 2p 3/2  L=1 -> linit=1
!     5+     etc.
      if (ihole.le. 2 .or. ihole.eq. 5 .or. ihole.eq.10 .or.            &
     &    ihole.eq.17 .or. ihole.eq.24 .or. ihole.eq.29)  then
!        hole in s state
         linit = 0
         kinit = -1
      elseif (ihole.eq. 3 .or. ihole.eq. 6 .or. ihole.eq.11 .or.        &
     &        ihole.eq.18 .or. ihole.eq.25 .or. ihole.eq.30)  then
!        hole in p 1/2 state
         linit = 1
         kinit = 1
      elseif (ihole.eq. 4 .or. ihole.eq. 7 .or. ihole.eq.12 .or.        &
     &        ihole.eq.19 .or. ihole.eq.26)  then
!        hole in p 3/2 state
         linit = 1
         kinit = -2
      elseif (ihole.eq. 8 .or. ihole.eq.13 .or.                         &
     &        ihole.eq.20 .or. ihole.eq.27)  then
!        hole in d 3/2 state
         linit = 2
         kinit = 2
      elseif (ihole.eq. 9 .or. ihole.eq.14 .or.                         &
     &        ihole.eq.21 .or. ihole.eq.28)  then
!        hole in d 5/2 state
         linit = 2
         kinit = -3
      elseif (ihole.eq.15 .or. ihole.eq.22)  then
!        hole in  f 5/2 state
         linit = 3
         kinit = 3
      elseif (ihole.eq.16 .or. ihole.eq.23)  then
!        hole in  f 7/2 state
         linit = 3
         kinit = -4
      else
!        some unknown hole
         call par_stop('invalid hole number in setkap')
      endif

      return
      end
