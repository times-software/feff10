!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: setgam.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setgam (iz, ihole, gamach)

!     Sets gamach, core hole lifetime.  Data comes from graphs in
!     K. Rahkonen and K. Krause,
!     Atomic Data and Nuclear Data Tables, Vol 14, Number 2, 1974.
!     output gamach is in eV

      implicit double precision (a-h, o-z)

      dimension zh(8,16), gamh(8,16)

      dimension zk(8), gamkp(8)
      parameter (ryd  = 13.605698d0)
      parameter (hart = 2*ryd)
      character*512 slog


!     Note that 0.99 replaces 1.0, 95.1 replaces 95.0 to avoid roundoff
!     trouble.
!     Gam arrays contain the gamma values.
!     We will take log10 of the gamma values so we can do linear
!     interpolation from a log plot.

      data  zh   / 0.99,  10.0, 20.0, 40.0, 50.0, 60.0, 80.0, 95.1,     &
     &              0.99, 18.0, 22.0, 35.0, 50.0, 52.0, 75.0,  95.1,    &
     &              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,   &
     &              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,   &
     &              0.99,  20.0, 28.0, 30.0, 36.0, 53.0,  80.0, 95.1,   &
     &              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,   &
     &              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,   &
     &              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,   &
     &              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,   &
     &              0.99,  30.0, 40.0, 47.0, 50.0, 63.0,  80.0, 95.1,   &
     &              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,   &
     &              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,   &
     &              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,   &
     &              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,   &
     &              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0,   &
     &              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0/

      data  gamh / 0.02,  0.28, 0.75,  4.8, 10.5, 21.0, 60.0, 105.0,    &
     &              0.07,  3.9,  3.8,  7.0,  6.0,  3.7,  8.0,  19.0,    &
     &              0.001, 0.12,  1.4,  0.8,  2.6,  4.1,   6.3, 10.5,   &
     &              0.001, 0.12, 0.55,  0.7,  2.1,  3.5,   5.4,  9.0,   &
     &              0.001,  1.0,  2.9,  2.2,  5.5, 10.0,  22.0, 22.0,   &
     &              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  15.0, 16.0,   &
     &              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  10.0, 10.0,   &
     &              0.0006,0.09, 0.07, 0.48,  1.0,  4.0,   2.7,  4.7,   &
     &              0.0006,0.09, 0.07, 0.48, 0.87,  2.2,   2.5,  4.3,   &
     &              0.001,0.001,  6.2,  7.0,  3.2, 12.0,  16.0, 13.0,   &
     &              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,   &
     &              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,   &
     &              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,   &
     &              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,   &
     &              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9,   &
     &              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9/

!     Since feff8 can be called any number of times . ALA

      if (ihole .le. 0)  then
         gamach = 0
         write(slog,'(a,1pe13.5)') ' No hole in SETGAM, gamach = ',     &
     &                             gamach
         call wlog(slog)
         return
      endif
      if (ihole .gt. 16)  then
         call wlog(' This version of FEFF will set gamach = 0.1 eV ' // &
     &             ' for O1 and higher hole')
         call wlog(' You can use CORRECTIONS card  to set ' //          &
     &   ' gamach = 0.1 + 2*vicorr ')
!        stop 'SETGAM-2'
      endif

      zz = iz
      if (ihole .le. 16)  then
         do 10  i = 1, 8
            gamkp(i) = log10 (gamh(i,ihole))
            zk(i) = zh(i,ihole)
   10    continue
         call terp (zk, gamkp, 8, 1, zz, gamach)
      else
!     include data from the tables later.
!     Now gamach=0.1eV for any O-hole for any element.
         gamach = -1.0
      endif

!     Change from log10 (gamma) to gamma
      gamach = 10.0 ** gamach


      return
      end
