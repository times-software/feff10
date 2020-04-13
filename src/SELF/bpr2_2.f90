!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: bpr2_2.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Written by Josh Kas.
!     For explanation see bpr1.ps
      FUNCTION bpr2(q,dppar,cpar)
!     PROGRAM bpr2
      COMPLEX*16 q, cpar(2)
      DOUBLE PRECISION dppar(4)
      
      COMPLEX*16 En, emk, k
      COMPLEX*16 rq(0:12), rq00, rq0L1, rq0L2, rq0L3
      COMPLEX*16 Amp, a1, a2, a3, a4, L1, L2, L3, L4
      COMPLEX*16 gam, wq
      INTEGER iqexp, iqmax
      LOGICAL OnShll

      DOUBLE PRECISION pi, wp2
      PARAMETER (pi = 3.1415926535897932384626433d0)
      
      COMPLEX*16 bpr2
            
      COMPLEX*16 Logi
      DOUBLE PRECISION Omegaq, Gamq
      EXTERNAL Logi, Omegaq, Gamq

      IF(dppar(3).gt.0.5d0) THEN
         OnShll = .true.
      ELSE
         OnShll = .false.
      END IF
      
      k   = cpar(1)
      En  = cpar(2)
      wp2 = dppar(1)**2
      wq = Omegaq(dppar(1), DBLE(q))
!      wq = SQRT(dppar(1)**2 + 2.d0*q**2 - 0.8*q**3 + q**4)
!      wq  = SQRT(wp2*theta((2.d0-q)/0.2d0) +
!     &     theta((q-2.d0)/0.2d0)*(4.d0*q**2/3.d0 + q**4 + wp2 - wq2Of2))
!     Small q dependence in gamma
      gam = Gamq(dppar(2),DBLE(q))
!     High q limit for gammma
!      gam = SQRT(dppar(2)**2 + (2.4d0*q)**2)
!     high q limit with low q cutoff
!      gam = SQRT(dppar(2)**2 +
!     &     (0.5 + 2.0d0**2*
!     &     (ATAN(DBLE((q**2+2*q-wq)/(dppar(1)*0.1d0))) + pi)/(2*pi))*
!     &     q**2)
      

      IF(OnShll) THEN
         emk = 0.d0
      ELSE
         emk = En - k**2
      END IF

      IF(DBLE(wq-gam).gt.0.d0) THEN
!      IF(.true.) THEN
         iqmax = 12

         Amp = 1.d0/(96.d0*gam**4*q*(gam**2 + 7.d0*wq**2))
         
         a1 = (emk + gam) + (2.d0*k*q - q**2 - wq)
         a2 = (emk + gam) - (2.d0*k*q + q**2 + wq)
         a3 = (emk - gam) + (2.d0*k*q - q**2 - wq)
         a4 = (emk - gam) - (2.d0*k*q + q**2 + wq)
      
!     La1M - La2M - La3M + La4M;
         L1 = LOG(ABS((a1/a3)*(a4/a2))) + Logi(-a1,-1) - Logi(-a2,-1) - &
     &        Logi(-a3,-1) + Logi(-a4,-1)
!     La1 - La2 - La3 + La4;
         L2 = LOG(ABS((a1/a3)*(a4/a2))) + Logi(a1,1) - Logi(a2,1) -     &
     &        Logi(a3,1) + Logi(a4,1) 
!     La1M + La2M - La3M - La4M;
         L3 = LOG(ABS((a1/a3)*(a2/a4))) + Logi(-a1,-1) + Logi(-a2,-1) - &
     &        Logi(-a3,-1) - Logi(-a4,-1)
!     La1 - La2 + La3 - La4;
         L4 = LOG(ABS((a1/a2)*(a3/a4))) + Logi(a1,1) - Logi(a2,1) +     &
     &        Logi(a3,1) - Logi(a4,1)
                  
         rq(12) = (21.*(-5.*L1*q**12 - 24.*L1*q**10*wq - 45.*L1*q**8*   &
     &        wq**2 - 40.*L1*q**6*wq**3 - 15.*L1*q**4*wq**4 + L2*wq**6))&
     &        /gam
         
         rq(11) = (1260.*k*L3*q**3*(q**2 + wq)**4)/gam
         
         rq(10) = (630.*L1*q**2*(q**2 + wq)**3*                         &
     &        (emk*(q**2 + wq) - 2.*k**2*(5.*q**2 + wq)))/gam
         
         rq(9) = (-84.*k*q*(q**2 + wq)**2*                              &
     &        (15.*emk*L3*(q**2 + wq)*(5.*q**2 + wq) -                  &
     &        2.*gam*(q**2 + wq)*(25.*q**2 + wq) -                      &
     &        40.*k**2*L3*q**2*(5.*q**2 + 2.*wq)))/gam
         
         rq(8) = (105.*(-3.*emk**2*L1*(q**2 + wq)**3*(5.*q**2 + wq) +   &
     &        48.*emk*k**2*L1*q**2*(q**2 + wq)**2*(5.*q**2 + 2.*wq) -   &
     &        48.*k**4*L1*q**4*(q**2 + wq)*(5.*q**2 + 3.*wq) +          &
     &        gam**2*(3.*L1*q**8 + 8.*L1*q**6*wq + 6.*L1*q**4*wq**2 -   &
     &        1.*L2*wq**4)))/gam
         
         rq(7) = (168.*k*q*(15.*emk**2*L3*(q**2 + wq)**2*(5.*q**2 + 2.  &
     &        *wq)- 4.*emk*(q**2 + wq)*(15.*k**2*L3*q**2*(5.*q**2 + 3.  &
     &        *wq) +gam*(q**2 + wq)*(25.*q**2 + 7.*wq)) +q**2*(-15.*gam &
     &        **2*L3*(q**2 + wq)**2 +24.*k**4*L3*q**2*(5.*q**2 + 4.*wq) &
     &        +8.*gam*k**2*(q**2 + wq)*(25.*q**2 + 13.*wq))))/gam
         
         rq(6) = (420.*L1*(emk**3*(q**2 + wq)**2*(5.*q**2 + 2.*wq) -    &
     &        18.*emk**2*k**2*q**2*(q**2 + wq)*(5.*q**2 + 3.*wq) +      &
     &        2.*k**2*q**2*(-8.*k**4*q**4 +                             &
     &        3.*gam**2*(q**2 + wq)*(3.*q**2 + wq)) +                   &
     &        3.*emk*q**2*(-1.*gam**2*(q**2 + wq)**2 +                  &
     &        8.*k**4*(5.*q**4 + 4.*q**2*wq))))/gam
      
         rq(5) = (-168.*k*q*(15.*emk**3*L3*(q**2 + wq)*(5.*q**2 + 3.*wq)&
     &        +4.*gam*(-20.*k**4*q**4 +gam**2*(q**2 + wq)*(10.*q**2 + wq&
     &        ) +5.*gam*k**2*L3*q**2*(3.*q**2 + 2.*wq)) -6.*emk**2*(10. &
     &        *k**2*L3*q**2*(5.*q**2 + 4.*wq) +gam*(q**2 + wq)*(25.*q**2&
     &        + 13.*wq)) +emk*(120.*k**4*L3*q**4 -15.*gam**2*L3*(q**2 + &
     &        wq)*(3.*q**2 + wq) +16.*gam*k**2*q**2*(25.*q**2 + 19.*wq))&
     &        ))/gam
         
         rq(4) = (315.*(-1.*(gam**4 -                                   &
     &        2.*gam**2*(3.*emk**2 - 24.*emk*k**2 + 8.*k**4) +          &
     &        5.*emk**2*(emk**2 - 16.*emk*k**2 + 16.*k**4))*L1*q**4 -   &
     &        8.*emk*(emk**2*(emk - 8.*k**2) - 1.*gam**2*(emk - 4.*k**2)&
     &        )*L1*q**2*wq + (-3.*emk**4*L1 + 2.*emk**2*gam**2*L1 +gam  &
     &        **4*L2)*wq**2))/gam
         
         rq(3) = (28.*k*q*(-5.*gam**3*(64.*k**2 - 9.*gam*L3)*q**2 +     &
     &        45.*emk**4*L3*(5.*q**2 + 4.*wq) +                         &
     &        24.*emk*gam**2*(5.*(4.*gam + 3.*k**2*L3)*q**2 +           &
     &        11.*gam*wq) -                                             &
     &        24.*emk**3*(25.*(gam + k**2*L3)*q**2 + 19.*gam*wq) +      &
     &        30.*emk**2*gam*(40.*k**2*q**2 - 3.*gam*L3*(3.*q**2 + 2.*wq&
     &        ))))/gam
         
         rq(2) = (42.*(15.*(emk - 1.*gam)*(emk + gam)*                  &
     &        (emk**2*(emk - 10.*k**2) - 1.*gam**2*(emk - 2.*k**2))*L1* &
     &        q**2 + 4.*(3.*emk**5*L1 - 5.*emk**3*gam**2*L1 -           &
     &        2.*gam**5*L4)*wq))/gam
         
         rq(1) = (84.*k*(50.*emk**4*gam - 80.*emk**2*gam**3 + 22.*gam**5&
     &        -15.*emk**5*L3 + 30.*emk**3*gam**2*L3 - 15.*emk*gam**4*L3)&
     &        *q)/gam
         
         rq(0) = (105.*(-1.*emk**6*L1 + 3.*emk**4*gam**2*L1 -           &
     &        3.*emk**2*gam**4*L1 + gam**6*L2))/gam

      ELSE
         iqmax = 6

         Amp = 7.d0/(120.d0*q*(gam + wq)**5*(8.d0*gam**2 - 5.d0*gam*wq +&
     &        wq**2))
!        + idelta
         a1 = emk - 1.*gam - 2.*k*q - 1.*q**2 - 1.*wq
!        + idelta         
         a2 = emk - 1.*gam + 2.*k*q - 1.*q**2 - 1.*wq
!        + idelta
         a3 = emk + 2.*k*q - 1.*q**2
!        + idelta
         a4 = emk - 2.*k*q - 1.*q**2

!        L1 = Log[a1] - Log[a2];
         L1 = LOG(ABS((a1/a2))) + Logi(a1,1) - Logi(a2,1)
!        L2 = Log[-a1] + Log[-a2] - Log[-a3] - Log[-a4];
         L2 = LOG(ABS((a1/a4)*(a2/a3))) + Logi(-a1,-1) + Logi(-a2,-1) - &
     &        Logi(-a3,-1) - Logi(-a4,-1)
!        L3 = Log[-a1] - Log[-a2] + Log[-a3] - Log[-a4];
         L3 = LOG(ABS((a1/a2)*(a3/a4))) + Logi(-a1,-1) - Logi(-a2,-1) + &
     &        Logi(-a3,-1) - Logi(-a4,-1)

         rq(6) = 300.*gam**6*L1
         
         rq(5) = 240.*gam**5*(11.*k*q + 4.*L1*wq)
         
         rq(4) = 60.*gam**4*(-15.*L3*(emk**2 - 2.*(emk - 2.*k**2)*q**2 +&
     &        q**4) + 15.*L1*wq**2 +2.*k*q*(5.*(-5. + 6.*L2)*(emk - 1.*q&
     &        **2) + 38.*wq))
         
         rq(3) = -320.*gam**3*k*q*(30.*emk**2 + 40.*k**2*q**2 + 3.*(q**2&
     &        + wq)*(10.*q**2 + wq) - 3.*emk*(20.*q**2 + 11.*wq))
         
         rq(2) = 60.*gam**2*(-5.*L1*wq**4 + 4.*k*q*(10.*(-1. + 3.*L2)*( &
     &        -1.*emk + q**2)*(emk**2 - 2.*(emk - 2.*k**2)*q**2 + q**4) &
     &        +2.*(-11. + 10.*L2)*(3.*emk**2 - 6.*emk*q**2 + 4.*k**2*q  &
     &        **2 + 3.*q**4)*wq +(-77. + 30.*L2)*(-1.*emk + q**2)*wq**2 &
     &        - 16.*wq**3) +5.*L3*(3.*(emk**4 + q**8) - 8.*emk**3*wq + 6.&
     &        *emk**2*wq**2 +6.*q**4*(3.*emk**2 + 8.*k**4 + 16.*k**2*wq&
     &        + wq**2 - 4.*emk*(6.*k**2 + wq)) -12.*q**2*(emk - 1.*wq)  &
     &        *(emk**2 + 2.*k**2*wq - 1.*emk*(6.*k**2 + wq)) -4.*q**6*(3.&
     &        *emk - 2.*(9.*k**2 + wq))))
         
         rq(1) = 240.*gam*k*q*(25.*emk**4 + 80.*k**4*q**4 + (q**2 + wq) &
     &        **3*(25.*q**2 + wq) +8.*k**2*q**2*(q**2 + wq)*(25.*q**2 + &
     &        13.*wq) - 4.*emk**3*(25.*q**2 + 19.*wq) +2.*emk**2*(100.*k&
     &        **2*q**2 + 3.*(q**2 + wq)*(25.*q**2 + 13.*wq)) -4.*emk*((q&
     &        **2 + wq)**2*(25.*q**2 + 7.*wq) + 4.*k**2*(25.*q**4 + 19. &
     &        *q**2*wq)))

!        Pieces of rq(0)
         rq00  = 40.*k*q*(30.*(5.*emk**4 - 20.*emk**2*(emk - 2.*k**2)*q &
     &        **2 +2.*(15.*emk**2 - 40.*emk*k**2 + 8.*k**4)*q**4 - 20.  &
     &        *(emk - 2.*k**2)*q**6 + 5.*q**8)*wq -516.*(emk - 1.*q**2) &
     &        *(emk**2 - 2.*(emk - 2.*k**2)*q**2 + q**4)*wq**2 +208.*(3.&
     &        *emk**2 - 6.*emk*q**2 + 4.*k**2*q**2 + 3.*q**4)*wq**3 -   &
     &        291.*(emk - 1.*q**2)*wq**4 +30.*wq**5)

         rq0L1 = 60.*wq**6*L1
         
         rq0L2 = -240.*k*q*(-15.*emk**5 + 15.*emk**4*(5.*q**2 + 4.*wq)  &
     &        -10.*emk**3*(20.*k**2*q**2 + 3.*(q**2 + wq)*(5.*q**2 + 3. &
     &        *wq)) -15.*emk*(16.*k**4*q**4 + (q**2 + wq)**3*(5.*q**2 + &
     &        wq) +8.*k**2*q**2*(q**2 + wq)*(5.*q**2 + 3.*wq)) +30.*emk &
     &        **2*((q**2 + wq)**2*(5.*q**2 + 2.*wq) + 4.*k**2*(5.*q**4 +&
     &        4.*q**2*wq)) +q**2*(15.*(q**2 + wq)**4 + 40.*k**2*(q**2 + &
     &        wq)**2*(5.*q**2 + 2.*wq) +48.*k**4*(5.*q**4 + 4.*q**2*wq))&
     &        )*L2
         
         rq0L3 = 60.*(-5.*(emk**2 - 2.*(emk - 2.*k**2)*q**2 + q**4)*    &
     &        (emk**4 - 4.*emk**2*(emk - 14.*k**2)*q**2 + 2.*(3.*emk**2 &
     &        - 56.*emk*k**2 + 8.*k**4)*q**4 -4.*(emk - 14.*k**2)*q**6 +&
     &        q**8) +24.*(emk - 1.*q**2)*(emk**4 - 4.*emk**2*(emk - 10. &
     &        *k**2)*q**2 +2.*(3.*emk**2 - 40.*emk*k**2 + 40.*k**4)*q**4&
     &        - 4.*(emk - 10.*k**2)*q**6 + q**8)*wq -45.*(emk**4 - 4.   &
     &        *emk**2*(emk - 6.*k**2)*q**2 + 2.*(3.*emk**2 - 24.*emk*k  &
     &        **2 + 8.*k**4)*q**4 -4.*(emk - 6.*k**2)*q**6 + q**8)*wq**2&
     &        +40.*(emk - 1.*q**2)*(emk**2 - 2.*(emk - 6.*k**2)*q**2 + q&
     &        **4)*wq**3 -15.*(emk**2 - 2.*(emk - 2.*k**2)*q**2 + q**4) &
     &        *wq**4)*L3
         
         rq(0) = rq00 + rq0L1 + rq0L2 + rq0L3
         
      END IF
      
      bpr2 = 0
      DO iqexp = 0, iqmax
         bpr2 = bpr2 + rq(iqexp)
      END DO
      bpr2 = Amp*bpr2

      RETURN
      END
