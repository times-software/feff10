!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: bpr1_2.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Written by Josh Kas.
!     For explanation see bpse.ps
      FUNCTION bpr1(q,dppar,cpar)
!     PROGRAM bpr1
      COMPLEX*16 q, cpar(2)
      DOUBLE PRECISION dppar(4)
      
      COMPLEX*16 En, emk, k, gam, wq
      COMPLEX*16 rq(0:12), rq00, rq0L1, rq0L2, rq0L3, rq0L4, rq0L5
      COMPLEX*16 Amp, a1, a2, a3, a4, b1, b2, b3, b4, L1, L2, L3, L4,   &
     &     L5, L6, L7, L8
      DOUBLE PRECISION delta, Zero, MM1
      COMPLEX*16 bpr1
      INTEGER iqexp, iqmax
      LOGICAL OnShll

      COMPLEX*16 Logi
      DOUBLE PRECISION Omegaq, Gamq
      EXTERNAL Logi, Omegaq, Gamq
      
      DOUBLE PRECISION pi
      PARAMETER (pi = 3.1415926535897932384626433d0)
      PARAMETER(delta = 1.d-10, Zero = 1.d-10)
      IF(dppar(3).gt.0.5d0) THEN
         OnShll = .true.
      ELSE
         OnShll = .false.
      END IF
      k = cpar(1)
      En = cpar(2)

      wq = Omegaq(dppar(1), DBLE(q))
!      wq = SQRT(dppar(1)**2 + 2.d0*q**2 - 0.8*q**3 + q**4)
!     Small q dependence in gamma
      gam = Gamq(dppar(2),DBLE(q))
!     High q limit for gammma
!      gam = SQRT(dppar(2)**2  + (2.4d0*q)**2)
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
         
         a1 = (emk - gam) + (2.d0*k*q + wq - q**2) 
         a2 = (emk + gam) - (2.d0*k*q + wq + q**2) 
         a3 = (emk + gam) + (2.d0*k*q + wq - q**2) 
         a4 = (emk - gam) - (2.d0*k*q + wq + q**2) 
         
         b1 = (1.d0 - gam) + (wq - En) 
         b2 = (1.d0 - gam) - (wq + En) 
         b3 = (1.d0 + gam) + (wq - En) 
         b4 = (1.d0 + gam) - (wq + En) 
         
!        La1 - La2M - La3 + La4M
         L1 = LOG(ABS((a1/a3)*(a4/a2))) + Logi(a1,-1) - Logi(-a2,-1) -  &
     &        Logi(a3,-1) + Logi(-a4,-1)
!        La1 + La2M - La3 - La4M
         L2 = LOG(ABS((a1/a3)*(a2/a4))) +  Logi(a1,-1) + Logi(-a2,-1) - &
     &        Logi(a3,-1) - Logi(-a4,-1) 
!        La1 - La2 - La3 + La4
         L3 = LOG(ABS((a1/a3)*(a4/a2))) + Logi(a1,-1) - Logi(a2,1) -    &
     &        Logi(a3,-1) + Logi(a4,1)   
!        La1 - La2 + La3 - La4
         L4 = LOG(ABS((a1/a4)*(a3/a2))) + Logi(a1,-1) - Logi(a2,1) +    &
     &        Logi(a3,-1) - Logi(a4,1)   
!        Lb1M + Lb2M - Lb3M - Lb4M
         L5 = LOG(ABS((b1/b3)*(b2/b4))) + Logi(-b1,1) + Logi(-b2,-1) -  &
     &        Logi(-b3,1) - Logi(-b4,-1) 
!        Lb1 - Lb2M - Lb3 + Lb4M
         L6 = LOG(ABS((b1/b3)*(b4/b2))) + Logi(b1,-1) - Logi(-b2,-1) -  &
     &        Logi(b3,-1) + Logi(-b4,-1) 
!        Lb1 + Lb2M - Lb3 - Lb4M
         L7 = LOG(ABS((b1/b3)*(b2/b4))) + Logi(b1,-1) + Logi(-b2,-1) -  &
     &        Logi(b3,-1) - Logi(-b4,-1) 
!        Lb1M - Lb2M + Lb3M - Lb4M
         L8 = LOG(ABS((b1/b2)*(b3/b4))) + Logi(-b1,1) - Logi(-b2,-1) +  &
     &        Logi(-b3,1) - Logi(-b4,-1)
      
         Amp = 7.d0/(480.d0*gam**4*q*(gam**2 + 7.d0*wq**2))

         rq(12) = (15.*((L3 + L5)*wq**6 +                               &
     &        8.*L2*q**6*wq*(3.*q**4 + 5.*wq**2) -                      &
     &        5.*L1*q**4*(q**8 + 9.*q**4*wq**2 + 3.*wq**4)))/gam
         
         rq(11) = (900.*k*q**3*(-4.*L1*q**2*wq*(q**4 + wq**2) +         &
     &        L2*(q**8 + 6.*q**4*wq**2 + wq**4)))/gam
      
         rq(10) = (30.*q**2*(5.*(2.*gam + 3.*emk*L1 - 30.*k**2*L1)*q**8 &
     &        -60.*(emk - 8.*k**2)*L2*q**6*wq +2.*(26.*gam + 45.*(emk - &
     &        6.*k**2)*L1)*q**4*wq**2 -60.*(emk - 4.*k**2)*L2*q**2*wq**3&
     &        +(2.*gam + 15.*(emk - 2.*k**2)*L1)*wq**4))/gam
         
         rq(9) = (60.*k*q*(-25.*(3.*emk - 8.*k**2)*L2*q**8 +            &
     &        8.*(19.*gam + 30.*(emk - 2.*k**2)*L1)*q**6*wq -           &
     &        90.*(3.*emk - 4.*k**2)*L2*q**4*wq**2 +                    &
     &        8.*(7.*gam + 15.*emk*L1 - 10.*k**2*L1)*q**2*wq**3 -       &
     &        15.*emk*L2*wq**4))/gam
         
         rq(8) = (15.*(-5.*(20.*gam*(emk - 8.*k**2) - 3.*gam**2*L1 +    &
     &        15.*(emk**2 - 16.*emk*k**2 + 16.*k**4)*L1)*q**8 +         &
     &        40.*(6.*emk**2 - 1.*gam**2 - 72.*emk*k**2 + 48.*k**4)*L2* &
     &        q**6*wq - 6.*(52.*gam*(emk - 4.*k**2) - 5.*gam**2*L1 +    &
     &        15.*(3.*emk**2 - 24.*emk*k**2 + 8.*k**4)*L1)*q**4*wq**2   &
     &        + 120.*emk*(emk - 4.*k**2)*L2*q**2*wq**3 -                &
     &        1.*(4.*emk*gam + 15.*emk**2*L1 +                          &
     &        gam*(4. - 4.*En + 5.*gam*(L3 + L5)) +                     &
     &        15.*(-1. + En)**2*L7)*wq**4))/gam
         
         rq(7) = (120.*k*q*(15.*(5.*emk**2 - 1.*gam**2 -                &
     &        20.*emk*k**2 + 8.*k**4)*L2*q**6 -                         &
     &        2.*(38.*gam*(3.*emk - 4.*k**2) - 15.*gam**2*L1 +          &
     &        6.*(15.*emk**2 - 40.*emk*k**2 + 8.*k**4)*L1)*q**4*wq -    &
     &        15.*(-9.*emk**2 + gam**2 + 12.*emk*k**2)*L2*q**2*wq**2 -  &
     &        2.*emk*(14.*gam + 15.*emk*L1)*wq**3))/gam
         
         rq(6) = (20.*(5.*(-8.*gam**3 +                                 &
     &        30.*gam*(emk**2 - 12.*emk*k**2 + 8.*k**4) -               &
     &        9.*gam**2*(emk - 6.*k**2)*L1 +                            &
     &        3.*(5.*emk**3 - 90.*emk**2*k**2 + 120.*emk*k**4 -         &
     &        16.*k**6)*L1)*q**6 -                                      &
     &        90.*(-1.*gam**2*(emk - 4.*k**2) +                         &
     &        2.*emk*(emk**2 - 12.*emk*k**2 + 8.*k**4))*L2*q**4*wq +    &
     &        3.*(-4.*gam**3 + 78.*emk*gam*(emk - 4.*k**2) +            &
     &        45.*emk**2*(emk - 6.*k**2)*L1 -                           &
     &        15.*gam**2*(emk - 2.*k**2)*L1)*q**2*wq**2 -               &
     &        30.*(emk**3*L2 - 1.*(-1. + En)**3*L6)*wq**3))/gam
         
         rq(5) = (120.*k*q*(-15.*(5.*emk**3 - 3.*emk*gam**2 -           &
     &        20.*emk**2*k**2 +                                         &
     &        4.*gam**2*k**2 + 8.*emk*k**4)*L2*q**4 +                   &
     &        4.*(-11.*gam**3 + 19.*emk*gam*(3.*emk - 4.*k**2) +        &
     &        30.*emk**2*(emk - 2.*k**2)*L1 -                           &
     &        5.*gam**2*(3.*emk - 2.*k**2)*L1)*q**2*wq -                &
     &        15.*emk*(3.*emk**2 - 1.*gam**2)*L2*wq**2))/gam
         
         rq(4) = 15.*((5.*(32.*gam**3*(emk - 4.*k**2) -                 &
     &        40.*emk*gam*(emk**2 - 12.*emk*k**2 + 8.*k**4) -           &
     &        3.*gam**4*L1 +                                            &
     &        6.*gam**2*(3.*emk**2 - 24.*emk*k**2 + 8.*k**4)*L1 -       &
     &        15.*emk**2*(emk**2 - 16.*emk*k**2 + 16.*k**4)*L1)*q**4)/  &
     &        gam + (120.*emk*(emk**2*(emk - 8.*k**2) -                 &
     &        1.*gam**2*(emk - 4.*k**2))*L2*q**2*wq)/gam +              &
     &        (-104.*emk**3 + 104.*(-1. + En)**3 +                      &
     &        16.*(1. + emk - 1.*En)*gam**2 + 15.*gam**3*(L3 + L5) +    &
     &        30.*gam*(emk**2*L1 + (-1. + En)**2*L7) -                  &
     &        (45.*(emk**4*L1 + (-1. + En)**4*L7))/gam)*wq**2)
         
         rq(3) = (60.*k*q*(5.*(3.*gam**4 + 5.*emk**3*(3.*emk - 8.*k**2) &
     &        -6.*emk*gam**2*(3.*emk - 4.*k**2))*L2*q**2 -4.*emk*(38.   &
     &        *emk**2*gam - 22.*gam**3 + 15.*emk**3*L1 -15.*emk*gam**2  &
     &        *L1)*wq))/gam
         
         rq(2) = (30.*((22.*gam**5 + 50.*emk**3*gam*(emk - 8.*k**2) -   &
     &        80.*emk*gam**3*(emk - 4.*k**2) +                          &
     &        15.*emk**4*(emk - 10.*k**2)*L1 -                          &
     &        30.*emk**2*gam**2*(emk - 6.*k**2)*L1 +                    &
     &        15.*gam**4*(emk - 2.*k**2)*L1)*q**2 -                     &
     &        4.*(3.*emk**5*L2 - 5.*emk**3*gam**2*L2 +                  &
     &        (-1. + En)**3*(-3.*(-1. + En)**2 + 5.*gam**2)*L6 +        &
     &        2.*gam**5*(L4 + L8))*wq))/gam
         
         rq(1) = (-900.*emk*(emk**2 - 1.*gam**2)**2*k*L2*q)/gam
         
         rq(0) = -300.*emk**5 + 300.*(-1. + En)**5 +                    &
     &        800.*(emk**3 - 1.*(-1. + En)**3)*gam**2 -                 &
     &        660.*(1. + emk - 1.*En)*gam**4 + 75.*gam**5*(L3 + L5) -   &
     &        225.*gam**3*(emk**2*L1 + (-1. + En)**2*L7) +              &
     &        225.*gam*(emk**4*L1 + (-1. + En)**4*L7) +                 &
     &        (-75.*emk**6*L1 - 75.*(-1. + En)**6*L7)/gam
         
      ELSE
         iqmax = 6

         Amp = 7.d0/(120.d0*q*(gam + wq)**5*(8.d0*gam**2 - 5.d0*gam*wq +&
     &        wq**2))
!        + i*delta
         a1 = emk - (2.d0*k*q + q**2 + wq) - gam 
!        - i*delta         
         a2 = emk + (2.d0*k*q - q**2 + wq) + gam 
!        - i*delta
         a3 = emk + 2.d0*k*q - q**2 
!        + i*delta
         a4 = emk - 2.d0*k*q - q**2 

!        + i*delta
         b1 = -1.d0 + (En - wq) - gam 
!        - i*delta
         b2 = -1.d0 + (En + wq) + gam 
!        + i*delta
         b3 = En - 1.d0 
!        - i*delta
         b4 = En - 1.d0 

!        L1 = Log[a1] - Log[a2] - Log[b1] + Log[b2];
         L1 = LOG(ABS((a1/a2)*(b2/b1))) + Logi(a1,1) - Logi(a2,-1) -    &
     &        Logi(b1,1) + Logi(b2,-1)
!        L2 = Log[-a1] - Log[a2] + Log[a3] - Log[-a4];
         L2 = LOG(ABS((a1/a2)*(a3/a4))) + Logi(-a1,-1) - Logi(a2,-1) +  &
     &        Logi(a3,-1) - Logi(-a4,-1)
!        L3 = Log[-a1] + Log[a2] - Log[a3] - Log[-a4];
         L3 = LOG(ABS((a1/a4)*(a2/a3))) + Logi(-a1,-1) + Logi(a2,-1) -  &
     &        Logi(a3,-1) - Logi(-a4,-1)
!        L4 = -Log[-b1] + Log[b2] + Log[-b3] - Log[b4]
         L4 = LOG(ABS(b2/b1)) - Logi(-b1,-1) + Logi(b2,-1) +            &
     &        Logi(-b3,-1) - Logi(b4,-1)
!        L5 = Log[-b1] + Log[b2] - Log[-b3] - Log[b4];
         IF(ABS(b3).lt.Zero) THEN
            L5 = 0.d0
         ELSE
            L5 = LOG(ABS((b1/b3)*(b2/b4))) + Logi(-b1,-1) + Logi(b2,-1) &
     &           - Logi(-b3,-1) - Logi(b4,-1)
         END IF
         
         rq(6) = 300.*gam**6*L1

         rq(5) = 120.*gam**5*(-11.*emk + 11.*(-1. + En + q**2) +        &
     &        8.*L1*wq)

         rq(4) = 60.*gam**4*(-38.*(1. + emk - 1.*En)*wq +               &
     &        2.*q**2*(15.*(emk - 2.*k**2)*L2 + 19.*wq) -               &
     &        5.*(3.*(-2. + En)*En*L4 + 2.*k*(-5. + 6.*L3)*q*           &
     &        (-1.*emk + q**2) + 3.*(L4 + L2*(emk**2 + q**4) -          &
     &        1.*L1*wq**2)))
         
         rq(3) = 160.*gam**3*(-10.*(-1. + En)**3 + 10.*(emk - 1.*q**2)* &
     &        (emk**2 - 2.*(emk - 6.*k**2)*q**2 + q**4) +               &
     &        66.*k*q*(emk - 1.*q**2)*wq +                              &
     &        3.*(1. + emk - 1.*En - 1.*q**2)*wq**2)

         rq(2) = 60.*gam**2*(15.*(-2. + En)*En*(2. + (-2. + En)*En)*L4 +&
     &        15.*(L4 + L2*(emk**4 + q**8)) + (emk**3*(44. - 40.*L3) +  &
     &        4.*(-1. + En)**3*(-11. + 10.*L5))*wq + 30.*(emk**2*L2 +   &
     &        (-1. + En)**2*L4)*wq**2 + 32.*(1. + emk - 1.*En)*wq**3 -  &
     &        5.*L1*wq**4 + 40.*k*q**5*(emk*(3. - 9.*L3) +              &
     &        4.*k**2*(-1. + 3.*L3) + 6.*L2*wq) - 4.*q**6*(15.*(emk -   &
     &        6.*k**2)*L2 + (11. - 10.*L3)*wq) + 6.*q**4*(5.*           &
     &        (3.*emk**2 - 24.*emk*k**2 + 8.*k**4)*L2 -                 &
     &        2.*(emk - 4.*k**2)*(-11. + 10.*L3)*wq + 5.*L2*wq**2) +    &
     &        4.*k*q**3*(10.*emk*(3.*emk - 4.*k**2)*(-1. + 3.*L3) -     &
     &        40.*(3.*emk - 2.*k**2)*L2*wq + (-77. + 30.*L3)*wq**2) +   &
     &        4.*k*q*(10.*(-1. + 3.*L3)*(-1.*emk**3 + q**6) +           &
     &        60.*emk**2*L2*wq - 1.*emk*(-77. + 30.*L3)*wq**2) +        &
     &        4.*q**2*(-15.*emk**2*(emk - 6.*k**2)*L2 + 3.*emk*(emk -   &
     &        4.*k**2)*(-11. + 10.*L3)*wq - 15.*(emk - 2.*k**2)*        &
     &        L2*wq**2 - 8.*wq**3))

         rq(1) = 120.*gam*(5.*(-1. + En)**5 - 5.*(emk - 1.*q**2)*       &
     &        (emk**4 - 4.*emk**2*(emk - 10.*k**2)*q**2 +               &
     &        2.*(3.*emk**2 - 40.*emk*k**2 + 40.*k**4)*q**4 -           &
     &        4.*(emk - 10.*k**2)*q**6 + q**8) + 152.*k*q*              &
     &        (-1.*emk + q**2)*(emk**2 - 2.*(emk - 2.*k**2)*q**2 +      &
     &        q**4)*wq - 26.*(-1.*(-1. + En)**3 + (emk - 1.*q**2)*      &
     &        (emk**2 - 2.*(emk - 6.*k**2)*q**2 + q**4))*wq**2 +        &
     &        56.*k*q*(-1.*emk + q**2)*wq**3 - 1.*(1. + emk - 1.*En -   &
     &        1.*q**2)*wq**4)

!        Pieces of rq(0)
         
         rq00  = 40.*wq*(15.*(-1. + En)**5 - 15.*(emk - 1.*q**2)*(emk**4&
     &        - 4.*emk**2*(emk - 10.*k**2)*q**2 + 2.*(3.*emk**2 - 40.   &
     &        *emk*k**2 + 40.*k**4)*q**4 - 4.*(emk - 10.*k**2)*q**6 + q &
     &        **8) + 516.*k*q*(-1.*emk + q**2)*(emk**2 - 2.*(emk - 2.*k &
     &        **2)*q**2 + q**4)*wq - 104.*(-1.*(-1. + En)**3 + (emk - 1.&
     &        *q**2)*(emk**2 - 2.*(emk - 6.*k**2)*q**2 + q**4))*wq**2 + &
     &        291.*k*q*(-1.*emk + q**2)*wq**3 - 15.*(1. + emk - 1.*En - &
     &        1.*q**2)*wq**4)
         
         rq0L1 = 60.*L1*wq**6

         rq0L2 = 60.*L2*(-5.*(emk**2 - 2.*(emk - 2.*k**2)*q**2 + q**4)* &
     &        (emk**4 - 4.*emk**2*(emk - 14.*k**2)*q**2 + 2.*(3.*emk**2 &
     &        - 56.*emk*k**2 + 8.*k**4)*q**4 -4.*(emk - 14.*k**2)*q**6 +&
     &        q**8) -48.*k*q*(5.*emk**4 - 20.*emk**2*(emk - 2.*k**2)*q  &
     &        **2 +2.*(15.*emk**2 - 40.*emk*k**2 + 8.*k**4)*q**4 - 20.  &
     &        *(emk - 2.*k**2)*q**6 + 5.*q**8)*wq -45.*(emk**4 - 4.*emk &
     &        **2*(emk - 6.*k**2)*q**2 + 2.*(3.*emk**2 - 24.*emk*k**2 + &
     &        8.*k**4)*q**4 -4.*(emk - 6.*k**2)*q**6 + q**8)*wq**2 -80. &
     &        *k*q*(3.*emk**2 - 6.*emk*q**2 + 4.*k**2*q**2 + 3.*q**4)*wq&
     &        **3 -15.*(emk**2 - 2.*(emk - 2.*k**2)*q**2 + q**4)*wq**4)
         
         rq0L3 = 240.*L3*(emk - 1.*q**2)*(5.*k*q*(emk**2 - 2.*(emk - 6. &
     &        *k**2)*q**2 + q**4)*(3.*emk**2 - 6.*emk*q**2 + 4.*k**2*q  &
     &        **2 + 3.*q**4) +6.*(emk**4 - 4.*emk**2*(emk - 10.*k**2)*q &
     &        **2 + 2.*(3.*emk**2 - 40.*emk*k**2 + 40.*k**4)*q**4 -4.   &
     &        *(emk - 10.*k**2)*q**6 + q**8)*wq + 90.*k*q*(emk**2 - 2.  &
     &        *(emk - 2.*k**2)*q**2 + q**4)*wq**2 +10.*(emk**2 - 2.*(emk&
     &        - 6.*k**2)*q**2 + q**4)*wq**3 + 15.*k*q*wq**4)

         rq0L4 = 300.*(-1. + En)**2*L4*(-1.*(-1. + En)**4 - 9.*(-1. + En&
     &        )**2*wq**2 - 3.*wq**4)

         rq0L5 = -480.*(-1. + En)**3*L5*wq*(3.*(-1. + En)**2 + 5.*wq**2)
         
         rq(0) = rq00 + rq0L1 + rq0L2 + rq0L3 + rq0L4 + rq0L5

      END IF
      bpr1 = 0.d0
      DO iqexp = 0, iqmax
         bpr1 = bpr1 + rq(iqexp)
      END DO
      bpr1 = Amp*bpr1

      RETURN
      END
