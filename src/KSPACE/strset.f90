!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: strset.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/02/17 07:39:12 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE STRSET(KX,KY,KZ,GNR,TAUKINV,P)
!
!   ********************************************************************
!   *                                                                  *
!   *   SETS UP NON-RELATIVISTIC STRUCTURE CONSTANTS  GNR              *
!   *   FOR  REAL SPHERICAL HARMONICS AND CONVERTS TO  RELAT.  G'S     *
!   *                                                                  *
!   *   NOTE:  - GRNR  IS STORED IN  SUPPLIED WORKSPACE                *
!   *          - ONLY THE ELEMENTS  G(I,J) FOR I<=J ARE EVALUATED      *
!   *            THOSE FOR J>I ARE OBTAINED FROM THE FACT THAT         *
!   *            G(I,J) IS HERMITIAN:        G(I,J) = G(J,I)*          *
!   *          - FOR IREL<2:                 G(I,J) = GNR(I,J)         *
!   *          - G(I,J) IS STORED AS -G IN TAUKINV  !!!!!              *
!   *                                                                  *
!   *  15/12/94  HE  arg. list modified                                *
!   ********************************************************************
!

      use workstrfacs
      use boundaries
      use controls,only : cplxylm,irel,singleprec
      IMPLICIT NONE


! PARAMETER definitions
      COMPLEX*16 C0,CI
      PARAMETER (C0=(0.0D0,0.0D0),CI=(0.0D0,1.0D0))
!
! Dummy arguments
      REAL*8, intent(in)     :: KX,KY,KZ
      COMPLEX*16, intent(in) :: P
      COMPLEX*16,intent(out) :: GNR(NKKRMAX,NKKRMAX),TAUKINV(NKKRMAX,NKKRMAX)
!
! Local variables
      COMPLEX*16 CIL1ML2,CSUM,CSUM1,CSUM2,DLLMMKE(NLLMMMAX,NQQPMAX)
      INTEGER I,I1,I2,IG123,IKM1,IKM2,ILM1,ILM2,IOFF1,IOFF2,IQ1,IQ2,    &
     &        IQQP,IS,J1,J2,JKM1,JKM10,JKM2,JLM1,JLM10,JLM2,JQ1,JQ2,    &
     &        JQQP,LM,LM1,LM1LM2,LM2,LM3,N1,N2
      INTEGER lm2max  !KJ added for complex spherical harmonics
      COMPLEX DLLMMKEb(NLLMMMAX,NQQPMAX)
!
	  gnr=c0
	  taukinv=c0

      if(singleprec) then
         CALL STRBBDD2(DLLMMKEb,real(KX),real(KY),real(KZ))
         DLLMMKE=dcmplx(DLLMMKEb)
      else
         CALL STRBBDD(DLLMMKE,KX,KY,KZ)
      endif


!	if(kx+ky+kz.gt.0.1) stop

!
      IF ( IREL.LT.2 ) THEN
!
!     NON-RELATIVISTIC CALCULATION
!     ============================
!
         DO IQQP = 1,NQQP
            IQ1 = IJQ(1,IQQP)/100
            IQ2 = IJQ(1,IQQP) - 100*IQ1
            IOFF2 = IND0Q(IQ2)
            IOFF1 = IND0Q(IQ1)
            !write(*,*) 'iqqp,iq1,iq2,ioff1,ioff2',iqqp,iq1,iq2,ioff1,ioff2
!
            IG123 = 0
            LM1LM2 = 0
            DO LM1 = 1,NL**2
                 if (cplxylm) then !this if-block added !KJ
                    lm2max=nl**2
                 else
                    lm2max=lm1  !as in the original statement
                 endif
               DO LM2 = 1,lm2max  !KJ   LM1
                  LM1LM2 = LM1LM2 + 1
                  CSUM = C0
!				  write(*,*) 'lm1,lm2',lm1,lm2
                  DO I = 1,NGNT(LM1LM2)
                     IG123 = IG123 + 1
                     LM3 = IGNT(IG123)
                     CSUM = CSUM + GNT(IG123)*DLLMMKE(LM3,IQQP)
!					 if((ioff1+lm1).eq.1 .and. (ioff2+lm2).eq.1) write(*,*) 'csum',gnt(ig123),lm3,dllmmke(lm3,iqqp),csum
                  END DO
                  CIL1ML2 = CIPWL(LM1)/CIPWL(LM2)
                  TAUKINV(IOFF1+LM1,IOFF2+LM2) = -CSUM*CIL1ML2
                  TAUKINV(IOFF1+LM2,IOFF2+LM1) = -CSUM/CIL1ML2
!				  if((ioff1+lm2).eq.1 .and. (ioff2+lm1).eq.1) write(*,*) 'ataukinv(1,1)',csum,cil1ml2,-CSUM/CIL1ML2,TAUKINV(IOFF1+LM2,IOFF2+LM1)
!				  if((ioff1+lm1).eq.1 .and. (ioff2+lm2).eq.1) write(*,*) 'btaukinv(1,1)',csum,cil1ml2,-CSUM/CIL1ML2,TAUKINV(IOFF1+LM1,IOFF2+LM2)
				  
               END DO
            END DO
!
            IF ( IQQP.EQ.1 ) THEN
               DO LM = 1,NLM
                  TAUKINV(IOFF1+LM,IOFF2+LM) = TAUKINV(IOFF1+LM,IOFF2+LM) - CI*P
!				  if((ioff1+lm).eq.1 .and. (ioff2+lm).eq.1) write(*,*) 'ctaukinv(1,1)',ci,p,- CI*P,TAUKINV(IOFF1+LM,IOFF2+LM)
               END DO
            END IF
!
!     ---------------------------------------------------
!                          COPY
!     ---------------------------------------------------
!
            DO JQQP = 2,NIJQ(IQQP)
               JQ1 = IJQ(JQQP,IQQP)/100
               JQ2 = IJQ(JQQP,IQQP) - 100*JQ1
               JLM2 = IND0Q(JQ2)
               JLM10 = IND0Q(JQ1)
               DO ILM2 = IND0Q(IQ2) + 1,IND0Q(IQ2) + NKMQ(IQ2)
                  JLM2 = JLM2 + 1
                  JLM1 = JLM10
                  DO ILM1 = IND0Q(IQ1) + 1,IND0Q(IQ1) + NKMQ(IQ1)
                     JLM1 = JLM1 + 1
                     TAUKINV(JLM1,JLM2) = TAUKINV(ILM1,ILM2)
!					 if(jlm1.eq.1 .and. jlm2.eq.1) write(*,*) 'dtaukinv(1,1)',taukinv(ilm1,ilm2),taukinv(jlm1,jlm2)
                  END DO
               END DO
            END DO
!stop
         END DO
!	do i=1,nqqp
!	write(*,*) i,nijq(i),((ijq(i1,i)/100,ijq(i1,i)-100*(ijq(i1,i)/100 )),i1=1,nijq(i))
!	enddo
!	stop
!call writematrixdblefree(taukinv,nkkrmax,'taukinv.txt')
!stop
         RETURN
!
      ELSE
!
!     RELATIVISTIC CALCULATION
!     ========================
!
!     ---------------------------------------------------
!     set up non-relativistic str.const.  GNR  for every
!     representative q-qp-block  IQQP, then transfrom to
!     relativ. representation and finally copy result
!     to the  NIJQ(IQQP)-1  equivalent q-qp-blocks  JQQP
!     ---------------------------------------------------
!
         DO IQQP = 1,NQQP
            IQ1 = IJQ(1,IQQP)/100
            IQ2 = IJQ(1,IQQP) - 100*IQ1
!
            IG123 = 0
            LM1LM2 = 0
            DO LM1 = 1,NL**2
                 if (cplxylm) then !this if-block added !KJ
                    lm2max=nl**2
                 else
                    lm2max=lm1  !as in the original statement
                 endif
               DO LM2 = 1,lm2max  !KJ   LM1
                  LM1LM2 = LM1LM2 + 1
                  CSUM = 0.0D0
!
                  DO I = 1,NGNT(LM1LM2)
                     IG123 = IG123 + 1
                     LM3 = IGNT(IG123)
                     CSUM = CSUM + GNT(IG123)*DLLMMKE(LM3,IQQP)
                  END DO
!
                  CIL1ML2 = CIPWL(LM1)/CIPWL(LM2)
                  GNR(LM1,LM2) = CSUM*CIL1ML2
                  GNR(LM1,LM2) = CSUM
                  GNR(LM2,LM1) = CSUM
!30    GNR(LM2,LM1) = CSUM / CIL1ML2
!k
                  GNR(LM1,LM2) = CSUM*CIL1ML2
                  GNR(LM2,LM1) = CSUM/CIL1ML2
!k
!
               END DO
            END DO
!
!                                 SITE-DIAGONAL BLOCK
            IF ( IQQP.EQ.1 ) THEN
               DO LM = 1,NLM
                  GNR(LM,LM) = GNR(LM,LM) + CI*P
               END DO
            END IF
!
!
!     ---------------------------------------------------
!                          TRANSFORM
!     ---------------------------------------------------
!
            IOFF2 = IND0Q(IQ2)
            IOFF1 = IND0Q(IQ1)
            DO IKM2 = 1,NKMQ(IQ2)
               DO IKM1 = 1,NKMQ(IQ1)
!
                  CSUM1 = C0
                  DO IS = 1,2
                     N1 = NRREL(IS,IKM1)
                     N2 = NRREL(IS,IKM2)
                     DO I1 = 1,N1
                        J1 = IRREL(I1,IS,IKM1)
!
                        CSUM2 = C0
                        DO I2 = 1,N2
                           J2 = IRREL(I2,IS,IKM2)
                           CSUM2 = CSUM2 + GNR(J1,J2)*SRREL(I2,IS,IKM2)
                        END DO
!
                        CSUM1 = CSUM1 + DCONJG(SRREL(I1,IS,IKM1))*CSUM2
                     END DO
                  END DO
                  TAUKINV(IOFF1+IKM1,IOFF2+IKM2) = -CSUM1
               END DO
            END DO
!
!
!     ---------------------------------------------------
!                          COPY
!     ---------------------------------------------------
!
            DO JQQP = 2,NIJQ(IQQP)
               JQ1 = IJQ(JQQP,IQQP)/100
               JQ2 = IJQ(JQQP,IQQP) - 100*JQ1
               JKM2 = IND0Q(JQ2)
               JKM10 = IND0Q(JQ1)
               DO IKM2 = IND0Q(IQ2) + 1,IND0Q(IQ2) + NKMQ(IQ2)
                  JKM2 = JKM2 + 1
                  JKM1 = JKM10
                  DO IKM1 = IND0Q(IQ1) + 1,IND0Q(IQ1) + NKMQ(IQ1)
                     JKM1 = JKM1 + 1
                     TAUKINV(JKM1,JKM2) = TAUKINV(IKM1,IKM2)
                  END DO
               END DO
!
            END DO
!
         END DO
!
      END IF
! --------------------------------------- NON.-REL. / REL. CALCULATION
!

!call writematrixdblefree(taukinv,nkkrmax,'taukinv.txt')
      END
