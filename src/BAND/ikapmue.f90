      FUNCTION IKAPMUE(KAPPA,MUEM05)
!   *  INDEXING OF MATRIX-ELEMENTS:                                    *
!   *  I = 2*L*(J+1/2) + J + MUE + 1                                   *

      IMPLICIT NONE

      INTEGER KAPPA,MUEM05
      INTEGER IKAPMUE
      INTEGER IABS
      INTEGER JP05,L

      JP05 = IABS(KAPPA)
      IF ( KAPPA.LT.0 ) THEN
         L = -KAPPA - 1
      ELSE
         L = KAPPA
      END IF
      IKAPMUE = 2*L*JP05 + JP05 + MUEM05 + 1
      END
