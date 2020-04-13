! Written by J. Kas 8/1/2014
MODULE xmuio
  USE IOMOD
  REAL(8) xmuio_efermi, xmuio_gamch, xmuio_ne
  ! read columns of xmu and return whichever columns are requested, along with a few other bits of information.

  CONTAINS

    SUBROUTINE ReadXMU(file,ne,omega,xe,xk,mu,mu0,chi)
      CHARACTER*(*), INTENT(IN) :: file
      REAL(8),OPTIONAL :: omega(:), xe(:), xk(:), mu(:), mu0(:), chi(:)
      INTEGER,INTENT(OUT) :: ne
      INTEGER ihead, iE
      REAL(8) omegatmp, xetmp, xktmp, mutmp, mu0tmp, chitmp
      CHARACTER(LEN=300) word1, word2
      
      CALL OpenFl(file)
      ne = NumberOfLines(file)
      DO iE = 1, ne
         CALL ReadData(file, Double1 = omegatmp, Double2 = xetmp, &
              & Double3 = xktmp, Double4 = mutmp, Double5 = mu0tmp, &
              & Double6 = chitmp)
         IF(PRESENT(omega)) omega(iE) = omegatmp
         IF(PRESENT(xe)) xe(iE) = xetmp
         IF(PRESENT(xk)) xk(iE) = xktmp
         IF(PRESENT(mu)) mu(iE) = mutmp
         IF(PRESENT(mu0)) mu0(iE) = mu0tmp
         IF(PRESENT(chi)) chi(iE) = chitmp
      END DO
    END SUBROUTINE ReadXMU
  END MODULE xmuio
         
            
