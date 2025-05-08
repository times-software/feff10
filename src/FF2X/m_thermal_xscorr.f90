!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $Revision:  1.0 $
! $Author:  tts $
! $Date: 2018/11/1 $
!
! Thermal version of xscorr.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE m_thermal_xscorr
  USE constants
  USE dimsmod, only: nex
  IMPLICIT NONE
  PUBLIC :: thermal_xscorr, sommerfeld_xscorr, fermilorenT, sommerfeld_xscorr0
  PRIVATE :: cauchy, reg_integration, reciprocal_integration, fermiDirac, &
    &         fermilorenT_terp, cchi_eff, cchi_fine, lorenz0, fermilorenH, fermilorenTH
CONTAINS

  SUBROUTINE thermal_xscorr(ispec, emxs, ne1, ne, ik0, xsec, xsnorm, chia, &
    &                       vrcorr, vicorr, cchi, electronic_temperature)
    IMPLICIT none

    INTEGER, INTENT(IN) :: ispec, ne1, ne, ik0
    REAL(8), INTENT(IN) :: electronic_temperature, vrcorr, xsnorm(nex), vicorr
    COMPLEX*16, INTENT(IN) :: xsec(nex), chia(nex)
    COMPLEX*16, INTENT(INOUT) :: emxs(nex)
    COMPLEX*16, INTENT(OUT) :: cchi(nex)

    INTEGER :: polex, idx, ic0, ip, nc
    INTEGER :: ne2, ne4, npole
    REAL(8) :: tk, xloss, efermi
    REAL(8) :: omega(nex)
    COMPLEX*16 :: corr, corr1, corr2, residue, leg3
    COMPLEX*16 :: ec(nex), fc(nex), xmu(nex), pole
    COMPLEX*16 :: xmu1(nex), xmu2(nex), xmu3(nex), xmu4(nex)
    REAL(8), PARAMETER :: eps4 = 1.0d-4
    LOGICAL, PARAMETER :: print_out = .TRUE.
    EXTERNAL :: lorenz, astep
    COMPLEX*16 :: lorenz
    REAL(8) :: astep

    INTEGER :: ie, ic
    REAL(8) :: dele, w1, w2 ,w3
    REAL(8) :: ra, rb, cheight
    COMPLEX*16 :: xmu0, ff(nex)
    COMPLEX*16 :: bb, c1, f1, f2, e1, e2, z1, z2
    COMPLEX*16 :: dummy, dummy2, aa

    REAL(8) :: minee, maxee, dee
    COMPLEX*16 :: ee, f4
    INTEGER :: extra
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! tk = xloss/pi ! Critical temp in Hatree
    ! tk = 1e-6/hart
    tk = electronic_temperature/hart
    ne2 = 10
    extra = 1
    npole = ne-2*ne1-ne2-extra
    ne4 = ne1

    efermi = DBLE(emxs(ne-extra))
    xloss = DIMAG(emxs(ne4+1))
    cheight = DIMAG(emxs(1))

    ! IF (ABS(cheight-xloss).LT.eps4) STOP "ERROR: contour height <= xloss"

    DO ie = 1, ne
      xmu(ie) = xsec(ie) + xsnorm(ie)*chia(ie)
    ENDDO


    OPEN(UNIT=155, file='curve.dat', status='REPLACE')
    OPEN(UNIT=150, file="prexmu.dat", status='REPLACE')
    OPEN(UNIT=149, file="raw.dat", status='REPLACE')
    OPEN(UNIT=153, file="leg1.dat", status='REPLACE')
    OPEN(UNIT=151, file='residue.dat', status='REPLACE')
    OPEN(UNIT=152, file='contour.dat', status='REPLACE')
    OPEN(UNIT=1543, file='ratio.dat', status='REPLACE')
    WRITE(149,*) "Temperature (Hatree)", tk
    WRITE(149,*) "Electronic Temperature (eV)", electronic_temperature
    WRITE(149,*) "xloss = ", xloss, " Hatree"
    WRITE(149,*) "Chemical potential = ", efermi*hart, " eV with shift", vrcorr*hart
    WRITE(149,*) "Number of poles = ", npole
    WRITE(149, *) 'Omega(Hart)  \t  Re CCHI  \t   Im CCHI \t  1-Fermi \t  Re xmu0  \t  Im xmu0'

    DO ie =1, ne
      WRITE(155, '(I4,20E20.10E3)') ie, DBLE(emxs(ie)), DIMAG(emxs(ie))
    ENDDO
    CLOSE(155)


    !---------- Construct the integration contour
    nc = 0
    DO ie = 1, ne2
        nc = nc + 1
        ec(nc) = emxs(2*ne1+ie)
        fc(nc) = xmu(2*ne1+ie)
    ENDDO

    DO ie=1, ne4
      nc = nc + 1
      ec(nc) = emxs(ie)
      fc(nc) = xmu(ie)
    ENDDO
    !--------- End of countour construction

    DO ie = 1, ne4  !---------- Real frequencies from ecv to edge + 5
      omega(ie) = DBLE(emxs(ie))
      xmu2(ie) = xmu(ie)     !--- xmu(omega+i*height)
      xmu3(ie) = xmu(ne4+ie) !--- xmu(omega+i*xloss)
    ENDDO           !---------- End of Real frequencies from ecv to edge + 5

    DO ie = 1, npole
      xmu4(ie) = xmu(2*ne1+ne2+ie) !--- xmu(efermi+i*En)
    ENDDO


    if (abs(vrcorr).gt.eps4) then
      ! account for the fermi level shift
      efermi = efermi - vrcorr
    endif

    DO ie = 1, ne4   !--------- Perform convolution
      bb = 1

      xmu0 = xmu3(ie)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Analytic Part (no pole omega+i*xloss)!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cchi(ie) = 0.d0
      !z1 = 7.662230927071878E-003
      CALL fermilorenT(omega(ie), tk, efermi, xloss, cchi(ie))
      ! CALL fermilorenT_terp(omega(ie), tk, efermi, xloss, cchi(ie))

      dummy = cchi(ie)
      cchi(ie) = xmu0*cchi(ie)
      IF (ispec.ne.2) cchi(ie) = xmu0 - cchi(ie) ! Absorption
      dummy = 1-dummy
      z1 = omega(ie)
      WRITE(149, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie)), DBLE(dummy), DBLE(xmu0), DIMAG(xmu0), DBLE(fermiDirac(z1,tk,efermi))

      !!!!!!!!!!!!!!!!!!!!!!
      !! RESIDUE of Fermi !!
      !!!!!!!!!!!!!!!!!!!!!!
      ! At zero temperature this will become last leg
      residue = 0.d0
      dele = omega(ie) - efermi
      DO ic = 1, npole
        z1 = emxs(2*ne1+ne2+ic) ! pole
        if (abs(vrcorr).gt.eps4) z1 = z1 - vrcorr ! Account for fermi shift
        f1 = xmu4(ic)-xmu0
          residue = residue + f1 * xloss/((z1-omega(ie))**2 + xloss**2)
      ENDDO
      residue = -residue*(2.d0*coni*tk)

      IF (ispec.ne.2) residue = - residue ! Absorption
      WRITE(151, '(20E20.10E3)') DBLE(omega(ie)), DBLE(residue), DIMAG(residue), DBLE(efermi)

      !!!!!!!!!!!!!!!
      !! First leg !!
      !!!!!!!!!!!!!!!
      ! Integrate from 0 to cheight with x in R-space.

      corr1 = 0.d0
      WRITE(153, '(20E20.10E3)') DBLE(omega(ie)), DBLE(corr1), DIMAG(corr1)

      !!!!!!!!!!!!!!!!!
      !! Second  leg !!
      !!!!!!!!!!!!!!!!!
      corr2 = 0.d0
      CALL cchi_eff(emxs(1:ne4), xmu2, ne4, efermi, tk, xmu0, xloss, omega(ie), ispec, corr2)
      ! CALL fermilorenH(emxs(1:ne4), ne4, efermi, tk, xloss, omega(ie), ispec, dummy2)
      ! CALL fermilorenTH(omega(ie),tk,efermi,xloss,cheight,dummy2)
      ! WRITE(152, '(20E20.10E3)') DBLE(omega(ie)), DBLE(corr2), DIMAG(corr2), DBLE(xmu2(ie)-xmu0), DIMAG(xmu2(ie)-xmu0)
      ! WRITE(1543, '(20E20.10E3)') DBLE(omega(ie)), DBLE(corr2), DIMAG(corr2), &
      !                            & DBLE(cchi(ie)), DIMAG(cchi(ie)), &
      !                            & DBLE(dummy), DBLE(xmu0), DIMAG(xmu0), DBLE(dummy2), DIMAG(dummy2)
      WRITE(1543, '(20E20.10E3)') DBLE(omega(ie)), DBLE(1-fermiDirac(omega(ie)+coni*0.d0,tk,efermi)), &
                                 &  DBLE(dummy),&
                                 &  DBLE(cchi(ie)), DIMAG(cchi(ie)), &
                                 &  DBLE(corr2), DIMAG(corr2), &
                                 &  DBLE(cauchy(omega(ie)+coni*0.0d0,efermi,xloss)), &
                                 & DIMAG(xmu2(ie)-xmu0)
      !!!!!!!!!!!!!!!!!!!!!
      !! Combine results !!
      !!!!!!!!!!!!!!!!!!!!!
      corr = corr1 + corr2
      cchi(ie) = cchi(ie) + corr + residue
      WRITE(150, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie))

      ! Return the result of convolution minus bare value
      cchi(ie) = cchi(ie) - xmu2(ie)
    ENDDO !------- End of convolution

    CLOSE(149)
    CLOSE(150)
    CLOSE(151)
    CLOSE(152)
    CLOSE(153)
    CLOSE(1543)

    RETURN
  END SUBROUTINE

  SUBROUTINE fermilorenT(omega, tk, efermi, xloss, result)
    ! Perform convolution between fermi function
    ! and lorentzian.
    USE constants
    IMPLICIT none
    ! Input:
    !   omega  - energy in Hatree
    !   tk     - electronic temperature in Hatree
    !   efermi - chemical potential in Hatree
    !   xloss  - broadening in Hatree
    REAL(8), INTENT(IN) :: tk, efermi, xloss, omega

    ! Output:
    !   result
    COMPLEX*16, INTENT(OUT) :: result

    ! Local varaibles:
    !   z1        - integration variable
    !   z2        - another name for omega
    !   dummy     - temporary variable
    !   Nanalytic - integration grid size
    INTEGER :: ic, Nanalytic
    REAL(8) :: domega, omegaMax, omegaMin, z1, z2, dummy
    !REAL(8), PARAMETER :: pi = 3.141592653589793d0

    omegaMax = 1.0
    omegaMin = 0.d0
    Nanalytic = int(20000)
    domega = (omegaMax-omegaMin)/DBLE(Nanalytic-1)

    result = 0.d0
    DO ic = 1, Nanalytic-1
      ! From -Inf(Avoid) to -1 [Change of variable 1/x]
      z1 = 1/DBLE(-omegaMax+DBLE(ic-1)*domega)
      z2 = omega
      IF (DBLE((z1-efermi)/tk).lt.100) THEN
        dummy = 1/(EXP((z1-efermi)/tk)+1)
        result = result+dummy*(xloss/pi)/((z1-z2)**2+xloss**2)*(z1**2)
      ENDIF
    ENDDO

    DO ic = 1, Nanalytic-1
      ! From -1 to 0
      z1 = -omegaMax+DBLE(ic-1)*domega
      z2 = omega
      dummy = 1/(EXP((z1-efermi)/tk)+1)
      result = result+dummy*(xloss/pi)/((z1-z2)**2 + xloss**2)
      ! From 0 to 1
      z1 = omegaMin+DBLE(ic-1)*domega
      z2 = omega
      dummy = 1/(EXP((z1-efermi)/tk)+1)
      result = result+dummy*(xloss/pi)/((z1-z2)**2 + xloss**2)
    ENDDO


    DO ic = Nanalytic, 2, -1
      ! From 1 to Inf(Avoid) [Change of variable 1/x]
      z1 = 1/DBLE(omegaMax-DBLE(Nanalytic-ic)*domega)
      z2 = omega
      IF (DBLE((z1-efermi)/tk).lt.100) THEN
        dummy = 1/(EXP((z1-efermi)/tk)+1)
        result = result+dummy*(xloss/pi)/((z1-z2)**2+xloss**2)*(z1**2)
      ENDIF
    ENDDO
    result = result*domega
    RETURN
  END SUBROUTINE fermilorenT

  SUBROUTINE fermilorenT_terp(omega, tk, efermi, xloss, cchi)
    ! Perform convolution between fermi function and lorentzian.
    ! This gets very close to the zero temperature limit but there is a kink
    ! at very low omega. [NEED TO FIX THIS]
    IMPLICIT none
    ! Input:
    !   omega  - energy in Hatree
    !   tk     - electronic temperature in Hatree
    !   efermi - chemical potential in Hatree
    !   xloss  - broadening in Hatree
    REAL(8), INTENT(IN) :: tk, efermi, xloss, omega

    ! Output:
    !   result
    COMPLEX*16, INTENT(OUT) :: cchi
    ! COMPLEX*16, INTENT(INOUT) :: result
    ! COMPLEX*16, INTENT(OUT) :: cchi
    ! Local varaibles:
    !   z1        - integration variable
    !   z1p        - integration variable
    !   z2        - another name for omega
    !   dummy     - temporary variable
    !   dummy2    - temporary variable
    !   Nanalytic - integration grid size
    !   Nxmu      -
    !   xmu_p3
    !   xmu_m3
    !   window_size
    !   dxmu
    INTEGER :: ic
    COMPLEX*16 :: z1, z1p
    REAL(8) :: z2, dummy, dummy2, left, center, right
    REAL(8) :: xmu_m3, xmu_p3, dxmu
    REAL(8), PARAMETER :: window_size = 12.d0
    INTEGER, PARAMETER :: Nxmu = INT(1E5)

    xmu_m3 = efermi - window_size*tk
    xmu_p3 = efermi + window_size*tk
    dxmu = (2.d0*window_size*tk)/DBLE(Nxmu-1)
    left = 0.d0
    right = 0.d0
    center = 0.d0
    cchi = 0.d0

    ! -inf to xmu_m3
    left = left + reciprocal_integration(1.d0/MIN(-1.d0, xmu_m3), Nxmu, -1, omega, xloss, tk, efermi) ! -inf to xmu_m3 or -1
    IF (xmu_m3.GT.-1) THEN
      left = left +  reg_integration(-1.d0, xmu_m3, Nxmu, omega, xloss, tk, efermi) ! -1 to xmu_m3
    ENDIF

    ! xmu_m3 to xmu_p3
    DO ic = 1, Nxmu-1
      z1 = xmu_m3 + DBLE(ic)*dxmu
      z1p = xmu_m3 + DBLE(ic-1)*dxmu
      z2 = omega
      ! dummy = 1/(EXP((z1-efermi)/tk)+1)
      ! dummy2 = 1/(EXP((z1p-efermi)/tk)+1)
      dummy = fermiDirac(z1, tk, efermi)
      dummy2 = fermiDirac(z1p, tk, efermi)
      center = center + 0.5d0*(dummy*cauchy(z1,z2,xloss) + dummy2*cauchy(z1p,z2,xloss))
    ENDDO
    center = center*dxmu

    ! xmu_p3 to inf
    right = right + reciprocal_integration(1.d0/MAX(1.d0, xmu_p3), Nxmu, 1, omega, xloss, tk, efermi) ! xmu_m3 or 1 to inf
    IF (xmu_m3.GT.1) THEN
      right = right +  reg_integration(xmu_p3, 1.d0, Nxmu, omega, xloss, tk, efermi) ! xmu_p3 to 1
    ENDIF

    cchi = left + right + center
    RETURN
  END SUBROUTINE fermilorenT_terp

  SUBROUTINE fermilorenH(emxs, ne4, efermi, tk, xloss, omega, ispec, result)
    ! Interpolation of contour integral for cchi
    ! We interpolation fermi function near the efermi
    REAL(8), INTENT(IN) :: efermi, tk, xloss, omega
    COMPLEX*16, INTENT(IN) :: emxs(:)
    INTEGER, INTENT(IN) :: ne4, ispec

    COMPLEX*16, INTENT(OUT) :: result

    ! Local variables
    COMPLEX*16 ::  efermi_m3, efermi_p3, z1, z2, f1, f2, x1, x2, k1, k2
    COMPLEX*16 ::  de, e1, e2, aa
    INTEGER :: ind_m3, ind_p3, ie, ne4_terp
    REAL(8) :: defermi
    COMPLEX*16, ALLOCATABLE :: emxs_terp(:)
    INTEGER, PARAMETER :: nxmu = 10000
    REAL(8), PARAMETER :: window_size_terp = 10, eps4=1e-4
    LOGICAL :: BUG=.FALSE.

    defermi = (2*window_size_terp*tk)/DBLE(nxmu)
    efermi_m3 = efermi - window_size_terp*tk + coni*DIMAG(emxs(1))
    efermi_p3 = efermi + window_size_terp*tk + coni*DIMAG(emxs(1))

    ! Find where does xmu +/- window_terp lies on the contour
    ind_m3 = binarysearch(DBLE(efermi_m3), DBLE(emxs), ne4)-1
    ind_p3 = binarysearch(DBLE(efermi_p3), DBLE(emxs), ne4)

    IF (ind_m3.EQ.0) BUG=.TRUE.
    IF (BUG)  PRINT*, "ERROR: xmu_m3 less than minimum. ", &
                      & DBLE(efermi_m3*hart), DBLE(emxs(1)*hart)


    ne4_terp = ind_m3 + nxmu + (ne4-ind_p3+1)
    ALLOCATE(emxs_terp(ne4_terp))
    emxs_terp = 0.d0

    DO ie=1, ind_m3
      ! Copy vertical points and all points less than xmu-3kT
      emxs_terp(ie) = emxs(ie)
    ENDDO

    DO ie=ind_m3+1, ind_m3+nxmu
      ! Remove xmu +/- 3kT from original grid and replaced by interpolated points
      emxs_terp(ie) = efermi_m3 + (ie-ind_m3-1) * defermi
    ENDDO

    ! Insert xmu_p3 into the thing
    emxs_terp(ind_m3+nxmu+1) = efermi_p3

    DO ie=1, (ne4-ind_p3)
      !  Going back to original grid at ind_p3
      emxs_terp(ind_m3+nxmu+1+ie) = emxs(ie+ind_p3)
    ENDDO

    ! OPEN(UNIT=154, file='lorentzian.dat', status='OLD')
    result = 0.d0
    e1 = omega + coni*xloss
    e2 = omega - coni*xloss
    DO ie = 1, ne4_terp-1
      z1 = emxs_terp(ie)
      x1 = DBLE(z1)
      f1 = fermiDirac(x1,tk,efermi)
      IF (ispec.ne.2) f1 = 1.d0-f1
      z2 = emxs_terp(ie+1)
      x2 = DBLE(z2)
      f2 = fermiDirac(x2,tk,efermi)
      IF (ispec.ne.2) f2 = 1.d0-f2
      de = z2-z1
      !! result = result + f1*f2*cauchy(z1, omega, xloss)*de
      IF ((DBLE(z1).GT.DBLE(efermi_m3)).AND.(DBLE(z1).LT.DBLE(efermi_p3))) THEN
        !! Trapezoidal rule
        result = result + 0.5d0*de*(f1*cauchy(z1, omega, xloss) + &
                                  & f2*cauchy(z2,omega,xloss))
      ELSE
        ! Use linear interpolation of xmu*fermi between z1 and z2
        aa = 0
        IF (abs(z1-e1).gt.eps4 .and. abs(z2-e1).gt.eps4) THEN
          aa = - log((z2-e1)/(z1-e1)) *(f1*(z2-e1)+f2*(e1-z1))
          !! z1 or z2 equal to e1; in this case contribution to corr is exactly zero
        ENDIF
        !! second pole
        aa = aa + log((z2-e2)/(z1-e2)) *(f1*(z2-e2)+f2*(e2-z1))
        result = result - aa/ (z2-z1) /2/pi/coni
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE fermilorenH

  SUBROUTINE fermilorenTH(omega, tk, efermi, xloss, height, result)
    ! Perform convolution between fermi function
    ! and lorentzian.
    USE constants
    IMPLICIT none
    ! Input:
    !   omega  - energy in Hatree
    !   tk     - electronic temperature in Hatree
    !   efermi - chemical potential in Hatree
    !   xloss  - broadening in Hatree
    REAL(8), INTENT(IN) :: tk, efermi, xloss, omega, height

    ! Output:
    !   result
    COMPLEX*16, INTENT(OUT) :: result

    ! Local varaibles:
    !   z1        - integration variable
    !   z2        - another name for omega
    !   dummy     - temporary variable
    !   Nanalytic - integration grid size
    INTEGER :: ic, Nanalytic
    REAL(8) :: domega, omegaMax, omegaMin,dummy
    COMPLEX*16 :: z1, z2
    !REAL(8), PARAMETER :: pi = 3.141592653589793d0

    omegaMax = 1.0
    omegaMin = 0.d0
    Nanalytic = int(20000)
    domega = (omegaMax-omegaMin)/DBLE(Nanalytic-1)

    result = 0.d0
    DO ic = 1, Nanalytic-1
      ! From -Inf(Avoid) to -1 [Change of variable 1/x]
      z1 = 1/DBLE(-omegaMax+DBLE(ic-1)*domega+coni*height)
      z2 = omega
      IF (DBLE((z1-efermi)/tk).lt.100) THEN
        dummy = 1/(EXP((z1-efermi)/tk)+1)
        result = result+dummy*(xloss/pi)/((z1-z2)**2+xloss**2)*(z1**2)
      ENDIF
    ENDDO

    DO ic = 1, Nanalytic-1
      ! From -1 to 0
      z1 = -omegaMax+DBLE(ic-1)*domega+coni*height
      z2 = omega
      dummy = 1/(EXP((z1-efermi)/tk)+1)
      result = result+dummy*(xloss/pi)/((z1-z2)**2 + xloss**2)
      ! From 0 to 1
      z1 = omegaMin+DBLE(ic-1)*domega+coni*height
      z2 = omega
      dummy = 1/(EXP((z1-efermi)/tk)+1)
      result = result+dummy*(xloss/pi)/((z1-z2)**2 + xloss**2)
    ENDDO


    DO ic = Nanalytic, 2, -1
      ! From 1 to Inf(Avoid) [Change of variable 1/x]
      z1 = 1/DBLE(omegaMax-DBLE(Nanalytic-ic)*domega+coni*height)
      z2 = omega
      IF (DBLE((z1-efermi)/tk).lt.100) THEN
        dummy = 1/(EXP((z1-efermi)/tk)+1)
        result = result+dummy*(xloss/pi)/((z1-z2)**2+xloss**2)*(z1**2)
      ENDIF
    ENDDO
    result = result*domega
    RETURN
  END SUBROUTINE fermilorenTH

  SUBROUTINE cchi_eff(emxs, xmu2, ne4, efermi, tk, xmu0, xloss, omega, ispec, result)
    ! Interpolation of contour integral for cchi
    ! We interpolation fermi function near the efermi
    REAL(8), INTENT(IN) :: efermi, tk, xloss, omega
    COMPLEX*16, INTENT(IN) :: emxs(:), xmu2(:), xmu0
    INTEGER, INTENT(IN) :: ne4, ispec

    COMPLEX*16, INTENT(OUT) :: result

    ! Local variables
    COMPLEX*16 ::  efermi_m3, efermi_p3, z1, z2, f1, f2, x1, x2, k1, k2
    COMPLEX*16 ::  de, e1, e2, aa
    INTEGER :: ind_m3, ind_p3, ie, ne4_terp
    REAL(8) :: defermi
    COMPLEX*16, ALLOCATABLE :: xmu2_terp(:), emxs_terp(:)
    INTEGER, PARAMETER :: nxmu = 10000
    REAL(8), PARAMETER :: eps4=1e-4
    REAL(8) :: window_size_terp
    LOGICAL :: BUG=.FALSE., TERP=.TRUE.
    CHARACTER*512 message

    window_size_terp = 10

    IF (efermi.GT.DBLE(emxs(1))) THEN
      window_size_terp = MIN(FLOOR((efermi-DBLE(emxs(1)))/tk), 10)
    ELSE
      PRINT*, "WARNING: efermi less than minimum.", efermi*hart, DBLE(emxs(1)*hart)
    ! ELSE
    !   IF ((tk*hart).GT.0.01) window_size_terp = 3.0
    ENDIF

    defermi = (2*window_size_terp*tk)/DBLE(nxmu)
    efermi_m3 = efermi - window_size_terp*tk + coni*DIMAG(emxs(1))
    efermi_p3 = efermi + window_size_terp*tk + coni*DIMAG(emxs(1))
    IF(DBLE(efermi_m3).LT.DBLE(emxs(1))) THEN
       call wlog("WARNING: xmu_m3 less than first point in energy grid. Check results.")
    ELSEIF(DBLE(efermi_p3).GT.DBLE(emxs(ne4))) THEN
       call wlog("WARNING: xmu_p3 greater than last point in energy grid. Check results.")
    END IF
    ! Find where does xmu +/- window_terp lies on the contour
    ! JK - Adding MAX/MIN to be sure no out of bounds arrays.
    ind_m3 = MAX(binarysearch(DBLE(efermi_m3), DBLE(emxs), ne4)-1,1)
    ind_p3 = MIN(binarysearch(DBLE(efermi_p3), DBLE(emxs), ne4),ne4)

    IF (ind_m3.EQ.0) BUG=.TRUE.
    IF (BUG) THEN
      ! This shouldn't happen any more - JK.
      ! Need to be careful regarding the window_size_terp
      ! efermi_m3 must be greater than emxs(1)
      ! TODO: implement optimum window_size_terp
      write(message,'(a,f10.5)') "ERROR: xmu_m3 less than minimum. ", &
          & DBLE(efermi_m3*hart), DBLE(emxs(1)*hart)
      call wlog(message)
      STOP 
    ENDIF

    ! OPEN(UNIT=12312, FILE='petak.dat', STATUS='REPLACE')
    ! IF (DBLE(efermi_m3).LT.DBLE(emxs(1))) THEN
    !   PRINT*, "YES"
    !   TERP = .FALSE.
    !   ALLOCATE(emxs_terp(nxmu))
    !   DO ie = 1, nxmu
    !     emxs_terp(ie) = efermi_m3 + DBLE(ie-1)*defermi
    !   ENDDO
    !   ind_m3 = binarysearch(DBLE(emxs(1)), DBLE(emxs_terp), nxmu)-1
    !   DEALLOCATE(emxs_terp)
    !
    !   ne4_terp = ind_m3 + nxmu + (ne4-ind_p3+1)
    !   ALLOCATE(xmu2_terp(ne4_terp))
    !   ALLOCATE(emxs_terp(ne4_terp))
    !
    !   emxs_terp(1) = emxs(1)
    !   xmu2_terp(1) = xmu2(1)
    !   WRITE(12312,'(20E20.10E3)') DBLE(emxs_terp(1)), DIMAG(emxs_terp(1)), &
    !                  & DBLE(xmu2_terp(1)), DIMAG(xmu2_terp(1))
    !   DO ie = ind_m3, nxmu
    !     ! Copy from the point where its more than emxs(1) to xmu+a*kT
    !     emxs_terp(ie-ind_m3+2) = efermi_m3 + DBLE(ie-1)*defermi
    !     xmu2_terp(ie-ind_m3+2) = interp1d(DBLE(emxs_terp(ie-ind_m3+2)),  &
    !                                      & DBLE(emxs), xmu2)
    !     WRITE(12312,'(20E20.10E3)') DBLE(emxs_terp(ie-ind_m3+2)), DIMAG(emxs_terp(ie-ind_m3+2)),&
    !                    & DBLE(xmu2_terp(ie-ind_m3+2)), DIMAG(xmu2_terp(ie-ind_m3+2))
    !   ENDDO
    !
    !   ! Insert xmu_p3 into the thing
    !   emxs_terp((nxmu-ind_m3)+2) = efermi_p3
    !   xmu2_terp((nxmu-ind_m3)+2) = interp1d(DBLE(efermi_p3),  DBLE(emxs), xmu2)
    !   WRITE(12312,'(20E20.10E3)') DBLE(emxs_terp((nxmu-ind_m3)+2)), DIMAG(emxs_terp((nxmu-ind_m3)+2)), &
    !                 & DBLE(xmu2_terp((nxmu-ind_m3)+2)), DIMAG(xmu2_terp((nxmu-ind_m3)+2))
    !   DO ie=1, (ne4-ind_p3)
    !     !  Going back to original grid at ind_p3
    !     emxs_terp((nxmu-ind_m3)+2+ie) = emxs(ie+ind_p3)
    !     xmu2_terp((nxmu-ind_m3)+2+ie) = xmu2(ie+ind_p3)
    !     WRITE(12312,'(20E20.10E3)') DBLE(emxs_terp((nxmu-ind_m3)+2)+ie), DIMAG(emxs_terp((nxmu-ind_m3)+2)+ie), &
    !                   & DBLE(xmu2_terp((nxmu-ind_m3)+2)+ie), DIMAG(xmu2_terp((nxmu-ind_m3)+2)+ie)
    !   ENDDO
    !   STOP
    ! ENDIF
    ! CLOSE(12312)
    IF (TERP) THEN
      ne4_terp = ind_m3 + nxmu + (ne4-ind_p3+1)
      ALLOCATE(xmu2_terp(ne4_terp))
      ALLOCATE(emxs_terp(ne4_terp))
      emxs_terp = 0.d0
      xmu2_terp =  0.d0

      DO ie=1, ind_m3
        ! Copy all points less than xmu-3kT
        emxs_terp(ie) = emxs(ie)
        xmu2_terp(ie) = xmu2(ie)
      ENDDO

      DO ie=ind_m3+1, ind_m3+nxmu
        ! Remove xmu +/- 3kT from original grid and replaced by interpolated points
        emxs_terp(ie) = efermi_m3 + (ie-ind_m3-1) * defermi
        xmu2_terp(ie) = interp1d(DBLE(emxs_terp(ie)),  DBLE(emxs), xmu2)
      ENDDO

      ! Insert xmu_p3 into the thing
      emxs_terp(ind_m3+nxmu+1) = efermi_p3
      xmu2_terp(ind_m3+nxmu+1) = interp1d(DBLE(efermi_p3),  DBLE(emxs), xmu2)

      DO ie=1, (ne4-ind_p3)
        !  Going back to original grid at ind_p3
        emxs_terp(ind_m3+nxmu+1+ie) = emxs(ie+ind_p3)
        xmu2_terp(ind_m3+nxmu+1+ie) = xmu2(ie+ind_p3)
      ENDDO
    ! ELSE
    !   ! Do not do interpolation
    !   ne4_terp = ne4
    !   ALLOCATE(xmu2_terp(ne4_terp))
    !   ALLOCATE(emxs_terp(ne4_terp))
    !   emxs_terp = emxs
    !   xmu2_terp = xmu2
    ENDIF

    ! OPEN(UNIT=154, file='lorentzian.dat', status='OLD')
    result = 0.d0
    e1 = omega + coni*xloss
    e2 = omega - coni*xloss
    DO ie = 1, ne4_terp-1
      z1 = emxs_terp(ie)
      k1 = xmu2_terp(ie)-xmu0
      x1 = DBLE(z1)
      f1 = fermiDirac(x1,tk,efermi)
      IF (ispec.ne.2) f1 = 1.d0-f1
      z2 = emxs_terp(ie+1)
      k2 = xmu2_terp(ie+1)-xmu0
      x2 = DBLE(z2)
      f2 = fermiDirac(x2,tk,efermi)
      IF (ispec.ne.2) f2 = 1.d0-f2
      de = z2-z1
      ! result = result + k1*f1*cauchy(z1, omega, xloss)*de
      IF ((DBLE(z1).GT.DBLE(efermi_m3)).AND.(DBLE(z1).LT.DBLE(efermi_p3))) THEN
        !! Trapezoidal rule
        result = result + 0.5d0*de*(f1*k1*cauchy(z1, omega, xloss) + &
                                  & f2*k2*cauchy(z2,omega,xloss))
      ELSE
        ! Use linear interpolation of xmu*fermi between z1 and z2
        aa = 0
        IF (abs(z1-e1).gt.eps4 .and. abs(z2-e1).gt.eps4) THEN
          aa = - log((z2-e1)/(z1-e1)) *(k1*f1*(z2-e1)+k2*f2*(e1-z1))
          !! z1 or z2 equal to e1; in this case contribution to corr is exactly zero
        ENDIF
        !! second pole
        aa = aa + log((z2-e2)/(z1-e2)) *(k1*f1*(z2-e2)+k2*f2*(e2-z1))
        result = result - aa/ (z2-z1) /2/pi/coni
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE cchi_eff

  COMPLEX*16 FUNCTION cchi_fine(omegaMin, omegaMax, n, omega, xloss, tk, efermi, &
                      & emxs, xmu2m0, ne4, ispec)
    IMPLICIT none
    INTEGER, INTENT(IN) :: n, ne4, ispec
    REAL(8), INTENT(IN) :: omegaMin, omegaMax, omega, xloss, tk, efermi
    COMPLEX*16, INTENT(IN) :: emxs(:), xmu2m0(:)

    ! Local variables
    COMPLEX*16 :: zp, zm, f1, f2, result, a1, a2, z1, corr2, zb
    REAL(8) :: domega, z2, dummy, dummy2
    INTEGER :: ic

    corr2 = 0.d0
    OPEN(UNIT=9991, FILE='Source.dat', STATUS='REPLACE')
    DO ic = 1, ne4-1
      z1 = emxs(ic)
      f1 = xmu2m0(ic)
      zb = DBLE(z1)
      f2 = fermiDirac(zb,tk,efermi)
      IF (ispec.ne.2) f2 = 1.d0-f2
      corr2 = corr2 + f1*f2*cauchy(z1,omega,xloss)
      WRITE(9991,'(20E20.10E3)') DBLE(z1), DIMAG(z1), DBLE(f1), DBLE(f2), DBLE(f1*f2*cauchy(z1,omega,xloss)), DBLE(emxs(ic+1)-emxs(ic))
    ENDDO
    corr2 = corr2*DBLE(emxs(2)-emxs(1))
    CLOSE(9991)

    OPEN(UNIT=9991, FILE='Interp.dat', STATUS='REPLACE')
    domega = (omegaMax-omegaMin)/DBLE(n-1)
    result = 0.d0
    DO ic = 1, n-1
      zm = omegaMin+DBLE(ic-1)*domega + coni*DIMAG(emxs(1))
      zp = omegaMin+DBLE(ic)*domega + coni*DIMAG(emxs(1))
      z2 = omega
      ! dummy = 1/(EXP((zm-efermi)/tk)+1)
      ! dummy2 = 1/(EXP((zp-efermi)/tk)+1)
      a1 = DBLE(zm)
      a2 = DBLE(zp)
      dummy = fermiDirac(a1,tk,efermi)
      dummy2 = fermiDirac(a2,tk,efermi)
      f1 = interp1d(DBLE(zm),DBLE(emxs),xmu2m0)
      f2 = interp1d(DBLE(zp),DBLE(emxs),xmu2m0)
      IF (ispec.ne.2) dummy = 1.d0-dummy
      IF (ispec.ne.2) dummy2 = 1.d0-dummy2
      ! result = result + 0.5d0*(f1*dummy*cauchy(zm,z2,xloss) + f2*dummy2*cauchy(zp,z2,xloss))

      WRITE(9991,'(20E20.10E3)') DBLE(zm), DIMAG(zm), DBLE(f1), DBLE(dummy), DBLE(f1*dummy*cauchy(zm,z2,xloss)), DBLE(domega)
    ENDDO
    CLOSE(9991)
    cchi_fine = result*domega
    RETURN
  END FUNCTION

  SUBROUTINE sommerfeld_xscorr(ispec, emxs, ne1, ne, ik0, xsec, xsnorm, chia, &
    &                       vrcorr, vicorr, cchi, electronic_temperature)
    IMPLICIT none

    INTEGER, INTENT(IN) :: ispec, ne1, ne, ik0
    REAL(8), INTENT(IN) :: electronic_temperature, vrcorr, xsnorm(nex), vicorr
    COMPLEX*16, INTENT(IN) :: xsec(nex), chia(nex)
    COMPLEX*16, INTENT(INOUT) :: emxs(nex)
    COMPLEX*16, INTENT(OUT) :: cchi(nex)

    INTEGER :: polex, idx, ic0, ip, nc
    INTEGER :: ne2, ne4, npole
    REAL(8) :: tk, xloss, efermi
    REAL(8) :: omega(nex)
    COMPLEX*16 :: corr, corr1, corr2, residue, leg3, correction
    COMPLEX*16 :: ec(nex), fc(nex), xmu(nex), pole
    COMPLEX*16 :: xmu1(nex), xmu2(nex), xmu3(nex), xmu4(nex)
    REAL(8), PARAMETER :: eps4 = 1.0d-4
    LOGICAL, PARAMETER :: print_out = .TRUE.
    EXTERNAL :: lorenz, astep
    COMPLEX*16 :: lorenz, astep

    INTEGER :: ie, ic
    REAL(8) :: dele, w1, w2 ,w3
    REAL(8) :: ra, rb, cheight
    COMPLEX*16 :: xmu0, ff(nex)
    COMPLEX*16 :: bb, c1, f1, f2, e1, e2, z1, z2
    COMPLEX*16 :: dummy, dummy2, aa, fd1, fd2

    REAL(8) :: minee, maxee, dee
    COMPLEX*16 :: ee, f4
    INTEGER :: extra
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! tk = xloss/pi ! Critical temp in Hatree
    ! tk = 1e-6/hart
    tk = electronic_temperature/hart
    ne2 = 10
    extra = 1
    npole = ne-2*ne1-ne2-extra
    ne4 = ne1

    efermi = DBLE(emxs(ne-extra))
    xloss = DIMAG(emxs(ne4+1))
    cheight = DIMAG(emxs(1))

    ! IF (ABS(cheight-xloss).LT.eps4) STOP "ERROR: contour height <= xloss"

    DO ie = 1, ne
      xmu(ie) = xsec(ie) + xsnorm(ie)*chia(ie)
    ENDDO

    OPEN(UNIT=155, file='curve.dat', status='REPLACE')
    OPEN(UNIT=150, file="prexmu.dat", status='REPLACE')
    OPEN(UNIT=149, file="raw.dat", status='REPLACE')
    OPEN(UNIT=153, file="leg1.dat", status='REPLACE')
    OPEN(UNIT=151, file='residue.dat', status='REPLACE')
    OPEN(UNIT=152, file='contour.dat', status='REPLACE')
    ! OPEN(UNIT=154, file='lorentzian.dat', status='REPLACE')
    WRITE(149,*) "Temperature (Hatree)", tk
    WRITE(149,*) "Electronic Temperature (eV)", electronic_temperature
    WRITE(149,*) "xloss = ", xloss, " Hatree"
    WRITE(149,*) "Chemical potential = ", efermi*hart, " eV"
    WRITE(149,*) "Number of poles = ", npole
    WRITE(149, *) 'Omega(Hart)  \t  Re CCHI  \t   Im CCHI \t  1-Fermi \t  Re xmu0  \t  Im xmu0'

    DO ie =1, ne
      WRITE(155, '(20E20.10E3)') ie, DBLE(emxs(ie)), DIMAG(emxs(ie))
    ENDDO
    CLOSE(155)
    !---------- Construct the integration contour
    nc = 0
    ! DO ie = 1, ne2
    !     nc = nc + 1
    !     ec(nc) = emxs(2*ne1+ie)
    !     fc(nc) = xmu(2*ne1+ie)
    ! ENDDO

    DO ie=1, ne4
      IF (ispec.NE.2) THEN
        IF ((DBLE(emxs(ie))-efermi).GT.eps4) THEN
          nc = nc + 1
          ec(nc) = emxs(ie)
          fc(nc) = xmu(ie)
        ENDIF
      ELSE
        IF (efermi-DBLE(emxs(ie)).GT.eps4) THEN
          nc = nc + 1
          ec(nc) = emxs(ie)
          fc(nc) = xmu(ie)
        ENDIF
      ENDIF
    ENDDO
    !--------- End of countour construction

    DO ie = 1, ne4  !---------- Real frequencies from ecv to edge + 5
      omega(ie) = DBLE(emxs(ie))
      xmu2(ie) = xmu(ie)     !--- xmu(omega+i*height)
      xmu3(ie) = xmu(ne4+ie) !--- xmu(omega+i*xloss)
    ENDDO           !---------- End of Real frequencies from ecv to edge + 5

    DO ie = 1, npole
      xmu4(ie) = xmu(2*ne1+ne2+ie) !--- xmu(efermi+i*En)
    ENDDO

    DO ie = 1, ne4   !--------- Perform convolution
      bb = 1

      xmu0 = xmu3(ie)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Analytic Part (no pole omega+i*xloss)!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cchi(ie) = 0.d0
      CALL fermilorenT(omega(ie), tk, efermi, xloss, cchi(ie))
      ! CALL fermilorenT_terp(omega(ie), tk, efermi, xloss, cchi(ie))

      dummy = cchi(ie)
      cchi(ie) = xmu0*cchi(ie)
      IF (ispec.ne.2) cchi(ie) = xmu0 - cchi(ie) ! Absorption
      dummy = 1-dummy
      WRITE(149, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie)), DBLE(dummy), DBLE(xmu0), DIMAG(xmu0)

      !!!!!!!!!!!!!!!!!!!!!!
      !! RESIDUE of Fermi !!
      !!!!!!!!!!!!!!!!!!!!!!
      ! At zero temperature this will become last leg
      residue = 0.d0
      dele = omega(ie) - efermi
      DO ic = 1, npole
        z1 = emxs(2*ne1+ne2+ic) ! pole
        f1 = xmu4(ic)-xmu0
          residue = residue + f1 * xloss/((z1-omega(ie))**2 + xloss**2)
      ENDDO
      residue = -residue*(2.d0*coni*tk)

      IF (ispec.ne.2) residue = - residue ! Absorption
      WRITE(151, '(20E20.10E3)') DBLE(omega(ie)), DBLE(residue), DIMAG(residue), DBLE(efermi)

      !!!!!!!!!!!!!!!!
      !! Correction !!
      !!!!!!!!!!!!!!!!

      corr1 = 0.d0
      dele = efermi-omega(ie)
      aa = interp1d(efermi, DBLE(emxs(1:ne4)), xmu2)
      dummy = interp1d(efermi*(1.d0-1e-3), DBLE(emxs(1:ne4)), xmu2)  ! temporary
      dummy2 = interp1d(efermi*(1.d0+1e-3), DBLE(emxs(1:ne4)), xmu2) ! temporary
      fd2 = (dummy2-2.d0*aa+dummy)/DBLE(2d-5)**2
      fd1 = (dummy2-dummy)/DBLE(2e-5)

      ! from pole above real axis
      dummy = (fd1*(dele+coni*(cheight-xloss))+xmu0-aa)/(dele+coni*(cheight-xloss))**2
      ! from pole below real axis
      dummy2 = (fd1*(dele+coni*(cheight+xloss))+xmu0-aa)/(dele+coni*(cheight+xloss))**2

      corr1 = (pi*tk**2)*(dummy-dummy2)/(12*coni*xloss)
      corr1 = -corr1*pi*tk**2/(24.d0*xloss)
      IF (ispec.eq.2) corr1 = -corr1


      WRITE(153, '(20E20.10E3)') DBLE(omega(ie)), DBLE(corr1), DIMAG(corr1)

      !!!!!!!!!!!!!!!!!
      !! Second  leg !!
      !!!!!!!!!!!!!!!!!
      corr2 = 0.d0
      DO ic = 1, nc-1
        z1 = ec(ic)
        z2 = ec(ic+1)
        f1 = fc(ic)-xmu0
        f2 = fc(ic+1)-xmu0
        aa = 0
        if (abs(z1-e1).gt.eps4 .and. abs(z2-e1).gt.eps4) then
          aa = log((z2-e1)/(z1-e1)) *(f1*(z2-e1)+f2*(e1-z1))
          ! z1 or z2 equal to e1; in this case contribution to corr is exactly zero
        endif
        ! second pole
        aa = aa - log((z2-e2)/(z1-e2)) *(f1*(z2-e2)+f2*(e2-z1))
        corr2 = corr2 + aa/ (z2-z1) /2/pi/coni
      ENDDO
      IF (ispec.EQ.2) corr2 = -corr2

      WRITE(152, '(20E20.10E3)') DBLE(omega(ie)), DBLE(corr2), DIMAG(corr2)

      !!!!!!!!!!!!!!!!!!!!!
      !! Combine results !!
      !!!!!!!!!!!!!!!!!!!!!
      corr = corr1 + corr2
      cchi(ie) = cchi(ie) + corr + residue
      WRITE(150, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie))

      ! Return the result of convolution minus bare value
      cchi(ie) = cchi(ie) - xmu2(ie)
    ENDDO !------- End of convolution

    CLOSE(149)
    CLOSE(150)
    CLOSE(151)
    CLOSE(152)
    CLOSE(153)
    ! CLOSE(154)

    RETURN
  END SUBROUTINE


  SUBROUTINE sommerfeld_xscorr0(ispec, emxs, ne1, ne, ik0, xsec, xsnorm, chia, &
    &                       vrcorr, vicorr, cchi, electronic_temperature)

    IMPLICIT none

    INTEGER, INTENT(IN) :: ispec, ne1, ne, ik0
    REAL(8), INTENT(IN) :: electronic_temperature, vrcorr, xsnorm(nex), vicorr
    COMPLEX*16, INTENT(IN) :: xsec(nex), chia(nex)
    COMPLEX*16, INTENT(INOUT) :: emxs(nex)
    COMPLEX*16, INTENT(OUT) :: cchi(nex)

    INTEGER :: polex, idx, ic0, ip, nc
    INTEGER :: ne2, ne4, npole
    REAL(8) :: tk, xloss, efermi, dele
    REAL(8) :: omega(nex)
    COMPLEX*16 :: corr, corr1, corr2, correction, leg3, residue
    COMPLEX*16 :: ec(nex), fc(nex), xmu(nex), pole
    REAL(8), PARAMETER :: eps4 = 1.0d-4
    EXTERNAL :: lorenz, astep
    REAL(8) :: astep
    COMPLEX*16 :: lorenz

    INTEGER :: ie, ic
    REAL(8) :: x1, x2, cheight
    DOUBLE PRECISION :: w1
    COMPLEX*16 :: xmu0, ff(nex)
    COMPLEX*16 :: bb, c1, f1, f2, e1, e2, z1, z2
    COMPLEX*16 :: dummy, dummy2, aa

    tk = electronic_temperature/hart
    ne2 = ne-ne1-3
    npole = 0
    ne4 = 0

    efermi = DBLE(emxs(ne-3))
    xloss = DIMAG(emxs(1))
    cheight = DIMAG(emxs(1))

    ! Only happens at low temperature
    IF (MOD(xloss, tk*pi).LT.eps4) PRINT*, "WARNING: xloss close to pole"

    OPEN(UNIT=150, file="prexmu.dat", status='REPLACE')
    OPEN(UNIT=149, file="raw.dat", status='REPLACE')
    OPEN(UNIT=153, file="leg1.dat", status='REPLACE')
    OPEN(UNIT=151, file='residue.dat', status='REPLACE')
    OPEN(UNIT=152, file='contour.dat', status='REPLACE')
    OPEN(UNIT=155, file='curve.dat', status='REPLACE')
    OPEN(UNIT=156, file='correction.dat', status='REPLACE')
    ! OPEN(UNIT=154, file='lorentzian.dat', status='REPLACE')
    WRITE(149,*) "Temperature (Hatree)", tk
    WRITE(149,*) "Electronic Temperature (eV)", electronic_temperature
    WRITE(149,*) "xloss = ", xloss, " Hatree"
    WRITE(149,*) "Chemical potential = ", efermi*hart, " eV"
    WRITE(149,*) "Number of poles = ", npole
    WRITE(149, *) 'Omega(Hart)  \t  Re CCHI  \t   Im CCHI \t  1-Fermi \t  Re xmu0  \t  Im xmu0'

    DO ie = 1, ne
      xmu(ie) = xsec(ie) + xsnorm(ie)*chia(ie)
    ENDDO

    DO ie = 1, ne1  !---------- Real frequencies
      omega(ie) = DBLE(emxs(ie))
    ENDDO           !---------- End of Real frequencies

    if (abs(vrcorr).gt.eps4) then
      !       account for the fermi level shift
      bb = xmu(ik0)
      efermi = efermi - vrcorr
      call terpc(omega, xmu ,ne1, 1, efermi, bb)

      !       shift the vertical axis
      do ie = 1, ne2+1
        emxs(ne1+ie) = emxs(ne1+ie) - vrcorr
      enddo

      !       rescale values on vertical axis
      bb = bb/xmu(ik0)
      do ie = ne1+1, ne
        xmu(ie) = xmu (ie) * bb
      enddo
    else
      bb = 1
    endif
    !---------- Construct the integration contour
    nc = 0
    DO ie = 1, ne2
      IF (DIMAG(emxs(ne1+ie)).LT.xloss) THEN
        nc = nc + 1
        ec(nc) = emxs(ne1+ie)
        fc(nc) = xmu(ne1+ie)
        WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
      ENDIF
    ENDDO

    ! Add corner at efermi + xloss*i
    nc = nc+1
    ic0 = nc
    if (abs(vrcorr).gt.eps4) then
      ec(nc) = efermi + coni*xloss
      fc(nc) = bb * xmu(ik0)
      WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
    else
      ec(nc) = emxs(ik0)
      fc(nc) = xmu(ik0)
      WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
    endif

    IF (ispec.NE.2) THEN
      DO ie =1,ne1
        IF (DBLE(emxs(ie))-efermi.GT.eps4) THEN
          nc = nc+1
          ec(nc) = emxs(ie)
          fc(nc) = xmu(ie)
          WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
        ENDIF
      ENDDO
    ELSE
      DO ie =ne1,1,-1
        IF (efermi-DBLE(emxs(ie)).GT.eps4) THEN
          nc = nc+1
          ec(nc) = emxs(ie)
          fc(nc) = xmu(ie)
          WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
        ENDIF
      ENDDO
    ENDIF
    !--------- End of countour construction


    !--------- Perform convolution
    DO ie = 1, ne1
      bb = 1
      IF (omega(ie).ge.efermi) THEN
        xmu0 = xmu(ie)
        IF (ispec.eq.2) xmu0 = xmu(ik0)*bb
      ELSE
        xmu0 = xmu(ik0)*bb
        IF (ispec.eq.2) xmu0 = xmu(ie)
      ENDIF

      e1 = omega(ie) + coni*xloss
      e2 = omega(ie) - coni*xloss
      do ic = 1, nc
        ff(ic) = fc(ic) - xmu0
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Analytic Part (no pole omega+i*xloss) !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cchi(ie) = (1-astep(xloss, omega(ie)-efermi))
      dummy = cchi(ie)
      cchi(ie) = xmu0*cchi(ie)
      IF (ispec.ne.2) cchi(ie) = xmu0 - cchi(ie) ! Absorption
      dummy = 1-dummy
      WRITE(149, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie)), DBLE(dummy), DBLE(xmu0), DIMAG(xmu0)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! sommerfeld correction !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      x1 = 0
      x2 = 0
      dele = omega(ie)-efermi
      if (abs(dele).lt.eps4) dele = 0.d0
      w1 = DIMAG(emxs(ne-1))

      ! ! Top pole
      ! dummy = (xmu(ne)-2.d0*xmu(ne-1)+xmu(ne-2))/(DIMAG(emxs(ne)-emxs(ne-2))**2)
      ! dummy = -pi*tk**2 * dummy/(24.d0*coni*xloss)
      !
      ! ! Bottom pole
      ! dummy2 = (xmu(ne)-xmu(ne-2))/DIMAG(emxs(ne)-emxs(ne-2))
      ! dummy2 = - pi*(tk)**2 * dummy2/(24.d0*xloss**2)
      ! correction = dummy - dummy2
      ! IF (ispec.NE.2) correction = -correction
      ! WRITE(156, '(20E20.10E3)') DBLE(omega(ie)), DBLE(correction), DIMAG(correction), DBLE(lorenz0(xloss,w1,dele)), DIMAG(lorenz0(xloss,w1,dele))



      !!!!!!!!!!!!!!
      !!  Contour !!
      !!!!!!!!!!!!!!
      corr = 0.d0
      ! add half matsubara pole contribution
      ! equivalent to integral from efermi to efermi+i*w1
      dele = omega(ie)-efermi
      if (abs(dele).lt.eps4) dele = 0.0d0
      w1 = DIMAG(ec(1))
      corr = corr + lorenz0(xloss,w1,dele)*ff(1)*coni*w1
      leg3 = lorenz0(xloss,w1,dele)
      WRITE(151, '(20E20.10E3)') DBLE(omega(ie)), DBLE(leg3), DIMAG(leg3)

      ! Horizontal  Contour
      DO ic = 1, nc-1
        z1 = ec(ic)
        z2 = ec(ic+1)
        f1 = ff(ic)
        f2 = ff(ic+1)
        aa = 0
        IF (abs(z1-e1).gt.eps4 .and. abs(z2-e1).gt.eps4) THEN
          aa = log((z2-e1)/(z1-e1)) *(f1*(z2-e1)+f2*(e1-z1))
          ! z1 or z2 equal to e1; in this case contribution to corr is exactly zero
        ENDIF
        ! second pole
        aa = aa - log((z2-e2)/(z1-e2)) *(f1*(z2-e2)+f2*(e2-z1))
        corr = corr + aa/ (z2-z1) /2/pi/coni
      ENDDO
      IF (ispec.eq.2) corr = -corr

      cchi(ie) = cchi(ie) + corr + correction !+  correction + residue
      WRITE(152, '(20E20.10E3)') DBLE(omega(ie)), DBLE(corr), DIMAG(corr)
      WRITE(150, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie))




      ! Return the result of convolution minus bare value
      cchi(ie) = cchi(ie) - xmu(ie)
    ENDDO !------- End of convolution

    CLOSE(149)
    CLOSE(150)
    CLOSE(151)
    CLOSE(152)
    CLOSE(153)
    CLOSE(155)
    CLOSE(156)
    ! end of cycle over frequency points

    if (abs(vrcorr).gt.eps4) then
      do  ie = ne1+1, ne
        emxs(ie) = emxs(ie) + vrcorr
      enddo
    endif

    RETURN
  END SUBROUTINE

  REAL(8)FUNCTION reg_integration(omegaMin, omegaMax, n, omega, xloss, tk, efermi)
    ! Perform trapezoidal integration of fermi-loren function
    ! with finite end-points
    IMPLICIT none
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: omegaMin, omegaMax, omega, xloss, tk, efermi

    ! Local variables
    COMPLEX*16 :: zp, zm
    REAL(8) :: domega, z2, dummy, dummy2, result
    INTEGER :: ic

    domega = (omegaMax-omegaMin)/DBLE(n-1)
    result = 0.d0
    DO ic = 1, n-1
      zm = omegaMin+DBLE(ic-1)*domega
      zp = omegaMin+DBLE(ic)*domega
      z2 = omega
      dummy = fermiDirac(zm, tk, efermi)
      dummy2 = fermiDirac(zp, tk, efermi)
      result = result + 0.5d0*(dummy*cauchy(zm,z2,xloss) + dummy2*cauchy(zp,z2,xloss))
    ENDDO
    reg_integration = result*domega
    RETURN
  END FUNCTION

  REAL(8) FUNCTION reciprocal_integration(omegaMax, n, sign, omega, xloss, tk, efermi)
    ! Return the integral using 1/x substitution method with inf endpoint
    ! for fermi-loren function

    ! Input variables
    !   omega: Reciprocal of interval
    !   sign : +1 for  omegaMax to Inf
    !         -1 for -Inf to omegaMax
    !      n : Number of points
    IMPLICIT none
    REAL(8), INTENT(IN) :: omegaMax, omega, xloss, tk, efermi
    INTEGER, INTENT(IN) :: n, sign

    ! Local variable
    COMPLEX*16 :: zm, zp
    REAL(8) :: z2, dz, dummy, dummy2, result, domega
    INTEGER :: ic

    ! Assertion
    IF (ABS(sign).NE.1) STOP "Sign must be +/- 1"
    domega = ABS(omegaMax)/DBLE(n-1)
    result = 0.d0

    IF (sign.EQ.1) THEN
      DO ic = n, 2, -1
        ! From 1 to Inf(Avoid) [Change of variable 1/x]
        zm = 1/DBLE(ABS(omegaMax)-DBLE(n-ic)*domega)
        z2 = omega
        IF (DBLE((zm-efermi)/tk).lt.100) THEN
          ! dummy = 1/(EXP((zm-efermi)/tk)+1)
          dummy = fermiDirac(zm, tk, efermi)
          result = result + dummy*cauchy(zm, z2, xloss)*(zm**2)
        ENDIF
      ENDDO
      result = result
    ELSE
      DO ic = 1, n-1
        ! From -Inf(Avoid) to -1 [Change of variable 1/x]
        zm = 1/DBLE(-ABS(omegaMax)+DBLE(ic-1)*domega)
        z2 = omega
        ! IF (DBLE((z1-efermi)/tk).lt.100) THEN
          ! dummy = 1.d0/(EXP((zm-efermi)/tk)+1.d0)
          dummy = fermiDirac(zm, tk, efermi)
          result = result + dummy * cauchy(zm, z2, xloss) * (zm**2)
        ! ENDIF
      ENDDO
    ENDIF

    reciprocal_integration = result*domega
    RETURN
  END FUNCTION reciprocal_integration

  ! Binary Insert: Return index i such that a[i-1] < x0 <= a[i]
  INTEGER FUNCTION binarysearch(x0, x, n)
    IMPLICIT none
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x0, x(n)

    INTEGER :: L, R, M, P

    IF (n.LE.0) goto 1234

    L = 0
    R = n-1

    ! Put End
    IF (x0.GE.x(n)) THEN
      binarysearch = n+1
      GOTO 222
    ENDIF

    ! Put Front
    IF (x0.LE.x(1)) THEN
      binarysearch = 1
      GOTO 222
    ENDIF

    111 M = (R+L)/2

    IF (M.LE.0) GOTO 1234

    IF(M.EQ.L) THEN
      binarysearch = M + 1
      GOTO 222
    ENDIF

    IF (x0.LT.x(M)) THEN
      GOTO 123 !LEFT
    ELSE IF (x0.GT.x(M)) THEN
      GOTO 321 !RIGHT
    ELSE
      binarysearch = M ! FOUND
      GOTO 222
    ENDIF

    ! LEFT
    123 IF (x0.EQ.x(M)) THEN
          binarysearch = M
          GOTO 222
        ELSE
          R = M
          GOTO 111 ! Divide again
        ENDIF


    ! RIGHT
    321 IF (x0.EQ.x(R)) THEN
          binarysearch = R
          GOTO 222
        ELSE
          L = M
          GOTO 111 ! Divide again
        ENDIF


    1234 binarysearch = -1 ! Not found or too short
    PRINT*, "ERROR: Binary search failed"
    222 return
  END FUNCTION binarysearch

  ! Linear interpolation
  COMPLEX*16 FUNCTION interp1d(x0, x, y)
    IMPLICIT none
    INTEGER :: n, idx
    REAL(8), INTENT(IN) :: x0, x(:)
    COMPLEX*16, INTENT(IN) :: y(:)
    COMPLEX*16 :: y1, y2
    REAL(8) :: x1, x2
    n = size(x)
    if ((x(1).LT.x0).and.(x0.LT.x(n))) then
      idx = binarysearch(x0,x,n)
      y2 = y(idx)
      x2 = x(idx)
      y1 = y(idx-1)
      x1 = x(idx-1)
      interp1d = (y2-y1)*(x0-x1)/(x2-x1) + y1
    else
      IF (x(1).EQ.x0) THEN
        interp1d = y(1)
      ELSE IF (ABS(x(n)-x0).LT.1e-10) THEN
        interp1d = y(n)
      ELSE
        PRINT*, "ERROR: Out of range !",x(1)*hart,"<=",x0*hart,"<=",x(n)*hart
        STOP
      ENDIF
    endif
    RETURN
  END FUNCTION interp1d

  COMPLEX*16 FUNCTION fermiDirac(z, tk, efermi)
    IMPLICIT NONE
    COMPLEX*16, INTENT(IN) :: z
    REAL(8) :: tk, efermi
    fermiDirac = 0.d0
    IF (tk.GT.0.d0) THEN
      IF (DBLE((z-efermi)/tk).lt.100) THEN
        fermiDirac = 1.d0/(EXP((z-efermi)/tk)+1.d0)
      ENDIF
    ELSE
      IF (DBLE(z).GT.efermi) fermiDirac=1.d0
    ENDIF
    RETURN
  END FUNCTION

  COMPLEX*16 FUNCTION cauchy(x, w, xloss)
    IMPLICIT none
    COMPLEX*16, INTENT(IN) :: x
    REAL(8), INTENT(IN) :: w, xloss

    cauchy =  (xloss/pi)/((x-w)**2+xloss**2)
    RETURN
  END FUNCTION cauchy

  COMPLEX*16 FUNCTION lorenz0(xloss,w,dele)
    ! use constants
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: xloss, dele, w
    lorenz0 = xloss /pi / (xloss**2+(coni*w-dele)**2)

    return
  END FUNCTION lorenz0

  REAL(8) FUNCTION astep (xloss, dele)
    USE constants
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: xloss, dele
    !  1 - fermi*cauchy
    astep = 0.5d0 + atan(dele/xloss) /pi
    if (astep.lt.0.d0) astep = 0.d0
    if (astep.gt.1.d0) astep = 1.d0

    RETURN
  END FUNCTION astep
END MODULE m_thermal_xscorr
