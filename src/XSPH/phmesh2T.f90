!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finite temperature version of phmesh2
!
! $Revision: 1.0 $
! $Author: tts $
! $Date: 2018/6/1 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE phmesh2T (iprint, ispec, edge, emu, vi0, gamach,        &
  &     ecv, xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3, iGrid)
  !     This subroutine makes the energy mesh used for phases and cross sections,
  !     as well as for the fms routine, path, and genfmt. For EXAFS, the final output
  !     chi is on a different (usually finer) grid with mu0 interpolated.
  !     This will reproduce the old (FEFF84) grids, as well as any combination of user
  !     defined energy, k, exponential, or arbitrary (read from file) grids. The input
  !     for the user defined grids is read from grid.inp. Details of grid.inp are given
  !     in rdgrid.f
  USE constants
  USE DimsMod, ONLY: nex

  IMPLICIT NONE
  ! Input:
  !     iprint -
  !     ispec  - controls which grid to use (0=EXAFS,1=XANES,2=XES,3=DANES,4=FPRIME)
  !     edge   - This name is misleading and is not the x-ray edge energy.
  !              edge = xmu - vr0, where vr0 is given as an the first option in the
  !              EXCHANGE card.
  !     vi0    - Contant imaginary part added to the potential, second option in the
  !              EXCHANGE card.
  !     ecv    - Core-valence separation
  !     gamach - Core-Hole broadening. Gives an additional constant imaginary part to
  !              the potential.
  !     xkmax  - maximum k value for EXAFS/XANES calculations. holds emin for f'
  !              calculations.
  !     xkstep - k-grid spacing for XANES calculations. holds emax for f' calculations
  !     vixan  - energy step for FMS calculations (grid is even in energy near edge)

  INTEGER :: iprint, ispec, iGrid
  DOUBLE PRECISION:: edge, emu, vi0, ecv, gamach, xkmax, xkstep, vixan

  ! Output:
  !     ne     - Total number of energy points.
  !     ne1    - Number of energy points on the horizontal grid.
  !     ik0    - point where k=0
  !     em(ne) - energy array
  INTEGER ne, ne1, ne3, ik0
  COMPLEX*16 em(nex)

  ! Local Variables:
  !     xloss  - total constant imaginary part of em
  !     xim    - energy step near the fermi level
  !     deltak - k step
  !     emin   - minimum e for exponential grid used by DANES
  !     emax   - maximum e for exponential grid used by DANES
  !     del    - step for exponential grid
  !     ios    - i/o errors
  !     nemax  - temp variable to hold max # of energy points
  DOUBLE PRECISION xloss, xim, deltak, emin, emax, del, del2
  INTEGER ios, nemax
  ! User defined grid variables
  !     nGridMax  - max number of grids
  !     nGrid     - number of grids.
  !     iGridType - Type of grid (1 = energy, 2 = k, 3 = exp)
  !     GridMin   - minimum k or E of grid. k for k-grids, e for others
  !     GridMax   - Maximum k or E of grid
  !     GridStep  - step size.
  INTEGER nGridMax
  PARAMETER(nGridMax=10)
  INTEGER nGrid, iGridType(nGridMax)
  DOUBLE PRECISION GridMin(nGridMax), GridMax(nGridMax), GridStep(nGridMax)
  DOUBLE PRECISION estep0, expdel, de

  ! Loop Variables:
  INTEGER :: ii, nv, n1, NPts, ie
  DOUBLE PRECISION :: getxk
  LOGICAL :: normal_mesh
  CHARACTER(512) message
  EXTERNAL :: getxk

  ! Initialization
  normal_mesh = .TRUE.

  ! Set total imaginary part, must be >= 0.02 eV
  xloss = MAX(gamach/2.d0 + vi0, 0.02d0/hart)
  ! xloss = MAX(gamach/2.d0 + vi0, 0.01d0/hart)
  ! Set energy step to half of imaginary part, or vixan if vixan is set.
  ne = 0
  ne3 = 0

  ! Find nv, number of points in the verticle grid.
  estep0 = 0.01/hart

  ! Exponential grid em(ne+1*j) = emin*exp(j*del)
  ! del = 0.6 is ok for Cu K edge, but needs more testing
  del = 0.4d0

  ! n1 is the # of points in a grid defined by estep0*exp(j*del) that lie below xloss.
  n1 = NINT(LOG(xloss/estep0)/del - 0.5)
  if (n1.le.0) n1 = 1
  ! Now redefine the grid so that xloss is halfway between em(n) and em(n+1)
  ! Solving
  !    xloss = [em(n1) + em(n1+1)]/2 = emin*[exp(n1*del) + exp((n1+1)*del)]/2
  ! gives emin = 2*xloss/(1+exp(del))*exp(-n1*del)
  expdel = EXP(del)
  emin = 2*xloss /(1+expdel)/expdel**n1
  if (emin.le.estep0) emin = emin*expdel
  ! Now change grid so that endpoint is at emax.
  emax = MIN(50.d0/hart,20.d0*xloss)
  nv = NINT(log(emax/emin) / del ) + 3

  IF(vixan.gt.0.0001) THEN
    xim = vixan
  ELSE
    xim = xloss/2.d0
  END IF

  IF(ispec.eq.2) THEN
    xkmax = xkmax/bohr/hart
    xkmax = xkmax - edge
    xkstep = xkstep/bohr/hart
    xkstep = xkstep - edge
  END IF
  ik0 = 0

  IF (normal_mesh) THEN
    CALL mesh_normal_T(em, ne, ne1, ne3, nex, ik0, ecv, edge, xloss, iGrid, &
                        & ispec, estep0/2.d0)
  ELSE
    ! Not tested
    CALL mesh_sommerfeld(em,ne,ne1,ne3,nex,ik0,ecv,edge,xloss,iGrid,ispec,&
                        & estep0,xkmax,xkstep,nv)
  ENDIF

  OPEN (unit=44, file='emesh.dat', status='unknown')
  WRITE(44,'(a,3(1x,f12.5))') '# edge, bohr, edge*hart ', edge, bohr, edge*hart
  WRITE(44,'(a,2(1x,i5))') '# ispec, ik0 ', ispec, ik0
  WRITE(44,*) '# ie, em(ie)*hart, xk(ie)'
  DO ie = 1, ne
    WRITE (44,'(i5, 3f20.5)') ie, dble(em(ie))*hart, getxk(dble(em(ie))-edge)/bohr
  END DO
  CLOSE(unit=44)

  open(unit=44, file='emesh.bin',form='unformatted')
  write(44) ne, ne1, ne3
  write(44) em(:ne)
  close(44)
  RETURN
END

SUBROUTINE mesh_normal_T(em, ne, ne1, ne3, nex, ik0, ecv, edge, xloss, iGrid,&
                            & ispec, vemin)
  USE constants
  USE xsph_inp, ONLY: electronic_temperature
  IMPLICIT NONE
  ! Adapted from mk_rhorrp_grid
  ! Strange: i have problem passing gamach and vi0 from phmesh2T. Why ?
  ! Well, given that the zero temperature grid is divided into
  ! horizontal part and vertical part
  ! we must use the same convention

  ! Grid connects vertices:
  !   1) (ecv,0)
  !   2) (ecv,e1)
  !   3) (emax,e1)
  !   4) (emax,0) The fermi function vanishes on leg (3)->(4)[ not included ]
  !
  !   (leg 2)                         (leg 3)
  !  ecv+1j*e1---------------------> emax+1j*e1
  !     ^               X                |
  !     |               X                |  (vanishes)
  !     |               X                v
  !    ecv-------------edge-------------emax
  !  (leg 1)                           (leg 4)
  !
  ! This routine determines values for e1 and emax,
  ! and generates the energy grid and calculates the
  ! Matsubara poles (np) lying within the contour

  ! Input:
  ! em      - Empty energy grid
  ! nex     - Total number of points for horizontal and vertical grid
  ! ecv     - Core-valence separation
  ! edge    -  xmu - vr0, where vr0 from the EXCHANGE card
  ! xloss   - Broadening
  ! iGrid   - is there any user defined grid
  ! vemin   - The smallest step in vertical grid [Used with user-defined grid]
  INTEGER, INTENT(IN) :: nex, iGrid, ispec
  DOUBLE PRECISION, INTENT(IN) :: ecv, edge, vemin

  ! Output:
  ! em      - Filled energy grid
  ! ne      - Total number of energy points [ np = ne-10-2*ne1 ]
  ! ne1     - Number of points on horizontal grid (leg 2)
  ! ne3     - 0 since last leg vanishes (leg 3)
  ! ik0     - this is obselete
  COMPLEX*16, INTENT(INOUT) :: em(nex)
  INTEGER, INTENT(OUT) :: ne, ne1, ne3, ik0

  ! Local Variables:
  ! e1      - Maximum height for imaginary plane
  ! min_e1  - Minimum maximum Height for imaginary plane
  ! emax    - Maximum real energy for grid
  ! temp    - Electronic temperature in Hartree
  ! np      - Number of poles enclosed by contour
  ! n1      - Number of points for first leg  (vertical) [ Default: 10 ]
  ! n2      - Number of points for second leg (horizontal)
  ! de      - Separation between imaginary points
  ! NPts    - index of the last point added to the energy grid
  ! nGrid   - number of user defined grid
  REAL(8)  :: e1, emax, emin, de, temp, min_e1, xloss
  INTEGER :: n1, n2, np, ne4
  INTEGER :: NPts, nGrid, nemax
  INTEGER, PARAMETER :: nGridMax=10
  INTEGER :: iGridType(nGridMax)
  DOUBLE PRECISION :: GridMin(nGridMax),GridMax(nGridMax),GridStep(nGridMax)
  DOUBLE PRECISION :: del
  COMPLEX*16 :: emH(nex)
  ! Loop variables:
  INTEGER :: i
  CHARACTER(512) message

  ! Initalize
  temp = electronic_temperature / hart
  ! IF (temp < 0.001) temp = 0.001 ! ~ 315 K (Should not be here)

  ! Determine imaginary part based on temperature (this should be configurable / overridable)
  ! For now, use the smallest multiple of 2*pi*kT that is greater than 0.05 Hartree (~1.4 eV).
  ! At low temperatures, there is a tradeoff between the number of points needed on the
  ! horizontal leg (which decreases with e1), and the number of Matsubara poles (which increases
  ! with e1). We could be a bit more intelligent here...


  ! For a reference, xloss/np.pi = 0.011 hatree
  np = 1
  e1 =  np * 2*pi*temp
  ! min_e1 = xloss
  min_e1 = 0.05 !TODO: At low temperature this gives too many poles.
  ! min_e1 = 0.1
  IF (e1 < min_e1) THEN
    np = ceiling(min_e1 / (2*pi*temp))
    e1 = np * 2*pi*temp
  END IF

  WRITE(message,'(a,i3)') "Number of poles = ", np
  CALL wlog(message)
  ! The Fermi function is negligible (5e-5) at 10T past the edge.
  ! emax = edge + 10 * temp
  ! emax = edge + 5
  emax = 5.8 ! Should make it temperature dependent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Construct Horizontal legs !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (iGrid.eq.0) THEN ! Default grid
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Horizontal leg (2) -> (3) !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Imaginary part gives broadening of similar order,
    ! so set step size based on imaginary part
    de = e1 / 4.
    emin = ecv

    n2 =  min(ceiling((0 - emin) / de), 21)
    ik0 = n2 ! Determine the 0 position
    de = (0 - emin) / n2 ! Determine the spacing
    n2 = ceiling((emax-emin)/ de) ! Determine the total number of points

    WRITE(message,'(a,i3)') "Horizontal grid, n2 = ", n2
    CALL wlog(message)

    DO i = 1, n2
      em(i) = emin + i * de + e1 * coni
    END DO
    ne1 = n2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Horizontal along xloss    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO i = 1, n2
      em(n2+i) = emin + i * de + xloss * coni

    END DO
    ne1 = n2
  ELSE ! User defined grid
    ! Make sure there are enough points left over to make vertical grid etc.
    nemax = nex
    n2 = 0

    ! Read grid.inp
    CALL RdGrid(emH,n2,nGrid,iGridType,GridMin,GridMax,GridStep,nGridMax,nemax)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Horizontal leg (2) -> (3) !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO i = 1, nGrid
      IF(iGridType(i).eq.1) THEN
        ! grid is regular in energy
        n2 = n2 + 1
        CALL MkEMesh(emH, n2, GridMin(i), GridMax(i),GridStep(i), NPts, nex)
        n2 = MIN(n2 + NPts, nemax)
      ELSEIF(iGridType(i).eq.2) THEN
        ! grid is regular in k
        n2 = n2 + 1
        CALL MkKMesh(emH, n2, GridMin(i), GridMax(i),GridStep(i), NPts, nex)
        n2 = MIN(n2 + NPts, nemax)
      ELSEIF(iGridType(i).eq.3) THEN
        ! grid is exponential
        n2 = n2 + 1
        del = LOG(GridStep(i))
        CALL MkExpMesh(emH, n2, 1.d0, GridMax(i)-GridMin(i)+1,GridStep(i), NPts, nex)
        emH(n2:) = emH(n2:) + GridMin(i) - 1.d0
        n2 = MIN(n2 + NPts, nemax)
      END IF
    END DO

    ! Dont know why the sort does not remove the degenerate 0
    ! Add a point at E = 0 in case there is not one.
    ! If FPRIME, add a point at edge.
    ! IF(n2+1.lt.nex) THEN
    !   emH(n2+1) = 0.d0
    !   n2 = n2 + 1
    ! ELSE
    !   emH(1) = 0.d0
    ! END IF

    ! Now, sort energy grid and remove degenerate points.
    CALL SortE(emH,n2,ik0,nex)

    ! Add imaginary part
    DO i = 1, n2+n2
      IF (i.LT.n2+1) THEN
        em(i) = emH(i)  + edge+ e1*coni
      ELSE
        em(i) = emH(i-n2)  +edge + xloss*coni
      ENDIF
    ENDDO
    ne1 = n2
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Vertical leg (1) -> (2) !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Since we start from -Inf to Inf, n1 = 0 (WRONG !!!! we start from ecv)

  n1 = 10 ! TODO: Remove Hard coded into many files.
  IF (n1.GT.0) THEN
    de = e1 / n1**2
    DO i = 1, n1
      em(2*ne1+i) = ecv + de * i**2 * coni
    END DO
  ENDIF

  !!!!!!!!!!!!!!!!!!!
  ! Matsubara poles !
  !!!!!!!!!!!!!!!!!!!
  DO i = 1,np
    em(2*ne1 + n1 + i) = edge + coni * (2*i - 1) * pi * temp
  END DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fake pole From T=0 module !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  em(2*ne1 + n1 + np + 1) = edge + coni * vemin

  ne3 = 0 ! TODO: let ne3 be 1st leg instead of 3rd leg
  ne = 2*ne1 + n1 + np + 1
END SUBROUTINE mesh_normal_T

SUBROUTINE mesh_sommerfeld(em, ne, ne1, ne3, nex, ik0, ecv, edge, xloss, iGrid,&
                            & ispec, estep0, xkmax, xkstep, nv)
  ! This is the same phmesh2 subroutine but with an extra point
  ! for computing derivative at efermi
  USE constants
  USE xsph_inp, ONLY: electronic_temperature
  IMPLICIT NONE

  ! Input:
  ! em      - Empty energy grid
  ! nex     - Total number of points for horizontal and vertical grid
  ! ecv     - Core-valence separation
  ! edge    -  xmu - vr0, where vr0 from the EXCHANGE card
  ! xloss   - Broadening
  ! iGrid   - is there any user defined grid
  ! vemin   - The smallest step in vertical grid [Used with user-defined grid]
  INTEGER, INTENT(IN) :: nex, iGrid, ispec, nv
  DOUBLE PRECISION, INTENT(IN) :: ecv, edge, estep0, xkmax, xkstep

  ! Output:
  ! em      - Filled energy grid
  ! ne      - Total number of energy points [ np = ne-10-2*ne1 ]
  ! ne1     - Number of points on horizontal grid (leg 2)
  ! ne3     - 0 since last leg vanishes (leg 3)
  ! ik0     - this is obselete
  COMPLEX*16, INTENT(INOUT) :: em(nex)
  INTEGER, INTENT(OUT) :: ne, ne1, ne3, ik0

  ! Local Variables:
  ! e1      - Maximum height for imaginary plane
  ! min_e1  - Minimum maximum Height for imaginary plane
  ! emax    - Maximum real energy for grid
  ! temp    - Electronic temperature in Hartree
  ! np      - Number of poles enclosed by contour
  ! n1      - Number of points for first leg  (vertical) [ Default: 10 ]
  ! n2      - Number of points for second leg (horizontal)
  ! de      - Separation between imaginary points
  ! NPts    - index of the last point added to the energy grid
  ! nGrid   - number of user defined grid
  REAL(8)  :: e1, emax, emin, de, tk, xloss
  INTEGER :: n1, n2, np, ne4
  INTEGER :: NPts, nGrid, nemax
  INTEGER, PARAMETER :: nGridMax=10
  INTEGER :: iGridType(nGridMax)
  DOUBLE PRECISION :: GridMin(nGridMax),GridMax(nGridMax),GridStep(nGridMax)
  DOUBLE PRECISION :: del
  COMPLEX*16 :: emH(nex)
  ! Loop variables:
  INTEGER :: ii


  tk = electronic_temperature / hart

  !--- Construct horizontal grid ----
  IF (iGrid.EQ.0) THEN
    CALL XanesGrid84(em,xkmax,xkstep,estep0,ne,ik0,nex,ispec)
    ne1 = ne
  ELSE
    nemax = nex - nv - 1
    ne = 0
    CALL RdGrid(em, ne, nGrid, iGridType, GridMin, GridMax, GridStep, nGridMax,&
              & nemax)
    DO ii = 1, nGrid
      IF (iGridType(ii).EQ.1) THEN
        ne = ne+1
        CALL MkEMesh(em, ne, GridMin(ii), GridMax(ii), GridStep(ii), NPts, nex)
        ne = MIN(ne + NPts, nemax)
      ELSEIF (iGridType(ii).EQ.2) THEN
        ne = ne+1
        CALL MkKMesh(em, ne, GridMin(ii), GridMax(ii), GridStep(ii), NPts, nex)
        ne = MIN(ne + NPts, nemax)
      ELSEIF (iGridType(ii).EQ.3) THEN
        ne = ne+1
        del = LOG(GridStep(ii))
        CALL MkKMesh(em, ne, 1.d0, GridMax(ii)-GridMin(ii)+1, GridStep(ii), &
                     & NPts, nex)
        em(ne:) = em(ne:) + GridMin(ii) - 1.d0
        ne = MIN(ne + NPts, nemax)
      ENDIF
    ENDDO
  ENDIF

  ! Add a point at E = 0 in case there is not one.
  ! If FPRIME, then add a point at edge.
  IF (ne+1.LT.nex) THEN
    em(ne+1) = 0.d0
    ne=ne+1
  ELSE
    em(1) = 0.d0
  ENDIF

  ! Now sort energy grid and remove degenerate points
  CALL SortE(em, ne, ik0, nex)
  ne1 = ne
  ne3 = 0

  DO ii = 1, ne
    em(ii) = em(ii) + edge + coni*xloss
  ENDDO


  !--- Construct Vertical grid ----
  ne = ne + 1
  CALL MkVGrid84(em,ne,xloss,nex)
  DO ii = ne1+1, ne
    ! Shift vertical grid by edge
    em(ii) = em(ii) + edge
  ENDDO

  !--- Add points for gradient ----
  ne = ne + 3
  ! em(ne-2) = edge*(1.d0-1e-4) + coni*xloss
  ! em(ne-1) = edge + coni*xloss
  ! em(ne) = edge*(1.d0+1e-4) + coni*xloss
  em(ne-2) = edge + coni*estep0
  em(ne-1) = edge + coni*estep0/2.0
  em(ne) = edge
END SUBROUTINE mesh_sommerfeld
