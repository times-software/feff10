!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdgrid.f90,v $:
! $Revision: 1.7 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE RdGrid(em,ne,nGrid,iGridType,GridMin,GridMax,GridStep, &
     &     nGridMax,nex)
!     Read data from grid.inp.
!     file should have lines with the following format:
!     
!     Grid_Type    GridMin    GridMax    GridStep
!     
!     Grid_Type can be any of the following (case insensitive):
!     e_grid    - (regular in energy)
!     k_grid    - (regular in k)
!     exp_grid  - (exponential in energy)
!     user_grid - (read a grid from the file)
!
!     Energy, and k are given relative to the edge.
!     For e_grid, GridMin, GridMax, and GridStep are
!     given in (eV).
!     For exp_grid, E(i) = GridMin + exp(Ln(GridStep)*i) so that GridStep
!     is the size of the first energy step.
!     For k_grid, units are inverse angstroms.
!     user_grid is a special case and is followed by one energy
!     point per line, i.e.
!            user_grid
!            -1.01
!            -0.55
!            10.01
!              .!              .
!              .
!
!     More than one grid may be specified, and grids can
!     overlap, for example:
!
!     e_grid -10 10 0.1
!     k_grid  0  15 0.5
!
!     will make overlapping grids going from E = -10 eV to
!     k = 15 Angstrom**(-1). Up to 10 different grids can be
!     defined.
!     If the 'last' keyword is used in the GridMin field, i.e.
!        expgrid  last  100
!     the specified grid will start where the last grid ended.
!     This is usefull when defining non-overlapping k/e grids.
!     Comments lines have #,!,c, or * at the beginning.
      use constants

!     Input:
!     nGridMax - max number of grids that can be defined.
!     nex      - max number of energy points
      INTEGER nGridMax, nex

!     Output:
!     nGrid     - number of grids defined in file
!     iGridType - Type of grid. (0 = user, 1 = energy, 2 = k, 3 = exponential)
!     GridMin   - Minimum value of grid.
!     GridMax   - Maximum value of grid.
!     GridStep  - Step size.
!     ne        - number of energy points
!     em(nex)        - energy grid
      INTEGER nGrid, iGridType(nGridMax), ne, iErr
      DOUBLE PRECISION GridMin(nGridMax), GridMax(nGridMax),            &
     &     GridStep(nGridMax)
      COMPLEX*16 em(nex)

!     Local Variables:
!     ios      - i/o error flag
!     iUGrid   - unit number for grid.inp
!     RealE    - real part of energy
!     ImagE    - imaginary part of energy
!     nWords   - number of words in the line
!     Words(4) - array of words
!     line     - string to hold line
!     ieMin    - index of minimum of user defined grid
!     ieMax    - index of max of user defined frid
      INTEGER ios, iUGrid, nWords, ieMin, ieMax
      DOUBLE PRECISION RealE, ImagE
      CHARACTER(20) Words(10)
      CHARACTER(100) line
      CHARACTER(512) message

!     Loop Variables:
      INTEGER i1, i2

!     Externals
      LOGICAL isnum
      EXTERNAL isnum

      iUGrid = 22
      OPEN(unit=iUGrid,file='grid.inp',status='old',iostat=ios)
      CALL CHOPEN(ios, 'grid.inp', 'xsph')
      
      DO nGrid = 1, nGridMax
!        if all points used exit
         IF(ne.eq.nex) THEN
            CALL wlog("WARNING:" // &
                & " Too many energy points defined in grid.inp.")
            WRITE(message,'(a,i3)') "Maximum is ", nex
            CALL wlog(message)
            CALL wlog("Grid will be truncated.")
            EXIT
         END IF
!        Read comment lines
         CALL rdcmt(iUGrid,'#!*C',iErr)
         IF(iErr.EQ.-1) GOTO 5 ! End of file found
!        Read data line into string variable "line" and change to
!        lowercase.
         READ(iUGrid,'(A)',END=5) line
!         CALL lower(line)
!        bwords breaks line into words which are then passed
!        back in Words array
         nWords = 4
         CALL untab(line)
         CALL bwords(line,nWords,Words)

!        Set iGridType
         IF(Words(1).eq.'user_grid') THEN
            iGridType(nGrid) = 0
         ELSEIF(Words(1).eq.'e_grid') THEN
            iGridType(nGrid) = 1
!            IF(nwords.ne.4) CALL GridError('Error in grid.inp', line)
         ELSEIF(Words(1).eq.'k_grid') THEN
            iGridType(nGrid) = 2
         ELSEIF(Words(1).eq.'exp_grid') THEN
            iGridType(nGrid) = 3
         ELSE
            PRINT*, 'Error in grid.inp'
            PRINT*, 'Bad line:'
            PRINT*, line
            STOP
         END IF
         
         IF(iGridType(nGrid).ne.0) THEN
            READ(Words(3),*) GridMax(nGrid)
            READ(Words(4),*) GridStep(nGrid)
            IF(Words(2).eq.'last') THEN
               ! Set the grid minimum to the max of the last grid.
               IF(nGrid.gt.1) THEN
                  CALL SetGridMin(GridMin,GridMax,GridStep,iGridType,   &
     &                 nGrid)
               ELSE
                  GridMin(1) = 0.d0
               END IF
            ELSE
               READ(Words(2),*) GridMin(nGrid)
            END IF
         END IF

         IF(iGridType(nGrid).eq.0) THEN
!        User defined points: read from file.
            DO i2 = 1, nex
               ! Read comments
               CALL rdcmt(iUGrid,'#!*C',iErr)
               IF(iErr.EQ.-1) GOTO 5 ! End of file found
               ! Read line
               READ(iUGrid,'(A)',END=5) line
               nwords = 2
               ! break line into words
               CALL untab(line)
               CALL bwords(line,nWords,Words)
               ! if first word is number, Real(em) = num
               IF(isnum(Words(1))) THEN
                  READ(Words(1),*) RealE
                  ! if second word exists and is a num, Im(em) = num
                  ImagE = 0.d0
                  IF((nWords.ge.2).and.isnum(Words(2)))                 &
     &                 READ(Words(2),*) ImagE                  
                  em(i2) = (RealE + coni*ImagE)
                  ne = ne + 1
               ! If first word is not a number, exit loop and read line again.   
               ELSE
                  ! Set GridMax and GridMin for reference
                  GridMin(nGrid) = DBLE(em(ne - i2 + 1))
                  GridMax(nGrid) = DBLE(em(ne))

                  BACKSPACE(iUGrid)
                  EXIT
               END IF
            END DO
         END IF
      END DO
 5    CONTINUE
      nGrid = nGrid - 1

      DO i1 = 1, nGrid
         IF(iGridType(i1).eq.2) THEN
!     k-Grid. Set units to bohr**(-1)
            GridMin(i1) = GridMin(i1)*bohr
            GridMax(i1) = GridMax(i1)*bohr
            GridStep(i1) = GridStep(i1)*bohr
         ELSE
!     e-grid. Set units to hartrees
            GridMin(i1) = GridMin(i1)/hart
            GridMax(i1) = GridMax(i1)/hart
            GridStep(i1) = GridStep(i1)/hart
         END IF
      END DO
      DO i1 = 1, ne
         em(i1) = em(i1)/hart
      END DO

      CLOSE(iUGrid)
      RETURN
      END

      SUBROUTINE SetGridMin(GridMin, GridMax, GridStep, iGridType,      &
     &     nGrid)
!     This sets the minimum of the current grid to the maximum of the last grid + GridStep
      use constants
!     Input:
!     GridMin   - array that holds grid minima
!     GridMax   - array that holds grid maxima
!     GridStep  - array of steps
!     iGridType - array of grid types
!     nGrid     - current grid
      INTEGER nGrid
      INTEGER iGridType(nGrid)
      DOUBLE PRECISION GridMin(nGrid), GridMax(nGrid), GridStep(nGrid)

!     Output: GridMin(nGrid) (minimum of current grid.

      IF((iGridType(nGrid).ne.2).and.(iGridType(nGrid-1).ne.2).or.      &
     &    (iGridType(nGrid).eq.iGridType(nGrid-1))) THEN
!     If neither grid is a k grid, or if both are k-grid, just set the minimum to the previous
!     maximum.
         GridMin(nGrid) = GridMax(nGrid-1) + GridStep(nGrid)
      ELSEIF(iGridType(nGrid).eq.2) THEN
!     If current grid is k, kmin = sqrt(2*emax)
         GridMin(nGrid) =                                               &
     &        SQRT(2*GridMax(nGrid-1)/hart)/bohr + GridStep(nGrid)
      ELSE
!     If current grid is e, emin = k**2/2
         GridMin(nGrid) =                                               &
     &        (GridMax(nGrid-1)*bohr)**2/2*hart + GridStep(nGrid)
      END IF

      RETURN
      END

      SUBROUTINE GridError(message, line)
      CHARACTER(300) message, line
      
      CALL wlog(message)
      CALL wlog(line)
      STOP

      RETURN
      END
