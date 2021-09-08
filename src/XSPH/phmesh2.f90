!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: phmesh2.f90,v $:
! $Revision: 1.17 $
! $Author: jorissen $
! $Date: 2012/11/20 02:16:17 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE phmesh2 (iprint, ispec, edge, emu, vi0, gamach,        &
     &     ecv, xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3, iGrid)
!     This subroutine makes the energy mesh used for phases and cross sections,
!     as well as for the fms routine, path, and genfmt. For EXAFS, the final output
!     chi is on a different (usually finer) grid with mu0 interpolated.
!     This will reproduce the old (FEFF84) grids, as well as any combination of user
!     defined energy, k, exponential, or arbitrary (read from file) grids. The input
!     for the user defined grids is read from grid.inp. Details of grid.inp are given
!     in rdgrid.f
      use constants
      use DimsMod, only: nex 
	  
      IMPLICIT NONE
!     Input:
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

      INTEGER iprint, ispec, iGrid
      DOUBLE PRECISION edge, emu, vi0, ecv, gamach, xkmax, xkstep, vixan
      
!     Output:
!     ne     - Total number of energy points.
!     ne1    - Number of energy points on the horizontal grid.
!     ik0    - point where k=0
!     em(ne) - energy array
      INTEGER ne, ne1, ne3, ik0
      COMPLEX*16 em(nex)

!     Local Variables:
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
!     User defined grid variables
!     nGridMax  - max number of grids
!     nGrid     - number of grids.
!     iGridType - Type of grid (1 = energy, 2 = k, 3 = exp)
!     GridMin   - minimum k or E of grid. k for k-grids, e for others
!     GridMax   - Maximum k or E of grid
!     GridStep  - step size.
      INTEGER nGridMax
      PARAMETER(nGridMax=10)
      INTEGER nGrid, iGridType(nGridMax)
      DOUBLE PRECISION GridMin(nGridMax), GridMax(nGridMax),            &
     &     GridStep(nGridMax)
      DOUBLE PRECISION estep0, expdel, de
      
!     Loop Variables:
      INTEGER ii, nv, n1, NPts, ie
      DOUBLE PRECISION getxk
      EXTERNAL getxk
!     Initialization


!	  write(*,*) 'in phmesh2'
!	  write(*,*) 'iprint, ispec, edge, emu, vi0, gamach',iprint, ispec, edge, emu, vi0, gamach
!	  write(*,*) 'xkmax, xkstep, vixan, ne, ne1,ne3, iGrid',xkmax, xkstep, vixan, ne, ne1,ne3, iGrid



!     Set total imaginary part, must be >= 0.02 eV
      xloss = MAX(gamach/2.d0 + vi0, 0.02d0/hart)
!     Set energy step to half of imaginary part, or
!     vixan if vixan is set.
      ne = 0
      ne3 = 0

!     Find nv, number of points in the verticle grid.
      estep0 = 0.01/hart
!     Exponential grid em(ne+1*j) = emin*exp(j*del)
!     del = 0.6 is ok for Cu K edge, but needs more testing
      del = 0.4d0
!     n1 is the # of points in a grid defined by estep0*exp(j*del) that lie below xloss.
      n1 = NINT(LOG(xloss/estep0)/del - 0.5)
      if (n1.le.0) n1 = 1
!     Now redefine the grid so that xloss is halfway between em(n) and em(n+1) 
!     Solving
!     xloss = [em(n1) + em(n1+1)]/2 = emin*[exp(n1*del) + exp((n1+1)*del)]/2
!     gives
!     emin = 2*xloss/(1+exp(del))*exp(-n1*del)
      expdel = EXP(del)
      emin = 2*xloss /(1+expdel)/expdel**n1
      if (emin.le.estep0) emin = emin*expdel
!     Now change grid so that endpoint is at emax.
      emax = MIN(50.d0/hart,20.d0*xloss)
      nv = NINT( log(emax/emin) / del ) + 3



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
      
      IF(iGrid.eq.0) THEN
!     Use FEFF84 grids
         IF(ispec.eq.0) THEN
            ne = 1
            CALL ExafsGrid84(em, xkmax, ne, nex)
            ne1 = ne
            ik0 = 1
         ELSEIF((ABS(ispec).gt.0).and.((ABS(ispec).lt.4).and.(ispec.ne.2))) THEN
!     Use same grid for XANES, DANES
            CALL XanesGrid84(em, xkmax, xkstep, xim, ne, ik0, nex, ispec)
            ne1 = ne
         ELSEIF(ispec.eq.2) THEN ! XES
            ne = 1
            CALL MkEMesh(em, ne, xkmax, xkstep,xim, NPts, nex-nv-1)
            ne = ne + NPts
!           Add a point at the fermi level.
!            em(ne+1) = 0.d0
!            ne = ne+1
            CALL SortE(em,ne,ik0,nex-nv)
            ne1 = ne
         ELSEIF(ispec.eq.4) THEN
!     FPRIME
            CALL FPrimeGrid84(em, xkmax, xkstep, vixan, emu, edge, ne, ne1, ne3, nex)
			!KJ 11-2012 ik0 is not initialized!  Leaving it to ik=0 causes array bounds problems in path.
			ik0=1
         ELSEIF(ispec.eq.5) THEN
           call mk_rhorrp_grid(em, ne, ne1, ne3, nex, ecv, edge)
         END IF
         
!     If ispec is negative, we are not running FMS. Make EXAFS grid
!     for points above the fermi level.
         IF(ispec.lt.0) THEN
            ne = 11
            CALL ExafsGrid84(em, xkmax, ne, nex)
            ne1 = ne
         END IF
      ELSE
!     User defined grids.
         ! Make sure there are enough points left over to make vertical grid etc.
         nemax = nex - nv - 1
         ne = 0
         CALL RdGrid(em,ne,nGrid,iGridType,GridMin,GridMax,GridStep,    &
     &        nGridMax,nemax)

         DO ii = 1, nGrid
            IF(iGridType(ii).eq.1) THEN
               ! grid is regular in energy
               ne = ne + 1
               CALL MkEMesh(em, ne, GridMin(ii), GridMax(ii),           &
     &              GridStep(ii), NPts, nex)
               ne = MIN(ne + NPts, nemax)
            ELSEIF(iGridType(ii).eq.2) THEN
               ! grid is regular in k
               ne = ne + 1
               CALL MkKMesh(em, ne, GridMin(ii), GridMax(ii),           &
     &              GridStep(ii), NPts, nex)
               ne = MIN(ne + NPts, nemax)
            ELSEIF(iGridType(ii).eq.3) THEN
               ! grid is exponential
               ne = ne + 1
               del = LOG(GridStep(ii))
               CALL MkExpMesh(em, ne, 1.d0, GridMax(ii)-GridMin(ii)+1,    &
     &              GridStep(ii), NPts, nex)
               em(ne:) = em(ne:) + GridMin(ii) - 1.d0
               ne = MIN(ne + NPts, nemax)
            END IF
         END DO
         
!        Add a point at E = 0 in case there is not one. If FPRIME
!        add a point at edge.
         IF(ne+1.lt.nex) THEN
            em(ne+1) = 0.d0         
            ne = ne + 1
         ELSE
            em(1) = 0.d0
         END IF
         
!        Now, sort energy grid and remove degenerate points.
         CALL SortE(em,ne,ik0,nex)
         ne1 = ne
         ne3 = 0
      END IF

!     If XES, flip grid about 0.0
!      IF(ABS(ispec).eq.2) THEN
!         CALL ReverseGrid(em,ne,0.d0)
!         ik0 = ne1 + 1 - ik0 
!      END IF
!     Shift horizontal grid by edge + coni*xloss.
      IF(abs(ispec).lt.4) THEN
         DO ii = 1, ne
            em(ii) = em(ii) + edge + coni*xloss
         END DO
!      ELSE
!     Now fill another grid for the KK-Transform in FPRIME
!         nemax = MIN(nex-ne,100)
!         del = 3.d0/hart
         
!     Set elimit = 20*emu, but make sure that 1.d3 .le. elimit .le. 2.d5
!         elimit = MAX(1.d3/hart,MIN(20.d0*emu,2.d5/hart))
!         elimit = elimit - emu
      
!         ne3 = nemax
!         em(ne1+1) = edge
!         DO ii = 1, ne3-1
!            del2 = 0.d0
!            IF(DBLE(em(ne1+ii)).gt.0.d0) del2 = em(ne1+ii)*             &
!     &           (EXP( LOG( elimit/em(ne1+ii) ) / (ne3-ii) ) -1)
!            em(ne1+ii+1) = em(ne1+ii) + MAX(del,del2)
!         END DO
!         ne = ne1 + ne3
      END IF

!     If not fprime or rhorrp calculation, make vertical grid
      IF(0.le.abs(ispec) .and. abs(ispec) .le. 3) THEN
         ne = ne + 1
         CALL MkVGrid84(em, ne, xloss, nex)
!     Shift vertical grid by edge.
         DO ii = ne1+1, ne
            em(ii) = em(ii) + edge
         END DO
      END IF
         
      IF(ABS(ispec).eq.3) THEN
!     DANES: add more points to horizontal grid.
         ne3  = MIN(nex,150) - ne
         emin = DBLE(2*em(ne1)-em(ne1-1))
         emax = 7.d4
         del  = LOG(emax/emin)/(ne3-1)
         ne = ne + 1
         CALL MkExpMesh(em, ne, emin, emax, del, ne3, nex)
         em(ne:nex) = em(ne:nex) + coni*1.d-8
         DO ii = 0, ne3 - 1
            em(ne+ii) = em(ne+ii) + coni*1.d-8
         END DO
         ne = ne + ne3 - 1
      END IF      

!      IF (iprint .ge. 1)  THEN
         OPEN (unit=44, file='emesh.dat', status='unknown')
         WRITE(44,'(a,3(1x,f12.5))') '# edge, bohr, edge*hart ', edge, bohr, edge*hart
         WRITE(44,'(a,2(1x,i5))') '# ispec, ik0 ', ispec, ik0
         WRITE(44,*) '# ie, em(ie)*hart, xk(ie)'
         DO ie = 1, ne
           WRITE (44,'(i5, 3f20.5)') ie, dble(em(ie))*hart, getxk(dble(em(ie))-edge)/bohr
         END DO
         CLOSE(unit=44)
!      END IF

      open(unit=44, file='emesh.bin',form='unformatted')
      write(44) ne, ne1, ne3
      write(44) em(:ne)
      close(44)

      RETURN
      END

      SUBROUTINE MkVGrid84(em, ne, xloss, nex)
!     make the vertical grid in energy plane
!     first point is at 0.005 ev, second at 0.01 ev and
!     exponential grid with step 0.4 after that up to 50 eV
      use constants
!     Input:
!     ne    - first energy point
!     nex   - length of em array
!     xloss - total imaginary part of horizontal grid
      INTEGER ne, nex
      DOUBLE PRECISION xloss

!     Output:
!     ne      - number of energy points
!     em(nex) - energy grid
      COMPLEX*16 em(nex)

!     Local Variables:
!     n1     - number of points in exponential grid
!     estep0 - first two points are at estep0/2 and estep0
!     del    - spacing: em(j) = emin*exp(j*del)
!     expdel - exp(del)
!     emin   - minimum energy in exponential grid.
!     emax   - max energy in exponential grid
      INTEGER n1, nv
      DOUBLE PRECISION estep0, del, expdel, emin, emax

!     Loop Variables:
      INTEGER ii

      estep0 = 0.01/hart
!     Exponential grid em(ne+1*j) = emin*exp(j*del)
!     del = 0.6 is ok for Cu K edge, but needs more testing
      del = 0.4d0

!     n1 is the # of points in a grid defined by estep0*exp(j*del) that lie below xloss.
      n1 = NINT(LOG(xloss/estep0)/del - 0.5)
      if (n1.le.0) n1 = 1

!     Now redefine the grid so that xloss is halfway between em(n) and em(n+1) 
!     Solving
!     xloss = [em(n1) + em(n1+1)]/2 = emin*[exp(n1*del) + exp((n1+1)*del)]/2
!     gives
!     emin = 2*xloss/(1+exp(del))*exp(-n1*del)
      expdel = EXP(del)
      emin = 2*xloss /(1+expdel)/expdel**n1
      if (emin.le.estep0) emin = emin*expdel

!     Now change grid so that endpoint is at emax.
      emax = MIN(50.d0/hart,20.d0*xloss)

!     Make sure there are enough points for the vertical grid.
!      nv = NINT( log(emax/emin) / del ) + 2
!      ne = MIN(ne, nex-nv)
      
      em(ne) = coni*estep0/2
      em(ne+1) = coni*estep0
      ne = ne + 2
!     Josh         if (emin.le.estep0 .or. emin.ge.xloss) 
!     Josh     .     call par_stop(' Bad mesh in phmesh')
!     delk = log (xloss/tempk) /(n1+0.5)

      CALL MkExpMesh(em, ne, emin, emax, del, n1, nex)
      DO ii = 0, n1
         em(ne+ii) = (0,1)*em(ne+ii)
      END DO
      ne = ne + n1

      RETURN
      END
      
      SUBROUTINE MkExpMesh(em, iStart, emin, emax, del, NPts, nex)
      
      INTEGER iStart, nex
      DOUBLE PRECISION emin, emax, del
      COMPLEX*16 em(nex)
      
      INTEGER NPts

      INTEGER ii

      NPts = NINT( log(emax/emin) / del )
      
!     Fill grid
      DO ii = 0, NPts
         em(iStart+ii) = emin*exp(del*ii)
      END DO
         
      RETURN
      END


      SUBROUTINE ExafsGrid84(em, xkmax, ne, nex)
!     Make old (FEFF8.4) grid for EXAFS calculations.
      use constants
!     Input: 
!     xkmax - maximum k for grid
!     nex   - length of array em
      INTEGER nex
      DOUBLE PRECISION xkmax

!     Output:
!     em(nex) - energy grid array
!     ne      - number of points in energy grid
      COMPLEX*16 em(nex)
      INTEGER ne

!     Local Variables:
!     NPts   - Number of points that have been added to grid after a call
!              to MkKMesh
!     nemax  - maximum number of energy points (100)
!     deltak - k step (used when calling MkKMesh)
!     xkmin  - minimum (k used when calling MkKMesh)
!     xkmax2 - maximum (k used when calling MkKMesh)
      INTEGER NPts, nemax
      DOUBLE PRECISION deltak, xkmin, xkmax2, eps
      PARAMETER (small = 1.d-20)
      nemax = (MIN(100,nex))

!     20 pts (0 le k le 1.9, delk=0.1 ang(-1) )
      deltak = bohr/10
      xkmin = 0.d0
      xkmax2  = bohr*1.9d0 + deltak*0.01d0
      CALL MkKMesh(em, ne, xkmin, xkmax2, deltak, NPts, nemax)
!     20 pts (2 le k le 5.8, delk=0.2 ang(-1) )
      ne = ne + NPts + 1
      deltak = bohr/5
      xkmin  = bohr*2.d0
      xkmax2 = bohr*5.8d0+deltak*0.01d0
      CALL MkKMesh(em, ne, xkmin, xkmax2, deltak, NPts, nemax)
!     9 pts (6 le k le 10., delk=0.5 ang(-1) )
      ne = ne + NPts + 1
      xkmin = bohr*6.d0
      xkmax2 = bohr*10.d0+deltak*0.01d0
      deltak = bohr*0.5d0
      CALL MkKMesh(em, ne, xkmin, xkmax2, deltak, NPts, nemax)
!     make the rest of the points pts with deltak = 1.0 ang(-1)
      ne = ne + NPts + 1
      deltak = bohr
      xkmin = SQRT(2*DBLE(em(ne-1))) + deltak
!     Fill to end of grid, or max # of points.
      NPts = MIN(nemax-ne,NINT((xkmax-xkmin)/deltak)+1)
      xkmax2 = xkmin + (NPts-1)*deltak+deltak*0.01d0
      CALL MkKMesh(em, ne, xkmin, xkmax2, deltak, NPts, nemax)
      ne = ne + NPts

      RETURN
      END

      SUBROUTINE XanesGrid84(em, xkmax, xkstep, estep, ne, ik0, nex, ispec)
!     Make old (FEFF8.4) grid for XANES calculations.
      use constants
!     Input: 
!     xkmax  - maximum k for grid
!     xkstep - kstep at high k
!     estep  - estep near the fermi level
!     nex    - length of array em
      INTEGER nex, ispec
      DOUBLE PRECISION xkmax, xkstep, estep

!     Output:
!     emin    - -xim*n1
!     em(nex) - energy grid array
!     ik0     - zero point for k grid
      INTEGER ne
      COMPLEX*16 em(nex)

!     Local Variables:
      INTEGER n1, n2, nk, nemax, NPts
      DOUBLE PRECISION emin, emax, xkmin, dk, xkmax2
!     Make 10 points below fermi level
      nemax = 10
!     double k step below fermi level 
      dk = 2*xkstep
!     Not sure why to pick this number of steps regular in e?
      n1 = INT(estep/2/dk**2)
!     n2 is starting point of k grid minus 1 (int(k(emax)/dk) + 1)
      n2 = INT(SQRT(n1*2*estep)/dk)

!     If we can fit one more point in the egrid, do it
      If( (dk*(n2+1))**2 .gt. (n1+1)*2*estep ) n1 = n1+1

!     Make sure we don't use more than nemax points
      n1 = MIN(n1,nemax)
!     nk is number of points in k grid
      nk = nemax - n1

!     Fill k grid below edge
      xkmin = -dk*(n2+nk)
      xkmax2 = -dk*(n2+1)
      ne = 1
      CALL MkKMesh(em, ne, xkmin, xkmax2, dk, nk, nex)  
      
!     Fill e grid near edge
      ne = ne + nk + 1
       
      emin = -estep*n1
      emax = 0.d0
      CALL MkEMesh(em, ne, emin, emax, estep, NPts, nex)
      ne = ne + NPts + 1
      ik0 = ne - 1
!     Fill grid above the fermi level.
!     Same grid as before except that k spacing is xkstep, and 90 points
      nemax = nex - ne
!     Not sure why to pick this number of steps regular in e?
      n1 = INT(estep/2/xkstep**2)
!     n2 is starting point of k grid minus 1
      n2 = INT(SQRT(n1*2*estep)/xkstep)
      n1 = n1 + 1
!     If we can fit one more point in the egrid, do it
      If( (xkstep*(n2+1))**2 .gt. (n1+1)*2*estep ) n1 = n1+1
!     Make sure we don't use more than nemax points
      n1 = MIN(n1,nemax)
!     For XES, xkmax is energy
      IF (ispec.ne.2) then
         nb = int (xkmax**2 /estep/2) + 1
      ELSE
         nb = int (abs(xkmax) /estep) + 1
      END IF
      IF(nb.le.n1) n1 = nb
!     nk is number of points in k grid
      nk = nemax - n1
      
!     This time fill e grid first
      emin = estep
      emax = (n1-1)*estep
!     If k(emax) > xkmax set emax = e(xkmax) and nk = 0
      IF(SQRT(2*emax).gt.xkmax.and.ispec.ne.2) THEN
         emax = xkmax**2/2
         nk = 0
      END IF
      CALL MkEMesh(em, ne, emin, emax, estep, NPts, nex)
 
!     Now fill k grid
      ne = ne + NPts + 1
      xkmin = xkstep*(n2+1)
      xkmax2 = xkstep*(n2+nk)
!     if xkmax2 > xkmax, set xkmax2 = xkmax
      IF(xkmax2.gt.xkmax) xkmax2 = xkmax
      CALL MkKMesh(em, ne, xkmin, xkmax2, xkstep, NPts, nex)  
      ne = ne + NPts
      
      RETURN
      END
      
      SUBROUTINE FPrimeGrid84(em,emin,emax,estep,emu,edge,ne,           &
     &     ne1,ne3,nex)
!     Make old (FEFF84) grid for FPRIME calculation
      use constants
!     Input:
!     emin  - minimum energy
!     emax  - maximum energy
!     estep - energy step
!     emu   - x-ray edge energy
!     edge  - fermi level (xmu-vr0)
!     nex   - size of em array
      DOUBLE PRECISION emin, emax, estep, emu, edge
      INTEGER nex

!     Output:
!     ne      - total number of points in energy grid
!     ne1     - number of energy points in regular grid
!     em(nex) - energy grid
      INTEGER ne, ne1, ne3
      COMPLEX*16 em(nex)
      
!     Local variables:
!     nemax - maximum number of points in constant energy grid
!     del   - step for exponential grid.
!     
      INTEGER nemax
      DOUBLE PRECISION del, del2, elimit
!     Loop Variables:
      INTEGER ii
!     Initialization
      nemax = 100

!      write(*,*) 'emin,emax,estep',emin,emax,estep
!	  write(*,*) 'emu,edge',emu,edge
!	  write(*,*) 'bohr, hart, 1/a0Ha ',bohr, hart, 1.d0/(bohr*hart)

      emin  = emin/bohr/hart - emu 
      emax  = emax/bohr/hart - emu
	  !rexsph added a spurious bohr to emin,emax.  It put estep (vixan) already in Ha.
	  !as a result, emin emax estep are now all in Ha (as opposed to eV in feff.inp).
	  
!	  write(*,*) 'emin,emax',emin,emax

!     Fill a grid from emin to emax taking steps estep.
      em(1) = emin
      ne = 1
      IF(emin.lt.emax) THEN
         IF(estep.le.0.d0) estep = (emax-emin)/(nemax-1)
         ne = MIN(nemax,NINT((emax-emin)/estep)+1) !KJ added +1   
!		 write(*,*) 'estep,ne,nemax',estep,ne,nemax   
         DO ii = 1, ne  
            em(ii) = emin + (ii-1)*estep !KJ added -1
         END DO
      END IF
      ne1 = ne

!     Now fill another grid for the KK-Transform
      nemax = MIN(nex-ne,100)
      del = 3.d0/hart

!     Set elimit = 20*emu, but make sure that 1.d3 .le. elimit .le. 2.d5
      elimit = MAX(1.d3/hart,MIN(20*emu,2.d5/hart))
      elimit = elimit - emu
      
!	  write(*,*) 'elimit',elimit
	  
      ne3 = nemax
      em(ne1+1) = edge
      DO ii = 1, ne3-1
         del2 = 0
         IF(DBLE(em(ne1+ii)).gt.0.d0) del2 = em(ne1+ii)*                &
     &        (EXP( LOG( elimit/em(ne1+ii) ) / (ne3-ii) ) -1)
         em(ne1+ii+1) = em(ne1+ii) + MAX(del,del2)
      END DO
      ne = ne1 + ne3
      
      RETURN
      END

      SUBROUTINE ReverseGrid(em,ne,ZeroPoint)
!     Flips a grid about ZeroPoint.
!     Input:
!     em(ne)    - array to flip
!     ne        - number of elements
!     ZeroPoint - point to flip about
      INTEGER ne
      COMPLEX*16 em(ne), eTmp
      DOUBLE PRECISION ZeroPoint

!     Loop Variables
      INTEGER ii, i2, np
      np = ne/2
      DO ii = 1, ne
         em(ii) = ZeroPoint - em(ii)
      END DO

      DO ii = 1, np
         eTmp = em(ii)
         em(ii) = em(ne+1-ii)
         em(ne+1-ii) = eTmp
      END DO            
      
      RETURN
      END

      SUBROUTINE MkEMesh(em, iStart, emin, emax, estep, NPts, nex)
!     Make a grid even in k-space from xkmin to xkmax with grid spacing
!     deltak. If xkmin > xkmax, do nothing
      IMPLICIT NONE
!     Input
!     em(nex) - energy grid
!     iStart  - index of em to start at.
!     emin    - starting k
!     emax    - ending k
!     estep   - k spacing
!     nex     - lenth of em array
      INTEGER iStart, nex
      COMPLEX*16 em(nex)
      DOUBLE PRECISION estep, emin, emax

!     Output:
!     NPts    - index of the last point added to the energy grid.
      INTEGER NPts

!     Loop variables
      INTEGER ii
      NPts = NINT((emax - emin)/estep) 
      IF(NPts.lt.0) THEN
         NPts = -1
         RETURN
      END IF
      DO ii=0, NPts
         IF(ii.lt.nex) THEN 
            em(iStart + ii) = emin + estep*ii
         ELSE
!            CALL wlog('Energy grid is too large: truncating.')
            NPts = ii-1
            EXIT
         END IF
      END DO

      RETURN
      END

      
      SUBROUTINE MkKMesh(em, iStart, xkmin, xkmax, deltak, NPts, nex)
!     Make a grid even in k-space from xkmin to xkmax with grid spacing
!     deltak.
      IMPLICIT NONE
!     Input
!     em(nex) - energy grid
!     iStart  - index of em to start at.
!     xkmin   - starting k
!     xkmax   - ending k
!     deltak  - k spacing
!     nex     - lenth of em array
      INTEGER iStart, nex
      COMPLEX*16 em(nex)
      DOUBLE PRECISION deltak, xkmin, xkmax

!     Output:
!     NPts    - index of the last point added to energy grid.
      INTEGER NPts

!     Loop variables
      INTEGER ii, isgn

      NPts = NINT((xkmax - xkmin)/deltak) 
      IF(NPts.le.0) THEN
         NPts = -1
         RETURN
      END IF
      isgn = 1
      IF(xkmin.lt.0.d0) isgn = -1
      DO ii=0, NPts
         IF(ii.lt.nex) THEN 
            em(iStart + ii) = isgn*(xkmin + deltak*(ii))**2/2
         ELSE
!            CALL wlog('Energy grid is too large: truncating.')
            EXIT
         END IF
      END DO

      RETURN
      END
 

      SUBROUTINE WrtE(em, ne, fl)
!     WrtE made for debugging.
      INTEGER ne, iU
      COMPLEX*16 em(ne)
      CHARACTER*(*) fl
      CHARACTER(300) fl2
      INTEGER ii

      fl2 = 'DEBUG/' // fl
      iU = 23
      OPEN(unit=iU,file=fl2,status='replace')
      DO ii = 1, ne
         WRITE(iU,*) ii, em(ii)
      END DO
      CLOSE(iU)

      RETURN
      END

      SUBROUTINE SortE(em,ne,ik0,nex)
!     Sorts energy array em, eliminating degenerate points.
!     Also, set ik0.
!     Input:
!     ne     - number of energy points
!     em(ne) - energy grid
      INTEGER ne
      COMPLEX*16 em(ne)

!     Output: sorted array of energies, and number of unique energy points.
!     Also ik0
      INTEGER ik0

!     Local Variables:
!     RealE(nex)  - Re[em]
!     iOrder(nex) - Holds ordering for em.
!     nUE         - number of unique energy points
!     tol         - tolerence for degeneracy of energy points (in eV)
      DOUBLE PRECISION RealE(nex), E0, tol
      INTEGER iOrder(nex), nUE
      
!     Loop Variables:
      INTEGER ii, i2
      
      PARAMETER(tol = 0.001d0)
      ik0 = -1
      
      DO ii = 1, ne
         RealE(ii) = DBLE(em(ii))
      END DO

!     Do sorting of RealE.
      CALL qsortd(iOrder,ne,RealE)
      
!     Replace em with sorted values and remove degeneracy.
      nUE   = 1
      IF((ABS(RealE(iOrder(1))).lt.tol)) THEN
         em(1) = 0.d0
         ik0 = 1
      ELSE
         em(1) = RealE(iOrder(1))
      END IF

!     Remove degenerate points
      nDeg = 0
      DO ii = 2, ne
         
!        find next point that is not degenerate and set next em to the
!        value of the non-degenerate point.
         DO i2 = ii + nDeg, ne
            IF(ABS(RealE(iOrder(i2))-DBLE(em(ii-1))).gt.tol) THEN
               nUE = nUE + 1
               em(nUE) = RealE(iOrder(i2))
               EXIT
            ELSE
               nDeg = nDeg + 1
            END IF
         END DO
      END DO
      ne = nUE
      
!     Set ik0
      ik0 = 1
      E0 = ABS(DBLE(em(1)))
      DO ii = 1, nUE
         IF(ABS(DBLE(em(ii))).lt.E0) THEN
            E0 = ABS(DBLE(em(ii)))
            ik0 = ii
         END IF
      END DO
      em(ik0) = 0.d0
      
      RETURN
      END

      subroutine mk_rhorrp_grid(em, ne, ne1, ne3, nex, ecv, edge)
        use constants
        use potential_inp, only: scf_temperature
        implicit none

        ! Grid connects vertices:
        !    1) (ecv,0)
        !    2) (ecv,e1)
        !    3) (emax,e1)
        !    [4) (emax,0)] ! The fermi function vanishes on leg (3)->(4), so it is not included
        !
        ! This routine determines values for e1 and emax, generates the energy grid and calculates the
        ! Matsubara poles (np) lying within the contour
        !
        ! ne1 is the number of points in legs 1->2->3
        ! ne2 is the number of Matsubara poles
        ! ne3 = 0
        ! ne = ne1 + ne2 + ne3  (implicity defines ne2)

        complex*16, intent(inout) :: em(nex)
        integer, intent(out) :: ne, ne1, ne3
        integer, intent(in) :: nex
        double precision, intent(in) :: ecv, edge

        real*8  :: e1, emax, de, temp, min_e1
        integer :: n1, n2, np, i

        temp = scf_temperature / hart
        if (temp < 0.001) temp = 0.001 ! ~ 315 K

        ! Determine imaginary part based on temperature (this should be configurable / overridable)
        ! For now, use the smallest multiple of 2*pi*kT that is greater than 0.05 Hartree (~1.4 eV).
        ! At low temperatures, there is a tradeoff between the number of points needed on the
        ! horizontal leg (which decreases with e1), and the number of Matsubara poles (which increases
        ! with e1). We could be a bit more intelligent here...
        e1 = 2 * pi * temp
        np = 1
        min_e1 = 0.05
        if (e1 < min_e1) then
          np = ceiling(min_e1 / (2*pi*temp))
          e1 = np * 2*pi*temp
        end if

        ! The Fermi function is negligible (5e-5) at 10T past the edge.
        emax = edge + 10 * temp

        ! vertical leg (1) -> (2)
        n1 = 10
        de = e1 / n1**2
        do i = 1,n1
          em(i) = ecv + de * i**2 * coni
        end do

        ! horizontal leg (2) -> (3)
        ! Imaginary part gives broadening of similar order, so set step size based on imaginary part
        de = e1 / 4.
        n2 = min(ceiling((emax - ecv) / de), 101)
        de = (emax - ecv) / n2
        do i = 1, n2
          em(n1+i) = ecv + i * de + e1 * coni
        end do

        ne1 = n1 + n2

        ! now include Matsubara poles
        do i = 1,np
          em(ne1 + i) = edge + coni * (2*i - 1) * pi * temp
        end do

        ne3 = 0
        ne = ne1 + np
      end subroutine
