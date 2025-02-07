SUBROUTINE fkgk(norb,dgc,dpc,rhoval_l,dr,inrm,iz)
  IMPLICIT NONE
  ! JK - Calculate Slater-Condon parameters Fk and Gk between core-orbitals and self-consistent valence
  INTEGER kap(41), norb, inrm, ilast, iz
  REAL(8) dgc(251,41), dpc(251,41), dr(251), rhoval_l(251), rho0(251), norm
  INTEGER nnum(41) ! Principal quantum number
  INTEGER li, lj ! Orbital angular momentum
  INTEGER ji,jj ! Total angular momentu *2
  INTEGER i,j,k,ki,kj,kmi
  REAL(8) xji, xjj ! Total angular momentum
  REAL(8) Fk, Gk, Fk0, Gk0, alpha
  REAL(8), EXTERNAL :: fkgk_int
  REAL(8), PARAMETER :: hart = 27.2114d0
  CHARACTER orb_ang_mom(0:4)
  CHARACTER(3) tot_ang_mom(0:4)
  CHARACTER(100) frmt
  !INTEGER, EXTERNAL, DIMENSION(41) :: principle_qn
!     principal quantum number (energy eigenvalue)
!  data nnum  /1,2,2,2,3,  3,3,3,3,4,  4,4,4,4,4, &
!       &            4,5,5,5,5,  5,5,5,6,6,  6,6,6,7,7, &
!       &            7,8,8,8,7,  7,6,6,5,5/
  data orb_ang_mom /'s','p','d','f','g'/
  data tot_ang_mom /'1/2','3/2','5/2','7/2','9/2'/
!  data kap /-1,-1, 1,-2,-1,  1,-2, 2,-3,-1,  1,-2, 2,-3, 3, &
!     &            -4,-1, 1,-2, 2, -3, 3,-4,-1, 1, -2, 2,-3,-1, 1, &
!     &            -2,-1, 1,-2, 2, -3, 3,-4, 4,-5/

  call get_qns(iz,nnum,kap)
  ilast = inrm
  DO i = 1, norb
     IF(nnum(i).LT.0) CYCLE
     norm = NORM2(dgc(:,i))
     IF(norm.lt.1.d-12) THEN
        CYCLE
     END IF
>>>>>>> 8094982 (Fixed slater condon printout)
     li = abs(kap(i))
     IF(kap(i).LT.0) THEN
        li = li - 1
        ji = 2*li + 1
        xji = li + 1.d0/2.d0
     ELSE
        xji = li - 1.d0/2.d0
        ji = 2*li - 1
     END IF
     DO j = 1, i
<<<<<<< HEAD
=======
        IF(nnum(j).LT.0) CYCLE
        norm = NORM2(dgc(:,j))
        IF(norm.LT.1.d-12) CYCLE
        IF(norm.lt.1.d-12) THEN
           CYCLE
        END IF
>>>>>>> 8094982 (Fixed slater condon printout)
        ki = abs(kap(i)) - 1
        kj = abs(kap(j)) - 1
        lj = abs(kap(j))
        IF(kap(j).LT.0) THEN
           lj = lj - 1
           jj = 2*lj + 1
           xjj = lj + 1.d0/2.d0
        ELSE
           xjj = lj - 1.d0/2.d0
           jj = 2*lj - 1
        END IF
        !PRINT*, nnum(i), orb_ang_mom(li), tot_ang_mom((ji-1)/2)
        !PRINT*, nnum(j), orb_ang_mom(lj), tot_ang_mom((jj-1)/2)
        kmi = 2*min(ki,kj)
        k = 0
        !PRINT*, 'k, Fk'
        frmt = '(I1,A1,A3,X,I1,A1,A3,X,A1,I1,X,3F10.5)'
        DO WHILE (k.LE.kmi) 
           Fk = fkgk_int(i,i,j,j,k,dgc,dpc,rhoval_l,dr,ilast,alpha)
           rho0(:) = dgc(:,i)*dgc(:,i)+dpc(:,i)*dpc(:,i)
           Fk0 = fkgk_int(i,i,j,j,k,dgc,dpc,rho0,dr,ilast,alpha)
           IF(Fk0.GT.0.d0) THEN
              WRITE(33,fmt = frmt) nnum(i), orb_ang_mom(li), tot_ang_mom((ji-1)/2), nnum(j), orb_ang_mom(lj),tot_ang_mom((jj-1)/2),'F', k, Fk*hart, Fk0*hart, Fk/Fk0
<<<<<<< HEAD
           ELSE
              RETURN
=======
           !ELSE
           !   RETURN
>>>>>>> 8094982 (Fixed slater condon printout)
           END IF
           !PRINT*, nnum(i), orb_ang_mom(li), tot_ang_mom((ji-1)/2)
           !PRINT*, nnum(j), orb_ang_mom(lj), tot_ang_mom((jj-1)/2)
           !PRINT*, k, Fk*hart
           k = k + 2
        END DO

        IF(i.EQ.j) CYCLE
        ! Gk integrals
        ki = abs(kap(i))
        kj = abs(kap(j))
        k = abs(ki-kj)
        !PRINT*, 'k, Gk'
        IF ((kap(i)*kap(j)).LT.0) k = k + 1
        kmi = ki+kj-1
        DO WHILE (k.LE.kmi)
           Gk = fkgk_int(i,j,i,j,k,dgc,dpc,rhoval_l,dr,ilast, alpha)
           rho0(:) = dgc(:,i)*dgc(:,i)+dpc(:,i)*dpc(:,i)
           Gk0 = fkgk_int(i,j,i,j,k,dgc,dpc,rho0,dr,ilast, alpha)
<<<<<<< HEAD
           WRITE(33,fmt = frmt) nnum(i), orb_ang_mom(li), tot_ang_mom((ji-1)/2), nnum(j), orb_ang_mom(lj), tot_ang_mom((jj-1)/2), 'G', k, Gk*hart, Gk0*hart, alpha
=======
           IF(Fk0.GT.0.d0) THEN
              WRITE(33,fmt = frmt) nnum(i), orb_ang_mom(li), tot_ang_mom((ji-1)/2), nnum(j), orb_ang_mom(lj), tot_ang_mom((jj-1)/2), 'G', k, Gk*hart, Gk0*hart, Gk/Gk0
           END IF
>>>>>>> 8094982 (Fixed slater condon printout)
           !PRINT*, k, Gk*hart
           k=k+2
        END DO
     END DO
  END DO
END SUBROUTINE fkgk

REAL(8) FUNCTION fkgk_int(i,j,l,m,k,dgc,dpc,rhoval_l,dr,ilast,alpha)
  !INPUT
  INTEGER i,j,l,m,ir1, ir2
  REAL(8) dgc(251,41), dpc(251,41), dr(251), dg1(251),dg2(251), rhoval_l(251)
  REAL(8) f(251), g(251), ri(251), zk(251), yk(251), alpha
  LOGICAL, PARAMETER :: UseProjection = .FALSE.

  alpha = 1.d0
  ! Double integral of phi_i(r)*phi_j(r)*r<**k/r>**(k+1)phi_l(r')\phi_m(r')
  fkgk_int = 0.d0
  ! Check that orbital indices are reasonable
  IF(i.LE.0.OR.j.le.0.OR.l.LE.0.OR.m.LE.0) THEN
     PRINT*, "In fkgk: Orbitals indices must be positive:", i, j, l, m
     STOP
  END IF
  IF(i.EQ.8.OR.i.EQ.9) THEN ! Only works for 3d states.
     !Integrate density to get total charge.
     dg1(1:251) = MAX(rhoval_l(1:251),0.d0)
     CALL trap(dr(1:ilast),dg1(1:ilast),ilast,yk(1))
     dg1(1:251) = sqrt(max(rhoval_l(1:251),0.d0)/yk(1))
  ELSE
     dg1(1:251) = sqrt(dgc(1:251,i)*dgc(1:251,i) + dpc(1:251,i)*dpc(1:251,i))
  END IF
  
  IF(j.EQ.8.OR.j.EQ.9) THEN
     !Integrate density to get total charge.
     dg2(1:251) = MAX(rhoval_l(1:251),0.d0)
     CALL trap(dr(1:ilast),dg2(1:ilast),ilast,yk(1))
     dg2(1:251) = sqrt(max(rhoval_l(1:251),0.d0)/yk(1))
  ELSE
     dg2(1:251) = sqrt(dgc(1:251,j)*dgc(1:251,j) + dpc(1:251,j)*dpc(1:251,j))
  END IF

  
  !f(1:251) = (dgc(1:251,i)*dgc(1:251,j) + dpc(1:251,i)*dpc(1:251,j))*dr(1:251)**k
  f(1:251) = dg1(1:251)*dg2(1:251)*dr(1:251)**k
  !g(1:251) = (dgc(1:251,i)*dgc(1:251,j) + dpc(1:251,i)*dpc(1:251,j))*dr(1:251)**(-k-1)
  g(1:251) = dg1(1:251)*dg2(1:251)*dr(1:251)**(-k-1)
  yk = 0.d0
  DO ir1 = 1, ilast
     IF(ir1.EQ.1) THEN
        yk(ir1) = 0.d0
     ELSE
        !PRINT*, 'Calculating yk', ir1
        CALL trap(dr(1:ir1),f(1:ir1),ir1,yk(ir1))
        yk(ir1) = yk(ir1)*dr(ir1)**(-k)
     END IF
     IF(ir1.EQ.ilast) THEN
        zk(ir1) = 0.d0
     ELSE
        !PRINT*, 'Calculating zk', ir1
        CALL trap(dr(ir1:ilast),g(ir1:ilast),ilast-ir1+1,zk(ir1))
        zk(ir1) = zk(ir1)*dr(ir1)**(k+1)
     END IF
  END DO
  yk(:) = yk(:) + zk(:)

  IF(l.NE.9) THEN
  !IF(.TRUE.) THEN
     dg1(1:251) = sqrt(dgc(1:251,l)*dgc(1:251,l) + dpc(1:251,l)*dpc(1:251,l))
  ELSE
     !Integrate density to get total charge.
     dg1(1:251) = MAX(rhoval_l(1:251),0.d0)
     CALL trap(dr(1:ilast),dg1(1:ilast),ilast,f(1))
     dg1(1:251) = sqrt(max(rhoval_l(1:251),0.d0)/f(1))
     IF(UseProjection) THEN
        ! Project onto atomic state
        dg2(1:251) = dg2(1:251)*sqrt(dgc(1:251,l)*dgc(1:251,l) + dpc(1:251,l)*dpc(1:251,l))
        CALL trap(dr(1:ilast),dg2(1:ilast),ilast,alpha)
        dg2(1:251) = sqrt((dgc(1:251,l)*dgc(1:251,l) + dpc(1:251,l)*dpc(1:251,l))*alpha**2)
     END IF
  END IF
  
  IF(m.NE.9) THEN
  !IF(.TRUE.) THEN
     dg2(1:251) = sqrt(dgc(1:251,m)*dgc(1:251,m) + dpc(1:251,m)*dpc(1:251,m))
  ELSE
     !Integrate density to get total charge.
     dg2(1:251) = MAX(rhoval_l(1:251),0.d0)
     CALL trap(dr(1:ilast),dg2(1:ilast),ilast,f(1))
     dg2(1:251) = sqrt(max(rhoval_l(1:251),0.d0)/f(1))
     IF(UseProjection) THEN
        ! Project onto atomic state
        dg2(1:251) = dg2(1:251)*sqrt(dgc(1:251,m)*dgc(1:251,m) + dpc(1:251,m)*dpc(1:251,m))
        CALL trap(dr(1:ilast),dg2(1:ilast),ilast,alpha)
        dg2(1:251) = sqrt((dgc(1:251,m)*dgc(1:251,m) + dpc(1:251,m)*dpc(1:251,m))*alpha**2)
     END IF
  END IF

  !dg1(1:251) = sqrt(dgc(1:251,l)*dgc(1:251,l) + dpc(1:251,l)*dpc(1:251,l))
  !dg2(1:251) = sqrt(dgc(1:251,m)*dgc(1:251,m) + dpc(1:251,m)*dpc(1:251,m))
  !f(1:251) = (dgc(1:251,l)*dgc(1:251,m) + dpc(1:251,l)*dpc(1:251,m))*yk(:)/dr(:)
  f(1:251) = dg1(1:251)*dg2(1:251)*yk(1:251)/dr(1:251)
  CALL trap(dr,f,ilast,fkgk_int)
  return
end FUNCTION fkgk_int

