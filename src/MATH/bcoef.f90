!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: bcoef.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks,     &
     &                 kiind, lind, bmat)
!     written by alexei ankudinov; march 2000
!     calculate bmat: the energy independent sum over polarization and
!     angular momenta indices
!     bmat = \sum_{p,p', all m_j} <LS|J><J|R|J1><J1|\alpha_p exp(i kz)|I>
!                    ptz(p,p') 
!            <I|\alpha_p'^* exp(-i kz) J2><J2|R'|J'><J'|L'S'>
!     where R is rotation from spin vector to x-ray k-vector
!     and R' is rotation back
!     see Eq.10 and 11 in Ankudinov,Rehr, Phys.Rev.B (accepted),
!     Theory of solid state contribution to the x-ray elastic scattering
!     aditional rotation matrices are needed when x-ray k-vector
!     is not along the spin-axis (see rotations in rdinp)

!     more precisely it is
!     bmat(l1 l1' j l ml ms; l2 l2' j' l' ml' ms') =
!        (-)**(j-j'+l2'+1) i**(l'-l) \sum_{p,p',mi,m1,mj,m2,mj'}
!        <LS|J>   r^j_{m1,mj}(angks)   3j( j l1 i ; -m1 p mi)
!        (-p)**(l1+l1'+1) ptz(p,p') (-p')**(l2+l2'+1) 
!        3j( j' l2 i ; -m2  p' mi)   r^j'_{m2,mj'}(angks)   <J'|L'S'>
!     where l1 l1' are set by the multipole moment(E1-l1=1,l1'=0;
!     E2-l1=2,l1'=1; M1-l1=1,l1'=1; etc.;
!     j and l define quantum number kappa and for each multipole moment
!     Only few final kappa's are allowed and  it is convinient
!     to denote (l1 l1' j l) by one index 'k'
!     thus  k=1-8 to include both E1 and E2 transitions;
!     ml and ms are projections of orbital and spin moments.

!     bmat  is used to calculate absorption fine structure (chi) via
!       chi = \sum_{k ms,k' ms'}  rkk(k,ms)  rkk(k',ms')
!       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms')  G_(l' ml' ms';l ml ms)
!     where sum over spins can be moved from first sum to second for
!     spin independent systems. The above expression is suitable for FMS
!     and for MS expansion on can use Eq.15 in RA paper to obtain
!     expression for the termination   matrix
!     T_{lam1 ms,lamN ms'} = \sum_{k k'} rkk(k,ms) rkk(k',ms')
!       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms') gam(l,lam1,rho1,ms)
!        gamtl(l',lamN,rhoN,ms')
!     Notice that for spin-dependent systems the scattering F matrices
!     in RA paper also should have additional spin indices. In genfmt
!     we currently neglect spin-flip processes which simplifies
!     calculations with MS expansion. (T and F are diagonal in ms,ms')
       
!     This subroutine is written for general spin-dependent asymmetric
!     system and arbitrary polarization tenzor. The symmetry of the 
!     system and polarization tenzor can be used
!     to speed up FMS or MS calculations in appropriate subroutines.
!     (see comments in subroutines mpprmp, fmstot)

!     input:
!       kinit - kappa for initial orbital
!       ipol - polarization type measurement
!       ptz  - polarization tensor (needed only for ipol=1 case)
!       le2  - sets which multipole moments to include (see mkptz)
!       ltrace- .true. for xsect.f, where need to perform trace over ml
!       angks - angle between k-vector and spin-vector 

!     output
!       lind  - orb.mom.(kappa)  needed in fmstot only (for indexing)
!       bmat  - energy independent matrix to calculate absorption 
!       in many cases bmat is diagonal due to the choice of xyz frame,
!       but for general case full 16*(2*lx+1)*16*(2*lx+1) matrix is kept

      use DimsMod ,only: lx,nspx=>nspu

      complex*16 coni
      parameter (coni = (0,1))

!     need only parameter lx to set max orb momentum
!KJ      complex*16 ptz, bmat, pmat, tmat
      complex*16 ptz(-1:1,-1:1),  bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
!       to include all possible dipole and quadrupole transitions 
!       final kp, and kpp have 8 possibilities
      logical ltrace

!     local staff
!KJ      dimension  t3j( 8, 0:1, -lx:lx+1), x3j(8, -1:1, -lx:lx+1)
      real*8  t3j( 8, 0:1, -lx:lx+1), x3j(8, -1:1, -lx:lx+1) !KJ
!     qmat = <J2|R'|J'><J'|L'S'> - diagonal in kappa index
      real*8 qmat( -lx:lx+1, -lx:lx, 0:1, 8)
!     pmat = <J1|\alpha_j exp(i kz)|I> ptz <I|\alpha_k^* exp(-i kz)|J2>
      complex*16 pmat( -lx:lx+1, 8, -lx:lx+1, 8)
!     tmat = pmat*qmat ; bmat = qmat^T*tmat
      complex*16 tmat( -lx:lx+1, 8, -lx:lx, 0:1, 8)
!     total and orbital momenta for 8 possible final kappa
!KJ      dimension jind(8), lind(8), kind(8)
      integer jind(8), lind(8), kiind(8) !KJ

      integer i1,i2,i3,i4,i5,i6,k,kap,kinit,jkap,le2,ipol,ms,ml,lkap,mp,k1,j1,j2,j3,m1,m2,i,mj,jj,mmj,mmp,j,is,ms1,ml1,ms2,ml2  !KJ
	  integer ispin
	  real*8 value,angks
      real*8, external :: cwig3j, rotwig


       bmat = dcmplx(0,0)

!     3 dipole transitions
      do 20 k=-1,1
         kap=kinit+k
         if (k.eq.0) kap=-kap
         jkap = abs(kap)
         lkap = kap
         if (kap.le.0) lkap = abs(kap) -1
!        check that orbital momentum does not exceed max allowed
         if (lkap .gt. lx) then
!          set final j and l to unphysical values
           jkap = 0
           lkap = -1 
           kap = 0
         endif
         jind(k+2) = jkap
         lind(k+2) = lkap
         kiind(k+2) = kap
  20  continue

!     include 5 quadrupole or 3 mag.dipole  transitions
      do 120 k=-2,2
         jkap = abs(kinit) + k
         if (jkap.le.0) jkap = 0
         kap= jkap
         if (kinit.lt.0 .and. abs(k).ne.1) kap=-jkap
         if (kinit.gt.0 .and. abs(k).eq.1) kap=-jkap
         lkap = kap
         if(kap.le.0) lkap = - kap - 1
         if (lkap.gt.lx .or. le2.eq.0                                   &
     &                  .or. (le2.eq.1 .and. abs(k).eq.2)) then
!           set unphysical jkap and lkap to make shorter calculations
            jkap = 0
            lkap = -1
            kap = 0
         endif
         jind(k+6) = jkap
         lind(k+6) = lkap
         kiind(k+6) = kap
 120  continue

      if (ipol.eq.0) then
!       polarization average case; bmat is diagonal and simple
        do 100 k = 1, 8
        do 100 ms = 0 ,1
        do 100 ml = -lind(k), lind(k)
!         i2 = (2*l1+1) , where l1 is defined by multipole moment
          i2 = 3
          if (le2.eq.2 .and. k.gt.3) i2 = 5
          bmat(ml,ms,k, ml,ms,k) = 0.5d0 / (2*lind(k)+1.d0) / i2
          if (k.le.3) bmat(ml,ms,k, ml,ms,k) = - bmat(ml,ms,k, ml,ms,k)
 100    continue
      else
!       more complicated bmat for linear(ipol=1) and circular(ipol=2)
!       polarizations
!       Put 3j factors in x3j and t3j. t3j are multiplied by
!       sqrt(2*j'+1) for  further convinience.
        do 30  mp=-lx,lx+1
        do 30  ms=0,1
        do 30  k1=1,8
  30    t3j(k1,ms,mp) = 0.0d0
        do 40  mp=-lx,lx+1
        do 40  ms=-1,1
        do 40  k1=1,8
  40      x3j(k1,ms,mp) = 0.0d0

        do 70  k1 = 1,8
        do 70  mp = -jind(k1)+1,jind(k1)
          do 50 ms=0,1
            j1 = 2 * lind(k1)
            j2 = 1
            j3 = 2 * jind(k1) - 1
            m1 = 2*(mp-ms)
            m2 = 2*ms - 1
            t3j(k1,ms,mp)=sqrt(j3+1.0d0) * cwig3j(j1,j2,j3,m1,m2,2)
            if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0)                        &
     &          t3j(k1,ms,mp) = - t3j(k1,ms,mp)
!           t3j(m0,i)    are Clebsch-Gordon coefficients
  50      continue
          do 60 i=-1,1
            j1 = 2 * jind(k1) - 1
            j2 = 2
            if (k1.gt.3 .and. le2.eq.2) j2 = 4
            j3 = 2 * abs(kinit) - 1
            m1 = -2*mp + 1
            m2 = 2*i
            x3j(k1,i,mp)= cwig3j(j1,j2,j3,m1,m2,2)
  60      continue
  70    continue

!       calculate qmat
        do 220 i=1,8
        do 220 ms=0,1
        do 220 ml= -lind(i), lind(i)
        do 220 mj= -jind(i)+1, jind(i)
          mp = ml+ms
          jj = 2*jind(i) - 1
          mmj = 2*mj - 1
          mmp = 2*mp - 1
          value = rotwig(angks, jj, mmj, mmp, 2)
          qmat(mj,ml,ms,i) = value * t3j(i,ms,mp)
 220    continue

!       calculate pmat
        do 240 i2 = 1,8
        do 240 m2 = -jind(i2)+1, jind(i2)
        do 240 i1 = 1,8
        do 240 m1 = -jind(i1)+1, jind(i1)
          pmat(m1,i1,m2,i2) = 0
          if (abs(m2-m1).le.2) then
            do 230 j=-1,1
            do 230 i=-1,1
!             check that initial moment is the same
              if (m1-i.eq.m2-j) then
                is = 1
!               (-p) factors for M1 transitions
                if (le2.eq.1 .and. i.gt.0 .and. i1.gt.3) is = -is
                if (le2.eq.1 .and. j.gt.0 .and. i2.gt.3) is = -is
                pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) +                 &
     &          is * x3j(i1,i,m1) * ptz(i,j) * x3j(i2,j,m2)
              endif
 230        continue
!           multiply by (-)^(j-j'+l2'+1) i**(l'-l) factor
!           additional (-) is from Eq.10 (-2*ck)
            is = 1
            if (mod(jind(i1)-jind(i2), 2) .ne.0) is = -is
            if (i2.le.3) is = -is
            pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) * is                  &
     &           * coni**(lind(i2)-lind(i1))
          endif
 240    continue

!       calculate tmat = pmat*qmat
        do 270 i1=1,8
        do 270 ms=0,1
        do 270 ml=-lind(i1), lind(i1)
        do 270 i2=1,8
        do 270 mj=-jind(i2)+1, jind(i2)
          tmat(mj,i2, ml,ms,i1) = 0
          do 260 mp = -jind(i1)+1, jind(i1)
            tmat(mj,i2, ml,ms,i1) = tmat(mj,i2, ml,ms,i1)+              &
     &           pmat(mj,i2,mp,i1) * qmat(mp,ml,ms,i1)
 260      continue
 270    continue
         
!       calculate bmat = qmat^T * tmat
        do 300 i1=1,8
        do 300 ms1=0,1
        do 300 ml1=-lind(i1), lind(i1)
        do 300 i2=1,8
        do 300 ms2=0,1
        do 300 ml2=-lind(i2), lind(i2)
          bmat(ml2,ms2,i2, ml1,ms1,i1) = 0
          do 280 mj=-jind(i2)+1, jind(i2)
            bmat(ml2,ms2,i2, ml1,ms1,i1) = bmat(ml2,ms2,i2, ml1,ms1,i1)+&
     &      qmat(mj,ml2,ms2,i2) * tmat(mj,i2,ml1,ms1,i1) 
 280      continue
 300    continue
!       end of ipol=1,2 cases
      endif 

      if (ltrace) then
!       need to trace bmat over ml for xsect.f
        do 390 i1 = 1, 8
        do 390 ms1 = 0,1
        do 390 i2 = 1, 8
        do 390 ms2 = 0,1
          if (lind(i1).ne.lind(i2) .or. ms1.ne.ms2) then
               bmat(0,ms2,i2, 0,ms1,i1) = 0
          else
             do 360 ml = 1, lind(i1)
               bmat(0,ms1,i2, 0,ms1,i1) =  bmat(0,ms1,i2, 0,ms1,i1) +   &
     &         bmat(-ml,ms1,i2, -ml,ms1,i1) + bmat(ml,ms1,i2, ml,ms1,i1)
 360         continue
          endif
 390    continue
      endif

      if (ispin .eq. 0) then
!       G(Ls,L's') is spin diagonal; trace over spin
        do 480 i1 = 1, 8
        do 480 i2 = 1, 8
        do 480 ml1 = -lind(i1), lind(i1)
        do 480 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,0,i2, ml1,0,i1) +      &
     &                                  bmat(ml2,1,i2, ml1,1,i1)
 480    continue
      elseif (ispin.eq.2 .or. (ispin.eq.1 .and. nspx.eq.1)) then
!       move spin up part into the position of spin-down
        do 490 i1 = 1, 8
        do 490 i2 = 1, 8
        do 490 ml1 = -lind(i1), lind(i1)
        do 490 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,1,i2, ml1,1,i1)
 490    continue

      endif

      return
      end
