!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: vlda.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine vlda(ia, xnval,srho, srhovl,vtrho, ilast, idfock)
!    this program calculates xc potential, using core-vlaence separation
!    discussed in ankuodinov's thesis.  
!    written by alexei ankoudinov. 11.07.96
      use constants
      implicit double precision (a-h,o-z)
! Josh Kas - Changed array dimensions from 30 to 41 (and others) for high Z elements
! according to Pavlo Baranov's changes.

      dimension xnval(41), srho (251), srhovl(251), vtrho(251)
      common cg(251,41), cp(251,41), bg(10,41), bp(10,41),              &
     &        fl(41), fix(41), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),dv(251),av(10), &
     &              eg(251),ceg(10),ep(251),cep(10)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
 
      do 10 i = 1,251
        srhovl(i) = 0.0d0
 10     srho(i) = 0.0d0
 
!  find total and valence densities. Remove self-interaction if SIC
      do 50 j = 1, norb
        a = xnel(j)
        b = xnval(j)
!      use to test SIC
!       if (j .eq. ia) a=a-1.0d0
!       if (j .eq. ia) b=b-1.0d0
      do 50 i = 1,nmax(j)
        srho(i) = srho(i) + a * (cg(i,j)**2+cp(i,j)**2)
 50     srhovl(i) = srhovl(i) + b * (cg(i,j)**2+cp(i,j)**2)

!  constract lda potential. Put your favorite model into vbh.f.
!  exch=5,6 correspond to 2 ways of core-valence separation of V_xc.
      do 90 i = 1,251
        rho = srho(i) / (dr(i)**2)
        if (idfock.eq.5) then
!          for exch=5 valence density*4*pi
           rhoc = srhovl(i) / (dr(i)**2)
        elseif (idfock.eq.6) then
!          for exch=6 core density*4*pi
           rhoc = (srho(i)-srhovl(i)) / (dr(i)**2)
        elseif (idfock.eq.1) then
           rhoc = 0.0d0
        elseif (idfock.eq.2) then
           rhoc = srho(i) / (dr(i)**2)
        else
            call par_stop(' undefined idfock in subroutine vlda')
        endif

        if (rho .gt. 0.0 ) then
           rs = (rho/3)**(-third)
           rsc =101.0
           if (rhoc .gt.0.0) rsc = (rhoc/3)**(-third)
           xm = 1.0d0
!          vbh and edp in Hartrees
           if (idfock.eq.5 .or. idfock.eq.2) then
!             for exch=5, 2
              call vbh(rsc, xm, vxcvl)
           elseif (idfock.eq.6) then
!             for exch=6
              call vbh(rs, xm, vvbh)
                 xf = fa/rs
              call edp(rsc,xf,vdh)
              vxcvl = vvbh - vdh
           elseif (idfock.eq.1) then
!          for pure Dirac-Fock
              vxcvl = 0.0d0
           endif

!   contribution to the total energy from V_xc:=\int d^3 r V_xc * rho/2
           if (ilast.gt.0) vtrho (i) = vtrho(i) +                       &
     &         vxcvl * srho(i)
!    1         vxcvl * xnel(ia)*(cg(i,ia)**2+cp(i,ia)**2)
!           use to test SIC

!  add to the total potential and correct it's development coefficients
           if (i.eq.1) av(2) = av(2) +vxcvl/cl
           dv(i) = dv(i) +vxcvl/cl
        endif
 90   continue
 999  continue

      return
      end
