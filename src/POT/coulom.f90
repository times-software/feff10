!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: coulom.f90,v $:
! $Revision: 1.6 $
! $Author: hebhop $
! $Date: 2012/11/29 23:20:18 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine coulom( icoul, npot, ilast, rhoval, edenvl, edens,     &
     &     nat, rat, iatph, iphat, rnrm, dq, iz, vclap)
!     searches for fermi level in comlex energy plane 
!     Output:
!       rhoval - new valence density !KJ I think vclap is the only output!
!       vclap  - coulomb potential
!       qnrm   - charge inside each norman sphere
      use DimsMod, only: nphx=>nphu, natx
      use constants
      implicit double precision (a-h, o-z)

      dimension ilast(0:nphx)
      dimension ri05(251)
      dimension rhoval(251,0:nphx+1), edenvl(251,0:nphx)
      dimension edens(251,0:nphx), dq(0:nphx), iz(0:nphx)
      dimension rat(3,natx), iatph(0:nphx), iphat(natx), rnrm(0:nphx)
      dimension vclap(251,0:nphx)

!     work space
      dimension  drho(251), dvcl(251)
      external dist

!     make  radial grid with 0.05 step
      dx05=0.05d0
      do 10 i=1,251
         ri05(i) = exp(-8.8+dx05*(i-1))
  10  continue
  
      do 600 ip=0,npot
        do 550 ir=1, ilast(ip)
           drho(ir)= (rhoval(ir,ip)-edenvl(ir,ip))*ri05(ir)**2
  550   continue
        call potslw(dvcl,drho, ri05,dx05, ilast(ip))

        do 560 ir = ilast(ip)+1, 251
           dvcl(ir) = 0.0d0
  560   continue

        if (icoul.eq. 1) then
!         find the change of coulomb potential at norman radius for
!         each type of iph
          jnrm = (log(rnrm(ip)) + 8.8) / 0.05  +  2
          dvnrm = dq(ip) / rnrm(ip)
          iat0 = iatph(ip)
          do 570 iat=1,nat
             if (iat.ne.iat0) then
               rr = dist( rat(1,iat), rat(1,iat0))
               if (rr.lt.rnrm(ip)) rr=rnrm(ip)
               dvnrm = dvnrm + dq(iphat(iat)) / rr
             endif
  570     continue

!         transfer condition to r(jnrm) instead of r_nrm.
          dr = ri05(jnrm) - rnrm(ip)
!         xx = dr/rnrm(ip)
!         correction using linear expansion of drho
!         neglecting terms xx**4 and higher
          bb = (drho(jnrm)-drho(jnrm-1)) / (ri05(jnrm)-ri05(jnrm-1))
!         dvnrm = dvnrm - xx* (dq(ip)/ri05(jnrm) + xx* (drho(jnrm)*
!    1    (rnrm(ip)/ri05(jnrm)-0.5) + xx*(drho(jnrm)-bb*rnrm(ip))/3 ))
          dvnrm = dvnrm - dr / 2 * ( dq(ip) / rnrm(ip)**2 +             &
     &       (dq(ip)+drho(jnrm)*dr-bb/2*dr**2) / ri05(jnrm)**2 )


!         dvcl is calculated correct up to constant shift which is
!         fixed by the condition at R_nrm
          dvnrm = dvnrm - dvcl(jnrm)

        else
!         now this is default (icoul=0)
!         then do normalization based on norman picture
!         i.e. total density is approximated by a sum of densities
!          which are zero outside each norman sphere. use this
!         approximation only for the difference between 2 potentials
!         This is needed for infinite solid where the algorithm for
!         icoul=1 is unstable due to long range Coulomb potential
!         probably better fix will be to use Ewald summation to figure
!         out the Madelung constants (icoul=2 option to be done later).
          
          call frnrm (edens(1,ip), iz(ip), rnrm1)
          do 710 i = 1,251
  710        drho(i) = edens(i,ip) - edenvl (i,ip) +rhoval(i,ip)
          call frnrm (drho, iz(ip), rnrm2)
          rmin = min (rnrm1, rnrm2)
          inrm = (log(rmin) + 8.8) / 0.05  +  1
          r0 = ri05(inrm)

          delv = 0.d0
          if (rnrm2.gt.rnrm1) then
            aa = (drho(inrm+1)-drho(inrm)) / (ri05(inrm+1)-ri05(inrm))
            bb = drho(inrm) - aa * ri05(inrm)
            delv = delv - fab (aa, bb, r0, rnrm1, rnrm2)
          else
            aa = (edens(inrm,ip)-edens(inrm+1,ip))  / (ri05(inrm+1)-ri05(inrm))
            bb = - edens(inrm,ip) - aa * ri05(inrm)
            delv = delv - fab (aa, bb, r0, rnrm2, rnrm1)
          endif
          aa = (drho(inrm+1)-drho(inrm)+edens(inrm,ip)-edens(inrm+1,ip))  /  (ri05(inrm+1)-ri05(inrm))
          bb = drho(inrm) - edens(inrm,ip) - aa * ri05(inrm)
          delv = delv - fab (aa, bb, r0, r0, rmin)

          dvnrm = delv - dvcl(inrm)
        endif

        do 580 ir=1,ilast(ip)
           vclap(ir,ip) = vclap(ir,ip) + dvcl(ir) + dvnrm 
  580   continue
        do 590 ir=ilast(ip)+1,251
  590   vclap(ir,ip)=0.0d0
  600 continue

      return
      end

      double precision function fab (aa,bb,r0,r1,r2)
!     it is the \int_r1^r2 dr 4\pi\rho(r) r**2 (1/r0 - 1/r)
!     where 4\pi\rho(r) = aa*r + bb
!     you arrive to this integral as a result of norman picture
!     for normalization of coulomb potential just below the rmin
      implicit double precision (a-h, o-z)

      a2 = (r2**2-r1**2)/2.d0
      a3 = (r2**3-r1**3)/3.d0
      a4 = (r2**4-r1**4)/4.d0
      fab = aa*(a4/r0-a3) + bb*(a3/r0-a2)
      return
      end

