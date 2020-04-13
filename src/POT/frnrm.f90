!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: frnrm.f90,v $:
! $Revision: 1.7 $
! $Author: hebhop $
! $Date: 2012/11/29 23:20:18 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine frnrm (rho, iz, rnrm)
      use dimsmod, only: nrptx
      implicit double precision (a-h, o-z)

      !KJ dimension rho(nrptx)
      dimension rho(251) !KJ 11-2009: frnrm only called from coulom.f90 for 'rho' arguments of size 251 whereas
      ! nrptx is currently set to 1251 in m_dimsmod.f90 ...
      dimension xpc(251), ri(251)
      character*256 slog
      external rr
      
!     finds norman radius

!     Need overlapped densities.  We'll get them in the form
!     4*pi*density = rho.  Also need z of atom

!     Then integrate out to the point where the integral of
!     4*pi*density*r**2 is equal to iz
      sum= (9*rho(1)*rr(1)**3+28*rho(2)*rr(2)**3+23*rho(3)*rr(3)**3)/480
!     add initial point (r=0) correction (see subroutine somm2)
      dpas = 0.05
      d1 = 3.0
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=rr(1)/db
      dd=rr(1)*(1.0+1.0/(dd*(d1+1.0)))/d1
      sum = sum + dd*rho(1)*rr(1)**2 - db*rho(2)*rr(2)**2

      fl = rho(4) *rr(4)**3
      fr = rho(5) *rr(5)**3
      frr = rho(6) *rr(6)**3
      sum = sum + (25*fl + 12 *fr -frr)/480
      do 10  i = 7, nrptx
         fll = fl
         fl = fr
         fr = frr
		 if (i.le.251) then  !KJ to avoid array bounds being overstepped ...  3-2012
            frr = rho(i) * rr(i)**3
		 else
		    frr=0.d0
		 endif
         sumsav = sum
         sum = sum + (13*(fr+fl) -fll -frr)/480
         if (sum .ge. iz)  then
            inrm = i-2
            x= (iz-sumsav)/(sum-sumsav)
            goto 20
         endif
   10 continue
      call wlog(' FRNRM Could not integrate enough charge to reach required z.')
      write(slog,*) 'Atom type Z= ',iz,' charge found= ',sum,' r=',rr(251)
      call wlog(slog)
      call wlog('Cannot determine Norman radius -- quitting.')
      call par_stop('FRNRM-1')
   20 continue
      rnrm = rr(inrm)*(1 + x*0.05)
     
!     add next order correction ALA 3/97
        dx05 = 0.05
        x0 = 8.8
        jnrm =  (log(rnrm) + x0) / dx05  +  2
        i0=jnrm+1
        xirf = 2
        do 710 ir = 1, jnrm+2
           ri(ir) = rr(ir)
           xpc(ir) = rho(ir)*ri(ir)**2
  710   continue

        call somm2 (ri, xpc, dx05, xirf, rnrm,0,i0)
!       dq is how many new electrons are within norman sphere
        dn1 = xirf-iz
        x2 = x - dn1/((1-x)*xpc(inrm) + x*xpc(inrm+1))
        if (abs(x2-x).gt.0.0001) then
          xirf = 2
          rnrm = rr(inrm)*(1 + x2*0.05)
          call somm2 (ri, xpc, dx05, xirf, rnrm,0,i0)
          dn2 = xirf-iz
!         Newton-Raphson methof to find zeroes
          x = x2 - dn2 * (x2-x)/(dn2-dn1)
        endif
        rnrm = rr(inrm)*(1 + x*0.05)

      return
      end
