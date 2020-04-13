!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ellfun.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/06/29 01:05:24 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function ellpi(x)
      implicit double precision (a-h, o-z)

      en=-x
!     sign differs in Num.Rec. and Byrd/Friedman
      ellpi= rf(0.d0, 0.5d0, 1.d0) - en/3.d0*rj(0.d0,0.5d0,1.d0,1.d0+en)

      return
      end
      function rf(x,y,z)
      implicit double precision (a-h, o-z)
      parameter (errtol=0.0025, tiny=1.5e-38, big=3.e37, third=1./3.,   &
     &  c1=1./24., c2=0.1, c3=3./44., c4=1./14. )

      if (min(x,y,z).lt.0. .or. min(x+y,x+z,y+z).lt.tiny .or.           &
     &   max(x,y,z).gt.big) stop 'invalid arguments in rf'  !pause to stop  KJ 6-2012

      xt=x
      yt=y
      zt=z
  1   continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=third*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if (max(abs(delx), abs(dely), abs(delz)).gt.errtol) goto 1
      e2 = delx*dely-delz**2
      e3 = delx*dely*delz
      rf=(1+(c1*e2-c2-c3*e3)*e2+c4*e3)/sqrt(ave)
      
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .

      function rj(x,y,z,p)
      implicit double precision (a-h, o-z)
      parameter (errtol=0.0015, tiny=2.5e-13, big=9.e11, c1=3./14.,     &
     &  c2=1./3., c3=3./22., c4=3./26., c5=.75*c3, c6=1.5*c4, c7=.5*c2, &
     &  c8=c3+c3)

      if (min(x,y,z).lt.0. .or. min(x+y,x+z,y+z,abs(p)).lt.tiny .or.    &
     &   max(x,y,z,abs(p)).gt.big) stop 'invalid arguments in rj'   !pause to stop  KJ 6-2012

      sum=0
      fac=1
      if (p.gt.0) then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt= min(x,y,z)
        zt= max(x,y,z)
        yt = x+y+z -xt-zt
        a=1./(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif

  1   continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum = sum + fac*rc(alpha,beta)
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        pt=.25*(pt+alamb)
        ave=.2*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if (max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.errtol) goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3*ec
      ee=eb+2*delp*(ea-ec)
      rj=3*sum+fac*(1+ed*(-c1+c5*ed-c6*ee)+eb*(c7+delp*(-c8+delp*c4))   &
     &   +delp*ea*(c2-delp*c3)-c2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0) rj=a*(b*rj+3*(rcx-rf(xt,yt,zt)))
      
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .

      function rc(x,y)
      implicit double precision (a-h, o-z)
      parameter (errtol=0.0012, tiny=1.69e-38, sqrtny=1.3e-19,big=3.e37,&
     &  tnbg=tiny*big, comp1=2.236/sqrtny, comp2=tnbg*tnbg/25.,         &
     &  third=1./3., c1=.3, c2=1./3., c3=.375, c4=9./22.)

      if (x.lt.0 .or. y.eq.0 .or. (x+abs(y)).lt.tiny .or.               &
     &  (x+abs(y)).gt.big .or. (y.lt.-comp1.and.x.gt.0..and.x.lt.comp2))&
     &  stop 'invalid argument in rc'    !pause to stop  KJ 6-2012
      if (y.gt.0) then
        xt=x
        yt=y
        w=1.
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
  1   continue
        alamb=2.*sqrt(xt)*sqrt(yt)+yt
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        ave=third*(xt+yt+yt)
        s=(yt-ave)/ave
      if (abs(s).gt.errtol) goto 1
      rc = w*(1+s*s*(c1+s*(c2+s*(c3+s*c4)))) /sqrt(ave)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .
