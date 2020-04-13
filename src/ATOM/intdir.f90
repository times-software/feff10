!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: intdir.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine intdir(gg,gp,ag,ap,ggmat,gpmat,en,fl,agi,api,ainf,max0)
!            solution of the inhomogenios dirac equation
! gg gp initially exch. terms, at the time of return are wave functions
! ag and ap development coefficients of  gg and gp
! ggmat gpmat  values at the matching point for the inward integration
! en one-electron energy
! fl power of the first development term at the origin
! agi (api) initial values of the first development coefficients
! at the origin of a large (small) component
! ainf initial value for large component at point dr(max0) 
!   - at the end of tabulation of gg gp
 
      implicit double precision (a-h,o-z)
      common/comdir/cl,dz,bid1(522),dv(251),av(10),bid2(522)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      common/subdir/ell,fk,ccl,imm,nd,node,mat
      common/messag/dlabpr,numerr
      character*8 dlabpr
      dimension gg(251),gp(251),ag(10),ap(10),coc(5),cop(5),dg(5),dp(5)
      save
      data cop/2.51d+02,-1.274d+03,2.616d+03,-2.774d+03,1.901d+03/,     &
     &coc/-1.9d+01,1.06d+02,-2.64d+02,6.46d+02,2.51d+02/,               &
     &cmixn/4.73d+02/,cmixd/5.02d+02/,hxd/7.2d+02/,npi/5/,icall/0/
 
! numerical method is a 5-point predictor-corrector method
! predicted value    p(n) = y(n-1) + c * somme de i=1,5 cop(i)*y'(n-i)
! corrected value    c(n) = y(n-1) + c * somme de i=1,4 coc(i)*y'(n-i)
!                                  + coc(5)*p'(n)
! final value        y(n) = cmix*c(n) + (1.-cmix)*p(n)
!                           cmix=cmixn/cmixd
 
      if (icall.eq.0) then
         icall=1
         c=cmixn/cmixd
         a=1.0d00-c
         cmc=c*coc(5)
         f=coc(1)
         do 1 j=2,npi
            g=coc(j)
            coc(j)=c*f+a*cop(j)
 1          f=g
         coc(1)=c*cop(1)
      endif
      c=hx/hxd
      ec=en/cl
      ag(1)=agi
      ap(1)=api
      if (imm) 81,15,26
!      search for the second sign change point
 15   mat=npi
      j=1
 16   mat=mat+2
         if (mat.ge.np) then
!   i had trouble with screened k-hole for la, for f-electrons.
!   below i still define matching point if one electron energy not less
!   than -1ev.
            if (ec .gt. -0.0003) then
              mat = np - 12
              go to 25
            endif
            numerr=56011
!          * fail to find matching point
            return
         endif
         f=dv(mat)+ell/(dr(mat)*dr(mat))
         f=(f-ec)*j
         if (f) 25,25,16
 25      j=-j
      if (j.lt.0) go to 16
      if (mat .ge. np-npi) mat=np-12
 
!     initial values for the outward integration
 26   do 35 j=2,ndor
         k=j-1
         a=fl+fk+k
         b=fl-fk+k
         ep=a*b+av(1)*av(1)
         f=(ec+ccl)*ap(k)+ap(j)
         g=ec*ag(k)+ag(j)
         do 31 i=1,k
            f=f-av(i+1)*ap(j-i)
 31         g=g-av(i+1)*ag(j-i)
 
         ag(j)=(b*f+av(1)*g)/ep
 35      ap(j)=(av(1)*f-a*g)/ep
      do 41 i=1,npi
         gg(i)=0.0d00
         gp(i)=0.0d00
         dg(i)=0.0d00
         dp(i)=0.0d00
         do 41 j=1,ndor
            a=fl+j-1
            b=dr(i)**a
            a=a*b*c
            gg(i)=gg(i)+b*ag(j)
            gp(i)=gp(i)+b*ap(j)
            dg(i)=dg(i)+a*ag(j)
 41         dp(i)=dp(i)+a*ap(j)
      i=npi
      k=1
      ggmat=gg(mat)
      gpmat=gp(mat)
 
!     integration of the inhomogenious system
 51   cmcc=cmc*c

 55   continue
         a=gg(i)+dg(1)*cop(1)
         b=gp(i)+dp(1)*cop(1)
         i=i+k
         ep=gp(i)
         eg=gg(i)
         gg(i)=a-dg(1)*coc(1)
         gp(i)=b-dp(1)*coc(1)
         do 61 j=2,npi
            a=a+dg(j)*cop(j)
            b=b+dp(j)*cop(j)
            gg(i)=gg(i)+dg(j)*coc(j)
            gp(i)=gp(i)+dp(j)*coc(j)
            dg(j-1)=dg(j)
 61         dp(j-1)=dp(j)
         f=(ec-dv(i))*dr(i)
         g=f+ccl*dr(i)
         gg(i)=gg(i)+cmcc*(g*b-fk*a+ep)
         gp(i)=gp(i)+cmcc*(fk*b-f*a-eg)
         dg(npi)=c*(g*gp(i)-fk*gg(i)+ep)
         dp(npi)=c*(fk*gp(i)-f*gg(i)-eg)
      if (i.ne.mat) go to 55

      if (k.lt.0) go to 999
      a=ggmat
      ggmat=gg(mat)
      gg(mat)=a
      a=gpmat
      gpmat=gp(mat)
      gp(mat)=a
      if (imm.ne.0) go to 81
 
!     initial values for inward integration
      a=test1* abs(ggmat)
      if (ainf.gt.a) ainf=a
      max0=np+2
 73   a=7.0d+02/cl
 75   max0=max0-2
         if ((max0+1).le.(mat+npi)) then
            numerr=138021
!          *the last tabulation point is too close to the matching point
            return
         endif
      if (((dv(max0)-ec)*dr(max0)*dr(max0)).gt.a) go to 75

 81   c=-c
      a=- sqrt(-ec*(ccl+ec))
      if ((a*dr(max0)).lt.-1.7d+02) go to 73
      b=a/(ccl+ec)
      f=ainf/ exp(a*dr(max0))
      if (f.eq.0.0d00) f=1.0d00
      do 91 i=1,npi
         j=max0+1-i
         gg(j)=f* exp(a*dr(j))
         gp(j)=b*gg(j)
         dg(i)=a*dr(j)*gg(j)*c
 91      dp(i)=b*dg(i)
      i=max0-npi+1
      k=-1
      go to 51

 999  return
      end
