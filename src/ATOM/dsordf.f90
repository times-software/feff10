!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: dsordf.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function dsordf (i,j,n,jnd,a)

         USE ErrorMod
!              * calculation of diff. integrals*
!        integration by simpson method of the   hg*(r**n)
!        hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)  if jnd=1
!        hg=expression above multiplied by  dg  if jnd=-1
!        hg(l)=cg(l,i)*cp(l,j)                  if jnd=2
!        hg=expression above multiplied by  dg  if jnd=-2
!        hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)      if jnd=3
!        hg(l)=dg(l)*dg(l)+dp(l)*dp(l)          if jnd=4
!        hg is constructed by calling program   if jnd>=5
!                  cg(l,i)  large component of the orbital i
!                  cp(l,j)  small component of the orbital j
!        a is such that dg,dp or hg following the case
!        behave at the origin as cte*r**a
!        the integration is made as far as dr(j) for jnd>3
!
!        the development limits at the origin (used for calculation
!        of integral from 0 to dr(1) ) of functions dg,dp and hg are
!        supposed to be in blocks ag,ap and chg respectively
!        this program uses  aprdev
!
! Josh Kas - Changed array dimensions from 30 to 41 for high Z elements
      implicit double precision (a-h,o-z)
      common cg(251,41), cp(251,41), bg(10,41), bp(10,41),              &
     &         fl(41), fix(41), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
      dimension hg(251),chg(10)
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),            &
     &nq(41),kap(41),nmax(41)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      dimension bgi(10),bgj(10),bpi(10),bpj(10)
      external aprdev

      ! construction of the array hg
      SELECT CASE(ABS(jnd))
         CASE(1) 

            max0= min(nmax(i),nmax(j))
            bgi(:) = bg(:,i)
            bgj(:) = bg(:,j)
            bpi(:) = bp(:,i)
            bpj(:) = bp(:,j)
            
            DO l=1,max0
               hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)
            END DO
         
            DO l=1,ndor
               chg(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)
            END DO

            b=fl(i)+fl(j)

            IF(jnd.lt.0) THEN
               do l=1,max0
                  hg(l)=hg(l)*dg(l)
               end do
               
               do l=1,ndor
                  ap(l)=chg(l)
               end do
               
               b=b+a
               do l=1,ndor
                  chg(l)=aprdev(ap,ag,l)
               end do
            END IF

         CASE(2) ! hg(l)= cg(l,i)*cp(l,j)
            max0= min(nmax(i),nmax(j))
            bgi(:) = bg(:,i)
            bgj(:) = bg(:,j)
            bpi(:) = bp(:,i)
            bpj(:) = bp(:,j)

            do l=1,max0
               hg(l)=cg(l,i)*cp(l,j)
            end do
            
            do l=1,ndor
               chg(l)=aprdev(bgi,bpj,l)
            end do
         
            b=fl(i)+fl(j)

            IF(jnd.lt.0) THEN
               do l=1,max0
                  hg(l)=hg(l)*dg(l)
               end do
               
               do l=1,ndor
                  ap(l)=chg(l)
               end do
               
               b=b+a
               do l=1,ndor
                  chg(l)=aprdev(ap,ag,l)
               end do
            END IF

         CASE(3)
            max0= min(nmax(i),nmax(j))
            bgi(:) = bg(:,i)
            bgj(:) = bg(:,j)
            bpi(:) = bp(:,i)
            bpj(:) = bp(:,j)

            do l=1,max0
               hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)
            end do

            b=a+fl(i)
            do l=1,ndor
               chg(l)=aprdev(bgi,ag,l)+aprdev(bpj,ap,l)
            end do
         CASE(4)
            max0=j
            b=a

            do l=1,max0
               hg(l)=dg(l)*dg(l)+dp(l)*dp(l)
            end do
            
            b=b+b
            do l=1,ndor
               chg(l)=aprdev(ag,ag,l)+aprdev(ap,ap,l)
            end do

         CASE (5:)
            max0=j
            b=a

         CASE DEFAULT
            CALL Error('Illegal input to dsordf. jnl = ' // ACHAR(jnl), StopProgram = .FALSE.)
            CALL Error('jnl must be {-2, -1, 1, 2, 3, 4, 5, ...}')
      END SELECT

      dsordf=0.0d00
      io=n+1
      do l=1,max0
         hg(l)=hg(l)*(dr(l)**io)
      end do

      do l=2,max0,2
         dsordf=dsordf+hg(l)+hg(l)+hg(l+1)
      end do

      dsordf=hx*(dsordf+dsordf+hg(1)-hg(max0))/3.0d00

!        integral from 0 to dr(1)
      b=b+n

      do l=1,ndor
         b=b+1.0d00
         dsordf=dsordf+chg(l)*(dr(1)**b)/b
      end do
      return
      end
