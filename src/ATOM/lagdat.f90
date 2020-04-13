!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: lagdat.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lagdat (ia,iex)
!        * non diagonal lagrange parameteres *
! lagrange parameters involving orbital ia if ia is positive
! all lagrange parameters are calculated if ia is negative or zero
! contribution of the exchange terms is omitted if iex=0
!        this program uses akeato(bkeato) fdrirk multrk

      implicit double precision (a-h,o-z)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &nq(30),kap(30),nmax(30)
      common/scrhf1/eps(435),nre(30),ipl
      dimension ni(2),nj(2)
!#mn
       external akeato, bkeato, fdrirk
 
      i1= max(ia,1)
      idep=1
      if (ia.gt.0) go to 15
 11   idep=i1+1
 15   ni(1)=i1
      nj(2)=i1
      ji1=2* abs(kap(i1))-1
      do 201 i2=idep,norbsc
         if (i2.eq.i1.or.kap(i2).ne.kap(i1)) go to 201
         if (nre(i1).lt.0.and.nre(i2).lt.0) go to 201
! the following line was included to handle the case of 1 electron in
! 2 s-shells.
! Probably need to use schmidt orthogonalization in this case
         if (xnel(i1).eq.xnel(i2)) go to 201
         ni(2)=i2
         nj(1)=i2
         d=0.0d00
         do 101 l=1,norbsc
            k=0
            jjl=2* abs(kap(l))-1
            kma= min(ji1,jjl)
 41         a=akeato(l,i1,k)/xnel(i1)
            b=a-akeato(l,i2,k)/xnel(i2)
            c=b
            if (a.ne.0.0d00) c=c/a
            if ( abs(c).lt.1.0d-07) go to 51
            d=d+b*fdrirk(l,l,i1,i2,k)
 51         k=k+2
            if (k.le.kma) go to 41
            if (iex.eq.0) go to 101
            kma=(ji1+jjl)/2
            k= abs(jjl-kma)
            if ((kap(i1)*kap(l)).lt.0) k=k+1
 61         a=bkeato(l,i2,k)/xnel(i2)
            b=a-bkeato(l,i1,k)/xnel(i1)
            c=b
            if (a.ne.0.0d00) c=c/a
            if ( abs(c).lt.1.0d-07) go to 71
            d=d+b*fdrirk(i1,l,i2,l,k)
 71         k=k+2
            if (k.le.kma) go to 61
 101     continue
         i= min(i1,i2)
         j= max(i1,i2)
         eps(i+((j-1)*(j-2))/2)=d/(xnel(i2)-xnel(i1))
 201  continue
      if (ia.gt.0) go to 999
      i1=i1+1
      if (i1.lt.norbsc) go to 11
 999  return
      end
