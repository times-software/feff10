!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: dsordc.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      complex*16 function dsordc(j,a,dg,dp,ag,ap)
!              * calculation of overlap integrals*
!        integration by simpson method of the   hg*(r**0)
!        hg(l)=dg(l)*cg(l,j)+dp(l)*cp(l,j)
!                cg,cp(l,j)  orbital j
!        a is such that dg,dp or hg following the case
!        behave at the origin as cte*r**a
!        the development limits at the origin (used for calculation
!        of integral form 0 to dr(1) ) of functions dg,dp and hg are
!        supposed to be in blocks ag,ap and chg respectively
!        this program uses   aprdec
      use dimsmod, only: nrptx
      implicit double precision (a-h,o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      complex*16 aprdec
      common/dff/ cg(nrptx,41), cp(nrptx,41), bg(10,41), bp(10,41),     &
     &              fl(41), fix(41), ibgp
      complex*16 dg(nrptx),ag(10),dp(nrptx),ap(10)
      complex*16 hg(nrptx),chg(10)
!     common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),
!    1   nq(41),kap(41),nmax(41)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)

!        construction of the array hg
      do  15 l= 1,ibgp
        bgj(l) = bg(l,j)
 15     bpj(l) = bp(l,j)

      do 221 l=1,idim
 221  hg(l)=dg(l)*cg(l,j)+dp(l)*cp(l,j)
      b=a+fl(j)
      do 241 l=1,ndor
 241     chg(l) = aprdec(ag,bgj,l) + aprdec(ap,bpj,l)
 
!        integration of the hg
      dsordc = (0.0d0, 0.0d0)
      do 305 l=1,idim
 305     hg(l)=hg(l)*dr(l)
      do 311 l=2,idim,2
 311     dsordc=dsordc+hg(l)+hg(l)+hg(l+1)
      dsordc=hx*(dsordc+dsordc+hg(1)-hg(idim))/3.0d0
!        integral from 0 to dr(1)
      do 331 l=1,ndor
         b=b+1.0d00
 331     dsordc=dsordc+chg(l)*(dr(1)**b)/b
      return
      end
