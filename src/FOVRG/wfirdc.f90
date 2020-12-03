!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: wfirdc.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wfirdc (eph,kap,nmax,vxc,ps,qs,aps,aqs,irr,ic3,vm,     &
     &                   rmt,jri, iwkb) 
!     calculate photoelectron orbital using lda in dirac equation
!     cg (cp) large (small) radial components
!     bg (bp) development coefficients at the origin of cg (cp)
!     eph one-electron energy of photoelectron
!     fl power of the first term of development at the origin
!     kap quantum number "kappa"
!     nmax number of tabulation points for the orbitals
!     vxc  is initial lda potential for photoelectron
!     ibgp first dimension of the arrays bg and bp
!        this programmes utilises nucdec,dentfa,soldir et messer
      use dimsmod, only: nrptx
      implicit double precision (a-h,o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      common/dff/cg(nrptx,41),cp(nrptx,41),bg(10,41),bp(10,41),         &
     &             fl(41), fix(41), ibgp
      dimension kap(41),nmax(41)
!    for photoelectron potential and wavefunction will be complex
      complex*16 eph,dg,ag,dp,ap,dv,av,eg,ceg,ep,cep,vxc(nrptx)
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10),vm(nrptx)
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),            &
     &dv(nrptx),av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/messag/dlabpr,numerr
      character*8 dlabpr
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
 
! Debug: Fer
! Check the input:
!     print *, 'eph:  ', eph 
!     print *, 'kap:  ', kap
!     print *, 'nmax: ', nmax
!     print *, 'vxc:  ', vxc
!     print *, 'ps:   ', ps
!     print *, 'qs:   ', qs
!     print *, 'aps:  ', aps
!     print *, 'aqs:  ', aqs
!     print *, 'irr:  ', irr
!     print *, 'ic3:  ', ic3
!     print *, 'vm:   ', vm
!     print *, 'rmt:  ', rmt
!     print *, 'jri:  ', jri
!     print *, 'iwkb: ', iwkb

      cl=1.370373d+02
!     speed of light in atomic units
      dz = nz
!     make r-mesh and calculate nuclear potential
!     hx exponential step
!     dr1 first tabulation point multiplied by nz
      dr1= nz*exp(-8.8)
      call nucdec (anoy,dr,dvn,dz,hx,nuc,idim,10,dr1)
!     notice that here nuc=1, 
!     unless you specify nuclear mass and thickness in nucdec.f


      a=(dz/cl)**2
      if (nuc.gt.1) a=0.0d00
      do 11 j=1,norb
         b=kap(j)*kap(j)-a
         if (j.eq.norb) b=b+(kap(j)+1)*ic3
         fl(j)= sqrt(b)
 11      fix(j) = dr(1)**(fl(j)-abs(kap(j)))
!     if irregular solution
      if (irr.gt.0) then
         fl(norb) = -fl(norb)
         fix(norb) = 1.0/fix(norb)
      endif

!     use lda potential to calculate initial w.f.
      do 21 i=1,jri-1
 21   dv(i)= vxc(i)/cl
      do  i=jri,idim
        dv(i)= vxc(jri+1)/cl
      enddo
      if (numerr.ne.0) return
      do 51 i=1,idim
         eg(i)=0.0d00
 51      ep(i)=0.0d00
      do 61 i=1,ibgp
         ceg(i)=0.0d00
 61      cep(i)=0.0d00
      call potdvp
      av(2)=av(2)+(vxc(nuc)-dvn(nuc))/cl

!     resolution of the dirac equation to get initial orbital
      if (irr.lt.0) then
         if (a .gt. 0.0d0) then 
            aps(1) = 1.0
            if (kap(norb) .lt. 0) then
               aqs(1)=aps(1)*dz/(cl*(kap(norb)-fl(norb)))
            else
               aqs(1)=aps(1)*cl*(kap(norb)+fl(norb))/dz
            endif
         else
            if (kap(norb).lt.0)then
               aps(1)=1.0d00
               aqs(1)=0.0d00
            else
               aps(1)=0.0d00
               aqs(1)=1.0d00
            endif
         endif
      endif

 211  np=1+(8.8 + log(10.0))/hx
!     exp(-8.8+(np-1)*hx) = 10.0 bohrs - max distance
      if (idim .lt. np) np=idim
      if (nmax(norb) .gt. np) nmax(norb)=np
         
      if (irr.lt.0) then
         call solout( eph, fl(norb), aps(1), aqs(1), kap(norb), rmt,    &
     &              jri, nmax(norb), ic3, vm, iwkb)
      else
         call solin( eph, fl(norb), aps(1), aqs(1), kap(norb), rmt,     &
     &              jri, nmax(norb), ic3, vm, iwkb)
      endif
         
      do 261 i=1,10
         aps(i)=ag(i)
 261     aqs(i)=ap(i)
      do 271 i=1,idim
         ps(i)=dg(i)
 271     qs(i)=dp(i)
      return
      end
