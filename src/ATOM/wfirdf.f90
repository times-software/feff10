!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: wfirdf.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine wfirdf (en,ch,nq,kap,nmax,ido)
!     calculate initial orbiatls from integration of dirac equation
! cg (cp) large (small) radial components
! bg (bp) development coefficients at the origin of cg (cp)
! en one-electron energies 
! fl power of the first term of development at the origin
! ch ionicity (nuclear charge - number of electrons)
! nq principal quantum number
! kap quantum number "kappa"
! nmax number of tabulation points for the orbitals
! ibgp first dimension of the arrays bg and bp
!        this programmes utilises nucdev,dentfa,soldir et messer
 
      implicit double precision (a-h,o-z)
      common cg(251,30), cp(251,30), bg(10,30), bp(10,30),              &
     &         fl(30), fix(30), ibgp
      dimension en(30),nq(30),kap(30),nmax(30)
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),                &
     &dv(251),av(10),eg(251),ceg(10),ep(251),cep(10)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/inelma/nem
      common/messag/dlabpr,numerr
      character*8 dlabpr
      character*512 slog
      common/snoyau/dvn(251),anoy(10),nuc
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
!#mn
       external dentfa

      cl=1.370373d+02
!    speed of light in atomic units
      dz = nz
! make r-mesh and calculate nuclear potential
! hx exponential step
! dr1 first tabulation point multiplied by nz
      hx=5.0d-02
      dr1= nz*exp(-8.8)
      call nucdev (anoy,dr,dvn,dz,hx,nuc,idim,ndor,dr1)
! notice that here nuc=1, 
! unless you specify nuclear mass and thickness in nucdev.f

      a=(dz/cl)**2
      if (nuc.gt.1) a=0.0d00
      do 11 j=1,norb
         b=kap(j)*kap(j)-a
         fl(j)= sqrt(b)
!        quick fix of development coefficients. ala
 11      fix(j) = dr(1)**(fl(j)-abs(kap(j)))
! calculate potential from thomas-fermi model
      do 21 i=1,idim
 21   dv(i)=(dentfa(dr(i),dz,ch)+dvn(i))/cl
      if (numerr.ne.0) return
      do 51 i=1,idim
         eg(i)=0.0d00
 51      ep(i)=0.0d00
      do 61 i=1,ibgp
         ceg(i)=0.0d00
         cep(i)=0.0d00
 61      av(i)=anoy(i)/cl
      av(2)=av(2)+dentfa(dr(nuc),dz,ch)/cl
      test1=testy/rap(1)
      b=test1

! resolution of the dirac equation to get initial orbitals
      if (ido.ne.1) then
         call wlog('only option ido=1 left')
         ido = 1
      endif
!  here was a piece to read orbitals from cards
      do 281 j=1,norb
         bg(1,j)=1.0d00
         i=nq(j)- abs(kap(j))
         if (kap(j).lt.0) i=i-1
         if (mod(i,2).eq.0) bg(1,j)=-bg(1,j)
         if (kap(j).lt.0) go to 201
         bp(1,j)=bg(1,j)*cl*(kap(j)+fl(j))/dz
         if (nuc.gt.1) bg(1,j)=0.0d00
         go to 211

 201     bp(1,j)=bg(1,j)*dz/(cl*(kap(j)-fl(j)))
         if (nuc.gt.1) bp(1,j)=0.0d00
 211     np=idim
         en(j)=-dz*dz/nq(j)*nq(j)
         method=0
         ifail = 0
         call soldir                                                    &
     &     (en(j),fl(j),bg(1,j),bp(1,j),b,nq(j),kap(j),nmax(j),ifail)

         if (numerr.eq.0) go to 251
         call messer
         write(slog,'(a,2i3)')                                          &
     &   'soldir failed in wfirdf for orbital nq,kappa ',nq(j),kap(j)
         call wlog(slog)
         go to 281

 251     do 261 i=1,ibgp
            bg(i,j)=ag(i)
 261        bp(i,j)=ap(i)
         do 271 i=1,np
            cg(i,j)=dg(i)
 271        cp(i,j)=dp(i)
 281  continue
      nem=0
      return
      end
