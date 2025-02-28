!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: inmuat.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2012/09/11 22:52:14 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inmuat (ihole, xionin, iunf, xnval, iholep, xmag, iorb, iph, nq2)  !KJ 12-2010 added iph, nq2
      implicit double precision (a-h,o-z)
! Josh Kas - Changed array dimensions from 30 to 41 for high Z elements
! according to Pavlo Baranov's changes.
      dimension xnval(41), xmag(41), iorb(-5:4)
	  integer nq2(41)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
! the meaning of common variables is described below
      common/ratom1/xnel(41),en(41),scc(41),scw(41),sce(41),nq(41),kap(41),nmax(41)
! en one-electron energies
! scc factors for acceleration of convergence
! scw precisions of wave functions
! sce precisions of one-electron energies
! nmax number of tabulation points for orbitals
      ! JK - 435 below may need to be changed to 820 for high Z
      ! elements.
      common/scrhf1/eps(820),nre(41),ipl
! eps non diagonal lagrange parameters
! nre distingue: - the shell is closed (nre <0)
!                  the shell is open (nre>0)
!                - the orbitals in the integral rk if abs(nre) > or =2
! ipl define the existence of lagrange parameters (ipl>0)
      common/snoyau/dvn(251),anoy(10),nuc
! dvn nuclear potential
! anoy development coefficients at the origin of nuclear potential
! this development is supposed to be written anoy(i)*r**(i-1)
! nuc index of nuclear radius (nuc=1 for point charge)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      ! JK - 435 below may need to be changed to 820 for high Z
      ! elements.
      data nucm/11/,nesn/50/,ideps/820/ ! JK - changing nucm to 5 to see what happens.

      ndor=10

      ! testy precision for the wave functions
      testy=1.0d-05

      ! teste precision for the one-electron energies
      teste=5.0d-06

      ! rap tests of precision for soldir
      rap(1)=100.
      rap(2)=10.

      do i = 1, 41
         en(i) = 0.d0
         xmag(i) = 0
         xnval(i) = 0
      end do

!write(*,*) '*********** Calling getorb for :'
      call getorb (nz, ihole, xionin, iunf, norb, norbsc, iorb, iholep, nq, kap, xnel, xnval, xmag, iph)  !KJ 12-2010 added iph
	  nq2(:)=nq(:)
!      write(*,*) 'nz,ihole,xionin',nz,ihole,xionin
!	  write(*,*) 'iunf,norb,norbsc',iunf,norb,norbsc
!	  write(*,*) 'iorb',iorb
!	  write(*,*) 'iholep,nq',iholep,nq
!	  write(*,*) 'kap',kap
!	  write(*,*) 'xnel',xnel
!	  write(*,*) 'xnval',xnval
!	  write(*,*) 'xmag,iph',xmag,iph
!     stop
	  

      xk=0
      do i=1,norb
         xk=xk+xnel(i)
      end do
!write(*,*) 'norb,xnel',norb,xnel(1:norb)
!write(*,*) 'nz,xionin,xk',nz,xionin,xk
      !if ( abs(nz-xionin-xk) .gt. 0.001) call par_stop('check number of electrons in getorb.f')

      norbsc=norb
! nz atomic number     noi ionicity (nz-number of electrons)
! norb number of orbitals
! xnel(i) number of electrons on orbital i.
! first norbsc orbitals will be determined selfconsistently,
! the rest of orbitals are orthogonolized if iorth is non null,
! and their energies are those on cards if iene is non null
! or otherwise are the values obtained from solving dirac equation
      nes=nesn
! nes number of attempts in program soldir
      nuc=nucm
! nuc number of points inside nucleus (11 by default)
      do 171 i=1,ideps
 171  eps(i)=0.0d00

      idim = 251
      if (mod(idim,2) .eq. 0) idim=idim-1

      ipl=0
! if ipl non null, it permits a repartition of tabulation points and certain precision tests.
      do 401 i=1,norb
         nre(i)=-1
         llq= abs(kap(i))
         l=llq+llq
         if (kap(i).lt.0) llq=llq-1
         if (llq.lt.0.or.llq.ge.nq(i).or.llq.gt.4) call par_stop('kappa out of range, check getorb.f')
         nmax(i)=idim
         scc(i)=0.3
         if (xnel(i) .lt. l)  nre(i)=1
         if (xnel(i) .lt. 0.5)  scc(i)=1.0
         do 385 j=1,i-1
            if (kap(j).ne.kap(i)) go to 385
            if (nre(j).gt.0.or.nre(i).gt.0) ipl=ipl+1
 385     continue
 401  continue
 999  return
      end
