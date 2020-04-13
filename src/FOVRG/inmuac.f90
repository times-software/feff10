!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: inmuac.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inmuac (ihole, xionin, iunf, ikap, iph) !KJ 2-2011 added iph

      use dimsmod, only: nrptx
      implicit double precision (a-h,o-z)
	  integer,intent(in) :: iph !KJ
      common/dff/cg(nrptx,30),cp(nrptx,30),bg(10,30),bp(10,30),fl(30), fix(30), ibgp
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
! the meaning of common variables is described below
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30), nq(30),kap(30),nmax(30)
! en one-electron energies
! scc factors for acceleration of convergence
! scw precisions of wave functions
! sce precisions of one-electron energies
! nmax number of tabulation points for orbitals
      common/scrhf1/eps(435),nre(30),ipl
! eps non diagonal lagrange parameters
! nre distingue: - the shell is closed (nre <0)
!                  the shell is open (nre>0)
!                - the orbitals in the integral rk if abs(nre) > or =2
! ipl define the existence of lagrange parameters (ipl>0)
      common/snoyac/dvn(nrptx),anoy(10),nuc
! dvn nuclear potential
! anoy development coefficients at the origin of nuclear potential
! this development is supposed to be written anoy(i)*r**(i-1)
! nuc index of nuclear radius (nuc=1 for point charge)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension xnval(30), iorb(-4:3)
      data nucm/11/

      testy=10.**(-5)
! testy precision for the wave functions

      call getorb (nz, ihole, xionin, iunf, norb, norbsc, iorb, iholep, nq, kap, xnel, xnval, en, iph) !KJ 2-2011 added iph
!     don't need xmag here, so use en as a dummy

      ipl=0
      do 40 i=1,norb
         en(i) = 0.d0
         nre(i)=-1
         llq= abs(kap(i))
         l=llq+llq
!       find last tabulation point
         nmax(i)=0
         do 100  j = idim, 1, -1
            if ( abs(cg(j,i)) .ge. 1.0d-11 .or. abs(cp(j,i)) .ge. 1.0d-11 )  then
               nmax(i) = j
               goto 16
            endif
  100    continue
   16    continue

         scc(i)=0.3
         if (xnel(i) .lt. l)  nre(i)=1
         if (ikap.eq.kap(i)) ipl=ipl+1
  40  continue
      norbsc=norb
      norb = norb+1
      xnel(norb)=1
      kap(norb)=ikap
      nq(norb) =9
! nz atomic number     noi ionicity (nz-number of electrons)
! norb number of orbitals
! xnel(i) number of electrons on orbital i.
      nuc=nucm
! nuc number of points inside nucleus (11 by default)

      return
      end
