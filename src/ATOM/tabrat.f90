!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: tabrat.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine tabrat
!      tabulation of the results
! do identifications of orbitals
! nmax number of tabulation points for wave function
!      this programm uses dsordf

      implicit double precision (a-h,o-z)
      parameter (zero = 0)
      parameter (ryd  = 13.605698d0)
      parameter (hart = 2*ryd)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &nq(30),kap(30),nmax(30)
      character*2  titre(30)
      character*2  ttire(9)
      dimension at(8),mbi(8)
      logical open_16

      external dsordf
      data ttire /'s ', 'p*', 'p ', 'd*', 'd ', 'f*', 'f ','g*', 'g '/
 
      do 110 i=1,norb
         if (kap(i) .gt. 0) then
           j=2*kap(i)
         else
           j=-2*kap(i)-1
         endif
         titre(i)=ttire(j)
 110  continue

!     tabulation of number of points and of average values of
!                   r**n (n=6,4,2,1,-1,-2,-3)
      inquire(unit=16,opened=open_16)
      do 201 i=2,8
 201     mbi(i)=8-i-i/3-i/4+i/8
      if(open_16) then
        write(16,*)                                                     &
     &  'number of electrons nel and average values of r**n in a.u.'
        write(16,2061) (mbi(k),k=2,8)
 2061   format (5x,'nel     -E ','     n=',7(i2,8x))
      endif
      do 251 i=1,norb
         llq= abs(kap(i))-1
         j=8
         if (llq.le.0) j=7
         do 241 k=2,j
 241        at(k)=dsordf(i,i,mbi(k),1, zero)
 251     if (open_16)                                                   &
     &     write(16,2071) nq(i),titre(i),xnel(i),-en(i)*hart,           &
     &                  (at(k),k=2,j)
 2071 format(i1,a2,f6.3,8(1pe10.3))

!      overlap integrals
      if (norb.le.1) return
      if (open_16) write(16,321)
 321  format(10x,'overlap integrals')
      do 351 i=1,norb-1
         do 331 j=i+1,norb
            if (kap(j).ne.kap(i)) go to 331
            at(1)=dsordf(i,j,0,1, zero)
            if(open_16)                                                 &
     &        write(16,2091)  nq(i),titre(i),nq(j),titre(j),at(1)
 331     continue
 351  continue
 2091 format (4x,i3,a2,i3,a2,f14.7)
      return
      end
