!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: inipot.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inipot (dgc, dpc, edenvl, vvalgs, xnmues)
!     initialize values of arrays to zero
        use DimsMod, only: nphx=> nphu, lx
        
      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
      parameter (zero=0.0d0)

      dimension dgc(251,41,0:nphx+1), dpc(251,41,0:nphx+1)
      dimension edenvl(251,0:nphx), vvalgs (251,0:nphx)

      real*8, intent(inout) :: xnmues(0:lx,0:nphx)

      do 10 iph  = 0,nphx+1
      do 10 iorb = 1,41
      do 10 i = 1,251
   10    dgc(i,iorb,iph) = zero

      do 20 iph  = 0,nphx+1
      do 20 iorb = 1,41
      do 20 i = 1,251
   20    dpc(i,iorb,iph) = zero

      do 30 iph = 0, nphx
      do 30 i = 1, 251
   30    edenvl(i, iph) = zero

      do 40 iph = 0, nphx
      do 40 i = 1, 251
   40    vvalgs(i, iph) = zero

      do 50 iph = 0, nphx
      do 50 ll = 0, lx
   50    xnmues (ll, iph) = zero

      return
      end
