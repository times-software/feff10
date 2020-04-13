!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: dmscf.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dmscf(em, ne, nex, matsize, chi0, dipmat, xkmat,dipscf)

      implicit double precision (a-h, o-z)

      parameter (maxsize = 78)
      complex*16 coni
      parameter (coni = (0,1))
      dimension em(nex)
      complex*16 chi0(nex,maxsize,maxsize), dipscf(nex,maxsize)
      dimension dipmat(nex,maxsize)
      complex*16 xkmat(nex,maxsize,maxsize)
      dimension ipiv(matsize)
      complex amat(matsize,matsize)
      real xxr, xxim
      complex amatinv(matsize,matsize)
      complex rhs(matsize), aminv(matsize,matsize)
      character*1 trans

      do 200 ie = 1, ne
       do 10 i = 1, matsize
       do 10 j = 1, matsize
         amat(i,j) =0

!        Mscf = 1 - K*chi0
         if (i .eq. j) then
           delta = 1
         else 
           delta = 0
         endif
    
!        chi0 is diagonal, so it has only one index j
!        K  = 1 - K*chi0   
         xxr = real(delta)
         xxim = 0
         do 5 k=1,matsize
           xxr = xxr - real( dble(xkmat(ie,i,k) * chi0(ie,k,j)))
           xxim=xxim - real(dimag(xkmat(ie,i,k) * chi0(ie,k,j)))
   5     continue
         amat(i,j)=cmplx(xxr, xxim) 
   10  continue

       do 50 i = 1, matsize
         dipscf(ie,i) = 0
         do 40 j = 1, matsize
!          ddr =  dble(real(amatinv(i,j))) * dipmat(ie,j)
!          ddim = dble(aimag(amatinv(i,j))) * dipmat(ie,j)
           dipscf(ie,i) = dipscf(ie,i) + amatinv(i,j)*dipmat(ie,j)
           amatinv(i,j) = amat(i,j)
   40    continue
   50  continue
      
!      find screened matrix elements by precise matrix inversion 
       call cgetrf(matsize,matsize,amat,matsize,ipiv,info) 
       if (info .lt. 0) call wlog('  *** Error in cgetrf')
  
       do 60 i = 1, matsize
         rhs(i) = dipmat(ie,i)
   60  continue
  
       nrhs = 1
       trans = 'N'
       call cgetrs(trans,matsize,nrhs,amat,matsize,ipiv,                &
     &              rhs,matsize,info)   
   
       if (info .lt. 0) call wlog('  *** Error in cgetrc')

       do 70 i = 1, matsize
         dipscf(ie,i) = rhs(i)
   70  continue

  200 continue
      
      return
      end
