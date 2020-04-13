!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: writematrix.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************************************************

      subroutine writematrix(a,n,s)
      
      implicit none
      integer n,i,j
      complex a(n,n)
      character*(*) s
      
      open(93,file=s,form='formatted',status='unknown')
      do i=1,n
        write(93,2) (a(i,j),j=1,n)
      enddo
      close(93)
1      format(2000(f14.8,1x,f14.8,3x))      
2      format(2000(f10.4,1x,f10.4,2x))      
      
      return
      end
!**********************************************************************************************************

      subroutine writematrixfree(a,n,s)
      
      implicit none
      integer n,i,j
      complex a(n,n)
      character*(*) s
      
      open(93,file=s,form='formatted',status='unknown')
      do i=1,n
        write(93,2) (a(i,j),j=1,n)
      enddo
      close(93)
1      format(2000(e14.8,1x,e14.8,3x))      
2      format(2000(e10.4,1x,e10.4,2x))      
      
      return
      end
!**********************************************************************************************************

      subroutine writematrixdble(a,n,s)
      
      implicit none
      integer n,i,j
      complex*16 a(n,n)
      character*(*) s
      
      open(93,file=s,form='formatted',status='unknown')
      do i=1,n
        write(93,2) (a(i,j),j=1,n)
      enddo
      close(93)
1      format(2000(f14.8,1x,f14.8,3x))      
2      format(2000(f10.4,1x,f10.4,2x))      
      
      return
      end
!**********************************************************************************************************

      subroutine writematrixdblefree(a,n,s)
      
      implicit none
      integer n,i,j
      complex*16 a(n,n)
      character*(*) s
      
      open(93,file=s,form='formatted',status='unknown')
      do i=1,n
        write(93,2) (a(i,j),j=1,n)
      enddo
      close(93)
1      format(2000(e14.8,1x,e14.8,3x))      
2      format(2000(e10.4,1x,e10.4,2x))      
      
      return
      end
