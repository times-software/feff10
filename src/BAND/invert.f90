      subroutine invert(ndim,nadp,smat,suse,sinv)
      implicit none
      integer i,j,k,ii,ndim,nadp
      double precision smat(nadp,nadp),suse(nadp,nadp),sinv(nadp,nadp)
      double precision ratio,swap,rc1,rc2
!
!  copy, set up identity
!
      sinv=dble(0)
      do i=1,ndim
        do j=1,ndim
          suse(j,i)=smat(j,i)
        end do
        sinv(i,i)=1
      end do

!  do inversion by pivoted Gaussian elimination

      do i=1,ndim

        ii=i
        do j=i+1,ndim
          rc1=dabs(suse(j,i))
          rc2=dabs(suse(ii,i))
          if (rc1.gt.rc2) ii=j
        end do
        if (ii.gt.i) then
          do j=i,ndim
            swap=suse(i,j)
            suse(i,j)=suse(ii,j)
            suse(ii,j)=swap
          end do
          do j=1,ndim
            swap=sinv(i,j)
            sinv(i,j)=sinv(ii,j)
            sinv(ii,j)=swap
          end do
        end if
        if (suse(i,i).eq.dcmplx(0.d0,0.d0)) then
          write (6,*) 'ZERO DETERMINANT...'
          write (98,*) 'ZERO DETERMINANT...'
          stop
        end if
        do j=1,ndim
          if (j.ne.i) then
            ratio=-suse(j,i)/suse(i,i)
          else
            ratio=dcmplx(1.d0,0.d0)/suse(i,i)-dcmplx(1.d0,0.d0)
          endif
          do k=i,ndim
            suse(j,k)=suse(j,k)+ratio*suse(i,k)
          end do
          do k=1,ndim
            sinv(j,k)=sinv(j,k)+ratio*sinv(i,k)
          end do
        end do

      end do

      return
      end
