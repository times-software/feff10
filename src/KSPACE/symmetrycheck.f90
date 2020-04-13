!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: symmetrycheck.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2012/01/30 22:04:08 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine symmetrycheck(ierr,ntrans,mtrx,tnp,mult)

      use constants,only: pi2
      implicit none
      integer,intent(out) :: ierr
        integer,intent(in) :: mtrx(48,3,3),ntrans
      real*8,intent(in) :: tnp(48,3)

      integer mult(48,48),nerr(48),mtest(3,3),i,j,k,l,m,itest,maxerr
      real*8 ttest
        integer yes(ntrans),nspacegroup
      logical equiv
        real*8,parameter :: prec=0.00000001d0


      ierr=0
      if (ntrans .eq. 1) return

!  check for duplicate operations
      do i=2,ntrans
        do j=1,i-1
          equiv = .true.
          do k=1,3
            if (dabs(tnp(i,k)-tnp(j,k)).gt.prec) equiv= .false.
            do l=1,3
              if (mtrx(i,k,l).ne.mtrx(j,k,l)) equiv=.false.
            end do
          end do
          if ( equiv ) then
            write ( 6, '(/,1a11,2x,2i5)' ) ' equiv ops ', i, j
            write ( 6, * ) 'symmetry check fault'
          end if
        end do
      end do

!  make mul table
!    Check that the product of two sym ops is again a sym op
      do i=1,ntrans
        nerr(i) = 0
        do j=1,ntrans
          mult(i,j) = 0
            mtest=0
!  multiply i and j
          do k=1,3
            do l=1,3
              do m=1,3
                mtest(k,l) = mtest(k,l) + mtrx(i,k,m)*mtrx(j,m,l)
              end do
            end do
          end do
!  check for match
          do k=1,ntrans
            equiv = .true.
            l = 1
            m = 1
            do while ( equiv .and. ( l .le. 3 ) )
              if ( mtest( l, m ) .ne. mtrx( k, l, m ) ) equiv = .false.
              if ( m .lt. 3 ) then
                m = m + 1
              else
                m = 1
                l = l + 1
              end if
            end do
            if ( equiv ) mult( i, j ) = k
          end do
          if (mult(i,j).eq.0) stop 'no deal -- not a group'
        end do
      end do

!  if translations not correct set mult(i,j) to -1
      do i=1,ntrans
        do j=1,ntrans
          k = mult(i,j)
          l = 1
          do while ( ( mult( i, j ) .ne. -1 ) .and. ( l .lt. 4 ) )
            ttest = tnp(j,l)
            do m=1,3
              ttest = ttest + dble(mtrx(i,m,l))*(tnp(i,m)-tnp(k,m))
            end do
            ttest = dabs(ttest)/pi2
            itest = ttest * 1.001
            if (dabs(ttest-dble(itest)) .ge. 0.0001d0) mult(i,j) = -1
            l = l + 1
          end do
        end do
      end do


      yes=0
      nspacegroup=0
        do i=1,ntrans
!        Check multiplication table
         equiv=.true.
           do j=1,ntrans
              if(mult(i,j).le.0) equiv=.false.
              if(mult(j,i).le.0) equiv=.false.
         enddo
           if(equiv) then
              nspacegroup=nspacegroup+1
              yes(i)=1
           endif
        enddo

!        write(6,*) 'nspacegroup',nspacegroup
!        write(6,*) 'yes ',yes



!  check multiplication table
      do i=1,ntrans
      do j=1,ntrans
      if (mult(i,j) .le. 0) then
        nerr(i) = nerr(i) + 1
        nerr(j) = nerr(j) + 1
      end if
      end do
      end do
!  find element with max error
      ierr   = 0
      maxerr = 0
      do i=1,ntrans
        if (nerr(i) .gt. maxerr) then
          maxerr = nerr(i)
          ierr = i
        end if
      end do
!      write(6,'(21h1multiplication table,/)')
!
!      do i=1,ntrans
!        write(6,'(1x,48i3)') (mult(i,j),j=1,ntrans)
!      end do

!      write(6,*) 'Ierr = ',ierr
!        write(6,*) 'maxerr = ',maxerr


      return
      end







