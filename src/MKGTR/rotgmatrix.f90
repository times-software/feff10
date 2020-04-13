      subroutine rotgmatrix(nq,elpty,pha,beta,qweights,nsp,nspx,lx,gg,ggrot)  
      implicit none

      integer nq
      double precision elpty
      double precision qcst,qsnt,qcsf,qsnf !KJ,qweights
	  complex*16 qweights(nq)
      integer nsp,nspx,lx
      complex gg(nspx*(lx+1)**2, nspx*(lx+1)**2)
      complex ggrot(nspx*(lx+1)**2, nspx*(lx+1)**2)

      double precision rot
      complex rotm(-lx:lx,-lx:lx,0:lx) 
      integer ixx1,ixx2
      integer ixtot
      complex pham,pha
      double precision beta
      integer l1,m1,l2,m2
      integer mp1,mp2
      integer is1,is2
      integer ig1,ig2
      integer igp1,igp2
      double precision  rotwig
      external rotwig
      

      ixtot=nspx*(lx+1)**2

!     if nq.eq.0 .or. elpty.le.0 then just copy the original 
!     g-matrix and weight it. Any other case rotate it. 
      

      if (nq.eq.0 .or. elpty.lt.0.0d0) then 
!         write(6,*) "doing spherical in rotgmatrix"
         do ixx1=1,ixtot
            do ixx2=1,ixtot
               ggrot(ixx2,ixx1)=gg(ixx2,ixx1)
            end do
         end do
      else

!     use 
!     alpha=phi, beta=theta, gamma=0

!         pha=cmplx(qcsf,qsnf)
!         pha=conjg(pha)
!         beta=atan2(qcst,qsnt)

!     zero out and fill out the rotation matrix

		 rotm(:,:,:)=dcmplx(0.0d0,0.0d0)

         do l1=0,lx
            do m1=-l1,l1
               do m2=-l1,l1
!                 the inverse rotation
                  pham=conjg(pha)**m2
                  rot=rotwig(-beta,l1,m2,m1,1)
! oops
!                 rotm(m2,m1,l1)=pham*rot*sqrt(qweights)
                  rotm(m2,m1,l1)=pham*rot
               end do
            end do
         end do


         do is1=1,nsp
            do is2=1,nsp
               do l1=0,lx
                  do m1=-l1,l1
                     ig1=nsp*(l1*l1+l1)+nsp*m1+is1
                     do l2=0,lx
                        do m2=-l2,l2
                           ig2=nsp*(l2*l2+l2)+nsp*m2+is2
                           ggrot(ig2,ig1)=cmplx(0.0d0,0.0d0)
                           do mp1=-l1,l1
                              igp1=nsp*(l1*l1+l1)+nsp*mp1+is1
						      do mp2=-l2,l2
                                 igp2=nsp*(l2*l2+l2)+nsp*mp2+is2
								 ggrot(ig2,ig1)=ggrot(ig2,ig1)+rotm(m1,mp1,l1)*gg(igp1,igp2)* conjg(rotm(m2,mp2,l2))
                              end do                          
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         
      end if


!KJ : I moved the weight out of here, because I want to do weight(q) * weight(q') in MDFF calculations
      return
	  
!!     add the weight
!       do ixx1=1,ixtot
!        do ixx2=1,ixtot
!           ggrot(ixx2,ixx1)=ggrot(ixx2,ixx1)*qweights
!        end do
!      end do
!
!      return
      end
