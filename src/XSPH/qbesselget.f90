!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: qbesselget.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2010/12/17 01:12:33 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     simple subroutine for getting the 
!     j_l(q*x) for inelastic x-ray scattering 
!     matrix elements
!
!

      subroutine qbesselget(qtrans,nr,ljmax,ri,jbesq)
      implicit none 
      double precision qtrans
      integer nr
      integer ljmax
      double precision ri(nr)
      double precision jbesq(nr,0:ljmax)
      
      double complex nl(1:(ljmax+2)),jl(1:(ljmax+2))
      double complex qr

      integer ir,lj
      
      do ir=1,nr
         qr=dcmplx(qtrans*ri(ir),0.d0)

         if (dble(qr).lt.1.e+8) then 
            call besjnjas(qr,jl,nl,ljmax)
            do lj=0,ljmax
               jbesq(ir,lj)=dble(jl(lj+1))
            end do
         else
            do lj=0,ljmax
               jbesq(ir,lj)=0.0d0
            end do
         end if
!         write(6,*) ir,nr,qr,ri(ir),jbesq(ir,0)
      end do

      return
      end 
