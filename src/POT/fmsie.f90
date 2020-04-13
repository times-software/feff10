!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmsie.f90,v $:
! $Revision: 1.11 $
! $Author: jorissen $
! $Date: 2012/02/04 00:38:51 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fmsie( iph0, nph, lipotx, ie, em, eref, ph, iz, rfms, lfms, nat, iphat, rath, gtr,lverb)

!     full multiple scattering code for single energy point
!     written by a.ankudinov 06.1997 using earlier written subroutines
!     coded by b.ravel
!     modified by a.ankudinov 2001 for new matrix inversion algorithms
!     Feb. 2002, a.ankudinov: fixed logic for MPI calculations
!       lfms=0  - extended system calculations (e.g. crystal)
!       lfms=1  - small system calculations (e.g. molecule)
!       lfms=2  - same as 1 for MPI run (forces call yprep)

!KJ for k-space code
      use controls,only : ispace
      use DimsMod, only: natx, nphx=> nphu, lx, nspx=>nspu

      implicit double precision (a-h, o-z)

!     input
      dimension iphat(natx), rath(3,natx)
      real rat(3,natx), rfms, rdirec, toler1, toler2
      real rpart,aipart
      integer nph
      dimension iz(0:nphx)
      complex*16, intent(in) :: ph(lx+1, 0:nphx)
      complex, intent(inout) :: gtr(0:lx, 0:nphx)
      logical,intent(in) :: lverb

!     work space
      integer iph0
      complex*16 em, eref
      character*512 slog

!     fms staff
      integer lipotx(0:nphx)
      complex*16 dck
      complex conis
      parameter (conis = (0,1))
      real  temper, thetax, sig2
      complex ck(nspx)
      save

      logical, allocatable :: lcalc(:)
      complex, allocatable :: gg(:,:,:)
      complex, allocatable :: xphase(:,:,:)


      if (rfms .le. 0.0.and.(ispace.ne.0)) goto 900 !KJ added ispace

      allocate(lcalc(0:lx))
      allocate(gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx))
      allocate(xphase(nspx, -lx:lx, 0:nphx))
      gg(:,:,:)=cmplx(0)
	  xphase(:,:,:)=cmplx(0)
	  lcalc(:)=.true.

!     set default (LU) inv method
      minv = 0
      rdirec = 2*rfms
      toler1 = 0.e0
      toler2 = 0.e0

!     transform to single precision
      rat(:,:) = real (rath(:,:))
      temper =0.0e0
      thetax =0.0e0
      sig2  = 0.0e0

!     it will be nice to call yprep once for all energy points,
!     fix later, and now call it every time
      if (ie.eq.1 .or. lfms.eq.0 .or. lfms.eq.2) then
          if (ispace.eq.0) then
              continue
!	          call kprep(em,ne,nex)  !moved to xsphsub.f90 / reapot.f90 !KJ
          else
              call yprep(iph0,nat,inclus,nph,iphat,rfms,rat)
          endif
      endif

      if (inclus.gt.1.or.ispace.eq.0) then !KJ added ispace

!      call fms for a cluster around central atom
       if (ie.eq.1 .and. lverb .and. ispace.eq.1) then
	      if(ispace.eq.1) write (slog,35) inclus, iph0
  35      format ('FMS for a cluster of ',i4, ' atoms around atom type ',i3)
          if(ispace.eq.1) call wlog (slog)
       endif

       ck(1) = cmplx(sqrt(2*(em-eref)))
       do ipp = 0,nph
         do ill = -lipotx(ipp), lipotx(ipp)
           xphase(1, ill, ipp) = cmplx(ph( 1+abs(ill), ipp))
         enddo
       enddo
       iverb=0
       if (ie.eq.1) iverb = 1
       nsp = 1
       ispin = 0
!      Here at last real space and reciprocal space calculations separate :
       if (ispace.eq.0) then
         call fmskspace(ispin,ck,xphase,ie,em-eref,gg,iverb,dble(0))
       else
         call fms(lfms,nsp,ispin,inclus,nph,ck,lipotx,xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc,gg)
       endif


!      make ck= i, since coni is c*16
       do ip=0,nph
         if (lfms.ne.0 .or. ip.eq.iph0) then
           do il=0,lipotx(ip)
             ix = il**2
             do im=1,2*il+1
               gtr(il, ip) = gtr(il, ip) + gg(ix+im,ix+im,ip)
             enddo
             gtr(il,ip)= gtr(il,ip) * exp(2*conis*xphase(1,il,ip))/(2*il+1)
           enddo
         endif
       enddo
      endif


      ! Deallocate local variables
      deallocate(lcalc,gg,xphase)
 900  continue

      return
      end
