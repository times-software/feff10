!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: fmssz.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2010/05/27 23:14:52 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fmssz( iph0, ie, em, eref, ph, iz, nph,                &
     &           rfms, lfms, nat, iphat, rath, amat, lipotx, gctr, gtr)
!     uses Bruce Ravel subroutine to do FMS in self-consistency loop
!     written by alexei ankudinov 06.1997

      use DimsMod, only: natx,nphx=>nphu, lx, nspx=>nspu
      use controls,only : ispace

      implicit double precision (a-h, o-z)

!     input
      dimension iphat(natx), rath(3,natx)
      real rat(3,natx), rfms
      real rpart,aipart
      integer nph
      dimension iz(0:nphx)

      complex*16, intent(inout) :: ph(lx+1, 0:nphx)
      complex, intent(inout) :: gtr(2,2, 3,0:lx, 0:nphx)
      real, intent(inout) :: amat(-lx:lx,2,2, 3,0:lx)
      real, intent(inout) ::gctr(2,2, 3,0:lx, 0:nphx)

!     work space
      integer iph
      complex*16 em, eref
      character*512 slog
!     fms staff
      integer lipotx(0:nphx)
      complex*16 dck
      real  rdirec, toler1, toler2
      complex ck(nspx)

      complex, allocatable :: xphase(:,:,:), gg(:,:,:)
      logical, allocatable :: lcalc(:)

      save
      
      ! Allocate local variables
      allocate(xphase(nspx, -lx:lx, 0:nphx))
      allocate(gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx))
      allocate(lcalc(0:lx))

      if (rfms .gt. 0.or.ispace.eq.0) then  !KJ added ispace
        do 25 iat=1,nat
        do 25 j=1,3
   25   rat(j,iat) = real (rath(j,iat))

!       transform to single precision
        minv = 0
        rdirec = 2*rfms
        toler1 = 0.e0
        toler2 = 0.e0


!       it will be nice to call yprep once for all energy points,
!       fix later, and now call it every time
!KJoriginal lines :        if (ie.eq.1 .or. lfms.eq.0) 
!     1    call yprep(iph0, nat, inclus, nph, iphat, rfms, rat,
!     2       iz, rdirec)
!KJ now mine :
        if (ie.eq.1 .or. lfms.eq.0) then
             if (ispace.eq.0) then
                 continue
!	    !    call kprep(em,ne,nex)  !moved to xsph.f !KJ
           else
              call yprep(iph0,nat,inclus,nph,iphat,rfms,rat)
           endif
        endif
!KJ

        if (inclus.gt.1.or.ispace.eq.0) then !KJ added ispace
!c        call fms for a cluster around central atom
          if (ie.eq.1) then
             write (slog,35) inclus, iph0
  35         format ('        Doing FMS for a cluster of ',i3,          &
     &       ' atoms around iph = ',i2)
             call wlog (slog)
          endif

          dck=sqrt(2*(em-eref))
          rpart  = real(dble(dck))
          aipart = real(dimag(dck))
          ck(1) = cmplx(rpart, aipart)
          do 50 ipp = 0,nph
            do 40 ill = -lipotx(ipp), lipotx(ipp)
              rpart  = real(dble( ph(abs(ill)+1,ipp)))
              aipart = real(dimag(ph(abs(ill)+1,ipp)))
              xphase(1, ill, ipp) = cmplx(rpart, aipart)
  40        continue
  50      continue
          iverb=0
          if (ie.eq.1) iverb = 1
!         neglect spin-flip processes (fix later for ispin=1)
          nsp = 1
          ispin = 0
          do 55 ill = 0, lx
  55      lcalc(ill) = .true.
!KJ original call to fms :
!KJ          call fms(lfms, nsp, ispin, inclus, nph, ck, lipotx, xphase,ie,
!KJ     1     iverb, minv, rdirec, toler1, toler2, lcalc,gg)

!KJ now my code :
! Here at last real space and reciprocal space calculations separate :
          if (ispace.eq.0) then ! KJ
            call fmskspace(ispin,ck,xphase,ie,em-eref,gg,iverb,dble(0))  ! KJ
          else           ! KJ
!            call fms(lfms,nsp,ispin,inclus,npot,ck,lmaxph,xphase,ie, iverb,minv,rdirec,toler1,toler2,lcalc,gg)
            call fms(lfms,nsp,ispin,inclus,npot,ck,lipotx,xphase,ie, iverb,minv,rdirec,toler1,toler2,lcalc,gg)
          endif      !  KJ
!KJ
        endif
      endif

      do 200 ip=0,nph

        if (lfms.ne.0 .or. ip.eq.iph0) then
          do 190 lpp =0,lipotx(ip)
             ix1 = lpp**2 
             do 170 im=1,2*lpp+1
!              now cycle over gtr dimensions
               do 100 iop = 1,3
               do 100 i2 = 1,2
               do 100 i1 = 1,2
                 if (rfms.gt.0 .and. inclus.gt.0) gtr(i1,i2,iop,lpp,ip)= &
     &             gtr(i1,i2,iop,lpp,ip) + amat(im-lpp-1,i1,i2,iop,lpp) &
     &             * gg(ix1+im,ix1+im,ip)
                 gctr(i1, i2, iop,lpp,ip)= gctr(i1, i2, iop,lpp,ip)     &
     &             + amat(im-lpp-1,i1,i2,iop,lpp)
 100           continue
 170         continue
 190      continue
        endif
 200  continue

      ! Deallocate local variables
      deallocate(xphase,gg,lcalc)

      return
      end
