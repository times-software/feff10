!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ldossub.f90,v $:
! $Revision: 1.10 $
! $Author: bmattern $
! $Date: 2013/01/15 05:20:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ldos ( nph, edens, edenvl, dmag, vtot, vvalgs,      &
     &            rmt, rnrm, ixc, rhoint, vint, xmu, jumprm,   &
     &            x0, dx, rgrd, xion, iunf, iz,                &
     &            xnval, adgc, adpc, dgc,dpc,                  &
     &            ihole, qnrm, xnmues,                         &
     &            emin, emax, eimag, rfms2, lfms2, lmaxph,     &
     &            nat, iphat, rat, minv, rdirec,               &
     &            toler1, toler2 )

  !     print out LDOS in files rholNN.dat, if requested
  use controls,only: ispace
  use DimsMod, only: nphx=>nphu, nrptx, lx, nex, ltot, natx
  use constants
  use par
  use ldos_inp,only: neldos,ldostype
  implicit none
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
  ! Input
  integer, intent(in) :: nph,nat,ihole,ixc,iunf
  real*8,  intent(in) :: xnval(41,0:nphx)
  real*8,  intent(in) :: xnmues(0:lx,0:nphx)
  integer, intent(in) :: iphat(natx)
  real*8,  intent(in) :: rat(3,natx)
  real,    intent(in) :: rfms2, rdirec, toler1, toler2
  real*8,  intent(in) :: emin,emax,rhoint,vint,x0,dx,rgrd

  integer, intent(in), dimension(0:nphx) :: iz,lmaxph
  real*8,  intent(in), dimension(0:nphx) :: rmt, rnrm,qnrm,xion
  real*8,  intent(in), dimension(251,0:nphx) :: dmag,vtot,vvalgs


  ! Input and Output
  integer, intent(inout) :: minv,lfms2,jumprm
  real*8,  intent(inout) :: xmu, eimag

  real*8,  intent(inout), dimension(251,0:nphx) :: edens,edenvl
  real*8,  intent(inout), dimension(251,41,0:nphx) :: dgc, dpc   !KJ 12-2011 changed "nphx+1" to "nphx" ; extra field isn't used anywhere, and conflicts with declaration in other routines
  real*8,  intent(inout), dimension(10,41,0:nphx) :: adgc, adpc  !KJ 12-2011 changed "nphx+1" to "nphx" ; extra field isn't used anywhere, and conflicts with declaration in other routines

  ! Local variables
  integer iph
  character*30  fname
  character*512 slog
  
  ! Added to satisfy implicit none
  integer :: msapp,lmaxsc,ne,ik0,jri,jri1,i,itmp,ie,ne1,ne3,idwopt,iph0
  real*8  :: de,enext,edge,vjump,tk,thetad,sig2g,critcw


  real*8,  dimension(nrptx) :: dmagx,ri
  integer, dimension(0:nphx) :: lmax

  ! work space
  real*8, allocatable, dimension(:) :: dum, vtotph, vvalph
  real*8, allocatable, dimension(:,:) :: dgcn, dpcn
 
  real*8, allocatable :: xrhoce(:,:,:)
  ! JPR - variable gtr is not used and its size can conflict with the 
  ! size in the fmsdos code.  Therefore, it has been eliminated.
  !  complex, allocatable ::  gtr(:,:,:)
  complex*16, allocatable, dimension(:,:,:) :: xrhole, ph
  integer, allocatable :: inclus(:)
  complex*16, allocatable :: em(:),eref(:)

  ! Allocate local variables
  !-JPR-See above ! allocate(gtr(0:lx,0:nphx,nex)  
  allocate(dum(nrptx), vtotph(nrptx), vvalph(nrptx),           &
       &   dgcn(nrptx,41), dpcn(nrptx,41),                     &
       &   xrhoce(0:lx,nex,0:nphx),                            &
       &   xrhole(0:lx,nex,0:nphx),   &
       &   ph(nex, ltot+1, 0:nphx), inclus(0:nphx),            & 
       &   em(nex), eref(nex)) 

  !  msapp=2 - G_c + FMS only (no paths added)
  msapp=2

  if (emin.lt.emax) then
     !c    plot DOS between emin and emax with eimag above real axis
     !if(master)call wlog('              LDOS calculation for specified grid')
     lmaxsc = lx

     ne = neldos
	 if (ne > nex) then
	   ne = nex
	   write(slog,*) "Warning: neldos = ", neldos, "is larger than nex = ", nex, "."
	   call wlog(slog)
	 endif

     de = (emax-emin)/(ne-1)
     if (eimag.lt.0) eimag=3*de
     enext=emin
     do i=1,ne
        em(i) = enext + coni*eimag
        enext= enext + de
     enddo

     !       ik0 is a starting point for path filters in energy
     ik0 = ne-45
     edge = xmu


     if(ispace.eq.0) call kprep(em,ne,.true.)  !KJ for rec space code

     do iph = 0, nph
        write(slog,30) iph
30      format('     potential type ', i2)
        !if(master)call wlog(slog)
        lmax(iph) = lx

        call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph),  &
             &       vint, rhoint, dx, rgrd, jumprm,                 &
             &       vjump, ri, vtotph, dum, dmagx)
        if (mod(ixc,10) .ge.5) then
           if (jumprm .gt. 0) jumprm = 2
           call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph),        &
                &       dmag(1,iph), vint, rhoint, dx, rgrd , jumprm,  &
                &       vjump, ri, vvalph, dum, dmagx)
           if (jumprm .gt. 0) jumprm = 1
        endif
        call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)

        jri = (log(rmt(iph)) + x0) / rgrd + 2
        jri1 = jri+1
        eref(1) = vtotph(jri1)
        do i = 1, jri1
           vtotph(i) = vtotph(i) - eref(1)
        enddo
        if (ixc.ge.5) then
           do i = 1, jri1
              vvalph(i) = vvalph(i) - eref(1)
           enddo
        else
           do i = 1, jri1
              vvalph(i) = vtotph(i)
           enddo
        endif

        itmp = 0
        call rhol( rgrd, x0, ri, ne, em,                             &
             &     ixc, rmt(iph), rnrm(iph),                         &
             &     vtotph, vvalph, dgcn, dpcn, eref(1),              &
             &     adgc(1,1,iph), adpc(1,1,iph), xrhole(0,1,iph),    &
             &     xrhoce(0,1,iph), ph(1,1,iph),                     &
             &     iz(iph), xion(iph), iunf, itmp, lmaxsc,           &
             &     xnval(1,iph),                    &
             &     iph) !KJ iph
     enddo


     !       Write out phases for fmsdos, paths and genfmt
     do ie =1, ne
        eref(ie) = eref(1)
     enddo
     ne1 = ne
     ne3 = 0

     !c      call fms for a cluster around central atom
     !JPR - gtr has been removed as the final arg to fmsdos
     if (lfms2.ne.0) then
        iph0 = 0
        call fmsdos(2, rfms2, lfms2, iph0, idwopt, tk,thetad,sig2g,   &
             &      lmaxph, nat, iphat, rat, inclus(0),               &
             &      ne, ne1, ne3, nph, em, eref, iz, ph,              &
             &      minv, rdirec, toler1, toler2, ldostype)
        do iph0 = 1, nph
           inclus(iph0) = inclus(0)
        enddo
     else
        do iph0 = 0, nph
           call fmsdos(2, rfms2, lfms2,iph0,idwopt,tk,thetad,sig2g,    &
                &      lmaxph, nat, iphat, rat, inclus(iph0),          &
                &      ne, ne1, ne3, nph, em, eref, iz, ph,            &
                &       minv, rdirec, toler1, toler2, ldostype)
        enddo
     endif

     do iph = 0,nph
        if(master)write (slog,'(a, i5)') 'Writing DOS for atom type ', iph
        call wlog(slog) 
        call ff2rho (critcw, ne, xrhoce(0,1,iph), xrhole(0,1,iph), iph, msapp, em, lfms2, qnrm, xnmues, xmu, inclus, ldostype)
     enddo

     call par_barrier
  endif

  ! Deallocate local variables
  deallocate(dum, vtotph, vvalph,dgcn, dpcn, xrhoce, xrhole,   &
       &     ph, inclus, em, eref)

  !-JPR-See above  deallocate(gtr)

  return
end subroutine ldos
