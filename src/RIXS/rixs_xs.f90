! Written by Josh Kas 9/1/2010
subroutine rixs_xs ( nph, edens, edenvl, dmag, vtot, vvalgs,      &
     &            rmt, rnrm, ixc, rhoint, vint, xmu, jumprm,   &
     &            x0, dx, rgrd, xion, iunf, iz,                &
     &            xnval, adgc, adpc, dgc,dpc,                  &
     &            ihole, qnrm, xnmues,                         &
     &            emin, emax, eimag, rfms2, lfms2, lmaxph,     &
     &            nat, iphat, rat, minv, rdirec,               &
     &            toler1, toler2 )

  !     print out LDOS in files rholNN.dat, if requested
  use controls,only: ispace
  use DimsMod, only: nphx=>nphu, natx, nrptx
  use constants
  use par
  implicit none

  ! Input
  integer, intent(in) :: nph,nat,ihole,ixc,iunf
  real*8,  intent(in) :: xnval(30,0:nphx)
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
  real*8,  intent(inout), dimension(251,30,0:nphx+1) :: dgc, dpc
  real*8,  intent(inout), dimension(10,30,0:nphx+1) :: adgc, adpc

  ! Local variables
  integer iph
  character*30  fname
  character*512 slog
  
  ! Added to satisfy implicit none
  integer :: msapp,lmaxsc,ne,ik0,jri,jri1,i,itmp,ie,ie2,ne1,ne3,idwopt,iph0
  real*8  :: de,enext,edge,vjump,tk,thetad,sig2g,critcw


  real*8,  dimension(nrptx) :: dmagx,ri
  integer, dimension(0:nphx) :: lmax

  ! work space
  real*8, allocatable, dimension(:) :: dum, vtotph, vvalph, vch1, vch2
  real*8, allocatable, dimension(:,:) :: dgcn, dpcn
 
  real*8, allocatable :: xrhoce(:,:,:)
  ! JPR - variable gtr is not used and it's size can conflict with the 
  ! size in the fmsdos code.  Therefore, it has been eliminated.
  !  complex, allocatable ::  gtr(:,:,:)
  complex*16, allocatable, dimension(:,:,:) :: xrhole, xbruce, ph
  integer, allocatable :: inclus(:)
  complex*16, allocatable :: em(:),eref(:)

  ! stuff from feff.f for rdinp, pathfinder and genfmt
  logical wnstar

  ! Following passed to pathfinder, which is single precision.
  ! Be careful to always declare these!
  integer, parameter:: necrit=9, nbeta=40
  real, allocatable, dimension(:,:,:) :: fbetac,fbeta
  real, allocatable, dimension(:) :: cksp,ckspc,xlam,xlamc
  real critpw, pcritk, pcrith
  complex*16 rl(nrptx), rlp(nrptx)
  
  ! Allocate local variables
  !-JPR-See above ! allocate(gtr(0:lx,0:nphx,nex)  
  allocate(dum(nrptx), vtotph(nrptx), vvalph(nrptx),           &
       &   dgcn(nrptx,30), dpcn(nrptx,30),                     &
       &   xrhoce(0:lx,nex,0:nphx),                            &
       &   xrhole(0:lx,nex,0:nphx), xbruce(0:lx,nex,0:nphx),   &
       &   ph(nex, ltot+1, 0:nphx), inclus(0:nphx),            & 
       &   em(nex), eref(nex), fbeta(-nbeta:nbeta,0:nphx,nex), &
       &   cksp(nex), xlamc(necrit), xlam(nex),vch1(nrptx),vch2(nrptx))



  !  msapp=2 - G_c + FMS only (no paths added)
  msapp=2

  if (emin.lt.emax) then
     !c    plot DOS between emin and emax with eimag above real axis
     call wlog('              LDOS calculation for specified grid')
     lmaxsc = lx
     ne = min (101, nex)
     !KJ debugging        ne = min (401, nex)
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


     if(ispace.eq.0) call kprep(em,ne,nex)  !KJ for rec space code


     write(slog,30) iph
30   format('     potential type ', i2)
     call wlog(slog)
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
     DO ie = 1, ne1
        call CalcRl( rgrd, x0, ri, ne, em,                           &
             &     ixc, rmt(iph), rnrm(iph),                         &
             &     vtotph, vvalph, dgcn, dpcn, eref(1),              &
             &     adgc(1,1,iph), adpc(1,1,iph), xrhole(0,1,iph),    &
             &     xrhoce(0,1,iph), ph(1,1,iph),                     &
             &     iz(iph), xion(iph), iunf, itmp, lmaxsc,           &
             &     xbruce(0,1,iph), xnval(1,iph), rl)

        DO ie2 = 1, ne1
           


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
             &      minv, rdirec, toler1, toler2)
        do iph0 = 1, nph
           inclus(iph0) = inclus(0)
        enddo
     else
        do iph0 = 0, nph
           call fmsdos(2, rfms2, lfms2,iph0,idwopt,tk,thetad,sig2g,    &
                &      lmaxph, nat, iphat, rat, inclus(iph0),          &
                &      ne, ne1, ne3, nph, em, eref, iz, ph,            &
                &       minv, rdirec, toler1, toler2)
        enddo
     endif

     do iph = 0,nph
        write (slog,'(a, i5)') ' Calculating chi and rho...', iph
        call wlog(slog) 
        call ff2rho (critcw, ne, xrhoce(0,1,iph), xrhole(0,1,iph), iph, msapp, em, lfms2, qnrm, xnmues, xmu, inclus)
     enddo


     call par_barrier
  endif

  ! Deallocate local variables
  deallocate(dum, vtotph, vvalph,dgcn, dpcn, xrhoce, xrhole,   &
       &     xbruce, ph, inclus, em, eref, fbeta, cksp, xlamc, xlam )
  !-JPR-See above  deallocate(gtr)

  return
end subroutine ldos
