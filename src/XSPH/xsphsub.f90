!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xsphsub.f90,v $:
! $Revision: 1.24 $
! $Author: jorissen $
! $Date: 2012/05/15 22:57:34 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xsph

  use controls,only : ispace
  use IOMOD
  use AtomicPotIO
  use DimsMod, only: nphx=>nphu, nex, nspx=>nspu, ltot, nrptx, lx
  use global_inp
  use atoms_inp
  use xsph_inp,  iPl=>IPlsmn
  use constants
  use potential_inp
  use ldos_inp
  use nrixs_inp
  use errorfile
  use hubbard_inp
! Cluster code -- multiple shell single scattering version of FEFF
! This program (or subroutine) calculates potentials and phase
! shifts for unique potentials specifed by atoms and overlap cards.
!
! Input files:  potph.inp    input data, atoms, overlaps, etc.
! Output:       phases.bin   phase shifts for use by the rest of the program.

  implicit none
  real*8 dx,x0,dxnew,rnrmav,xmu,vint,rhoint,emu_pot,s02,erelax,wp,rs,xf,qtotel,emu_apot,e_chsh,rint
  real*8 emu,tmpgam,edge,dum1,dum2,dum3,corr,xmunew,xmuvr,vjump,xnorm1,xnorm2,dum4,dum5,dum6,dum7
  integer i,i1,i2,ik0,iph,icenter,ios,ne,ne1,ne3,nsp,isp,idmag,ispinp,jnew,n,itmp,ie,kdif,iq
! Notes:
! nat    number of atoms in problem
! nph    number of unique potentials
! ihole  hole code of absorbing atom
! iph=0 for central atom

! Overlap calculation results
! edens(251,0:nphx)   -   overlapped density*4*pi
  real*8 edens(251,0:nphx)
! vtot (251,0:nphx)   -   overlapped total potential
  real*8 vtot (251,0:nphx), vstad(251), vclap (251,0:nphx)

! Muffin tin calculation results
! imt(0:nphx)  -  r mesh index just inside rmt
  integer imt(0:nphx), inrm(0:nphx)
! rmt(0:nphx)  -  muffin tin radius
  real*8 rmt(0:nphx)
! rnrm(0:nphx)  -  Norman radius
  real*8 rnrm(0:nphx), qnrm(0:nphx), folpx(0:nphx)

! PHASE output
! eref(nex, nspx)         -     interstitial energy ref
  complex*16 eref(nex, nspx)
! ph(nex,-ltot:ltot,nspx,0:nphx) - phase shifts
  complex*16 ph( nex, -ltot:ltot, nspx, 0:nphx)
! lmax(0:nphx)      -     number of ang mom levels
  integer lmax(0:nphx)

!  character*80 title(nheadx)

  complex*16  em(nex), xsec(nex,nspx)
  real*8 xsnorm(nex, nspx), dgc0(251), dpc0(251)

! additional data needed for relativistic version
  real*8 dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
  real*8 adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
  real*8 dgcn(nrptx,30), dpcn(nrptx,30)
  real*8 edenvl(251,0:nphx), eorb(30)
  integer kappa(30)
  real*8 vvalgs (251,0:nphx), xnval(30,0:nphx)
  integer iorb(-4:3,0:nphx)

! nrx = max number of r points for phase and xsect r grid
!KJ      parameter (nrx = nrptx)
  real*8 ri(nrptx), vtotph(nrptx), rhoph(nrptx)
  real*8  dmagx(nrptx), dmag(251,0:nphx)
  real*8 dgcx(nrptx), dpcx(nrptx), vvalph(nrptx), rhphvl(nrptx)
  real*8 vch (251), vchp(nrptx)
  real*8,allocatable :: vhelp(:,:) !KJ 12-2011 dummy

  logical UseStdm
  character*512 slog

! Josh - Added iPl for PLASMON card, and iint for interstitial index
  integer NPolesTmp, imtstd
  real*8 eorbTmp(30,0:nphx+1), rmtstd, a
!KJ 7-09 for NRIXS :
  real*8 dgcx0(nrptx), dpcx0(nrptx)
  complex*16,allocatable :: rkk(:,:,:,:) ! rkk(nex,nq,kfinmax,nspx)  !KJ kfinmax used to be 8 - set to 8 in nrixs_init if no nrixs
  !complex*16 rkk(nex,8,nspx) !KJ 12/2010 added 4th dimension "nq"
!KJ for HUBBARD:
  real*8   Vnlm(0:lx,(lx+1)**2,2,0:nphx)
  complex*16  aph(nex,lx+1,(lx+1)**2,2,0:nphx)

! Debug: Fer
! correorb input
  integer          sh_iz, sh_ihole, sh_jri, sh_kappa(30)
  real*8 sh_rmt, sh_dx, sh_p2f, sh_edge, sh_ri(nrptx), sh_dgcn(nrptx,30), sh_dpcn(nrptx,30), &
                       sh_adgc(10,30), sh_adpc(10,30), sh_eorb(30)
  complex*16       sh_vxc(nrptx)
! correorb output
  integer          sh_neg(30), sh_norbp
  real*8 sh_eng(nex,30), sh_rhoj(nex,30)
! Variables for testing
   complex*16 ts_eph, ts_vxc(nrptx), ts_ps(nrptx), ts_qs(nrptx), ts_aps(10), ts_aqs(10)
   integer    ts_kap(30), ts_nmax(30), ts_irr, ts_ic3, ts_vm(nrptx), ts_jri, ts_iwkb
   real*8     ts_rmt
! End of Fer

! NRIXS variables
!KJ      integer indmax
  integer iri
!
  real*8, allocatable :: xnmues(:,:)



  allocate(xnmues(0:lx,0:nphx))

  if(.not.(do_nrixs.eq.1)) nq=1
  allocate(rkk(nex,nq,kfinmax,nspx))
  rkk=dcmplx(0)
  corr = 0.d0

! Phase shift calculation
! Atom r grid
  UseStdm = .FALSE.
  dx = 0.05d0
  x0 = 8.8d0
! Phase r grid
  dxnew = rgrd

  call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,            &
                       emu_pot, s02, erelax, wp, ecv,rs,xf, qtotel,        &
                       imt, rmt, inrm, rnrm, folp, folpx, xnatph,      &
                       dgc0, dpc0, dgc, dpc, adgc, adpc,               &
                       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,&
                       eorb, kappa, iorb, qnrm, xnmues, nohole, ihole, &
                       inters, totvol, iafolp, xion, iunf, iz, jumprm)
!write(*,*) 'rnrm',rnrm
!write(*,*) 'nph,nphx',nph,nphx
if(.false.) then
  CALL ReadAtomicPots(nph = nph, iz = iz(0:nph), ihole = ihole, dmag = dmag(:,0:nph), &
            dgc0 = dgc0, dpc0 = dpc0, dgc = dgc(:,:,0:nph+1), dpc = dpc(:,:,0:nph+1),  &
            adgc = adgc(:,:,0:nph+1), adpc = adpc(:,:,0:nph+1), erelax = erelax,       &
            emu = emu_apot, xnval = xnval(:,0:nph+1), eorb = eorbTmp, rnrm = rnrm(0:nph))
else
            emu_apot=emu_pot
endif


    if (mldos_hubb .eq. 2) then
       Vnlm=0.d0
       open(21, file='v_hubbard.bin',form='unformatted',status='old')
       read(21) Vnlm
       close(21)
    endif

!write(*,*) 'rnrm',rnrm
!write(*,*) 'emu',emu_apot,emu_pot

      ! JK - change initial state here so that we can take matrix elements with one core
      !      state while using a different state for core-hole potential. For use with
      !      RIXS.
      IF(iCoreState.GT.0) ihole = iCoreState

! Get the occupation of this initial state. JK - find normalization for partially filled
! valence orbitals. For most cases this should be one.

! Default
  emu = emu_apot

! Use the one calculated in pot
  if ( ChSh_Type == 1 ) then
    emu = emu_pot
  end if

! Use the one calculated here
  if ( (ChSh_Type == 2) .or. (ChSh_Type == 3) ) then

! Set the target center
    iCenter = 0
    sh_iz            = iz(iCenter)
    sh_ihole         = ihole
    sh_rmt           = rmt(iCenter)
    sh_jri           = imt(iCenter)+1
    sh_dx            = dx
    sh_ri            = (/ (exp(-8.8+dx*(i-1)), i=1,nrptx) /)
    sh_p2f           = xmu-vint
    sh_edge          = xmu
    sh_vxc(1:251)    = vtot(1:251,iCenter) - vint
    sh_dgcn(1:251,:) = dgc(:,:,iCenter)
    sh_dpcn(1:251,:) = dpc(:,:,iCenter)
!   sh_adgc          = adgc(:,:,iCenter)
!   sh_adpc          = adpc(:,:,iCenter)
    sh_adgc          =  0.0d0
    sh_adpc          =  0.0d0
    sh_eorb          = eorbTmp(:,iCenter)
    sh_kappa         = kappa

    call correorb(sh_iz, sh_ihole, sh_rmt, sh_jri, sh_dx,sh_ri, sh_p2f, sh_edge, sh_vxc, &
                  sh_dgcn, sh_dpcn, sh_adgc, sh_adpc, sh_eorb, sh_neg, sh_eng, sh_rhoj, sh_kappa, sh_norbp, 0) !KJ iph=0
    e_chsh = -sh_eng(1,1)
    write(6,fmt='(a,f20.10)') 'e_chsh: ', e_chsh
    emu = e_chsh + erelax

  end if

  ! Write potentials for feffFP
  IF(.FALSE.) THEN
     OPEN(unit = 33, file = 'vtot.dat', status = 'replace')
     DO i1 = 0, nphx
        DO i2 = 1, 251
           WRITE(33, *) vtot(i2,i1)
        END DO
     END DO
  END IF

  vstad(:) = vtot(:,0)
  imtstd   = imt(0)+10
  rmtstd   = rmt(0)*2.d0 !exp(-x0+dx*(imtstd-1))
  rint = (rmt(0)+rmtstd)/2.d0
  a = pi/2.d0/(rmtstd - rint)
  !DO i = 1, imt(0)
  !   ri(i) = exp(-x0+dx*(i-1))
! !    write(48,*) ri(i), vstad(i), vtot(i,0) !KJ 7-09 commented out - looks like a debugging file
  !END DO
  !DO i = imt(0) + 1, 251
  !   ri(i) = exp(-x0+dx*(i-1))
  !   IF(ri(i).lt.rint) THEN
  !      vstad(i) = vint
  !   ELSEIF(ri(i).lt.rmtstd) THEN
  !      vstad(i) = vint*COS(a*(ri(i)-rint))**2
  !   ELSE
  !      vstad(i) = 0.d0
  !   END IF
! !    write(48,*) ri(i), vstad(i), vtot(i,0)  !KJ 7-09 commented out - looks like a debugging file
  !END DO
!  DO i = 1, 251
! dgc0(i) = dgc(i,ihole,0)
! dpc0(i) = dpc(i,ihole,0)
!  END DO
!  eorb(:) = eorbTmp(:,0)

!  lopt=true for the Rivas code of optical constants
!   lopt = .true.
!  lopt is now set by using the SETEDG card.
   if (lopt) then
     call getedg(ihole,iz(0), emu)
     ik0 = 1
     call wlog('Fixing edge energy from Elam table:')
     write(slog,fmt="('   emu = ',f10.3,' eV')") emu*hart
     call wlog(slog)
   endif


   novr(0:nph) = 0

!  update header, since e.g. one may use diff ixc for the same potential
   call sthead (ntitle, title, nph, iz, rmt, rnrm, xion, ihole, ixc, vr0, vi0, gamach, xmu, xf, vint, rs, nohole, lreal, rgrd)

   ! Josh - set losses to zero for now. This will be an option with CHBROAD card.
   tmpgam = gamach
   if(iGammaCH.eq.1)  gamach = 0.d0

! Make energy mesh
  edge = xmu - vr0
  if (.not.lopt) emu = emu - vr0

!c    manual input. Later make TDLDA and PMBSE cards
!c    TDLDA ifxc  (izstd=1 if TDLDA card is present)
! izstd = 0
! ifxc = 0
!c    PMBSE  ipmbse  nonlocal ifxc itdlda
! ipmbse = 2
! ibasis = 2
!c    ipmbse=0 (do not run); 1-LF only; 2-PM only; 3-combined;
!c           4-combined with s-function kernel splitting
! nonlocal = 0
!c    nonlocal = 0 (local fxc); 1-read W from pot.ch; 2-from yoshi.dat
! itdlda = 2
!c    itdlda = 0, 1, 2 should be run in sequence
!c    end of manual input

! check that logic is set up correctly
  if (ipmbse.le.0) itdlda = 0
  if (nohole.lt.0) then
!   core-hole potential is used already
    if (ifxc.ne.0) then
      call wlog(' Reset ifxc=0 since NOHOLE card is absent')
      ifxc = 0
      if (ipmbse.gt.0) nonlocal = 0
    endif
    if (ipmbse.eq.3 .and. izstd.eq.0) then
      call wlog(' Reset ipmbse=1 since NOHOLE card is absent')
      ipmbse = 1
    endif
  endif
  if (izstd.gt.0 .and. itdlda.gt.0) then
!   no need run PMBSE in this case
    call wlog(' Ignored PMBSE cards since TDLDA is present')
    itdlda = 0
  endif
  if (ipmbse.eq.2 .and. nonlocal.gt.0 .and. ifxc.gt.0) then
!   accounting for core-hole twice. reset ifxc=0
    call wlog(' Reset ifxc=0 since core-hole potential is used.')
    ifxc = 0
  endif
  if (ipmbse.eq.1 .and. nonlocal.gt.0) then
!   V_ch should be zero
    nonlocal = 0
  endif

! Josh - if nohole = 2, read wscrn.dat and add ch pot to vtot.
! Need to add file check and emesh check.
!  if (nohole.eq.2)  then
!     open (unit=13, file='wscrn.dat', status='old', iostat=ios)
!     call chopen (ios, 'wscrn.dat', 'ffmod2(xsph)')
!     read(13,*)  !KJ 12-2011 I added a header
!     open (unit=14, file='vtot.dat', status='replace',iostat=ios)
!     call chopen (ios, 'vtot.dat', 'ffmod2(xsph)')
!     do i = 1, 251
!        read(13,'(2e20.10)',end=20) dum1, dum2
!        dum3 = vtot(i,0)
!        vtot(i,0) = vtot(i,0) - dum2
!        write(14,'(3e20.10)') dum1, dum3, dum2
!     end do
! 20      continue
!
!     ! Add xc interaction. JK - 11/2010
!     DO i = 1, 251
!        ri(i) = exp(-x0+dx*(i-1))
!        ! Get rhotot
!         ! dum1 is rs here.
!        if(edens(i,0).le.0) then
!           dum1 = 10.d0
!        else
!           dum1 = (edens(i,0)/3)**(-third)
!        end if
!        CALL vbh(dum1,0.d0,dum2) ! dum2 is vxctot
!        ! Now do it for density with core-hole
!        ! dum1 = rho - rhoch
!        dum1 = edens(i,0) - (dgc0(i)**2 + dpc0(i)**2)/ri(i)**2
!        ! dum1 = rs
!        if(dum1.le.0.d0) then
!           dum1 = 10.d0
!        else
!           dum1 = (dum1/3)**(-third)
!        end if
!        CALL vbh(dum1,0.d0,dum3) ! dum3 is vxc with ch
!        dum1 = dum3 - dum2
!        !KJ WRITE(18,'(8(f10.4,x))') ri(i), dum1,dum2,dum3   !KJ 1-2012 This was going to "fort.18" so I'm guessing it's debugging output?
!        vtot(i,0) = vtot(i,0) + dum1
!     END DO
!
!     nohole = 0
!     close(13)
!     close(14)
!  end if
! Josh END

!  write(*,*) 'vtot(:,0)',vtot(:,0)

  if (itdlda.eq.0)  then
!KJ         IF(.FALSE.) THEN
!KJ            call phmesh (ipr2, ispec, edge, emu, vi0, gamach, ecv, xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3      &
!KJ                 &                 ,1,egrid3a,egrid3b,egrid3c,egridfile) !KJ added this line 01/07
!KJ         ELSE
!KJ from feffq:      if xkmax<0 call constant E step, otherwise call the 'normal' mesh code

        if ((xkmax.gt.0.0d0).or.(do_nrixs.ne.1)) then
          IF (electronic_temperature.GT.0) THEN
            !PRINT*, "USE PHMESH2T "
            CALL phmesh2T(ipr2, ispec, edge, emu, vi0, gamach, ecv, xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3, iGrid)
          ELSE
            !PRINT*, "USE PHMESH 0 "
            CALL phmesh2(ipr2, ispec, edge, emu, vi0, gamach, ecv, xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3, iGrid)
          ENDIF
!KJAleksi's code :                call phmesh (ipr2, ispec, edge, emu, vi0, gamach, ecv, xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3)
        else    !feff-q constant energy step grid
            xkmax = abs(xkmax)
            call phmeshjas (ipr2, ispec, edge, emu, vi0, gamach, ecv, xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3)
        end if
!KJ         END IF
  else
! nesvi TDLDA
        call meshlda (xkstep, ne, ne1, ne3, em, ik0)  !KJ added ne3
        corr = 1.0
  endif


!KJ for reciprocal space code :
!KJ kprep must run after e-mesh has been prepared.
!   Only do this if szlz is going to be called below (only place in xsph that actually goes to k-space by calling fmssz):
  if(ispace.eq.0) call kprep(em,ne,(ipr2.ge.3))
!KJ end


  if (itdlda.eq.1) then
!   to get the mesh only
    do 93  i = 1, ne1
      write(3,94) dble(em(i))*hart
 94       format (7f10.5)
 93     continue
    call WipeErrorfileAtFinish
    stop 'TDLDA energy mesh is written out'
!   end of itdlda=1 calculations
  endif

! Make old grid to report distances in xsect.bin
  do 95 i = 1, 251
 95   ri(i) = exp(-x0+dx*(i-1))

! open xsect.dat and write the header
  open (unit=1, file='xsect.dat', status='unknown', iostat=ios)
  call chopen (ios, 'xsect.dat', 'potph')
  call wthead (1, ntitle, title)
  write(1,45)
   45 format ('# ',1x, 71('-'))
  write(1,55) s02, erelax, wp, edge, emu
   55 format ('# ', 3e13.5, 2e15.7, ' method to calculate xsect')
  write(1,56) tmpgam*hart, ne1, ik0
   56 format ('# ',1p, e15.7, 3i7,' gamach in eV, # of points on horizontal axis')
  write(1,57)
   57 format ('# ','      em              xsnorm            xsec  ')
! end of the xsect.dat header

! nsp = 1 - spin average caclulations; 2 - spin up and down
  nsp = 1
  if (abs(ispin).eq.1) nsp = nspx
 !KJ new code 10-2014:
 !nsp = nspu
 !KJ changed again:
 ! nspx now linked to nspu through "use" statement
 !nsp = nspx ! JK - This was not correct for ispin = 2, which requires nsp = 1 even when nspx = 2.
 !PRINT*, 'nsp=', nsp
 if (nsp.gt.nspx) then
    call wlog('ERROR - FEFF wants to do a spin-polarized calculation in xsphsub but nspx=1 (arrays are too small).')
    stop
 endif

! scale spin density on each atom appropriately
  do iph = 0, nph
   do i = 1, 251
     dmag(i,iph) = dmag(i, iph) * spinph(iph)
   enddo
  enddo

  do 300 isp = 1, nsp
    if (ispin.ne.0) then
! make spin dependent potential if needed
! isp = 1-spin-down; 2-spin-up potentials
      idmag = (-1)**isp
      if (nsp.eq.1) then
         idmag = 1
         if(ispin.lt.0) idmag=-1
      endif
      call  istprm (nph, nat, iphat, rat, iatph, xnatph, novr, iphovr, nnovr, rovr, folp, folpx, iafolp,    &
     &               edens, edenvl, idmag, dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,   &
     &               ixc, rhoint,vint, rs, xf, xmu, xmunew, rnrmav, qtotel, inters, totvol)
      xmunew = xmu
      if (abs(ispin).eq.1 .and. nsp.eq.2) then
!   |ispinp| = |ispin|, but sign is determined by isp
         ispinp = abs(ispin)
         if (isp.eq.1) ispinp = -ispinp
      else
!    sign is determined by spin (always for ispin=-2,2)
         ispinp = ispin
      endif
    else
! spin-independent case
      ispinp = 0
    endif

! Josh - if nohole = 2, read wscrn.dat and add ch pot to vtot.
! Need to add file check and emesh check.
    if (nohole.eq.2)  then
       open (unit=13, file='wscrn.dat', status='old', iostat=ios)
       call chopen (ios, 'wscrn.dat', 'ffmod2(xsph)')
       read(13,*)  !KJ 12-2011 I added a header
       open (unit=14, file='vtot.dat', status='replace',iostat=ios)
       call chopen (ios, 'vtot.dat', 'ffmod2(xsph)')
       do i = 1, 251
          read(13,'(2e20.10)',end=20) dum1, dum2
          dum3 = vtot(i,0)
          vtot(i,0) = vtot(i,0) - dum2
          write(14,'(3e20.10)') dum1, dum3, dum2
       end do
20     continue

       ! Add xc interaction. JK - 11/2010
       DO i = 1, 251
          ri(i) = exp(-x0+dx*(i-1))
          ! Get rhotot
          ! dum1 is rs here.
          if(edens(i,0).le.0) then
             dum1 = 10.d0
          else
             dum1 = (edens(i,0)/3)**(-third)
          end if
          CALL vbh(dum1,0.d0,dum2) ! dum2 is vxctot
          ! Now do it for density with core-hole
          ! dum1 = rho - rhoch
          dum1 = edens(i,0) - (dgc0(i)**2 + dpc0(i)**2)/ri(i)**2
          ! dum1 = rs
          if(dum1.le.0.d0) then
             dum1 = 10.d0
          else
             dum1 = (dum1/3)**(-third)
          end if
          CALL vbh(dum1,0.d0,dum3) ! dum3 is vxc with ch
          dum1 = dum3 - dum2
          !KJ WRITE(18,'(8(f10.4,x))') ri(i), dum1,dum2,dum3   !KJ 1-2012 This was going to "fort.18" so I'm guessing it's debugging output?
          vtot(i,0) = vtot(i,0) + dum1
       END DO

       nohole = 0
       close(13)
       close(14)
    end if
! Josh END

!   calculate operators of interest (s_z, l_z, t_z)
    xmuvr = xmu - vr0
    if (ipr2.ge.3) call szlz(ispinp,ecv,nph,nat,rgrd,nohole,rfms2,  &
     &     lfms2, lmaxph, edens, edenvl, dmag, vtot, vvalgs, rmt, rnrm, &
     &     ixc, rhoint, vint, xmuvr, jumprm, xnval, iorb, x0, dx, xion, iunf, iz,    &
     &     adgc, adpc, dgc, dpc, ihole, rat, iphat, corr)
!    1                   em, ne1, ne, ik0 )

!   Cross section calculation, use phase mesh for now
!   Absorbing atom is iph=0
    write(slog,*) '   absorption cross section'
    call wlog(slog)
    iph = 0
    IF(UseStdm) THEN
       call fixvar (rmtstd, edens(1,0), vstad(1), dmag(1,0), 0.d0, 0.d0, dx, dxnew, jumprm, vjump, ri, vtotph, rhoph, dmagx)
    ELSE
       call fixvar (rmt(0), edens(1,0), vtot(1,0), dmag(1,0), vint, rhoint, dx, dxnew, jumprm, vjump, ri, vtotph, rhoph, dmagx)
    END IF
    call fixdsx (iph, dx, dxnew, dgc, dpc, dgcn, dpcn)
    if (mod(ixc,10) .ge. 5) then
       if (jumprm .gt. 0) jumprm = 2
       call fixvar (rmt(0), edenvl(1,0), vvalgs(1,0), dmag(1,0), vint, rhoint, dx, dxnew, jumprm, vjump, ri, vvalph, rhphvl, dmagx)
       if (jumprm .gt. 0) jumprm = 1
    endif
    call fixdsp (dx, dxnew, dgc0, dpc0, dgcx, dpcx, jnew)

    IF (do_nrixs .eq. 1) then
!c      Aleksi a.k.a. JAS wants to have the core state (dgcx0/dpcx0)
!c      calculated with the core hole present. It is read into dgc/dpc but
!c      needs to be stored into a new array that is cut
!c      the same way as dgc0/dpc0 (i.e. initial state core wf)
    call getholeorb0(dx,dxnew,ihole,jnew,iz,xion,iunf,dgc,dpc,dgcx0,dpcx0)

    call xsectjas (ipr2, dxnew, x0, ri, ne, ne1, ik0, em, edge, &
         ihole, emu, corr, dgcx, dpcx, jnew,dgcx0,dpcx0, & !dgcx0,dpcx0 are new for jas
         ixc0, lreal, rmt(0), rnrm(0), xmu, vi0, iPl, NPolesTmp, Eps0, EGap, & !iPl new JK, NPoles and Eps0 new JJK 3/9/2010 !KJ copied all of this stuff here
         gamach, vtotph, vvalph, rhoph, dmagx, rhphvl, &
         dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), xsec(1,isp), &
         xsnorm(1,isp), rkk(1,1,1,isp),iz(0), xion(0), iunf, &  !KJ 12/10 rkk(1,1,isp)
         xnval(1,iph), iorb(-4,iph), l2lp, &
         ipol, ispinp, abs(le2), angks,ptz, iph) !KJ iph      !, & ! abs is new
 !            qtrans,kfinmax,jmax,jinit,indmax,ldecmx) ! nrixs vars

    ELSE  ! no nrixs

    if (itdlda.eq.0) then
! Josh - added argument iPl to control many pole self energy
      NPolesTmp = NPoles
!      write(*,*) 'Breakdown at xsect'
      call xsect (ipr2, dxnew, x0, ri, ne, ne1, ik0, em, edge,      &
         ihole, emu, corr, dgcx, dpcx, jnew,                        &
         ixc0, lreal, rmt(0), rnrm(0), xmu, vi0, iPl, NPolesTmp, Eps0, EGap, & !iPl new JK, NPoles and Eps0 new JJK 3/9/2010
         gamach, vtotph, vvalph, rhoph, dmagx, rhphvl,              &
         dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), xsec(1,isp),     &
         xsnorm(1,isp), rkk(1:nex,1:nq,1:kfinmax,isp), iz(0), xion(0), iunf,         &  !KJ 12/10
         xnval(1,iph), izstd, ifxc, eorb, kappa, iorb(-4,iph), l2lp,& ! izstd,ifxc,eorb,kappa are new
         ipol, ispinp, le2, angks,ptz, iph) !KJ iph
 !    write(*,*) 'Not broken down yet'
    else
      if (nonlocal.gt.0) then
!   read potential with core-hole from a file
        if (nonlocal.eq.1) then
        !KJ 12-2011 To make sure there are no dimension mismatches, I have introduced a dummy array here.
          allocate(vhelp(251,0:nphx))
          vhelp=0.d0
          call rdpotp(vhelp)
          vch(:)=vhelp(:,0)
          deallocate(vhelp)
          !KJ call rdpotp(vch)
        elseif (nonlocal.eq.2) then
          open (unit=3, file='wscrn.dat', status='old')
          read(3,*)  !KJ 12-2011 I added a header

          n=0
 338      n = n+1
          read(3,337, end=339) dum1, dum2
!         use frac.ne.1  to mix bare and screened ch pot
!         frac = 0.80
          vch(n) = -1.d0*dum2
 337      format(6e20.10)
          goto 338
 339      continue
          close (unit=3)
        endif
        call fixvar (rmt(0), edens(1,0), vch, dmag(1,0), vint, rhoint, dx, dxnew, jumprm, vjump, ri, vchp, rhoph, dmagx)
        do 333 i = 1, nrptx
           if (ri(i).lt.rmt(0)) then
             if (nonlocal.eq.1) vchp(i) = vchp(i) - vtotph(i)
           elseif (ri(i).lt.40.d0) then
!             assume const/r behaviour
              vchp(i) = vchp(i-1) * ri(i-1) / ri(i)
           else
              vchp(i) = 0
           endif
!          testing: write core-hole potential in fort.17
           if (ri(i).lt.40.d0) write(17,'(2f30.5)') ri(i), vchp(i)
  333   continue

        close (unit=17)
      else
        vchp(1:nrptx)=0
      endif

! Josh - added argument iPl to control many pole self energy
      NPolesTmp = NPoles

      call xsectd (ipr2,dxnew, x0, ri, ne, ne1, ik0, em, edge,      &
     &       ihole, emu, corr, dgcx, dpcx, jnew,                        &
     &       ixc0, lreal, rmt(0), rnrm(0), xmu, vi0, iPl, NPolesTmp, Eps0, EGap, &
     &       gamach, vtotph, vvalph,vchp, rhoph, dmagx, rhphvl,         &
     &       dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), xsec(1,isp),     &
     &       xsnorm(1,isp), rkk(1,1,1,isp),iz(0), xion(0), iunf,          &
     &       xnval(1,iph), ipmbse, ifxc, ibasis, eorb, kappa,           &
     &       iorb(-4,iph), l2lp, ipol, ispinp, le2, angks,ptz, itdlda, iph) !KJ iph
    endif

  ENDIF !NRIXS?



    do iph = 0, nph
      write(slog,'(4x, a, i5)') 'phase shifts for unique potential', iph
      call wlog(slog)
! fix up variable for phase
      call fixvar (rmt(iph), edens(1,iph), vtot(1,iph), dmag(1,iph), &
        vint, rhoint, dx, dxnew, jumprm, vjump, ri, vtotph, rhoph, dmagx)
!  write(*,*) '##################### IPH ',iph
!  write(*,*) 'rmt(iph)',rmt(iph)
!  write(*,*) 'vtot(iph)',vtot(:,iph)
!  write(*,*) 'dmag(iph)',dmag(:,iph)
!  write(*,*) 'edens(iph)',edens(:,iph)
!  write(*,*)
!  write(*,*)

      if (mod(ixc,10) .ge.5) then
         if (jumprm .gt. 0) jumprm = 2
         call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph), dmag(1,iph), &
           vint, rhoint, dx, dxnew, jumprm, vjump, ri, vvalph, rhphvl, dmagx)
         if (jumprm .gt. 0) jumprm = 1
         call fixdsx (iph, dx, dxnew, dgc, dpc, dgcn, dpcn)
      endif
      if (iph .eq. 0)  then
         itmp = ihole
      else
         itmp = 0
      endif

      NPolesTmp = NPoles
      if (mldos_hubb .eq. 1) then
      call phase (iph, dxnew, x0, ri, ne, ne1, ne3, em, ixc, nsp,   &
     &            lreal, rmt(iph), xmu, vi0, iPl, NPolesTmp, Eps0, EGap,&
     &            gamach, vtotph, vvalph, rhoph, dmagx, rhphvl,         &
     &            dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), eref(1,isp),&
     &            ph(1,-ltot,isp,iph), lmax(iph), iz(iph), itmp,        &
     &            xion(iph), iunf, xnval(1,iph), ispinp)
     elseif (mldos_hubb .eq. 2) then
 !     write(*,*) 'Breakdown at phase_h'
      call phase_h (iph, dxnew, x0, ri, ne, ne1, ne3, em, ixc, nsp,   &
     &            lreal, rmt(iph), xmu, vi0, iPl, NPolesTmp, Eps0, EGap,&
     &            gamach, vtotph, vvalph, rhoph, dmagx, rhphvl,         &
     &            dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), eref(1,isp),&
     &            ph(1,-ltot,isp,iph), lmax(iph), iz(iph), itmp,        &
     &            xion(iph), iunf, xnval(1,iph), ispinp, isp, Vnlm(:,:,isp,iph), &
                  aph(:,:,:,isp,iph))
      endif
    enddo

 300  continue  ! end of big loop over spins

! write main output to xsect.dat
  340 format (e17.9, 4e13.5)
  if (abs(ispin).ne.1 .or. nspx.eq.1) then
    do ie = 1, ne
       write(1,340) dble(em(ie))*hart, dimag(em(ie))*hart, xsnorm(ie,1), dble(xsec(ie,1)), dimag(xsec(ie,1))
    enddo
  else
!   nspx = 2
    do ie = 1, ne
       write(1,340) dble(em(ie))*hart, dimag(em(ie))*hart,          &
     &             (xsnorm(ie,1) + xsnorm(ie,nspx)) / 2.d0 ,            &
     &           dble( (xsec(ie,1) + xsec(ie,nspx)) ),                  &
     &          dimag( (xsec(ie,1) + xsec(ie,nspx)) )
!  Normalize rkk to the average over up/down spin
!  nsp=2
       if (nq.eq.1) then  !KJ fix later this is "nq=0" in Aleksi's code ...
       xnorm1 = sqrt( 2*xsnorm(ie,1) / (xsnorm(ie,1) + xsnorm(ie,nspx)) )
       xnorm2 = sqrt( 2*xsnorm(ie,nspx) / (xsnorm(ie,1) + xsnorm(ie,nspx)) )
       do kdif = 1,8
       do iq=1,nq
         rkk (ie, iq, kdif, 1) = rkk (ie, iq, kdif, 1) * xnorm1
         rkk (ie, iq, kdif, nspx) = rkk (ie, iq, kdif, nspx) * xnorm2
       enddo
       enddo
       endif
    enddo
  endif
  close (unit=1)

! disable for now since dimensions are different
  if (ipr2 .ge. 2)  call wphase (nph, em, eref, lmax, ne, ph, ntitle, title)

! Write out phases for paths and genfmt
  if(do_nrixs.ne.1) then
     kfinmax=8
     indmax=8
  endif

  call wrxsph (nsp, ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge,ik0,&
               em, eref, lmax, iz, potlbl, ph, rkk,kfinmax,indmax,ipr2) !KJ 7-09 added ipr2

  if(mldos_hubb.eq.2) then
     open(27,file='aphase_hubbard.bin',form='unformatted')
     write(27) aph
     close(27)
  endif

  if (ipr2 .ge. 1) call axafs (em, emu, xsec(1,1), ne1, ik0)
!   calculate axafs.  axafs does not make sense for |ispin| = 1


  deallocate(xnmues)

  return
  end
