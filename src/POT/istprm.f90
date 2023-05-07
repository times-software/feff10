!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: istprm.f90,v $:
! $Revision: 1.15 $
! $Author: jorissen $
! $Date: 2011/11/30 22:57:15 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine istprm ( nph, nat, iphat, rat, iatph, xnatph,          &
  &                novr, iphovr, nnovr, rovr, folp, folpx, iafolp,   &
  &                edens, edenvl, idmag,                             &
  &                dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,  &
  &                ixc, rhoint, vint, rs, xf, xmu, xmunew,           &
  &                rnrmav, qtotel, inters, totvol)

  !     Finds interstitial parameters, rmt, vint, etc.
  use constants
  use DimsMod, only: natx, nphx=>nphu, novrx, novp
  use pdw_mod
  use potential_inp, only: iscfxc, scf_temperature, scf_thermal_vxc !-LC- iscfxc

  use pz_mod
  use ksdt_mod

  implicit double precision (a-h, o-z)

  integer iphat(natx)
  real*8 rat(3,natx)
  integer iatph(0:nphx)
  real*8 xnatph(0:nphx)
  integer novr(0:nphx)
  integer iphovr(novrx,0:nphx)
  integer nnovr(novrx,0:nphx)
  real*8 rovr(novrx,0:nphx)
  real*8 folp(0:nphx), folpx(0:nphx)
  real*8 edens(251,0:nphx), edenvl(251,0:nphx)
  real*8 dmag(251,0:nphx)
  real*8 vclap(251,0:nphx)
  real*8 vtot (251,0:nphx), vvalgs (251,0:nphx)
  integer imt(0:nphx)
  integer inrm(0:nphx)
  real*8 rmt(0:nphx)
  real*8 rnrm(0:nphx)
  real*8, parameter :: big = 5000
  character*512 slog
  logical,save,allocatable :: lnear(:) ! BUG fix: JJK - 2/20/2015
  !logical,lnear(0:nphx)
  integer  inn(0:nphx)
  real*8 rnnmin(0:nphx)
  real*8, external :: dist

  !     work space for linear algebra
  real*8 ri(251)
  !      parameter (novp=40)
  complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
  integer ipiv(novp*(nphx+1)+1)
  logical EmptyCell(0:nph)


  ! Debug: FDV
  !     write(6,fmt='(a)') 'Entering istprm'
  ! Find muffin tin radii.  We'll find rmt based on norman prescription,
  ! ie, rmt(i) = R * folp * rnrm(i) / (rnrm(i) + rnrm(j)),
  ! a simple average
  ! based on atoms i and j.  We average the rmt's from each pair of
  ! atoms, weighting by the volume of the lense shape formed by the
  ! overlap of the norman spheres.
  ! NB, if folp=1, muffin tins touch without overlap, folp>1 gives
  ! overlapping muffin tins.
  !
  ! rnn is distance between sphere centers
  ! rnrm is the radius of the norman sphere
  ! xl_i is the distance to the plane containing the circle of the
  !    intersection
  ! h_i  = rnrm_i - xl_i is the height of the ith atom's part of
  !    the lense
  ! vol_i = (pi/3)*(h_i**2 * (3*rnrm_i - h_i))
  !
  ! xl_i = (rnrm_i**2 - rnrm_j**2 + rnn**2) / (2*rnn)

  EmptyCell=.false.  !KJ 6-09

  !     find rmt from rnrm only on first call of istprm (rmt(0)=-1)
  if(.NOT.ALLOCATED(lnear)) THEN
     allocate(lnear(0:nphx))
     lnear = .FALSE. ! JJK - added initialization of lnear since some compilers initialize it to true.
  END IF
  if (rmt(0).le.0.0) then
    do 10 iph=0,nph
      10  lnear(iph)=.false.
      do 140  iph = 0, nph
        voltot = 0
        rmtavg = 0
        inrm(iph) = ii(rnrm(iph))
        ! skip empty cells. rmt will be found for these later.
        if(iph.gt.0) then  !KJ 7-09 added this check to keep from exceeding array bounds!!
          IF(EmptyCell(iph)) GOTO 140
        endif
        if (novr(iph) .gt. 0)  then
          !           Overlap explicitly defined by overlap card
          rnear = big
          inters = mod(inters,6)
          !           use Norman prescription only in this case

          do 124  iovr = 1, novr(iph)
            rnn  = rovr(iovr,iph)
            inph = iphovr(iovr,iph)
            if (rnn .le. rnear) then
              rnear = rnn
              rnnmin(iph) = rnn
              inn(iph) = inph
            endif
            !              Don't avg if norman spheres don't overlap
            if (rnrm(iph)+rnrm(inph) .le. rnn)  goto 124
            voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
            voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
            rmttmp = rnn * folp(iph) * rnrm(iph) / (rnrm(iph) + rnrm(inph))
            ntmp = nnovr(iovr,iph)
            rmtavg = rmtavg + rmttmp*voltmp*ntmp
            voltot = voltot + voltmp*ntmp
            124       continue
        else

          iat = iatph(iph)
          rnear = big
          rmt(iph) = big
          ! Debug: FDV
          !           write(6,fmt='(a,i4,2f16.10)') 'iat: ', iat, rnear, rmt(iph)
          do 130  inat = 1, nat
            if (inat .eq. iat)  goto 130
            rnn = dist (rat(1,inat), rat(1,iat))
            inph = iphat(inat)
            if (rnn .le. rnear) then
              rnear = rnn
              rnnmin(iph) = rnn
              inn(iph) = inph
            endif
            !              Don't avg if norman spheres don't overlap
            if (rnrm(iph)+rnrm(inph) .lt. rnn)  goto 130

            if (inters.lt.6) then
              !                Norman prescription
              voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
              voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
              rmttmp = rnn * folp(iph) * rnrm(iph) / (rnrm(iph) + rnrm(inph))
              rmtavg = rmtavg + rmttmp*voltmp
              voltot = voltot + voltmp
            else
              !                Matching point prescription
              do i=inrm(iph),1,-1
                j=ii(rnn-rnrm(iph))
                if (vclap(i,iph).le.vclap(j,inph)) then
                  d1 = (vclap(i+1,iph)-vclap(i,iph))/(rr(i+1)-rr(i))
                  d2 =(vclap(j,inph)-vclap(j-1,inph))/(rr(j)-rr(j-1))
                  rmtavg = rr(i) + (vclap(j,inph)+d2*(rnn-rr(i)-rr(j))-vclap(i,iph)) /(d1+d2)
                  goto 127
                  !                    exit from the loop
                endif
              enddo
              127            continue
              if (rmtavg.lt.rmt(iph)) rmt(iph) = rmtavg
            endif
            130       continue
        endif

    !        special situation if rnrm is too close or larger than the nearest neighbor distance
    if (rnrm(iph).ge.rnear) lnear(iph) = .true.

    ! Debug: FDV
    !        write(6,fmt='(a,i4,3f16.10)'), 'rmt: ', iph, rmt(iph), rmtavg
    if (rmtavg .le. 0)  then
      write(slog,132) iat, iph
      call wlog(slog)
      132       format (' WARNING: NO ATOMS CLOSE ENOUGH TO OVERLAP ATOM',  i5, ',  UNIQUE POT', i5, '!!  ',                    &
      'Rmt set to Rnorman.  May be error in input file.')
      rmt(iph) = rnrm(iph)
    elseif(inters.lt.6) then
      !           Norman prescription
      rmt(iph) = rmtavg / voltot
      if (rmt(iph) .ge. rnear)  then
        call wlog(' Rmt >= distance to nearest neighbor.  Not physically meaningful.')
        call wlog(' FEFF may crash.  Look for error in ATOM list or OVERLAP cards.')
      endif
      if (rnrm(iph) .ge. rnear) then
        imax = ii(rnear) - 1
        !             begin until loop
        133            if (vclap(imax,iph).lt.vclap(imax+1,iph)) goto 134
        imax = imax-1
        goto 133
        !             end of until loop
        134          continue
        rmt(iph) = exp(xx(imax)) - 0.0001
      endif
    endif
    
    140 continue

    ! Now set rmt for empty cells to touching.
    do iph = 1, nph !KJ changed 0 to 1 to prevent stepping outside array bounds
      if(EmptyCell(iph)) then
        ! Find distance to nearest mt boundary.
        rmtmin = 1000.d0
        rnnmin(iph) = 1000.d0
        do iat = 1, nat
          rmt(iph) = 0.d0
          if(iat.ne.iatph(iph)) then
            ! if two empty cells are being compared. And neither is defined, use the half distance.
            rnn = dist(rat(:,iatph(iph)),rat(:,iat))
            if(rnn.lt.rnnmin(iph)) rnnmin(iph) = rnn
            if(rmt(iphat(iat)).lt.0.d0) then
              rmt(iph) = rnn/2.d0
              ! Otherwise, use dist - rmt(iphat(iat))
            else
              rmt(iph) = rnn - rmt(iphat(iat))
            end if
            if(rmtmin.gt.rmt(iph)) rmtmin = rmt(iph)
          end if
        end do
        rmt(iph)  = rmtmin*folp(iph)
        if(rnrm(iph).lt.rmt(iph)) rnrm(iph) = rmt(iph)*1.06d0
      end if
    end do

    !     set maximum value for folp(iph) if AFOLP is in use
    !     LMTO lore says no more than 15% overlap
    !     do 144 iph = 0, nph
    ! 144 folpx(iph) = 1.15
    !     already done in pot.f

    do 145 iph = 0, nph
    if (iafolp .gt. 0 ) then
      temp = 0.2 + 0.8 * rnrm(iph) / rmt(iph)
    else
      temp = 0.3 + 0.7 * rnrm(iph) / rmt(iph)
    endif
    if (temp.lt.folpx(iph)) folpx(iph) = temp
    temp = rnnmin(iph)/rmt(iph)/1.06d0
    if (temp.lt.folpx(iph)) folpx(iph) = temp
    temp = exp( -(novp-3)*0.05d0)
    !      make sure that with given folpx(iph) the construction
    !      of the overlapping matrix in movrlp will not fail
    if (lnear(iph)) then
      !           lnear=.true. only when hydrogens are present in the system.
      !           want to scale both rmt for iph and inn, so that overlapping
      !           matrix calculations will not fail
      temp = rnnmin(iph) / (rmt(iph)*1.05d0 + temp*rmt(inn(iph)))
      if (temp.lt.folpx(iph)) folpx(iph) = temp
      if (temp.lt.folpx(inn(iph))) folpx(inn(iph)) = temp
    else
      temp = (rnnmin(iph) - rnrm(iph))/ (temp*rmt(inn(iph)))
      if (temp.lt.folpx(inn(iph))) folpx(inn(iph)) = temp
    endif
    145 continue

    ! Debug: FDV
    !       write(6,*) ' Entering afolp'
    !       do iph=0,nph
    !         write(6,fmt='(a,i4,3f16.10)'), &
    !         'rmt, folpx, folp: ', iph, rmt(iph)*bohr, folpx(iph), folp(iph)
    !       end do
    !       stop
  endif
  !     end of finding rmt from rnrm on first call of istprm.


  !     Need potential with ground state xc, put it into vtot
  do iph = 0, nph
    call sidx (edens(1,iph), 250, rmt(iph), rnrm(iph), imax, imt(iph), inrm(iph))
    do i = 1, imax
      if (edens(i,iph).le.0) then
        if(mod(i,10).eq.0) then
          write(slog, 149) 'negative density ', iph,edens(i,iph), &
          ' - usually harmless precision error, but check DOS if it persists'
          149          format (a, i3,f9.3, a)
          call wlog(slog)
        endif
        rs = 100
        xmag=1.0
      else
        rs = (edens(i,iph)/3)**(-third)
        !     spin dependent xc potential for ground state from Von Barth, Hedin
        !     J.Phys.C:Solid State Phys., 5, 1629 (1972).
        !     xmag/2 -fraction of spin up or down, depending on sign in renorm.f
        !     put xmag = 1.0 to calculate cmd with external potential difference
        xmag = 1.0 + idmag*dmag(i,iph)
      endif
      !           wrong for ferromagnets, need to overlap dmag(i)

      !            !-LC- checking if iscfxc has the right value (probably should go
      !            ! in the read input section, move it there)
      !            if ( iscfxc .lt. 1 .or. iscfxc .gt. 4) then
      !               if (scf_temperature .gt. 0 .and. scf_thermal_vxc .gt. 0) then
      !                 iscfxc = 3
      !               else
      !                 iscfxc = 1
      !            !   write(slog,149) "Error: iscfxc should take one of the values &
      !            !    1 (for vBH), 2 (for PZ), 3 (for PDW) or 4 (for KSDT) .. stopping"
      !                 iscfxc = 1
      !               endif
      !            endif


      ! if (scf_temperature .gt. 0 .and. scf_thermal_vxc .gt. 0) then
      select case(iscfxc) !-LC- two posibilities PDW or KSDT
        ! Modified by tts
        case(21)
          ! Finite T Vxc from Perrot, Dharma-Wardana 1984 paper
          ! N B: this isn't spin dependent!
          call pdw_vxc(rs, scf_temperature/hart, vvbh)
        case(22)
          !-LC- KSDT (No Spin)
          ! LDA XC free-energy parameterization from Monte Carlo data (unpol/pol)
          call ksdt_vxc(rs, scf_temperature/hart, vvbh)
        case(11)
          ! vvbh from Von Barth Hedin paper, 1971
          call vbh( rs, xmag, vvbh )
        case(12)
          call pz_vxc( rs, vvbh )
        case default
          write(slog, 911) 'Invalid choice of xc potential ', iscfxc
          911          format (a, i3)
          call wlog(slog)
      end select

      ! else

      ! select case(iscfxc) !-LC- two posibilities vBH or PZ
      ! case(11)
      ! vvbh from Von Barth Hedin paper, 1971
      ! PRINT*, "Using: VBH"
      ! call vbh( rs, xmag, vvbh )
      ! case(12)
      !-LC-  PZ (No Spin)
      ! PRINT*, "Using: PZ"
      ! call pz_vxc( rs, vvbh )
      ! end select
      ! end if

      vtot(i,iph) = vclap(i,iph) + vvbh
      ! if (mod(ixc,10).eq.5) then
      IF ((ixc.EQ.5).OR.(ixc.EQ.15)) THEN
        rsval = 10.0
        if (edenvl(i,iph) .gt. 0.00001) rsval = (edenvl(i,iph)/3)**(-third)
        if (rsval.gt.10.0) rsval = 10.0
        xmagvl = 1.0 + idmag * dmag(i,iph) * edens(i,iph) / edenvl(i,iph)
        call vbh(rsval,xmagvl,vvbhvl)
        vvalgs(i,iph) = vclap(i,iph) + vvbhvl
      ! elseif (mod(ixc,10) .ge. 6) then
      ELSEIF (ixc.EQ.9) THEN
        if (edens(i,iph).le.edenvl(i,iph)) then
          rscore =101.0
        else
          rscore = ((edens(i,iph)-edenvl(i,iph)) / 3)**(-third)
        endif
        rsmag = (edens(i,iph)*(1+idmag*dmag(i,iph)) / 3)**(-third)
        xfmag = fa/rsmag
        call edp(rscore,xfmag,vrdh)
        vvalgs(i,iph) = vclap(i,iph) + vvbh - vrdh
      else
        vvalgs(i,iph) = 0.d0
      endif
    enddo !i
  enddo !iph


  !     What to do about interstitial values?
  !     Calculate'em for all atoms, print'em out for all unique pots along
  !     with derivative quantities, like fermi energy, etc.
  !     Interstitial values will be average over all atoms in problem.

  !     rnrmav is averge norman radius,
  !     (4pi/3)rnrmav**3 = (sum((4pi/3)rnrm(i)**3)/n, sum over all atoms
  rnrmav = 0
  xn = 0
  !     volint is total interstitial volume
  volint = 0
  do iph = 0, nph
    rnrmav = rnrmav + xnatph(iph) * rnrm(iph)**3
    volint=volint-xnatph(iph) * rmt(iph)**3
    xn = xn + xnatph(iph)
  enddo
  if (totvol.le.0.0d0) then
    volint=4*pi/3 *(volint+rnrmav)
  else
    volint=4*pi/3 *volint + totvol
  endif
  !     volume of lenses from overlapping mt spheres is added in movrlp.
  rnrmav = (rnrmav/xn) ** third

  rs = 0
  vint   = 0
  rhoint = 0
  rsval = 0

  call movrlp(nph, nat, iphat, rat, iatph, xnatph, novr, iphovr, nnovr, rovr, &
  imt, rmt, rnrm, ri, lnear, cmovp,ipiv, volint,inters)

  !     If no contribution to interstitial from any atom, die.
  if (volint .le. 0)  then
    call wlog(' No interstitial density.  Check input file.')
    call par_stop('ISTPRM')
  endif
  !     find interstitial density
  call ovp2mt(nph, edens, 0, qtotel, ri, xnatph, lnear, inrm, imt, rnrm, rmt, cmovp, ipiv, rhoint,inters)
  rhoint = 4*pi * rhoint / volint

  ! if (ixc.ge.5) then
    IF ((ixc.EQ.5).OR.(ixc.EQ.9).OR.(ixc.EQ.10).OR.(ixc.EQ.13).OR.(ixc.EQ.15)) THEN ! Replace ixc.ge.5
    !        find valence potential inside mt sphere (vintvl -dummy)
    call ovp2mt(nph, vvalgs, 1, qtotel, ri, xnatph, lnear, inrm, imt, rnrm, rmt, cmovp, ipiv, vintvl,inters)
  endif

  !     find potential inside mt sphere and vint
  call ovp2mt(nph, vtot, 1, qtotel, ri, xnatph, lnear, inrm, imt, rnrm, rmt, cmovp, ipiv, vint,inters)

  if (vint.ge.xmu) then
    call wlog(' WARNING:interstitial level found above Fermi level.')
    call wlog(' Results may be unreliable. See manual for details.')
    vint = xmu - 0.05d0
    call ovp2mt(nph, vtot, 2, qtotel, ri, xnatph, lnear, inrm, imt, rnrm, rmt, cmovp, ipiv, vint,inters)
  endif
  call fermi (rhoint, vint, xmunew, rs, xf)

  return
end


double precision function calcvl (r1, r2, r)
  use constants, only: pi
  implicit double precision (a-h, o-z)
  xl = (r1**2 - r2**2 + r**2) / (2*r)
  h = r1 - xl
  calcvl = (pi/3) * h**2 * (3*r1 - h)
  return
end
