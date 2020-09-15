!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: correorb.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2012/06/29 01:05:24 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine correorb( iz, ihole, rmt, jri, dx, ri, p2f,edge,       &
     &       vxc, dgcn, dpcn, adgc, adpc, eorb,                         &
     &       neg, eng, rhoj, kappa, norbp, iph) !KJ iph
!     correct energy of the orbital for solid state potential
!     for deep core orbitals; output eng(1) (neg = 1, rhoj(1) = xx)
!     for valence band orbital create projected DOS on the 
!     orbital:  neg = number of energy points, eng - energy grid
!              rhoj(i) - projected DOS (sum_i rhoj(i) = xx)
!     coded by a.ankudinov 2004

!     input:
!        rmt     muffin-tin radius
!        jlast   last point for integration of Dirac eq.
!        jri     first interstitial grid point (imt + 1)
!        edge    shifted Fermi level
!        dx      dx in loucks' grid (usually .05)
!        ri(nr)  loucks' position grid, r = exp ((i-1)*dx - 8.8)
!        vxc(nr) coulomb+xc potential for total density
!        dgcn(dpcn) large(small) dirac components for 'iph' atom
!        adgc(adpc) their development coefficients
      use errorfile
      use dimsmod, only: nrptx, nex
	  use constants
      implicit double precision (a-h, o-z)
!Changed the dimensions to 40 to account for superheavy elements. Pavlo Baranov 07/2016
!     INPUT:
      integer,intent(in) :: iz, ihole, jri, iph
      real*8,intent(in) :: rmt,dx,ri(nrptx), p2f, edge
      complex*16,intent(inout) :: vxc(nrptx)
!     all atoms' dirac components and their development coefficients
!     orbitals(dgc,kappa), their energy(eorb) 
      real*8, intent(in) :: dgcn(nrptx,41), dpcn(nrptx,41)
      real*8, intent(in) :: adgc(10,41), adpc(10,41)
      real*8,intent(in) :: eorb(41)
      integer,intent(in) :: kappa(41)
!     OUTPUT: energy points and weights
      integer,intent(out) :: neg(41), norbp
      real*8,intent(out) :: rhoj(nex,41), eng(nex,41)
!     local variables
      complex*16 p2p, p2
      logical done


      do i = jri+1,nrptx
         vxc(i)=vxc(jri+1)
      enddo

!     calculate initial photoelectron orbital using lda
      ic3 = 0

      do iorb = 1, 41
        if ( eorb(iorb) .ge. 0.d0 ) cycle
        norbp = iorb
        kfin = kappa(iorb) 
        xx = 2.d0 * abs(kfin)
        neg(iorb) = 1
        eng(1,iorb) = eorb(iorb) 
        vint = edge - p2f
        eng(1,iorb) = eorb(iorb) - vint
        rhoj(1,iorb) = xx / 2 /abs(kfin)


        do 250 ib = 1,2
        p2p = eng(1,iorb)
!       return point with different p2p, if not done
        done = .false.

!       use broadening to remove possible resonant behaviour
!       when E_i - E_j = omega for occupied i and j
        gamb = max( 3.d0/hart, (p2f - dble(p2p))/5000)
        if (ib.eq.1) gamb = max( 3.d0/hart, (p2f - dble(p2p))/50)
!       gamb = max( 0.01d0/hart, (p2f - dble(p2p))/10000)
!       if (ib.eq.1) gamb = max( 0.01d0/hart, (p2f - dble(p2p))/100)


 100    continue
        p2p = dble(p2p) + coni*gamb
           
        x1 = -1
        x2 = 0
        x3 = 1
        jlast = jri+6

        itest = 0
        !if(iorb.eq.19) itest=1
        if (itest.eq.1 ) then
          do 120 ie = -100, 100
           p2 = p2p + ie*gamb/3
          call cdos( iorb,p2, kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos2, iph) !KJ iph
          write( 57, 130) dble(p2), dos2
 120      continue
 130      format (2f15.4)
          call WipeErrorfileAtFinish
          stop
        endif

        call cdos( iorb,p2p, kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos2, iph) !KJ iph
      
        p2 = p2p - gamb
        call cdos( iorb,p2,  kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos1, iph) !KJ iph
        p2 = p2p + gamb
        call cdos( iorb,p2,  kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos3, iph) !KJ iph
        ncdos = 3

        itry=0
        isituation=0

 200    continue
 itry=itry+1
 if(itry.gt.20) then
 !stop
    gamb=gamb*2.d0
    itry=0
    goto 100
  endif
        t1 = (x3-x2)/dos1
        t2 = (x1-x3)/dos2
        t3 = (x2-x1)/dos3
        aqd = t1 + t2 + t3
        bqd = (x3+x2)*t1 + (x2+x1)*t3 + (x1+x3)*t2

        if (aqd.eq.0) then
!         1/dos is linear or constant
!         increase broadening and try abain
          gamb = 2* gamb
          if (gamb.gt. abs(dble(p2p))) stop ' error in correrr.f'
    !      call wlog('going back to 100')
          goto 100
        endif

        x0 = bqd/aqd/2
        if ((dos1.le.0 .or. dos2.le.0 .or. dos3.le.0) .or. (1/dos1.gt.500 .and. 1/dos2.gt.500 .and. 1/dos3.gt.500) ) then
          p2p = p2p + gamb
          dos1 = dos2
          dos2 = dos3
          p2 = p2p + gamb
          call cdos(iorb, p2, kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos3, iph) !KJ iph
          ncdos = ncdos + 1
          isituation=1
        elseif (x0.le.x3 .and.x0.ge.x1) then
          done = .true.
          isituation=2
        elseif (x0.lt.x1) then
          p2p = p2p - gamb
          isituation=3
          dos3 = dos2
          dos2 = dos1
          p2 = p2p - gamb
          call cdos(iorb, p2, kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos1, iph) !KJ iph
          ncdos = ncdos + 1
          isituation=4
        elseif (x0.gt.x3) then
          p2p = p2p + gamb
          dos1 = dos2
          dos2 = dos3
          p2 = p2p + gamb
          isituation=5
          call cdos(iorb, p2, kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos3, iph) !KJ iph
          ncdos = ncdos + 1
        endif
 !       write(*,*) itry,isituation,x0,x0*gamb,p2p
        if (.not.done) goto 200

        eng(1,iorb) = dble(p2p) + x0*gamb

 250  continue
!       eigenenergy eng above is relative to vint; to get energy relative to absolute zero use line below
        eng(1,iorb) = eng(1,iorb) +vint
      end do

!     spread energy levels into bands above Ecut
      ecut = 0
      itest = 1
      do iorb = 1, norbp
      if (eng(1,iorb) .gt. ecut .and. itest.eq.0) then
        kfin = kappa(iorb) 
        xx = rhoj(1,iorb)
        gamb = max( 3.d0/hart, (p2f - dble(p2p))/50)
!       make energy grid 
        emin = ecut 
        emax = p2f
        de = gamb / 3
        ne = nint((emax-emin)/de)
        if (ne.le.1) ne = 2
        if (ne.gt.nex) ne = nex
        de =  (emax-emin) / (ne-1)
        neg(iorb) = ne
        sumx = 0
        do ie = 1, ne
          eng(ie,iorb) = emin + (ie-1)*de
          p2 = eng(ie,iorb) + coni*gamb
          call cdos(iorb, p2, kfin,  rmt, jri, jlast, dx, ri, vxc, iz, dgcn, dpcn, adgc, adpc, dos, iph) !KJ iph
          rhoj(ie,iorb) = dos
          sumx = sumx + dos
        enddo
        xx = xx/sumx
        sumx = 0
        do ie = 1,ne
          rhoj(ie,iorb) = rhoj(ie,iorb) * xx
          sumx = sumx + rhoj(ie,iorb)
        enddo
!       print*, 'checksum', sumx
!       fix later: don't need xnval to check sum?

      endif
      enddo

! Debug: Fer
      !write(6,fmt='(a,2f20.10)') 'eorb, eng: ', eorb(1), eng(1,1)

      return
      end
