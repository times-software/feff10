!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xscorr.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xscorr(ispec, emxs ,ne1, ne, ik0, xsec, xsnorm, chia,  &
     &                  vrcorr, vicorr, cchi)
!     convolute xmu(E)=xsec+xsnorm*chia with lorentzian using
!     calculations in the complex energy plane

!     Input: ispec - type of spectroscopy
!       emxs - complex energy grid
!       ne1 - number of points on horizonatal axis
!       ne - total number of points (ne-ne1) points on vertical axis
!       ik0 - Fermi level index on horizontal axis
!       xsec, xsnorm, chia - give function f in complex energy plain
!           xmu(ie) = xsec + xsnorm*chia
!       vrcorr = correction for the shift of the Fermi level
!       vicorr = 0 (disabled)
!     Output: cchi(w) - result of convolution for w = dble(emxs)
!       cchi(w) = \int_C dE xmu(E)*xloss/pi/((E-w)**2+xloss**2) =
!       xmu(w+i*xloss)* [1/2+atan(w-efermi/xloss)/pi] +
!       \int_C dE ff(E)*xloss/pi/((E-w)**2+xloss**2)
!       where ff(E)=xmu(E)-xmu(w+i*xloss) for w<efermi we use
!       xmu(efermi+i*xloss) instead of xmu(w+i*xloss);
!       contour C starts at efermi, goes vertically to efermi+i*xloss
!       and then goes horizontally to infinity + i*xloss
      use dimsmod, only: nex
	    use constants
      implicit double precision (a-h, o-z)

      dimension  xsnorm(nex), omega(nex)
      complex*16 emxs(nex), xsec(nex), chia(nex), cchi(nex)
      complex*16 xmu(nex), aa, bb, c1, f1, f2, ff(nex), xmu0
      parameter (eps4 = 1.0d-4)
      complex*16 ec(nex), fc(nex), e1,e2, z1,z2, corr, residue, corr2, leg3
      complex*16 lorenz
      external lorenz, astep
      real(8) dummy
      ne2 = ne-ne1
      efermi = dble(emxs(ne))
      xloss = dimag(emxs(1))

      OPEN(UNIT=150, file="prexmu.dat", status='REPLACE')
      OPEN(UNIT=151, file='residue.dat', status='REPLACE')
      OPEN(UNIT=152, file='contour.dat', status='REPLACE')
      OPEN(UNIT=155, file='curve.dat', status='REPLACE')
      OPEN(UNIT=149, file='raw.dat', status='REPLACE')
      WRITE(149,*) "Temperature (Hatree) = 0"
      WRITE(149,*) "Electronic Temperature (eV) = 0"
      WRITE(149,*) "xloss = ", xloss*hart, " eV"
      WRITE(149,*) "efermi = ", efermi*hart, " eV"
      WRITE(149,*) "Number of poles = 0"
      WRITE(149, *) 'Omega(Hart)    Re CCHI     Im CCHI   1-Fermi   Re xmu0    Im xmu0'

    ! WRITE(148, *) "IK0 = ", ik0
!     xmu - analytic function in complex energy plain
      !WRITE(149, *) 'Omega(Hart)   Re xmu0    Im xmu0'
      do  ie = 1,ne
        xmu (ie) = xsec(ie) + xsnorm(ie)*chia(ie)
      enddo

!     real frequencies
      do ie = 1, ne1
        omega(ie) = dble(emxs(ie))
      enddo

      if (abs(vrcorr).gt.eps4) then
        !       account for the fermi level shift
        bb = xmu(ik0)
        efermi = efermi - vrcorr
        call terpc(omega, xmu ,ne1, 1, efermi, bb)

        !       shift the vertical axis
        do ie = 1, ne2
          emxs(ne1+ie) = emxs(ne1+ie) - vrcorr
        enddo

        !       rescale values on vertical axis
        bb = bb/xmu(ik0)
        do ie = ne1+1, ne
          xmu(ie) = xmu (ie) * bb
        enddo
      else
        bb = 1
      endif

      ! construct the integration countur C
      nc = 0
      !  start with points on vertical axis below xloss
      do ie = 1,ne2
        if (dimag(emxs(ne1+ie)).lt.xloss) then
          nc = nc+1
          ec(nc) = emxs(ne1+ie)
          fc(nc) = xmu(ne1+ie)
          WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
        endif
      enddo
      !  add corner at efermi + xloss*i
      nc = nc+1
      ic0 = nc
      if (abs(vrcorr).gt.eps4) then
        ec(nc) = efermi + coni*xloss
        fc(nc) = bb * xmu(ik0)
        WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
      else
        ec(nc) = emxs(ik0)
        fc(nc) = xmu(ik0)
        WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
      endif
      ! add points on horizontal axis above efermi
      if (ispec.ne.2) then
        do ie = 1,ne1
          if (dble(emxs(ie))-efermi.gt.eps4) then
            nc = nc+1
            ec(nc) = emxs(ie)
            fc(nc) = xmu(ie)
            WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
          endif
        enddo
      else
        ! ispec=2 - emission calculations- need points below E_fermi
        do ie = ne1,1,-1
          if (efermi-dble(emxs(ie)).gt.eps4) then
            nc = nc+1
            ec(nc) = emxs(ie)
            fc(nc) = xmu(ie)
            WRITE(155, '(20E20.10E3)')  DBLE(ec(nc)), DIMAG(ec(nc)), DBLE(fc(nc)), DIMAG(fc(nc))
          endif
        enddo
      endif
!     endo of countour construction

!     cycle over frequency points
      do ie = 1, ne1
        if (omega(ie).ge.efermi) then
          xmu0 = xmu(ie)
          if (ispec.eq.2) xmu0 = xmu(ik0)*bb
        else
          xmu0 = xmu(ik0)*bb
          if (ispec.eq.2) xmu0 = xmu(ie)
        endif
        e1 = omega(ie) + coni*xloss
        e2 = omega(ie) - coni*xloss
        do ic = 1, nc
          ff(ic) = fc(ic) - xmu0
        enddo
        dele = omega(ie) - efermi
        cchi(ie) = xmu0 * astep( xloss, dele)
        dummy = astep( xloss, dele)
        if (ispec.eq.2) dummy = 1-dummy
        if (ispec.eq.2) cchi(ie) = xmu0 - cchi(ie)
        WRITE(149, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie)), DBLE(dummy), DBLE(xmu0), DIMAG(xmu0)
        !! WRITE(148, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie))

        corr = 0

        if (abs(dele).lt.eps4) dele = 0.0d0
        w1 = dimag(ec(1))
        w2 = dimag(ec(2))
        w3 = dimag(ec(3))
        ip =0

        ! add half matsubara pole contribution
        ! equivalent to integral from efermi to efermi+i*w1
        corr = corr + lorenz(ip,xloss,w1,dele)*ff(1) *coni*w1 !missing a factor of pi
        leg3 = lorenz(ip,xloss,w1,dele)*ff(1) *coni*w1
        WRITE(151, '(20E20.10E3)') DBLE(omega(ie)), DBLE(leg3), DIMAG(leg3)

        nc0 = 0 ! JK - set this to 0 for now. 9/2009
        if (nc0.gt.3) then
        !     add sommerfeld correction (correction for derivative)
        !       corr = corr + coni * w1**2 / 6   / (w3-w2) *
        !  2   (lorenz(ip,xloss,w3,dele)*ff(3)-lorenz(ip,xloss,w2,dele)*ff(2))
        endif

        ! cycle over contour points
        do ic = 1,nc-1
          ! perform integration over contour from efermi+i*2*w1 to efermi+i*xloss;
          ! linear interpolation of ff between  z1 and z2
          z1 = ec(ic)
          z2 = ec(ic+1)
          !         if (ic.eq.1) z1 = efermi+coni*2*w1
          f1 = ff(ic)
          f2 = ff(ic+1)
          ! if (ic.eq.1) f1 = (f1*(z2-z1) + f2*(z1-ec(ic))) / (z2-ec(ic))
          ! add correction from pole above real axis
          aa = 0
          if (abs(z1-e1).gt.eps4 .and. abs(z2-e1).gt.eps4) then
            aa = log((z2-e1)/(z1-e1)) *(f1*(z2-e1)+f2*(e1-z1))
            ! z1 or z2 equal to e1; in this case contribution to corr is exactly zero
          endif
          ! second pole
          aa = aa - log((z2-e2)/(z1-e2)) *(f1*(z2-e2)+f2*(e2-z1))
          corr = corr + aa/ (z2-z1) /2/pi/coni
        enddo

        ! end of cycle over contour points
        if (ispec.eq.2) corr = -corr
!       if (ispec.eq.2) corr = 0

        cchi(ie) = cchi(ie) +  corr
        WRITE(152, '(20E20.10E3)') DBLE(omega(ie)), DBLE(corr), DIMAG(corr)
        WRITE(150, '(20E20.10E3)') DBLE(omega(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie))
        ! return the result of convolution minus bare value
        cchi(ie) = cchi(ie) - xmu(ie)
      enddo
      CLOSE(151)
      CLOSE(152)
      CLOSE(155)
      CLOSE(150)
      CLOSE(149)
!     end of cycle over frequency points

!     restore the input energy mesh
      if (abs(vrcorr).gt.eps4) then
        do  ie = ne1+1, ne
          emxs(ie) = emxs(ie) + vrcorr
        enddo
      endif

      !open(unit=148, file="tunsheng.dat", status='REPLACE')
      !write(148, *) "Xloss = ", xloss, "EFermi = ", EFermi
      !write(148, *) "Ne1 = ", ne1, " and Ne = ", ne
      !write(148, *) "Re EXMS   Im EXMS   Re XMU    Im XMU     Re CCHI    Im CCHI"
      !do  ie =  ne
      !  write(148, '(20E20.10E3)') ie, DBLE(emxs(ie)), DIMAG(emxs(ie)), DBLE(xmu(ie)), &
      !            & DIMAG(xmu(ie)), DBLE(cchi(ie)), DIMAG(cchi(ie))
      !enddo
      !close(148)
      return
      end

      complex*16 function lorenz (ifp, xloss, w, dele)
      use constants
      implicit double precision (a-h, o-z)
!     ifp is dummy now. correspond to ifp=0 in old code
!     can remove it and change calls to lorenz in other routines

      lorenz = xloss /pi / (xloss**2+(coni*w-dele)**2)

      return
      end

      double precision function astep ( xloss, dele)
      use constants
      implicit double precision (a-h, o-z)
      !  1 - fermi*cauchy
      astep = 0.5d0 + atan(dele/xloss) /pi
      if (astep.lt.0.d0) astep = 0.d0
      if (astep.gt.1.d0) astep = 1.d0

      return
      end
