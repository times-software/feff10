!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sfconvsub.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2012/12/11 23:20:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sfconv(ekp,mu,gammach,npts2,wpts2,xchi,npts1,wpts1,    &
     &                spectf,weights,cchi,phase,iasym,icut,intout,omp)
! Convolutes the array xchi (signal) with the array spectf (spectral
! function or asymmetric broadening), where spectf
! has a delta function contribution of magnitude |weights(1)+i weights(2)|
! + weights(3).  Complex phase of weights(1) + i weights(2) is given 
! to overall phase of the otherwise real valued array spectf.
! Input: ekp - Photoelectron energy neglecting collective excitations.
!        mu - Chemical potential, position of edge.
!        gammach - Core hole lifetime.
!        npts2 - Dimension of signal array xchi.
!        wpts2 - Energy grid of signal.
!        xchi - Signal array.
!        npts1 - Dimension of spectral function array spectf.
!        wpts1 - Energy grid of spectral function.
!        weights - Array of weights of components of the spectral
!            function.  Only the delta function weights, weights(1)
!            and weights(3), and the imaginary part of the extrinsic
!            delta function, weights(2) are needed.
!        cchi - The magnitude of the convoluted signal.
!        phase - The phase of the convoluted signal.
!        iasym - set to 1 to include quasiparticle phase as an 
!            asymmetric 1/omega term to the etrinsic satellite
!            rather than as a complex spectral weight.  This is
!            necessary when convoluting with a real valued
!            function or one whose imaginary part is not known.
!        icut - set to 0 to avoid truncating the spectral function
!            at energies where there is insufficient energy to create
!            excitations.
!        intout - Flag to write diagnostic file, the running integration
!            of the convolution.
!        omp - Plasma frequency.
      implicit none
      integer npts1,npts2,i,j,intout,iasym,icut
      double precision wpts2(npts2),xchi(npts2),wpts1(npts1),           &
     &                 spectf(npts1),weights(8),cchi,phase,ekp,mu,      &
     &                 gammach,spectf2(npts1),wtmp,qpr,omp
!       spectf2 - The spectral function with a cutoff to prevent 
!           excitation energies higher than the total available 
!           energy in the system (ekp-mu).
!       wtmp - The new quasiparticle weight in going 
!           from spectf to spectf2.
!       qpr - Quasiparticle reduction, from renormalizing 
!           spectral function.
      double precision xrcchi,xicchi,xnorm,dw,w,www,xfact,              &
     &                 xrchi,am,phasez,efrac,pi,eV,store,               &
     &                 amp,del,lam
!       xrcchi - Real part of convoluted signal.
!       xicchi - Imaginary part of convoluted signal.
!       xnorm - Total weight of spectral function with cutoff,
!           used to normalize the spectral function.
!       dw - Energy interval.
!       w - Exitation energy (omega)
!       www - Quasiparticle energy, available energy minus excitation energy.
!       xfact - Used to store intermediate calculations.
!       xrchi - Convoluted signal before quasiparticle phase shift.
!       am - Magnitude of quasiparticle weight.
!       phasez - Quasiparticle phase.
!       efrac - Fractional distance from one energy grid point to
!           the next, used in interpolation.
!       pi - Ratio of circumference to diameter of a circle in
!           euclidian geometry.
!       eV - Unit conversion to electron volts.
!       store - Dummy variable used to store a value that will be 
!           overwritten but needs to be used again.
!       amp - "Amplitude", signal value at first data point,
!           used for extrapolating signal to lower energies to avoid
!           artifacts.
!       lam - "Lambda", decay constant for extrapolating signal below the
!           first data point.
!       del - "Delta", energy difference used in extrapolation of signal
!           below the first data point.
      parameter (eV=1.d0/27.21160d0)

      pi=dacos(-1.d0)
      xrcchi=0.d0
      xicchi=0.d0
!      iasym=0
! Extrinsic delta function (quasiparticle) weight.
      if (iasym.eq.1) then
        am=weights(1)
      else
        am=dsqrt(weights(1)**2+weights(2)**2)
      endif
! Quasiparticle phase.
      if (weights(1).ne.0.d0.and.iasym.ne.1) then
        phasez=datan(weights(2)/weights(1))
      else
        phasez=0.d0
      endif
! Reduce quasiparticle weight in accord with the same energy cutoff
! constraints used for the rest of the spectral function.
      if (icut.eq.0) then
        qpr=1.d0
      elseif (ekp-mu.ne.0.d0) then
        qpr=datan2(gammach,mu-ekp)/pi
      else
        qpr=0.5d0
      endif
      wtmp=qpr*(am+weights(3))
      xnorm=wtmp
! Cut off that portion of spectral function with more energy than is
! available in the system.
      do i=1,npts1
        if (i.eq.1) then
          dw=wpts1(2)-wpts1(1)
        elseif (i.eq.npts1) then
          dw=wpts1(npts1)-wpts1(npts1-1)
        else
          dw=(wpts1(i+1)-wpts1(i-1))/2.d0
        endif
        w=wpts1(i)
        www=ekp-w
! Construct weight function, cutoff at chemical potential, width 
! equal to core hole lifetime.
        if (icut.eq.0) then
          xfact=1.d0
        elseif (www-mu.ne.0.d0) then
          xfact=datan2(gammach,mu-www)/pi
        else
          xfact=0.5d0
        endif
! Multiply spectral function by cutoff weight function.
        if (icut.eq.0) then
          spectf2(i)=spectf(i)
        elseif (w.ge.0.d0) then
          spectf2(i)=spectf(i)*xfact
        else
          spectf2(i)=max(0.d0,spectf(i)*xfact)
        endif
        if (iasym.eq.1) then
          spectf2(i)=spectf2(i)-qpr*(weights(2)/(pi*am*dw))             &
     &      *log(((w+dw/2.d0)**2+(3.d0*dw)**2)                          &
     &         /((w-dw/2.d0)**2+(3.d0*dw)**2))                          &
     &         *exp(-((w)/(2*omp))**2)/2.d0
        endif
! Integration to find total spectral weight.
        xnorm=xnorm+spectf2(i)*dw
      enddo
! Main convolution loop.
      do i=1,npts1
        if (i.eq.1) then
          dw=wpts1(2)-wpts1(1)
        elseif (i.eq.npts1) then
          dw=wpts1(npts1)-wpts1(npts1-1)
        else
          dw=(wpts1(i+1)-wpts1(i-1))/2.d0
        endif
        w=wpts1(i)
        www=ekp-w
        xrchi=0
        if (www.gt.wpts2(npts2)) then
! Extrapolate signal to avoid artifacts.
          xrchi=xchi(npts2)
        elseif (www.le.wpts2(1)) then
! Extrapolate signal to avoid artifacts.
          amp=xchi(1)
          del=mu-wpts2(1)
          lam=del**2/(pi*dabs(amp)*(del**2+gammach**2))
          xrchi=amp*exp(lam*(www-wpts2(1)))
        else
! Interpolate signal onto spectral function energy.
          do j=1,npts2-1
            if (www.gt.wpts2(j).and.www.le.wpts2(j+1)) then
              efrac=(www-wpts2(j))/(wpts2(j+1)-wpts2(j))
              xrchi=xchi(j)+(xchi(j+1)-xchi(j))*efrac
              goto 50
            endif
          enddo
        endif
50      continue
		if (i.gt.1 .and. i.lt.npts1) then   ! KJ 12-2012 to prevent exceeding array bounds
        if ((w+wpts1(i-1))/2.d0.lt.0.d0.and.                            &
     &      (w+wpts1(i+1))/2.d0.ge.0.d0) then
! Add delta function contribution.
          xrcchi=xrcchi+wtmp*xrchi
        endif
		endif
! Convolution integration.
        xrcchi=xrcchi+spectf2(i)*dw*xrchi
! Print diagnostic files of running integration if requested.
        if (intout.ne.0) then
		  if (i.eq.1 .or. i.eq.npts1) then   ! KJ 12-2012 to prevent exceeding array bounds
            write(28,500) www,xrcchi/xnorm,xrchi,spectf(i),             &
     &            spectf2(i)/xnorm
          elseif ((w+wpts1(i-1))/2.d0.lt.0.d0.and.                          &
     &        (w+wpts1(i+1))/2.d0.ge.0.d0) then
            write(28,500) www,xrcchi/xnorm,xrchi,am+weights(3),         &
     &           spectf2(i)/xnorm
          else
            write(28,500) www,xrcchi/xnorm,xrchi,spectf(i),             &
     &            spectf2(i)/xnorm
          endif
        endif
      enddo
! Add overall phase equal to quasiparticle phase.
      store=xrcchi
      xrcchi=store*dcos(phasez)-xicchi*dsin(phasez)
      xicchi=xicchi*dcos(phasez)+store*dsin(phasez)
      xrcchi=xrcchi/xnorm
      xicchi=xicchi/xnorm
      cchi=dsqrt(xrcchi**2+xicchi**2)
      phase=datan2(xicchi,xrcchi)

 500  format(1x,5(e12.5,1x))
 501  format(1x,7(e10.3,1x))
      return
      end
