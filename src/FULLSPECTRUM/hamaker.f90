
      subroutine hamaker (omega,eps,nepts,epsi) 
! This subroutine takes a full epsilon spectrum pased in eps and
! outputs the  dielectric function on the energy grid going up the Im
! axis.  Now, omega is the energy grid, and nepts is the number of
! filled points in the arrays omega, eps, and epsi.

! This routine works with energy in any units; it returns
! the integral with respect to w' of (w'eps2(w')+w*eps1(w'))/(w'**2+w**2)
! over the energy interval covered by the energy grid omega.

! We should add 1/w**2 tail to the integral.

      use constants

      implicit none
      include 'HEADERS/params.h'

      double precision twopi,re,hc,e1(fullpts), n, K,e22,f22,           &
     &                 f1(fullpts), c, c2, wp2, tau, gam,               &
     &                 chi(fullpts), e2(fullpts), f2(fullpts),          &
     &                 einv1, einv2,x0(fullpts), chi22, dele1, icept

      real*8 omega(fullpts), klow
      real*8 eps1x0(fullpts),ysngl,x,khigh,x0real(fullpts)
      integer nn, i, j, nepts,ios
      double precision enrgy(fullpts), na, hc2re,slope,r,               &
     &             fpi,sigma, asig, rf,lfac,y,phi,e2drud,               &
     &             e2feff,theta,trnbeg, trnend,norm,ndrude
      double precision toobig
      complex*16 ref, ncmplx, epsi(fullpts)
      complex eps(fullpts)
      character line*75, slog*100
      logical drude,drudfl,flag

      parameter(na=6.0221367e-1)
      parameter(c=2.997925e+18)
      parameter(c2=c**2)
      parameter(hc = 12398.424468)
      parameter(re = 2.8179409187e-5)
      parameter(hc2re = .698760552714)
      parameter(twopi = 2.*pi)
      parameter(fpi = 4.*pi)
      parameter(drude=.false.) !true if drude term included
      parameter(drudfl=.false.) !true if drude term written to drude.dat
      parameter(tau=0.8e-14) !drude lifetime in s
      parameter(gam=hc/(twopi*tau*c)) !drude width in eV
      parameter(ndrude = 1.125) !number of electrons modeled in drude term
      parameter(trnbeg=0.066) !use drude below this energy in eV
      parameter(trnend=0.566) !use FEFF only above this energy in eV
      parameter(toobig=1d8) !use FEFF only above this energy in eV


      !double precision copies of the spectra
      do i=1, nepts
        enrgy(i)=omega(i)
        e2(i)=aimag(eps(i))
        e1(i)=dble(eps(i))
      enddo

      call wlog('Hamaker transform calculated up to...')

      !loop over the intervals between energy points. In each cycle we
      !compute the value of the KK integral at the frequency at the
      !midpoint x0(i) of the interval between enrgy(i) and enrgy(i+1).
      do 33 i = 1, nepts-1 

          !find midpoint of current interval. This is the frequency at
          !which we will calculate eps_2.
          x0(i) = (enrgy(i+1) + enrgy(i))/2.
          !save a single precision version of the grid of midpoints for
          !interpolation later.
          x0real(i) =dble(x0(i))
          !initialize final result
          epsi(i) = 0.

!         !debug
!         x0(i)=0.0

          !loop over the intervals. In each cycle we compute the
          !contribution to the KK integral from the interval between
          !enrgy(j) and enrgy(j+1).

          do 44 j = nepts-1, 1, -1
              if (.false.) then
!    $          (max(dabs(enrgy(j)-x0(i)),
!    $            dabs(enrgy(j+1)-x0(i))).lt.2.5d+1)
!    &        then
!             if (.true.) then !use analytic form everywhere for now.
              else
                  y = e2(j+1)*enrgy(j+1)!+e1(j+1)*x0(i)
                  y = y/(enrgy(j+1)**2+x0(i)**2)
                  y = y +                                               &
     &               (enrgy(j)*e2(j)) &  !+e1(j)*x0(i))                  &
     &               /(enrgy(j)**2+x0(i)**2)
                  y = y*(enrgy(j+1)-enrgy(j))/2.
                  epsi(i) = epsi(i) + y
              end if

   44     continue
          !add final factors
          epsi(i) =2*epsi(i)/pi
          !convert to single precision 
          eps1x0(i)=dble(epsi(i))

!         if(mod(i,5000).eq.0) then
          if(mod(i,500).eq.0) then
             write(slog,fmt="('         ',f15.2,' eV')") x0(i)*hart
             call wlog(slog)
          endif

   33 continue

      !interpolate back to the input grid... obviously this is not fast,
      !but it works. x0(i) is the midpoint of the ith interval between
      !points in enrgy -- calculating at the midpoint keeps the
      !singularity off grid points so, e.g. we never evaluate log(0).
      !We could use a better algorithm here or change the energy grid. I
      !don't like either of these options, so I'm just going to
      !interpolate back to my grid.
      flag=.true.
      do i=2,nepts
        x=omega(i)
        call lint(x0real,eps1x0,nepts,flag,khigh,klow,x,ysngl)
        epsi(i)=ysngl
      enddo
      !not really sure what to do with the first and last points of the
      !output grid since it is outside the interval covered by the x0
      !grid.
      epsi(1)=epsi(2)
      epsi(nepts)=epsi(nepts-1)
          

!     close(unit=88)
      end

