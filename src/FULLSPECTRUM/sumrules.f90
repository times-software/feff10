
!     This subroutine computes several sum rules for optical constants.
!     The optical constants should be stored in opcons.dat;
!     the results of the sum rule integrals will be written to
!     sumrules.dat. For formulas, see Shiles, et. al. phys rev B 22
!     1612.
!     MPP 7/8/03

!     sum1 -- eps_2 sumrule: \int omega*eps_2
!     sum2 -- mu sumrule : \int mu
!     sum3 -- epsinv sumrule \int omega*epsinv (eq 5.6 & 5.9 Shiles)
!     sum4 -- \int mu*(n-1)  
!     sum 5 -- n sumrule: \int (n-1)  (eq 5.1 Shiles)
!     stppwr --  \int omega*log(omega)*epsinv (numerator of RHS eq 5.13
!                                              Shiles)
!     iabsn --  \int abs(n-1) (denom. of RHS of eq 5.2 Shiles)


      subroutine sumrules(numden,infilnm,outfile)
      use constants
      implicit none
      include 'HEADERS/params.h'

      character*512 slog,infilnm,outfile
      real*8 emin, emax
      real*8 B, enrgy, omega, m, neff, n,                                 &
     &                 c, na, e, rho, am,                               &
     &                 mu, refl, eloss, eps1, e2curr, k,e2prev,         &
     &                 sum1, sum2, sum3, sum4, sum5, eprev,             &
     &                 factor,numden,elossl, stppwr, iabsn
      integer iepts, ii, ios
      character comment*100

!       numden=numden*bohr**3 !convert number density to a.u.
        B=1/(2*pi**2*numden) !for opcons.dat
!       write (slog, fmt="('constant in sumrule calc: ',e20.10)") B
!       call wlog(slog)
        sum1=0.0
        sum2=0.0
        sum3=0.0
        sum4=0.0
        sum5=0.0
        stppwr=0.0
        iabsn=0.0

        eprev=0
        e2prev=0

        open(unit=20,file=infilnm,status='old')
!       open(unit=20,file='opconsKK.dat',status='old')
!       open(unit=20,file='desy.dat',status='old')
!       open(unit=20,file='opcons.dat',status='old')
!       open(unit=20,file='opcons.gildardo',status='old')
!       open(unit=21,file='sumrules.dat',status='unknown')
        open(unit=21,file=outfile,status='unknown')
        do ii=1,fullpts
          read(20,*) comment
          if (comment(1:1).ne.'#') exit
        enddo

        ios=0
        ii=0
        do while(ios.eq.0)
          ii=ii+1

!         read in opcons.dat data
          read(20,*,iostat=ios) enrgy, eps1, e2curr                     &
     &       ,n, k, mu , refl, eloss
          !the following needed for Al desy data
          eloss=e2curr/(e2curr**2+(eps1+1.0)**2)


          enrgy=enrgy/hart !convert to a.u.
          mu=mu*bohr/1000 !convert to a.u.


!         add contribution from last read point to sumrules
          if (ios.eq.0) then
            sum1 = sum1 +(enrgy*e2curr+eprev*e2prev)/2                  &
     &                   *(enrgy-eprev)
            sum2 = sum2 + mu*(enrgy-eprev)
            sum3 = sum3 +(enrgy*eloss+eprev*elossl)/2                   &
     &                   *(enrgy-eprev)
            sum4 = sum4 + mu*n*(enrgy-eprev)
            sum5 = sum5+(enrgy-eprev)*n
            iabsn = iabsn+(enrgy-eprev)*abs(n)
            if (eprev.gt.0.0) then
              stppwr = stppwr +(log(enrgy)*enrgy*eloss+                 &
     &                            log(eprev)*eprev*elossl)/2            &
     &                            *(enrgy-eprev)
            else
              stppwr = stppwr +log(enrgy)*enrgy*eloss*(enrgy)
            endif
           
            eprev = enrgy
            e2prev=e2curr
            elossl=eloss
!           if(mod(ii,50).eq.1) then
              write(21,66) enrgy*hart, B*sum1, alpinv*B*sum2, B*sum3,   &
     &                     alpinv*B*sum4, sum5/iabsn, stppwr/sum3
!           endif
   66       format(7e24.10)


          endif
        enddo

!       if(mod(ii,50).ne.1) then !write final line if not already written
!         write(21,66) enrgy*hart, sum1, sum2, sum3, sum4, sum5
!       endif

!       call wlog('number of points read from opcons.dat') 
!       write (slog,fmt="('during computation of sumrules:',i9)") ii
!       call wlog (slog)
        write(slog,                                                     &
     &    fmt="('sumrules based on ',i5,' data points read from ',30a)")&
     &    ii-1,infilnm(1:30)
        call wlog (slog)
        write (slog, fmt="('epsilon_2 sumrule gives n_eff = ',f10.3)")  &
     &   sum1*B
        call wlog (slog)
        write (slog, fmt="('epsinv sumrule gives n_eff = ',f10.3)")     &
     &   sum3*B
        call wlog (slog)
!       write (slog, fmt="('number density used for sumrule: ', e10.3)")
!    &   numden
!       call wlog (slog)
        close(20)
        close(21)
      end
