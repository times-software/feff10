      subroutine rdst (path, omega, npts,                               &
     &                   sctfac, sctfc0, ehi, elo, ierr)
      !Sums results of XANES, DANES, and EXAFS calculations for a single
      !initial state.  Returns the x-ray scattering factor interpolated
      !onto the grid omega given that there are files
      !path/fms_im/xmu.dat, path/fms_re/xmu.dat, path/path_im/xmu.dat,
      !path/path_re/xmu.dat which contain the calculations to be
      !interpolated. The fms_im file contains photoabs. cross section
      !per site (XANES card); fms_re contains f'=Re(f) (DANES card). The
      !path files include the same information calculated via path
      !expansion. A smooth transition is made from FMS to path expansion
      !calculations between photoelectron momenta khi, klo.

      !inputs:
      !   path -- path to dir containing calculations for a given edge
      !   omega -- energy grid for epsilon
      !   npts -- number of points in omega
      !outputs:
      !   sctfac -- array holding complex (scalar) scattering factor
      !   ehi,elo -- frequency interval over which fine structure was
      !              calculated. 2 element arrays with first element for
      !              real part, second for im part.
      !   ierr -- error flag

      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu
      implicit none

      include 'HEADERS/params.h'

      real*8 eps1,eps1bk,eps2,eps2bk !these hold real or im part for one freq.
      complex sctfac (fullpts) !cmplx f for all omega
      complex sctfc0 (fullpts) !atomic background of sctfac
      real*8 fprime(fullpts), fdblp (fullpts) !real and im parts of sctfac
      real*8 fp0(fullpts), fpp0 (fullpts) !real and im parts of sctfc0
      integer npts,i,k,indexi,indexj,j, nfiles, ierr
      real*8 elo(2), ehi(2), omega(fullpts),y,y2,y0,y02,frac,numden
      character*512 path, slog,str, rfname(2), ifname(2), fname
      !rfcalc ifcalc hold real and im parts of f on the original grids
      !read from xmu.dat. rfcalc(i,1) is f' from fms calculation,
      !ifcalc(i,2) is f'' from path expansion, etc.
      real*8 rfcalc (maxpts,2), ifcalc (maxpts,2)
      !rf0, if0 are the backgrounds for rfcalc, ifcalc
      real*8 rf0 (maxpts,2), if0 (maxpts,2)
      !rtrint, itrint mark the energy intervals to transition from fms
      !calculations to path expansion for f' and f'' resp.
      real*8 rtrint(2), itrint(2)
      !energy grids from xmu.dat 
      real*8 emreal(maxpts,2), emim(maxpts,2)
      !number of points in emreal, rfcalc, emim, ifcalc
      integer renex(2) , imnex(2)

      logical conv,rdok,flag
      real*8 col1(maxpts), col2(maxpts),col3(maxpts),                     &
     &     col4(maxpts), col5(maxpts)
      real*8 scrtch(nex),scrtc2(nex), theta,x

!     call wlog ('starting rdst')
!     write(slog,fmt="('npts: ',i5)") npts
!     call wlog(slog)

      !iinitialize error flag
      ierr=0

      !construct filenames
      str='/fms_re/xmu.dat'
      call concat (path,str,slog,k)
      rfname(1)=slog
!     call wlog (slog)
      str='/fms_im/xmu.dat'
      call concat (path,str,slog,k)
      ifname(1)=slog
!     call wlog ('ifname(1): '//ifname(1))
      str='/path_re/xmu.dat'
      call concat (path,str,slog,k)
      rfname(2)=slog
!     call wlog (slog)
      str='/path_im/xmu.dat'
      call concat (path,str,slog,k)
      ifname(2)=slog
!     call wlog (slog)

!     goto 100
!     call wlog ('2 ifname(1): '//ifname(1))

      do i=1,2 !i=1 --> fms; i=2 --> path expansion

!Comment out reading of real part so we don't have to run DANES
!       !read the ith file for the real part
!       !read data from file
!       fname=rfname(i)
!       call rdxmunorm (fname,rdok,k,col1,col2,col3,col4,col5)
!       if(.not.rdok) ierr=ierr+1
!       renex(i)=k


        !save this calculation, convert to a.u.
        do j=1,renex(i)
          emreal(j,i)=col1(j)/hart
          !set interval to transition from fms to path exp.
          if (col3(j).le.khi) rtrint(2)=emreal(j,i)
!         write (unit=59,fmt="(4f20.10)") emreal(j,i)*hart,col3(j), khi,
!    &          rtrint(2)
          if (col3(j).le.klo) rtrint(1)=emreal(j,i)
          rfcalc(j,i)=col4(j)
          rf0(j,i)=col5(j)
        enddo


        !read the ith file for the im part
        !read data from file
        fname=ifname(i)
!       write (slog,fmt="(i1,': ',a)") i,ifname(i)
!       call wlog (slog)
!       call wlog ('fname is: ')
!       call wlog (fname)
        call rdxmu (fname,rdok,k,col1,col2,col3,col4,col5)
        if(.not.rdok) ierr=ierr+1
!       call wlog ('well, what is fname now??')
!       call wlog (fname)
        imnex(i)=k
!       call wlog('1here is the data from '//fname)
!       do j=1,k
!         write(slog,fmt="(20e20.10)") col1(j),col4(j),col5(j)
!         call wlog(slog)
!       enddo

        !save this calculation, convert to a.u.
        do j=1,imnex(i)
          emim(j,i)=col1(j)/hart
          !set interval to transition from fms to path exp. from the fms
          !file

          if(i.eq.1) then
            if (col3(j).le.khi) itrint(2)=emim(j,i)
            if (col3(j).le.klo) itrint(1)=emim(j,i)
!           write(slog,fmt="('khi, klo, k, E: ',20e20.10)") khi, klo,
!    &      col3(j), emim(j,i)*hart,itrint(1)*hart,itrint(2)*hart
!           call wlog (slog)
          endif
          !convert photoabs. cross section to f'' and save
          ifcalc(j,i)=col4(j)*alpinv*emim(j,i)*bohr**2
          if0(j,i)=col5(j)*alpinv*emim(j,i)*bohr**2
        enddo
!       call wlog('2here is the data from '//fname)
!       do j=1,imnex(i)
!         write(slog,fmt="(20e20.10)") 
!    &    emim(j,i)*hart,ifcalc(j,i),if0(j,i)
!         call wlog(slog)
!       enddo

      enddo !end loop over fms and path exp. 
      !all four files are now read and must be combined.

      !now interpolate real fms calculation onto output grid.
      indexi=0
      indexj=0
      do j=1,renex(1)
        col1(j)=emreal(j,1)
        col4(j)=rfcalc(j,1)
        col5(j)=rf0(j,1)
      enddo
!     write(slog,fmt="('trans. interval for fprime: ',4f20.10)")
!    &   rtrint(1)*hart,rtrint(2)*hart,rtrint(1),rtrint(2)
!     call wlog (slog)
!     open(unit=12,file='re_fms.dat')
      do j=1,npts !loop to read f'
        !if we need f' from fms calculation at this omega ...
        if (omega(j).ge.emreal(1,1).and.omega(j).le.rtrint(2)) then
          x=omega(j)
          call lint(col1, col4,nex,flag,indexi,indexj,                  &
     &              x,y)
          fprime(j)=y !store in output variable.
          call lint(col1, col5,nex,flag,indexi,indexj,                  &
     &              x,y)
          fp0(j)=y !store in output variable.
!         write(unit=12,fmt="(2e20.10)") x,y
        endif
      enddo
!     close (unit=12)

      !now interpolate real path exp calculation onto output grid.
      !we will smoothly transition into path exp.
      indexi=0
      indexj=0
      do j=1,renex(2)
        col1(j)=emreal(j,2)
        col4(j)=rfcalc(j,2)
        col5(j)=rf0(j,2)
      enddo
!     open(unit=12,file='re_path.dat')
      do j=1,npts !loop to read f'
        !if we need f' from path exp. calculation at this omega ...
        if                                                              &
     &    (omega(j).le.emreal(renex(2),2).and.omega(j).ge.rtrint(1))    &
     &    then

          !compute fraction of signal we want to come from fms at
          !this frequency. frac=0 above end of transition interval
          if (omega(j).le.rtrint(2).and.omega(j).ge.rtrint(1)) then
            frac=(rtrint(2)-omega(j))/(rtrint(2)-rtrint(1))
          else
            frac=0.0
          end if
          !compute mixing angle for smooth transition
          theta=frac*pi/2

          !interpolate onto output omega and store values
          x=omega(j)
          call lint(col1, col4,nex,flag,indexi,indexj,                  &
     &              x,y)
!         write(unit=12,fmt="(2e20.10)") x,y
          fprime(j)=y*cos(theta)**2+fprime(j)*sin(theta)**2 
          call lint(col1, col5,nex,flag,indexi,indexj,                  &
     &              x,y)
          fp0(j)=y*cos(theta)**2+fp0(j)*sin(theta)**2 
        endif
      enddo
!     close (unit=12)
      write(slog,fmt="('transition interval: ',20e20.10)")              &
     &    itrint(1)*hart,itrint(2)*hart
      call wlog (slog)

      !now interpolate im fms calculation onto output grid.
      indexi=0
      indexj=0
      do j=1,imnex(1)
        col1(j)=emim(j,1)
        col4(j)=ifcalc(j,1)
        col5(j)=if0(j,1)
      enddo
      open(file='fdplp.dat',unit=12)
      write(unit=12,fmt="('#npts: ',i7)") npts
      write(unit=12,fmt="('#emim(1,1): ',e20.10)") emim(1,1)
      write(unit=12,fmt="('#itrinit(1): ',e20.10)") itrint(1)*hart
      write(unit=12,fmt="('#itrinit(2): ',e20.10)") itrint(2)*hart
      write(unit=12,fmt="('#omega(1): ',e20.10)") omega(1)
      write(unit=12,fmt="('#omega(npts): ',e20.10)") omega(npts)
      do j=1,npts !loop to read f''
        !if we need f'' from fms calculation at this omega ...
        if (omega(j).gt.emim(1,1).and.omega(j).lt.itrint(2)) then
          x=omega(j)
          call lint(col1, col4,nex,flag,indexi,indexj,                  &
     &              x,y)
          fdblp(j)=max(y,0.0) !store in output variable.
          call lint(col1, col5,nex,flag,indexi,indexj,                  &
     &              x,y)
          fpp0(j)=max(y,0.0) !store in output variable.
          write(unit=12,fmt="(i4,2e20.10)") j,x,y
        endif
      enddo
      close(unit=12)

      !now interpolate im path exp calculation onto output grid.
      !we will smoothly transition into path exp.
      indexi=0
      indexj=0
      do j=1,imnex(2)
        col1(j)=emim(j,2)
        col4(j)=ifcalc(j,2)
        col5(j)=if0(j,2)
      enddo
      open(file='fdplp2.dat',unit=12)
      do j=1,npts !loop to read f'
        !if we need f' from path exp. calculation at this omega ...
        if (omega(j).le.emim(imnex(2),2).and.omega(j).ge.itrint(1)) then

          !compute fraction of signal we want to come from fms calc. at
          !this frequency. frac=0 above end of transition interval
          if (omega(j).le.itrint(2).and.omega(j).ge.itrint(1)) then
            frac=(itrint(2)-omega(j))/(itrint(2)-itrint(1))
          else
            frac=0.0
          end if

          !compute mixing angle for smooth transition, and store mixture
          theta=frac*pi/2

          !interpolate onto output omega
          x=omega(j)
          call lint(col1, col4,nex,flag,indexi,indexj,                  &
     &              x,y)
          fdblp(j)=max(y,0.0)*cos(theta)**2+fdblp(j)*sin(theta)**2 
          call lint(col1, col5,nex,flag,indexi,indexj,                  &
     &              x,y)
          fpp0(j)=max(y,0.0)*cos(theta)**2+fpp0(j)*sin(theta)**2 
        endif
        write(unit=12,fmt="(7e20.10)")frac,omega(j),fdblp(j),           &
     &        sin(theta)**2,                                            &
     &        cos(theta)**2,itrint(1),itrint(2)
      enddo
      close(unit=12)

      !record frequency interval over which calculations were performed
      ehi(1)=emreal(renex(2),2)
      ehi(2)=emim(imnex(2),2)
      elo(1)=emreal(1,1)
      elo(2)=emim(1,1)

      !combine real and im parts
!     open(unit=12,file='fprime_fine_st.dat')
      do i=1,npts
        sctfac(i)=fprime(i)+coni*fdblp(i)
        sctfc0(i)=fp0(i)+coni*fpp0(i)
!       write(unit=12,fmt="(3e20.10)") omega(i),sctfac(i)
      enddo
!     close(unit=12)

100   continue      
      end
