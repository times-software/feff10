      subroutine rdbkg (path, omega, npts,                              &
     &                   sctfac, neff, fp0, ierr)
      !Sums results of FPRIME calculations for a single initial state.
      !Returns the x-ray scattering factor interpolated onto the grid
      !omega given that there are files path/fprime1/xmu.dat,
      !path/fprime2/xmu.dat, etc. which contain the calculations to be
      !interpolated. Files numbered consecutively from 1 will be used.
      !For each energy point, the file with the smallest index will be
      !used (e.g. fprime1 takes precedence  over frpime2). Also computes
      !the eps_2 sumrule to find the effective number of electrons in
      !this edge, which is returned in neff.

      !inputs:
      !   path -- path to dir containing calculations for a given edge
      !   omega -- energy grid for epsilon
      !   npts -- number of points in omega
      !outputs:
      !   sctfac -- array holding complex (scalar) scattering factor
      !   neff -- number of effective electrons given by the eps_2
      !           sumrule.
      !   fp0 -- the value of f'(omega=0) needed by calling routines.
      !   ierr -- error flag

      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu
      implicit none

      include 'HEADERS/params.h'
      real*8 eps1,eps1bk,eps2,eps2bk !these hold real or im part for one freq.
      complex sctfac (fullpts) !cmplx f for all omega
      real fprime(fullpts), fdblp (fullpts) !real and im parts of sctfac
      integer npts,i,k,indexi,indexj,j, nfiles, ierr
      real*8 fulfdp(fullpts),efdp(fullpts) !f" on a big grid to neff integration
      real*8 elo, ehi, omega(fullpts),y,y2,y0,y02,frac,numden, eps_e
      character*512 path, slog,str, fname
      logical conv,rdok,flag
      real*8 col1(nex), col2(nex),col3(nex), col4(nex), col5(nex)
      real*8 co1(maxpts), co2(maxpts),co3(maxpts), co4(maxpts), co5(maxpts) !KJ fixing small bug
      real*8 scrtch(nex),scrtc2(nex), theta,x, neff,emin,emax
      !fp0  holds the real part of the scattering factor for the
      !lowest energy available; efp0 holds the energy of this point.
      real*8 fp0, efp0

!     call wlog ('starting rdbkgd')

      !set a small freq. to keep our calculation away from zero freq. to
      !avoid div. by 0 when we add 1/omega**2 factors
      eps_e=0.01/hart

      !initialize fp0, efp0, ierr
      fp0=0.0
      efp0=1e20
      ierr=0

      !first, count fprime files
      do i=1,1000
        !construct filename  
        str='/fprime'//char(i+48)//'/xmu.dat'
        call concat (path,str,fname,k)
        !see if it exists
        inquire(file=fname,exist=flag)
        if (flag) then
          !file exists, keep going
!         call wlog(fname)
        else
          !file doesn't exist, we'll assume we have found all files
          exit !escapes do loop over file names
        end if
      enddo
      nfiles=i-1
!     write(slog,fmt="('nfiles: ',i5)") nfiles
!     call wlog (slog)

      !form large grid for neff calculation
      emin=bigemin
      emax=bigemax
      call egrid_lin (fullpts, emin, emax, efdp)

      do i=nfiles,1,-1
        !read the ith file and interpolate it onto omega grid

        !construct filename
        str='/fprime'//char(i+48)//'/xmu.dat'
        call concat (path,str,fname,k)

        !read data from file
        call rdxmunorm (fname,rdok,k,co1,co2,co3,co4,co5)
        !KJ fixing small bug here:
        col1=co1(1:nex)
        col2=co2(1:nex)
        col3=co3(1:nex)
        col4=co4(1:nex)
        col5=co5(1:nex)
        if(.not.rdok) ierr=ierr+1

        !if  there is a point at 0, move it 
        if (col1(1).le.0.0) then
          eps_e=eps_e*hart !data is still in eV, so adjust delta
          !move lowest point to omega=eps_e
          indexi=0
          indexj=0
          call lint(col1,col2,nex,flag,indexi,indexj,eps_e,y)
          col2(1)=y
          indexi=0
          indexj=0
          call lint(col1,col3,nex,flag,indexi,indexj,eps_e,y)
          col3(1)=y
          indexi=0
          indexj=0
          call lint(col1,col4,nex,flag,indexi,indexj,eps_e,y)
          col4(1)=y
          indexi=0
          indexj=0
          call lint(col1,col5,nex,flag,indexi,indexj,eps_e,y)
          col5(1)=y
          !if there is more than one point at non-positive freq., complain to user and drop the edge.
          if(col1(2).le.0.0) then
            ierr=ierr+1
            call wlog ('Choking on non-positive energy points in')
            call wlog (path)
          end if
        end if

        !find boundaries of this calculation, convert to a.u.
        ehi=-1e20
        elo=1e20
        do j=1,k
          col1(j)=col1(j)/hart
          if (col1(j).gt.ehi) ehi=col1(j)
          if (col1(j).lt.elo) elo=col1(j)
        enddo

        !see if this calculation gives a lower energy than previous ones; if so, update fp0, efp0
        if (elo.lt.efp0) then
          efp0=elo
          fp0=col4(1)
        end if
 
!       call wlog ('starting interpolation')
        !now interpolate this calculation onto the output grid. note
        !that we will overwrite anything that is there already.
        !add 1/omega**2 for interpolation -- this amounts to usnig a
        !linear interpolation for epsilon instead of the scattering factor
        do j=1,k
          if (col1(j).ne.0) then
            scrtch(j)=(col4(j)-fp0)/col1(j)**2
          else
            scrtch(j)=0.0
          end if
        enddo
        indexi=0
        indexj=0
        do j=1,npts !loop to read f'
          if (omega(j).gt.elo.and.omega(j).lt.ehi) then
            x=omega(j)
            call lint(col1,scrtch,nex,flag,indexi,indexj,x,y)
            fprime(j)=y*omega(j)**2+fp0
          endif
        enddo
        indexi=0
        indexj=0
        do j=1,npts !loop to read f''
          if (omega(j).gt.elo.and.omega(j).lt.ehi) then
            x=omega(j)
            call lint(col1,col5,nex,flag,indexi,indexj,x,y)
!           fdblp(j)=max(y*omega(j)**2+fp0,0.0)
!           fdblp(j)=max(y/omega(j)**2,0.0)
            fdblp(j)=max(y,0.0)
          endif
        enddo
 
        !interpolate f'' onto big grid for calculation of neff     
        indexi=0
        indexj=0
        do j=1,fullpts !loop to read f''
          if (efdp(j).gt.elo.and.efdp(j).lt.ehi) then
            x=efdp(j)
            call lint(col1, col5,nex,flag,indexi,indexj,x,y)
            fulfdp(j)=max(y,0.0)
          endif
        enddo


      enddo

      !combine real and im parts
      do i=1,npts
        if(omega(i).le.efp0) then
          sctfac(i)=fp0
        else
          sctfac(i)=(fprime(i)+coni*fdblp(i))!/omega(i)**2
        end if
      enddo

      !add 1/omega**2 factor to fulfdp so we can use our sumrule routine which expects eps_2 as input
      do i=1,fullpts
        fulfdp(i)=fulfdp(i)/efdp(i)**2
      enddo
      !set number density argument so sumrule routine gives correct output
      x=1/(4*pi)
      !get neff for this edge
      call qsum(x,fulfdp,efdp,fullpts,neff)
!     write (slog,fmt="('neff used for f_0 estimate: ',f10.6)") neff
!     call wlog (slog)


      end
