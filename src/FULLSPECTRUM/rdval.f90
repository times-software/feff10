
      subroutine rdval (numden,iepts,omega,fname,eps2)
      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu
      implicit none

      include 'HEADERS/params.h'


      !number density of atoms
      real*8 numden
      !number of points in full spectrum
      integer iepts
      !spectrum and frequency grid
      real*8 eps2(fullpts),omega(fullpts)
      character*512 fname, slog
      integer ie
      logical flag,rdok
      !number of points in xmu.dat
      integer npts
      !columns of xmu.dat
      real*8 col1(maxpts),col2(maxpts),col3(maxpts)
      real*8 col4(maxpts),col5(maxpts)
      real*8 prefac, x,y
      integer indexi, indexj

      eps2(:)=0.0

      !read xmu.dat from valence calculation
!     fname="xmu.val"
      call rdxmu (fname,rdok,npts,col1,col2,col3,col4,col5)
      !check the return code of reading routine
      if (rdok) then
        !col5 now contains abs. cross-section per atom in square
        !angstroms, col1 has photon frequency in eV

        !prefac is such that eps_2=col4*prefac/col1; bohr**2 converts
        !square angstroms to atomic units.
        prefac=4.0*pi*alpinv*bohr**2*numden
        !convert cross-section to eps_2*omega
        do ie=1,npts
          col1(ie)=col1(ie)/hart
          col4(ie)=col4(ie)*prefac !/col1(ie)
        enddo !loop over lines read from xmu.dat

        open (unit=35,file="valence1.dat")
        do ie=1,npts
          write(unit=35,fmt="(20e20.10)") col1(ie)*hart, col4(ie)
        end do
        close (unit=35)

        !now we need to interpolate the current contribution onto the
        !output grid omega
        flag=.true.
        do ie=1,iepts
          if (omega(ie).lt.col1(npts)) then
            !x is the frequecny of the current point
            x=omega(ie)
            !find eps2 at this frequency by interpolation
            call lint(col1, col4,npts,flag,indexi,indexj,x,y)
            !store the interpolated value in the output array, adding 1/omega ommitted above
            eps2(ie)=y/x
          endif

        enddo !loop over frequency grid


      else
        !read failed; report this fact and zero our output
        eps2(:)=0.0
        call wlog ('Failed reading valence xmu.dat')
      end if

      !convert to eps2
      !interpolate and add

      end !subroutine rdval
