      subroutine opcons(omega,eps,npts,fname, path)
      !this subroutine takes the complex dielectric function and
      !computes the complex index of refraction, the absorption
      !coefficient, normal incidence reflectivity, and the energy loss
      !function. These optical constants are written into the file
      !opcons.dat.

      !note that eps should actually contain eps-1 so that the real part
      !is close to 0 instead of close to 1. Same comment applies to n.

      !inputs:
      !     omega -- frequency grid
      !     eps   -- complex scalar dielectric function
      !     npts  -- length of omega and eps
      !     fname -- filename for output file
      !     path  -- pathname to dir containing pot.bin. The xmu.dat
      !              header will be read from pot.bin and written at the
      !              top of the output file.

      !outputs: none, except the file.
      use dimsmod, only: nheadx,nex,nphx=>nphu
      use constants
      implicit none
      include 'HEADERS/params.h'
      complex eps(fullpts), refrac, cureps, curomg
      real*8 omega(fullpts), mu,reflct,  eloss, e1, e2
      integer npts, i, ios,j, istrln
      character*512 slog, str
      character*512 fname, path, path2
!     pot.bin stuff

      double precision   rnrm(0:nphx)
      integer  iz(0:nphx)
      double precision xnatph(0:nphx)
      character*80 title(nheadx)
      integer  ntitle,k, nph



      !open output file
      open(unit=66,file=fname,status='unknown',iostat=ios)
      call chopen (ios,fname,'opcons')


      !save path so opcons returns the same path it was passed
      path2=path
      !append filename to directory path
      slog='/fms_im/pot.bin'
      call concat (path,slog,str,k)
      path=str


      !read potentials info. rdpotp is the same as rdpot, except it takes
      !a pathname argument path to the file to read.
        call rdpotp_fs  (path, ntitle, title, xnatph, rnrm, iz, nph)

      !restore path
      path=path2
      do i=1,ntitle
        slog=title(i)
        j=istrln(slog)
        write(unit=66,fmt="(a)") '# '//slog(1:j)
      enddo
      slog="#   omega (eV)      epsilon_1       epsilon_2       n"//    &
     &"               kappa           mu (cm^(-1))    R"//              &
     &"               epsinv"
      write(unit=66,fmt="(a)") slog(1:125)


      !loop through the frequencies
      do i=1,npts

        !get other optical constants from epsilon
        call eps2opt(eps(i),omega(i),refrac,mu,reflct,eloss)

        !write them out to opconsKK.dat
        write (unit=66,fmt="(19e16.6)")                                 &
     &   omega(i)*hart,eps(i),refrac-1.0,mu,reflct,eloss

      enddo !end big loop over frequencies

      !close output file
      close(unit=66)

      end

      subroutine eps2opt(eps,omega,refrac,mu,reflct,eloss)
      !given complex dielectric contsant minus 1 in eps and frequency in
      !omega, computes the complex index of refraction (refrac), normal
      !incidence reflectivity (reflct), absorption coefficient in
      !inverse angstroms (mu), and energy loss function (eloss).
      use constants
      implicit none
      include 'HEADERS/params.h'

      complex refrac,eps
      real*8 eloss, reflct,mu, omega

      !index of refraction N=n+ik=(epsilon)^(1/2)
      refrac=sqrt(eps+1)
      !normal incidence reflectance
      reflct=abs((refrac-1)/(refrac+1))**2
      !absorption coefficient in inverse angstroms
      mu=2*omega*alphfs*aimag(refrac)/bohr*1000
      eloss=-1.0*aimag((eps+1)**(-1))

      end !subroutine eps2opt()
