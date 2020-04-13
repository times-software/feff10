!     written June 2004 Micah Prange

!     Reads absolute absorption cross section as a function of photon energy
!     from xmu.dat file at path fname.

!     input:
!     fname: filename of file to read from

!     output:
!     flag: false if read failed, true otherwise
!     npts: number of points read and returned in arrays
!           must be less than or equal to parameter maxpts
!     enrgy: photon energy from column one of xmu.dat (eV !!!)
!     eredge: photoelectron energy from column two of xmu.dat (eV !!!)
!     ck: photoelectron momentum from xmu.dat (1/angs !!!)
!     mu0: atomic background abs. cross section from xmu.dat (angs**2 !!!)
!     mu: abs. cross section from xmu.dat (angs**2 !!!)
 
      subroutine rdxmu (fname,flag,npts,enrgy,eredge,ck,mu,mu0)
      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu
      implicit none

      include 'HEADERS/params.h'

      character*512 fname, slog
      character  comment*36
      integer ios,npts,ihi,ilo,ie
      logical flag
      real*8 xnorm, enrgy(maxpts),eredge(maxpts),ck(maxpts)
      real*8 mu0(maxpts),mu(maxpts),xchi

!     call wlog (fname)

!     write(slog,fmt="('maxpts: ',i7)") maxpts
!     call wlog (slog)
  100 format(a36)
  101 format(a36,1pe20.4)
      flag=.false. 
      open(unit=42,file=fname,status='old',iostat=ios)
!     write(slog,fmt="('ios: ',i5)") ios
!     call wlog (slog)
      if (ios.ne.0) then
        do ie=1,maxpts
          mu(ie)=0
          mu0(ie)=0
          enrgy(ie)=0
        enddo
        npts=2
        call wlog('Read from file failed! ')
        call wlog('file:')
        slog='     '//fname
        call wlog(slog) 
        return
      endif

!     The mu data reported by Feff are normalized by the 
!     total absorbtion cross section at the edge energy
!     + 50 eV. This normalizing value (in square angs) is
!     written in each xmu.dat file on a line after one
!     containing 'paths used'. Here we find the normailization
!     to convert mu to an absolute value in angs^2. 

          do while (comment(14:23).ne.'paths used')
              read(42,100) comment
!             call wlog (comment)
          enddo
          read(42,101) comment, xnorm
          read(42,100) comment
          read(42,100) comment
!         write(slog,fmt="('xnorm: ',f20.10)") xnorm
!         call wlog (slog)

!     Read data. The variables are as follows
!     enrgy = abs photon energy
!     eredge = rel energy
!     ck     = photoelectron wave number
!     mu     = abs. cross section (Angs**2!!)
!     mu0    = atomic background cross section
!     xchi   =  chi, fine structure       

!     write(slog,fmt="('path: ',a60)") fname
!     call wlog(slog)

          npts = 0
          do while (ios.eq.0)
              npts = npts + 1
              read(42,*,iostat=ios) enrgy(npts), eredge(npts), ck(npts),&
     &                              mu(npts), mu0(npts), xchi
!             write(slog,fmt="(20e20.10)") enrgy(npts), mu(npts),
!    &        mu0(npts),mu(1)
!             write(slog,fmt="('point read at ',f20.10,' eV.')")
!    &          enrgy(npts)
!             call wlog (slog)
              mu(npts)=mu(npts)*xnorm !converts to angs^2
              mu0(npts)=mu0(npts)*xnorm !converts to angs^2
          enddo
!     write(slog,fmt="('mu(1): ',20e20.10)") mu(1)
!     call wlog(slog)
      if (npts.gt.1) then
        flag=.true. !normal exit
      else
        flag=.false. !at most one point was read, somthing is wrong
        do ie=1,maxpts
          mu(ie)=0
          mu0(ie)=0
          enrgy(ie)=0
        enddo
        npts=2
        call wlog('Read from file failed!')
        call wlog('file:')
        slog='     '//fname
        call wlog(slog) 
      endif
      !read of the (npts)th point failed, so do not report this point.
      npts=npts-1 

!     write(slog,fmt="('mu(1): ',20e20.10)") mu(1)
!     call wlog(slog)
!     write(slog,fmt="('number of points: ',i5)") npts
!     call wlog(slog)

      return
      end
