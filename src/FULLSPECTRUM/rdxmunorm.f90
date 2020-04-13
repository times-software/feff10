!     written June 2004 Micah Prange

!     Reads normalized absorption coef. as a function of photon energy
!     from xmu.dat file at path fname. See nearly identical subroutine
!     rdxmu.f for reading absolute abs. coef.

!     input:
!     fname: filename of file to read from

!     output:
!     flag: false if read failed, true otherwise
!     npts: number of points read and returned in arrays
!           must be less than or equal to parameter maxpts
!     enrgy: photon energy from column one of xmu.dat (eV !!!)
!     mu0: atomic background abs. from xmu.dat (1/angs !!!)
!     mu: abs. coef. from xmu.dat (1/angs !!!)
!     ihi: index of last point with momentum less than param khi
!     ilo: index of last point with momentum less than param klo
 
      subroutine rdxmunorm (fname,flag,npts,enrgy,eredge,ck,mu,mu0)
      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu
      implicit none

      include 'HEADERS/params.h'
      character comment*36
      character*512 fname, slog
      integer ios,npts,ihi,ilo,ie
      logical flag
      real*8 xnorm, enrgy(maxpts),mu0(maxpts),mu(maxpts)
      real*8 eredge(maxpts),ck(maxpts),xchi(maxpts)

  100 format(a36)
  101 format(a36,1pe20.4)

      flag=.false. 
      open(unit=42,file=fname,status='old',iostat=ios)
      if (ios.ne.0) then
        mu(:)=0
        mu0(:)=0
        enrgy(:)=0
        call wlog('Read from file failed!')
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
          enddo
          read(42,101) comment, xnorm
          read(42,100) comment
          read(42,100) comment

!     Read data. The variables are as follows
!     enrgy = abs photon energy
!     eredge = rel energy
!     ck     = photoelectron wave number
!     mu     = abs. cross section (Angs**2!!)
!     mu0    = atomic background cross section
!     xchi   =  chi, fine structure       


          npts = 0
          do while (ios.eq.0)
              npts = npts + 1
              read(unit=42,fmt=*,iostat=ios)                            &
     &          enrgy(npts), eredge(npts), ck(npts), mu(npts),          &
     &          mu0(npts), xchi(npts)
          enddo
      flag=.true. !normal exit
      npts=npts-1
      return
      end
