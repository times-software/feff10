!     written Nov. 2004 Micah Prange

!     Reads ldos from file at path fname (presumably ldosNN.dat).
!     input:
!     fname: filename of file to read from

!     output:
!     flag: false if read failed, true otherwise
!     npts: number of points read and returned in arrays
!           must be less than or equal to parameter maxpts
!     enrgy: photon energy from column one of xmu.dat (converted to hartrees)
!     ldos(l,i): the dos for angular momentum l listed on the ith line of fname
!                ldos is converted to states/hartree/atom
      subroutine rdldos (fname,flag,npts,efermi,enrgy,ldos)
      use constants
      implicit none
      include 'HEADERS/params.h'

      character*512 fname, slog
      character  comment*36
      integer ios,npts,ie,i,j,k
      logical flag
      real*8 ldos(0:lmax,maxpts), enrgy(maxpts), efermi
      real*8 col1, col2, col3, col4, col5, col6, col7, col8

  100 format(a36)
  101 format(a36,1pe20.4)
      flag=.false. 
      open(unit=42,file=fname,status='old',iostat=ios)
      if (ios.ne.0) then
        do ie=1,maxpts
          do j=0,lmax
            ldos(j,ie)=0.0
          enddo
        enddo
        npts=2
        call wlog('Read from file failed!')
        call wlog('file:')
        slog='     '//fname
        call wlog(slog) 
        return
      endif

!     This finds the fermi energy in the header and skips through
!     the rest of the header to the data.
      read(unit=42,fmt=55) comment, efermi
   55 format(a20,f12.3)
      write(slog,fmt="( 'fermi energy = ',f10.5, ' eV')") efermi
      call wlog(slog)
      efermi=efermi/hart
   11 continue
      read(unit=42,fmt=1936) comment
      k = index(comment, 'Lorentzian')
      if(k.eq.0) goto 11
 1936 format(a)
      read(unit=42,fmt=1936) comment
      read(unit=42,fmt=1936) comment

!     Read data. 
      npts = 0
      do while (ios.eq.0)
        npts = npts + 1
 1935   format ( f10.4, 4e13.6)
        read(unit=42,fmt=*,iostat=ios) col1, col2, col3, col4, col5
        enrgy(npts)=col1/hart
        ldos(0,npts)=col2*hart
        ldos(1,npts)=col3*hart
        ldos(2,npts)=col4*hart
        ldos(3,npts)=col5*hart
      enddo

      if (npts.gt.1) then
        flag=.true. !normal exit
      else
        flag=.false. !at most one point was read, somthing is wrong
        do ie=1,maxpts
          do j=0,lmax
            ldos(j,ie)=0.0
          enddo
        enddo
        npts=2
        call wlog('Read from file failed!')
        call wlog('file:')
        slog='     '//fname
        call wlog(slog) 
      endif
      !read of the (npts)th point failed, so do not report this point.
      npts=npts-1 
      

!     open (unit=12, file='testdos.dat')
!     do ie=1,npts
!       write(unit=12,fmt="(5e20.10)") enrgy(ie),
!    &     ldos(0,ie),ldos(1,ie),ldos(2,ie),ldos(3,ie)
!     enddo
!     close(unit=12)

      return
      end
