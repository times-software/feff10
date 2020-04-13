!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdeps.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdeps (omp,nplmax,npl,plengy,oscstr,plbrd)
      use dimsmod, only: nheadx
	  use constants
      implicit double precision (a-h, o-z)

      integer  npl,nplmax,ipl
      double precision  plengy(nplmax),oscstr(nplmax),plbrd(nplmax)

!     Local stuff
      character*512 slog
      character*80 head(nheadx),line
      dimension lhead(nheadx)

!     initialize
      do ipl=1,nplmax
        plengy(ipl)=0.d0
        oscstr(ipl)=0.d0
        plbrd(ipl)=0.d0
      enddo

!     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)

!     read s02.inp
      open (file='exc.dat', unit=3, status='old',iostat=ios)
      if (ios.eq.0) then
        npl=0
 100      read (3,'(a)',end=200)  line
          istart=1
 110      if (line(istart:istart).eq.' ') then
            istart=istart+1
            goto 110
          endif
          if (line(istart:istart).ne.'#') then
            npl=npl+1
            read (line,*) plengy(npl),plbrd(npl),oscstr(npl)
          endif
          goto 100
 200    continue
        do ipl=1,npl
          plengy(ipl)=plengy(ipl)/hart
          plbrd(ipl)=plbrd(ipl)/hart
        enddo
      else
        npl=1
        plengy(1)=omp
        plbrd(1)=0.001d0*omp
        oscstr(1)=1.d0
        open (file='exc.dat', unit=3, status='unknown')
        write(3,30) plengy(1)*hart,plbrd(1)*hart,oscstr(1)
      endif
      close(3)
      


!      write(6,*) '***',npl,'*',plengy(1),'*',plbrd(1),'*',oscstr(1),'***'


      return
      end
