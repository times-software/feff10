!     written June 2005 Micah Prange

!     Estimates number density of atoms of a specified type from 
!     a feff input file.

!     input: 
!     iz: atomic number of species for which we want to calculate number
!         density
!     path: filename of directory for some edge that has already been
!           calculated (so we can read info from geom.dat and pot.bin).

!     output:
!     numden: number density of species iz in inverse cubic bohrs

      subroutine rddens (path,iz,numden)
      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu
      implicit none

      include 'HEADERS/params.h'
      character*512 path, slog, str
      double precision  rnrm(0:nphx)

      integer  zz(0:nphx)
      double precision  xnatph(0:nphx),totn
      character*80 title(nheadx)
      double precision  dum(13)
      double precision totvol,  numden
      integer  i, ntitle,k,iz, nph


      !initialize variables
      do i=0,nphx
        xnatph(i)=0.0
        rnrm(i)=0.0
      enddo
      nph=0

      !append filename to directory path
      slog='/fms_im/pot.bin'
      call concat (path,slog,str,k)
      path=str
      !read potentials info to get totvol and xnatph
      !rdpotx is the same as rdpot, except it takes a pathname
      !argument path to the file to read.
        call rdpotp_fs  (path, ntitle, title, xnatph, rnrm, zz, nph)

      
      totn=0.0
      totvol=0.0
      do i=0,nph
        if (zz(i).eq.iz) then
          totn=totn+xnatph(i)
        end if
        totvol=totvol+xnatph(i)*rnrm(i)**3*4*pi/3
        write (slog,fmt=*) 'rnrm (Angstroms): ', rnrm(i)*bohr
        call wlog (slog)
        write (slog,fmt=*) 'totvol: ', totvol
        call wlog (slog)
        write (slog,fmt=*) 'totn: ', totn
        call wlog (slog)
        write (slog,fmt=*) 'xnatph(i): ', xnatph(i)
        call wlog (slog)
      enddo


      numden=totn/totvol

      return
      end
