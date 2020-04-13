!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ridxmu.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ridxmu(kinit, ne, ee, chil2, chil3, chil4, chil5, deltaso)
!c    reads several xmu.dat files and interpolates fine structure
!c    on energy grid given by array ee
!c    Input:
!c       kinit - initial state kappa
!c       ne - number of points in energy grid
!c       ee - energy grid
!c    Output:
!c       chil* - (1+chi) for various channels
!c       deltaso - spin-orbit splitting between 2 edges
!c    written by Alexei Ankudinov , Dec. 2004
      use dimsmod, only: nex 
	  use constants
      implicit double precision (a-h,o-z)
      dimension ee(nex), chil2(nex), chil3(nex), chil4(nex), chil5(nex)
      dimension x1(nex), x2(nex), x3(nex), x4(nex), x5(nex), x6(nex)
      dimension y1(nex), y2(nex), y3(nex), y4(nex), y5(nex), y6(nex)
      dimension s1(nex), s2(nex), s3(nex), s4(nex), s5(nex), s6(nex)
      dimension t1(nex), t2(nex), t3(nex), t4(nex), t5(nex), t6(nex)
      character*52 filename,a
      character*8,parameter :: xmuname="/xmu.dat"

!     initialize all (1+chi) to 1
        chil2(:) = 1.d0
        chil3(:) = 1.d0
        chil4(:) = 1.d0
        chil5(:) = 1.d0

      !KJ 2014
      open(5,file="listedges.pmbse",form="formatted",status="old")

!     get xmu.dat for the dominant l-->l+1 channel
      filename = 'Oddp1/xmu.dat'
      read(5,*) filename
      filename=filename(1:istrln(filename))//xmuname
      call xmudat(filename, x1, x2, x3, x4, x5, x6, n, e1)
      deltaso = 0.d0
      if (kinit.lt.-1) then
!       get second edge xmu.dat for l-->l+1 channel
        filename = 'Evenp1/xmu.dat'
        read(5,*) filename
        filename=filename(1:istrln(filename))//xmuname
        call xmudat(filename, y1, y2, y3, y4, y5, y6, m, e2)
        deltaso = (e2-e1) / hart
        do i = 1, m
          y2(i) = y2(i) + deltaso
        enddo
!       get  xmu.dat files  for l-->l-1 channels
        filename = 'Oddm1/xmu.dat'
        read(5,*) filename
        filename=filename(1:istrln(filename))//xmuname
        call xmudat(filename, s1, s2, s3, s4, s5, s6, n, ef)
        filename = 'Evenm1/xmu.dat'
        read(5,*) filename
        filename=filename(1:istrln(filename))//xmuname
        call xmudat(filename, t1, t2, t3, t4, t5, t6, m, ef)
      endif
      close(5)

!     construct energy grid (ee) from 2 grids (x2 and y2)
      if (deltaso.eq. 0.d0) then
!       special case if deltaso = 0, then ee=x2
        ne = n
        do ie = 1, ne
          ee(ie) = x2(ie)
        enddo
      else
!       need to combine 2 grids
!       start with x grid until first y point
        do i = 1, n
          if (x2(i).lt.y2(1)) then
            ix = i
            ee(ix) = x2(ix)
          endif
        enddo
        ne = ix
        iy = 1
!       continue with x grid until xstep is less than ystep
  20    continue
          if ( ix.eq.n) goto 30
          if ( y2(iy).ge.e2) goto 30

          xstep = x2(ix+1)-x2(ix)
          ystep = y2(iy+1)-y2(iy)
          if (ystep.le.xstep) goto 30

          ne = ne+1
          ix = ix+1
          ee(ne) = x2(ix)
!         next y-point if the end of y-interval beyond the end of x-interval
          if (x2(ix+1).gt.y2(iy+1)) iy = iy + 1
        goto 20

  30    continue
!       continue with y-grid; keep nex as a maximum number of points
        np = min(m, nex-ne+iy-1)
        do i = iy, np
          ne = ne + 1
          ee(ne) = y2(i)
        enddo
      endif
!     end of energy mesh construction

!     interpolate fine structure on the energy grid ee
      do ie = 1, ne
        z1 = ee(ie)
        xnew = 1
        if (z1.ge.x2(1).and.z1.le.x2(n))    call terp (x2,x6,n,3,z1, chil3(ie))
        if (kinit.lt.-1) then
          if (z1.ge.y2(1).and.z1.le.y2(m))      call terp (y2,y6,m,3,z1, chil2(ie))
          if (z1.ge.s2(1).and.z1.le.s2(n))      call terp (s2,s6,n,3,z1, chil5(ie))
          if (z1.ge.y2(1).and.z1.le.y2(m))      call terp (y2,t6,m,3,z1, chil4(ie))
        endif
      enddo

!     for each edge find energy point (ix) with first nonzero chi,
!     padd the energy points below ix with chi(ix) to avoid jumps in xmu
!     notice that at high energies chi=0 ,
!     since extrapolation has not been used
      ix = 0
      do ie = 1, ne
        if (ix.eq.0 .and. abs( chil3(ie)-1.d0 ).ne.0.d0) ix = ie
      enddo
      do ie = 1,ix-1
        chil3(ie) = chil3(ix)
      enddo
      if (kinit.lt.-1) then
!       do the same with the rest of the data
        do ie = 1,ix-1
          chil5(ie) = chil5(ix)
        enddo
        ix = 0
        do ie = 1, ne
          if (ix.eq.0 .and. abs( chil2(ie)-1.d0 ).ne.0.d0) ix = ie
        enddo
        do ie = 1,ix-1
          chil2(ie) = chil2(ix)
          chil4(ie) = chil4(ix)
        enddo
      endif

      do ie = 1, ne
      write(18,*) ee(ie)*hart, chil3(ie), chil2(ie)
      enddo

      return
      end

      subroutine xmudat (filename, x1, x2, x3, x4, x5, x6, n, ef)
	  use dimsmod, only: nex 
	  use constants
      implicit double precision (a-h,o-z)
      dimension x1(nex), x2(nex), x3(nex), x4(nex), x5(nex), x6(nex)
      character*512 string
      character*52 filename

      n = 0
      write(*,*) 'reading file ',filename
      open (file=filename, unit=3, status='unknown',iostat=ios)
 10   read(3,'(a)', end=20) string
      if (string(4:10).eq.'-------') goto 30
      goto 10
       
!     check if xmu.dat file does not have the data
 20   stop ' xmu.dat is empty'

 30   continue
      read(3,'(a)',end=20) string

 40   n=n+1
      read(3,*,end=50)    x1(n), x2(n), x3(n), x4(n), x5(n), x6(n)
      if (x3(n).gt.-0.01d0 .and. x3(n).lt.0.01d0) ef = x1(n)
      x6(n) = x6(n)/x5(n) + 1
      x2(n) = x2(n) / hart
      goto 40

 50   continue
      n = n -1

      close (unit=3)
      return
      end
