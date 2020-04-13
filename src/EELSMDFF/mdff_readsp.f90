!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mdff_readsp.f90,v $:
! $Revision: 1.3 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mdff_readsp

        ! This routine opens xmu.dat and saves its spectrum (col. 4) into the array s in column ind+1.
        ! Column 1 contains the energy mesh, which is assumed to be identical for all stored spectra.
        use eels_inp,only: ipmin,ipmax,ipstep,aver,cross,iinput,spcol
        use mdff_spectrum
        use constants,only: pi
        implicit none
! LOCALS
      real*8, allocatable :: xmufile(:,:,:)
      real*8 dummy
        character*100 a
        integer i,j,nheader,ios,ip,ipsteplocal,ncols
        character*14 f1
        character*6 ext


      !construct filenames
      do ip=ipmin,ipmax,ipstep

        if(ip.eq.1) then  !in this block we make the file extension
          ext(1:6)='.dat  '	  
        elseif(ip.eq.10) then
          ext(1:6)='10.dat' 
        elseif(ip.gt.1.and.ip.lt.10) then
          ext(1:1)='0'
          ext(2:2)= char(48+ip)
          ext(3:6)='.dat'
        else
          stop 'crazy ip in rdst'
        endif
        
        if     (iinput.eq.1) then !read from xmu.dat-file               !here we make the whole filename
          call mdff_concat('xmu',ext,f1,i)
          ncols=6  !This file contains 6 columns ; the spectrum is in the 4th
!	  spcol=4
        elseif (iinput.eq.2) then ! read from opconsKK.dat-file
          call mdff_concat('opconsKK',ext,f1,i)
          ncols=8  !This file contains 8 columns ; the spectrum is in the 3rd/8th
!	  spcol=3
        else
          stop 'crazy iinput in readsp'
        endif

!       Open xmu.dat and read its contents into the array s.

      open (unit=8, file=f1, status='unknown', iostat=ios)
      call chopen (ios, f1, 'eels')

!     count the number of header lines : check for the first non-blank character and see if it's #.
      do i=1,2000000
        read(8,*) a
          do j=1,100
            if(a(j:j).ne.' ') exit
          enddo
          if(a(j:j).ne.'#' .and. a(j:j).ne.' ') exit  !trying to cope with blank lines
        enddo
        nheader=i-1   !This is the number of header lines.
        i=0
        j=0
        do while(j.eq.0)
         read(8,*,end=9,err=9) dummy
         i=i+1
        enddo
9       continue
        ne=i+1

        if (ip.eq.ipmin) then
          allocate(xmufile(ne,ncols,10))
          xmufile(:,:,:)=dble(0)
        endif

!     Now skip header and start reading spectrum at first valid position.
      rewind(8)
        do j=1,nheader
          read(8,*)
        enddo

        do i=1,ne
          read(8,*,err=900,end=900) (xmufile(i,j,ip),j=1,ncols)
        enddo
      close(8)
      
      enddo  !loop over ip
      
      call allocate_spectrum_1(ne,9)
      
! Here it becomes a little tricky.  EELS thinks like this :
! ipmin,max,step tells it which files need to be read.
! aver and cross tell it which files need to be used.      
      
     
      s(:,1)=xmufile(:,1,ipmin)  !save energy grid
      ipsteplocal=ipstep

      if(aver.eq.0) then  ! orientation sensitive
         if(ipmin.eq.1.and.ipmax.eq.9) then
            if(cross.eq.1.and.ipstep.ne.1) then
               call wlog('Asked for                                     &
     &         cross-term spectrum, but necessary files are missing.')
               stop
            elseif(cross.eq.0.and.ipstep.eq.1) then
               call wlog('Cross-term files found but asked to neglect   &
     &           them.  Calculation proceeds without cross-terms.')
               ipsteplocal=4
            elseif(cross.eq.0.and.ipstep.eq.4) then
               call wlog('Calculation without cross-terms.')
            else
               call wlog('Calculation using cross-terms.')
            endif
           
            do ip=ipmin,ipmax,ipsteplocal
               s(:,ip+1)=xmufile(:,spcol,ip)  !save total spectrum into s
            enddo
         else
            call wlog ('Asked for orientation sensitive spectrum,       &
     &       but conflicting ipmin or ipmax in eels.inp.')
            stop
         endif
      elseif(aver.eq.1) then	
         if(cross.eq.1) call wlog('Asked for cross-terms                &
     &      in averaged calculation.  Ignoring this - cross terms       &
     &      disappear when averaging.')
         if(ipmin.eq.10.and.ipmax.eq.10) then
            call wlog('Averaged calculation - using averaged            &
     &           spectrum from xmu10.dat.')
            s(:,2)=xmufile(:,spcol,10)
            s(:,6)=xmufile(:,spcol,10)
            s(:,10)=xmufile(:,spcol,10)  !copy averaged spectrum onto diagonal


         elseif(ipmin.eq.1.and.ipmax.eq.9) then
             call wlog('Averaged calculation - averaging xx,yy and zz   &
     &          spectra.')
             s(:,2)=(xmufile(:,spcol,1)+xmufile(:,spcol,5)+             &
     &                     xmufile(:,spcol,9))/3 !*4*pi/3
             s(:,6)=s(:,2)
             s(:,10)=s(:,2)
         else
            call wlog('Readsp is confused.  Check ipmin&max.')
            write(*,*) 'ipmin= ',ipmin,'ipmax= ',ipmax
            stop 'crazy ipmin ipmax in readsp'
         endif
      endif


      return
900   write(*,*) 'Error reading xmu.dat in savspe for ip =',ip
      stop
      end



