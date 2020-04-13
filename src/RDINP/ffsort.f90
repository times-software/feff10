!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ffsort.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2010/11/30 19:41:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ffsort (iabs,nss,doptz)
! KJ 1-06 : I added second input argument doptz
! KJ 7-06 : I added third input argument nss      

!     finds iabs-th atom of 'iphabs' type in file atoms.dat and writes
!     a smaller list of all atoms within 'rclabs' of that particular
!     absorber into 'geom.dat' file.
!      first coded by a.l.ankudinov, 1998 for CFAVERAGE card
!      modified by a.l.ankudinov, march 2001 for new i/o structure

      use par
	  use constants
	  use dimsmod, only: natx, nphx=>nphu
	  use geometry_inp
	  use global_inp
	  use potential_inp,only: iz
      implicit double precision (a-h, o-z)

!c    INPUT
      integer iabs
      logical,intent(in) :: doptz  !KJ 1-06 : call mkptz or not?
      integer,intent(in):: nss !KJ 7-06 : for sanity check at bottom of file	
!c    OUTPUT: geom.dat
      integer  nat
      integer iatph(0:nphx), iphat(natx), index(natx)
      double precision  rat(3,natx)

!     Local stuff
      real*8,parameter :: big = 1.0e5
      character*512 slog
	  integer iz1,iz2

      external dist

!     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)
  45  format ( 2i8, f13.5)
  50  format ( 3i5, 2f12.4, i5) !KJ

!     Find the first absorber (iphabs type) in a long list (iabs.le.0),
!     or find iabs-th atom in the list of type iphabs (iabs.gt.0)
      iatabs = 0
      icount = 0
      ifound = 0
      do iat = 1, natt
         if (iphatx(iat) .eq. 0) iphatx(iat) = iphabs
         if (iphatx(iat) .eq. iphabs) icount = icount +1
         if (ifound.eq.0 .and. icount.gt.0 .and. (icount.eq.iabs .or.   &
     &                          (iabs.le.0 .and. icount.eq.1))) then
            iatabs = iat
            ifound =1
         endif
      enddo

!     Make several sanity checks
      if (iatabs.eq.0 .and. natt.gt.1) then
         call wlog(' No absorbing atom (unique pot 0 or iphabs in CFAVERAGE  card) was defined.')
         call par_stop('RDINP')
      endif
      if (iphabs.eq.0 .and. icount.gt.1) then
         call wlog(' More than one absorbing atom (potential 0)')
         call wlog(' Only one absorbing atom allowed')
         call par_stop('RDINP')
      endif
      if ((icount.gt.0 .and. icount.lt.nabs) .or. nabs.le.0) then
         nabs = icount
         call wlog(' Averaging over ALL atoms of iphabs type')
      endif

!     Make absorbing atom first in the short list
      if (iatabs .ne. 0) then
         rat(:,1) = 0
         iphat(1) = 0
         index(1) = iatabs
      endif
          
!     make a smaller list of atoms from a big one
      nat = 1
      do iat = 1,natt
         if (iat.ne.iatabs) then
            tmp = dist (ratx(1,iat), ratx(1,iatabs))
            if (tmp.gt.0.1 .and. tmp.le.rclabs) then
               nat = nat + 1
               if (nat.gt.natx) then
                 write (slog, 307) nat, natx
  307            format (' Number of atoms', i6, 'exceeds max allowed for the pathfinder =', i6)
                 call wlog (' Use or reduce rclabs in CFAVERAGE card')
                 call wlog (' Or increase parameter natx and recompile')
                 call par_stop('RDINP')
               endif
               rat(1,nat) = ratx(1,iat)-ratx(1,iatabs)
               rat(2,nat) = ratx(2,iat)-ratx(2,iatabs)
               rat(3,nat) = ratx(3,iat)-ratx(3,iatabs)
               iphat(nat) = iphatx(iat)
               index(nat) = iat
            endif
         endif
      enddo
!     sort atoms by distance
      do 315 iat = 1,nat-1
        r2min = rat(1,iat)**2 + rat(2,iat)**2 + rat(3,iat)**2
        imin = iat
        do 310 i = iat+1,nat
          r2 = rat(1,i)**2 + rat(2,i)**2 + rat(3,i)**2
          if (r2.lt.r2min) then
            r2min = r2
            imin = i
          endif
 310    continue
        if (imin.ne.iat) then
!         permute coordinates for atoms iat and imin
          do 311 i = 1,3
            r2 = rat(i,iat)
            rat(i,iat) = rat(i,imin)
            rat(i,imin) = r2
 311      continue
          i = iphat(iat)
          iphat(iat) = iphat(imin)
          iphat(imin) = i
          i = index(iat)
          index(iat) = index(imin)
          index(imin) = i
        endif
 315  enddo


!KJ 7-09 NOTE : Aleksi added following comment here in feff8q, which I don't understand :
!c     Bogus comment below. q-vector along z-axis and "polarization" along x+y vector. Not really used. 


!     rotate xyz frame for the most convenience and make polarization tensor
!     make polarization tensor when z-axis is along k-vector 
      if (doptz) call mkptz(nat,rat) !KJ I added the if-statement 1-06
!     rewrite global.inp for initial iteration to update 'ptz'
        call global_write(.false.) !Don't recalculate the norm of vectors - they're lost, since mkptz normalized everything ...

!     Find model atoms for unique pots that have them
!     Use atom closest to absorber for model
      iatph(:) = 0
!     By construction absorbing atom is first in the list
      iatph(0) = 1
      nph = 0
      do iph = 1, nphx
         rabs = big
         do iat = 2, nat
            if (iph .eq. iphat(iat))  then
               tmp = dist (rat(1,iat), rat(1,1))
               if (tmp .lt. rabs)  then
!                 this is the closest so far
                  rabs = tmp
                  iatph(iph) = iat
               endif
            endif
         enddo
         if (iatph(iph).gt.0) nph = iph
      enddo
!     if iatph > 0, a model atom has been found.

!     Check if 2 atoms are closer together than 1.75 bohr (~.93 Ang)
      ratmin = 1.0e20
      do iat = 1, nat
         do jat = iat+1, nat
            rtmp = dist(rat(1,iat),rat(1,jat))
            if (rtmp .lt. ratmin)  ratmin = rtmp
            if (rtmp .lt. 1.75 * bohr)  then
			
			   ! iat:  position in the list ordered by distance to absorber
			   ! iatx: position in the unordered list as entered in feff.inp
			   ! (in practice, many feff.inp have been ordered by another application and the 2 are the same)
               iatx = index(iat)
               jatx = index(jat)
			   iz1=iz(iphat(iat))
			   iz2=iz(iphat(jat))
			   if(iz1.ne.1 .or. iz2.ne.1 .or. rtmp.lt. 0.70) then
			      !KJ 11/2010:
				  ! added distance and atomic number to warning message (duh)
				  ! separate threshold for H-H bond (which is 0.74A = 1.4 bohr long)
                  call wlog(' :WARNING  TWO ATOMS VERY CLOSE TOGETHER. CHECK INPUT.')
                  write(slog,'(a,2i8,a,e13.5,a)') ' atoms ', iatx, jatx,' distance ',rtmp,' Angstrom'
                  call wlog(slog)
                  write(slog,'(i5,1p,3e13.5,a,i4)') iatx, (ratx(i,iatx),i=1,3), ' Z=',iz1
                  call wlog(slog)
                  write(slog,'(i5,1p,3e13.5,a,i4)') jatx, (ratx(i,jatx),i=1,3), ' Z=',iz2
                  call wlog(slog)
               endif
            endif
         enddo
      enddo

!     Write output geom.dat
      open (file='geom.dat', unit=3, status='unknown',iostat=ios)
        write (3,535) nat, nph
  535   format ('nat, nph = ', 2i5)
        write (3,516) (iatph(iph), iph=0,nph)
  516   format(16i5)
        write (3, 10) ' iat     x       y        z       iph  '
        write (3, 526)
  526   format (1x, 71('-'))
        ibounc = 1
        do 540  i = 1, nat
          write(3,536) i, rat(1,i), rat(2,i), rat(3,i), iphat(i), ibounc
  536     format(i4, 3f13.5, 2i4)
  540   continue
      close(3)

!     Atoms for the pathfinder
      if (iatabs.le.0 .and. nss.le.0 .and. nat.gt.0 )  then !KJ 7-06 added second and third condition
         call wlog(' Absorbing atom coords not specified.')
         call wlog(' Cannot find multiple scattering paths.')
         call par_stop('RDINP')
      endif

! 400 call par_barrier

      return
      end
