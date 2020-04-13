!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mpprmp.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2011/03/30 04:50:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mpprmp (npat, ipat, xp, yp, zp, ipol, ispin, evec, xivec,ica)   !KJ added ica 5/06

!     make path parameters,  xp, yp,zp for each atom for a given
!     path. The allowed symmetry operations are restricted by
!     polarization type ipol (evec, xivec) and spin type (ispin)
!     of calculations
!      ipol=0 - polarization average 
!      ipol=1 - linear (admixture of 2 linear for elpty.ne.0)
!      ipol=2 - circular dichroism
!      ispin=0 - spin-independent system V_up=V_dn=V_av
!      |ispin|=1 - V_up .ne. V_dn, sum over up and down calculations
!      ispin= 2 - V_up portion of |ispin|=1 
!      ispin=-2 - V_dn portion of |ispin|=1 
!    all possible cases fall into 7 categories of allowed symmetries
!    Spin-independent calculations 
!     1) IF ipol=0 
!        any path rotation, reflection and reversal are allowed
!     2) ELSEIF ispin=0 ipol=1  xivec.eq.0 (dipole transitions only)
!        any rotation around evec, reflections in planes normal
!        and parallel to evec, path reversal
!     3) ELSEIF ispin=0 ipol=1   xivec.ne.0
!        reflections in 2 planes (xivec, evec) and (xivec, B field)
!        reflection in (evec, B field) probably does not conserve
!        E1-E2 cross term (currently used; check and fix later)
!     4) ELSEIF ispin=0 ipol=2  
!        rotations around xivec, path reversal (? -check for XNCD)
!    Spin systems (only ispin.ne.0, ipol.ne.0 below)
!     5) ELSEIF  xivec(1)=xivec(2)=0 .and. ( ipol.eq.2 .or. 
!                 ( ipol.eq.1 .and.xivec(3)=evec(1)=evec(2)=0 ) )
!        rotations around spin axis
!     6) ELSEIF ipol=1 xivec(1)=xivec(2)=evec(3)=0 (XMLD)
!        only 180 degrees rotation around spin axis
!     7) ELSE   ipol=1,2 .and. (xivec(1).ne.0 or xivec(2).ne.0 )
!        NO symmetry operations
!    Disclaimer: the symmetry rules for might be
!    too restrictive and were checked for ipol=2 calculations (case=5)
!    One can always check the symmetry rules by comparing with case=7.

!    To exploit above symmetry, every path is recorded in a new frame of
!    reference, constructed for the calculations specified.

!     Input is list of atoms (npat, ipat(npat)), output are
!     x,y,z coord. of path in standard frame of reference
!     (see comments in timrep.f or here below)

      use dimsmod, only: natx, npatx
      double precision evec(3), xivec(3)
      double precision  ro2, norm, zvec, xvec, yvec, ri, xp1, yp1, zp1
      dimension ipat(npatx+1), zvec(3), xvec(3), yvec(3)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      dimension xp(npatx), yp(npatx), zp(npatx)
      dimension xp1(npatx), yp1(npatx), zp1(npatx)
      dimension ri(3,npatx)
      logical lkvec, lkz, lez
      integer ica !KJ new input parameter. overrides icase if positive.

      parameter (eps4 = 1.0E-4)

!     the atoms in this path are passed via the ipat and rat arrays




      if (ica.gt.0.and.ica.lt.8) then !KJ added this and next 2 lines
!KJ Purpose : for eels calculations, we want a result for many k and q.
!KJ This means we cannot use the symmetry operations.
         icase=ica
      else
      
!        determine which case we are dealing with (see above comments)
!        default is icase =7 - no symmetry at all (safe for untested cases)
!        logical lkvec answers whether xivec is a vector
         lkvec = .false.
         if (xivec(1)**2+xivec(2)**2+xivec(3)**2.gt.eps4) lkvec = .true.
!        logical lkz answers whether xivec is a vector along z at most
         lkz = .true.
         if (xivec(1)**2+xivec(2)**2.gt.eps4) lkz = .false.
!        logical lez answers whether evec is a vector along z at most
         lez = .true.
         if (evec(1)**2+evec(2)**2.gt.eps4) lez = .false.

         icase = 7
         if (ipol.eq.0) then
           icase = 1
         elseif (ispin.eq.0) then
           if (ipol.eq.1 .and. (.not.lkvec)) icase = 2
           if (ipol.eq.1 .and. lkvec) icase = 3
           if (ipol.eq.2) icase=4
         else
           if (ipol.eq.2 .and. lkz) icase = 5
           if (ipol.eq.1 .and. .not.lkvec  .and. lez) icase = 5
           if (ipol.eq.1 .and. lkz .and. evec(3)**2.lt.eps4) icase = 6
         endif
 
       endif  !KJ my block. 5/06

	
!     initialize staff
      xp = 0
      yp = 0
      zp = 0
      xp1 = 0
      yp1 = 0
      zp1 = 0
      nleg = npat + 1
      do 20  j = 1, npat
      do 20  i = 1, 3
         ri(i,j) = rat(i,ipat(j)) - rat(i,0)
   20 continue
      ri(:,nleg:npatx) = 0
	  xvec = 0
	  yvec = 0
	  zvec = 0

      if (icase.eq.1) then
!        z-axis along first leg
         norm = ri(1,1)*ri(1,1)+ri(2,1)*ri(2,1)+ri(3,1)*ri(3,1)
         norm = sqrt(norm)
         do 140 i = 1, 3
           zvec(i) = ri(i,1)/norm
  140    continue
      elseif (icase.eq.2 .or. icase.eq.3) then
!        z-axis in direction of polarization
         do 120 i = 1, 3
           zvec(i) = evec(i)
  120    continue
      else 
!        keep z-axis
         zvec(3) = 1.d0
      endif

      do 160 j = 1,npat
      do 160 i = 1, 3
        zp1(j) = zp1(j) + zvec(i)*ri(i,j)
  160 continue

!     if no symmetries, don't waste time
      if (icase.eq.7) then
         xvec(1) = 1.d0
         yvec(2) = 1.d0
         goto 390
      endif

      if (icase.eq.1 .or. icase.ge.4) goto 240
!     use z-->-z symmetry 
      num = 1
  200 continue
      if (abs(zp1(num)) .gt. eps4) then
         if (zp1(num) .lt. 0) then
!           inverse all z-coordinates and zvec, if 
!           first nonzero z-coordinate is negative 
            do 210 j = 1, 3
               zvec(j) = - zvec(j)
  210       continue
            do 220 j = 1, npat
               zp1(j) = - zp1(j)
  220       continue
         endif
         goto 240
      endif
      num = num +1
      if (num .lt. nleg) then
         goto 200
      endif
!     z--> -z symmetry has been used
  240 continue

!     use rotations around z and reflections containing z
      num = 1
  300 continue
      ro2 = 0
      do i =1, 3
         ro2 = ro2 + ri(i,num)*ri(i,num)
      enddo
!     looking for first atom which is not on z-axis
      ro2 = ro2 - zp1(num)*zp1(num)
      ro2 = sqrt(abs(ro2))
      if (ro2 .ge. eps4) then
!     if atom not on the z-axis then
         if (icase.eq.1.or.icase.eq.2.or.icase.eq.4.or.icase.eq.5) then
!           any rotation around z is allowed
!           choose x-axis so that x-coord. positive and y=0.
            do i = 1, 3
               xvec(i) = ri(i,num) - zvec(i)*zp1(num)
            enddo
            do i = 1, 3
               xvec(i) = xvec(i)/ro2
            enddo
         elseif (icase.eq.3) then
!           if elliptical polarization then
!           choose x-axis along incident beam
            do i =1, 3
               xvec(i) = xivec(i)
            enddo
         else
!           icase.eq.6 choose x-axis so that x-coord is positive
            xvec(1) = 1.d0
            if (ri(1,num).lt.0) xvec(1) = -1.d0
         endif
         yvec(1) = zvec(2)*xvec(3) - zvec(3)*xvec(2)
         yvec(2) = zvec(3)*xvec(1) - zvec(1)*xvec(3)
         yvec(3) = zvec(1)*xvec(2) - zvec(2)*xvec(1)
         goto 390
      endif
      num = num + 1
      if (num .lt. nleg) then
         goto 300
      endif

  390 continue

!     calculate x,y coord for each atom in chosen frame of reference
      do 400 j = 1, npat
      do 400 i =1,3
         xp1(j) = xp1(j) + xvec(i)*ri(i,j)
         yp1(j) = yp1(j) + yvec(i)*ri(i,j)
  400 continue

      if (icase.eq.3) then
!        check that first nonzero  x-coordinate is positive,
!        no need to check it in other cases.
         num = 1
  500    continue
         if (abs(xp1(num)) .ge. eps4) then
            if (xp1(num) .lt. 0) then
               do 510 j = 1, npat
                  xp1(j) = - xp1(j)
  510          continue
            endif
            goto 520
         endif
         num = num + 1
         if (num .lt. nleg) then
            goto 500
         endif
  520    continue
      endif

      if (icase.ge.4) goto 590
      num = 1
  570 continue
!     inverse all y-coordinates if first nonzero y-coord is negative
      if (abs(yp1(num)) .ge. eps4) then
         if (yp1(num) .lt. 0) then
            do 580 j = 1, npat
               yp1(j) = - yp1(j)
  580       continue
         endif
         goto 590
      endif
      num = num + 1
      if (num .lt. nleg) then
         goto 570
      endif
  590 continue

      do 595 j = 1, npat
        xp(j) = xp1(j)
        yp(j) = yp1(j)
        zp(j) = zp1(j)
  595 continue
!     now xp,yp,zp represent the path in standard order
      return
      end


!     icase = 7
!!     initialize staff
!     xp = 0
!	yp = 0
!	zp = 0
!	xp1 = 0
!	yp1 = 0
!	zp1 = 0
!     nleg = npat + 1
!     do 20  j = 1, npat
!     do 20  i = 1, 3
!        ri(i,j) = rat(i,ipat(j)) - rat(i,0)
!  20 continue
!        ri(:,nleg:npatx) = 0
!     do 40 i =1, 3
!        xvec(i) = 0
!         yvec(i) = 0
!        zvec(i) = 0
!  40 continue
!!
!!
!!
!!
!     if (elpty.ge.0.0d0) then 
!!
!!     evec is assumed along z-axis
!!
!      zvec(3) = 1.0d0
!     else
!!
!!     make the z-axis to be along the first leg
!!     
!        norm = ri(1,1)*ri(1,1)+ri(2,1)*ri(2,1)+ri(3,1)*ri(3,1)
!        norm = sqrt(norm)
!        do 140 i = 1, 3
!           zvec(i) = ri(i,1)/norm
!140     continue
!        
!     end if

!     do 160 j = 1,npat
!     do 160 i = 1, 3
!       zp1(j) = zp1(j) + zvec(i)*ri(i,j)
! 160 continue

!!
!! assuming here -q result equal to q, have to check
!!
!
!!     use z-->-z symmetry 
!     num = 1
! 200 continue
!     if (abs(zp1(num)) .gt. eps4) then
!        if (zp1(num) .lt. 0) then
!!           inverse all z-coordinates and zvec, if 
!!           first nonzero z-coordinate is negative 
!           do 210 j = 1, 3
!              zvec(j) = - zvec(j)
! 210       continue
!           do 220 j = 1, npat
!              zp1(j) = - zp1(j)
! 220       continue
!        endif
!        goto 240
!     endif
!     num = num +1
!     if (num .lt. nleg) then
!        goto 200
!     endif
!!     z--> -z symmetry has been used
! 240 continue

!!     use rotations around z and reflections containing z
!     num = 1
! 300 continue
!     ro2 = 0
!     do 310 i =1, 3
!        ro2 = ro2 + ri(i,num)*ri(i,num)
! 310 continue
!!     looking for first atom which is not on z-axis
!     ro2 = ro2 - zp1(num)*zp1(num)
!     ro2 = sqrt(abs(ro2))
!     if (ro2 .ge. eps4) then
!!     if atom not on the z-axis then
!
!!           any rotation around z is allowed
!!           choose x-axis so that x-coord. positive and y=0.
!           do 320 i = 1, 3
!              xvec(i) = ri(i,num) - zvec(i)*zp1(num)
! 320       continue
!           do 330 i = 1, 3
!              xvec(i) = xvec(i)/ro2
! 330       continue
!        yvec(1) = zvec(2)*xvec(3) - zvec(3)*xvec(2)
!        yvec(2) = zvec(3)*xvec(1) - zvec(1)*xvec(3)
!        yvec(3) = zvec(1)*xvec(2) - zvec(2)*xvec(1)
!        goto 390
!     endif
!     num = num + 1
!     if (num .lt. nleg) then
!        goto 300
!     endif

! 390 continue

!!     calculate x,y coord for each atom in chosen frame of reference
!     do 400 j = 1, npat
!     do 400 i =1,3
!        xp1(j) = xp1(j) + xvec(i)*ri(i,j)
!        yp1(j) = yp1(j) + yvec(i)*ri(i,j)
! 400 continue
!     if (icase.ge.4) goto 590
!     num = 1
! 570 continue
!!     inverse all y-coordinates if first nonzero y-coord is negative
!     if (abs(yp1(num)) .ge. eps4) then
!        if (yp1(num) .lt. 0) then
!           do 580 j = 1, npat
!              yp1(j) = - yp1(j)
! 580       continue
!        endif
!        goto 590
!     endif
!     num = num + 1
!     if (num .lt. nleg) then
!        goto 570
!     endif
! 590 continue

!     do 595 j = 1, npat
!       xp(j) = xp1(j)
!       yp(j) = yp1(j)
!       zp(j) = zp1(j)
! 595 continue
!!     now xp,yp,zp represent the path in standard order
!     return
!     end
