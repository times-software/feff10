!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mkptz.f90,v $:
! $Revision: 1.7 $
! $Author: jorissen $
! $Date: 2011/06/25 00:03:25 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine mkptz (nat,rat)
!     choose new right handed frame of reference with z along spvec,
!      y along (xivec cross spvec); simpler choice if one of them is 0.
!     get all vectors in new frame and
!     makes polarization tensor ptz when z is rotated along k-vector

!     input:
!     ipol = 0  random k-vector orientation in 3d; ptz(i,j)=\delta_{i,j}
!     ipol = 1 for polarizion vector eps and it's  complex conjugate epc
!        ptz(j,i) = 0.5 [(eps(-i))^* eps(-j) + (epc(-i))^* epc(-j)]
!        notice that complex conjugation and taking i-th component
!        are non commuting operations. (eps(-i))^* = (-)^i (epc(i))
!     ipol = 2 ptz(i,j)= i*\delta_{i,j}
!     elpty - ellipticity (optional for ipol=1)
!     xivec - direction of x-ray propagation
!     ispin - type of spin calculations
!        0 - spin independent
!        -1,1 - spin dependent potential
!        2 - calculations with spin-up potential
!       -2 - calculations with spin-down potential
!     spvec - direction of spin vector (along z at the output)
!     nat - number of atoms
!     rat - xyz cordinates of atoms (changed due to the rotations)

!     output:
!     angks - angle between k-vector and spin-vector (0-pi)
!     le2   - 0-only E1, 1-E1+M1, 2-E1+E2, 3-E1+E2+M1 transitions
!     ptz   - polarization tensor
      use constants
      use global_inp,only : ipol,elpty,evec,xivec,spvec,ispin,le2,angks,ptz,do_nrixs
      implicit double precision (a-h, o-z)

      integer,intent(in) :: nat
      double precision rat(3,nat)

!     additional local stuff to create polarization tensor ptz(i,j)
      dimension e2(3)
      complex*16  e(3),eps,epc
      dimension eps(-1:1),epc(-1:1)
      character*512 slog


!     make z axis along propagation (XIVEC).
!     le2=0 - only E1 transitions; le2=1 - E1+M1; le2=2 - E1+E2 
      rr = xivec(1)**2 + xivec(2)**2 + xivec(3)**2
      if (rr.eq.0) then
        angks = 0
!       special case when xivec is not specified
        if (ipol.eq.1) then
!         need to know xivec for E2 and M1 transitions
!         leave only E1 contribution
          if (le2.ne.0) call wlog('  Can do only E1 transitions. Specify k-vector for M1 or E2')
          le2 = 0
        else
!         for polarization average of circular dichroizm
          if (ispin.ne.0) then
!           spin-dependent case
            do 10 i = 1,3
  10        xivec(i) = spvec(i)
            rr = xivec(1)**2 + xivec(2)**2 + xivec(3)**2
          endif
        endif
      endif
            
              
!      if (rr.gt.0 .and. (.not.(do_nrixs.eq.1))) then
      if (rr.gt.0 ) then
         rsp = sqrt(rr)
         rr = xivec(1)**2 + xivec(2)**2
         if ( rr.ne.0 .or. xivec(3).lt.0) then
           if (rr.eq. 0) then
             cst = - 1
             snt = 0
             csf = 1
             snf = 0
           else
!            rotation is defined by angles theta and fi
             rr = sqrt(rr)
             cst = xivec(3) / rsp
             snt = rr / rsp
             csf = xivec(1) / rr
             snf = xivec(2) / rr
           endif
!          rotate all vectors
           do 20 i = 1, nat
 20        call rotate (rat(1,i), cst, snt, csf, snf)
           call rotate (evec, cst, snt, csf, snf)
           call rotate (xivec, cst, snt, csf, snf)
           call rotate (spvec, cst, snt, csf, snf)
         endif
      endif


!KJ 7-09 following code needs to run for feff8q
      if(do_nrixs.eq.1) then
!        Now the momentum transfer is along the z-axis 
!        Let us define a fake polarization vector that is perpendicular to it. 
         ipol=1
         evec(1)=0.707d0
         evec(2)=0.707d0
         evec(3)=0.0d0
      endif
!KJ

!     initialize ptz
      ptz(:,:) = 0

!     make ptz in the frame when z is along xivec, except ipol=0
      if (ipol .eq. 0) then
         do 40 i=-1,1
 40      ptz(i,i) = 1.d0 /3.d0
      elseif (ipol .eq. 2) then
         ptz( 1, 1) =  1.d0
         ptz(-1,-1) = -1.d0
      elseif (ipol .eq. 1) then
!       Normalize polarization vector
        x = sqrt (evec(1)**2 + evec(2)**2 + evec(3)**2)
        if (x .le. 0.000001) then
         call wlog(' STOP  Polarization vector of almost zero length.  Correct POLARIZATION card.')
         call par_stop('MKPTZ-1')
        endif
        do 50  i = 1, 3
         evec(i) = evec(i) / x
  50    continue
        x = sqrt (xivec(1)**2 + xivec(2)**2 + xivec(3)**2)
        if (x .gt. 0) then
!         run elliptical polarization code
          do 60  i = 1, 3
            xivec(i) = xivec(i) / x
  60      continue
          x = evec(1)*xivec(1)+evec(2)*xivec(2)+evec(3)*xivec(3)
          if (abs(x) .gt. 0.9) then
            call wlog(' polarization')
            write(slog,292)  (evec(i), i=1,3)
            call wlog(slog)
            call wlog(' incidence')
            write(slog,292) (xivec(i), i=1,3)
            call wlog(slog)
            call wlog(' dot product')
            write(slog,292)  x
            call wlog(slog)
  292       format (5x, 1p, 2e13.5)
            call wlog(' STOP polarization almost parallel to the incidence.')
            call wlog(' Correct ELLIPTICITY and POLARIZATION cards.')
            call par_stop('MKPTZ-2')
          endif
		  !KJ June 2013 : used to be a "if x .ne. 0.0" - holy fucking you-know-what!  Produced warnings for 10^-13 numerical noise ...
          if (dabs(x) .gt. 0.00001) then
!           if xivec not normal to evec then make in normal, keeping the
!           plane based on two vectors
            call wlog(' Changing polarization vector! Incidence is not normal to polarization.')
            call wlog(' Check your input for errors. Run continues.')
			write(*,*) 'evec . xivec=',x
			write(*,*) 'xivec=',xivec
			write(*,*) 'evec=',evec
            do 70  i = 1,3
              evec(i) = evec(i) - x*xivec(i)
  70        continue
            x = sqrt (evec(1)**2 + evec(2)**2 + evec(3)**2)
            do 80   i = 1, 3
               evec(i) = evec(i) / x
  80        continue
          endif
        else
!         elpty cannot be used with xivec=0
          elpty = 0.0
        endif 
     
        e2(1) = xivec(2)*evec(3)-xivec(3)*evec(2)
        e2(2) = xivec(3)*evec(1)-xivec(1)*evec(3)
        e2(3) = xivec(1)*evec(2)-xivec(2)*evec(1)
        do 90   i = 1,3
          e(i) = (evec(i)+elpty*e2(i)*coni)
  90    continue 
        eps(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
        eps(0)  =   e(3)
        eps(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
        do 100  i = 1,3
          e(i) = (evec(i)-elpty*e2(i)*coni)
  100   continue 
        epc(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
        epc(0)  =   e(3)
        epc(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
        do 110 i = -1,1
        do 110 j = -1,1
!         ptz(j,i) = (-1.0)**i * epc(i)*eps(-j) / (1+elpty**2)
!         above - true polarization tensor for given ellipticity, 
!         below - average over left and right in order to have
!         path reversal symmetry
          ptz(j,i) = ((-1.0)**i)*(epc(i)*eps(-j)+eps(i)*epc(-j))        &
     &               /(1+elpty**2)/2.0
  110   continue
      endif
!     end of making polarization tensor

      angks = 0


!     second rotate so that z parallel to spin
!     note that new y-axis is normal to spin AND incidence vector
!     which simplifies further expression for rotation matrix
      rr = spvec(1)**2 + spvec(2)**2 + spvec(3)**2
      if (rr.gt.0) then
         rsp = sqrt(rr)
         rr = spvec(1)**2 + spvec(2)**2
         if ( rr.ne.0 .or. spvec(3).lt.0) then
           if (rr.eq. 0) then
             cst = - 1
             snt = 0
             csf = 1
             snf = 0
             angks = pi
           else
!            rotation is defined by angles theta and fi
             rr = sqrt(rr)
             cst = spvec(3) / rsp
             snt = rr / rsp
             csf = spvec(1) / rr
             snf = spvec(2) / rr
             angks = acos( cst)
           endif
!          rotate all vectors
           do 120 i = 1, nat
 120       call rotate (rat(1,i), cst, snt, csf, snf)
           call rotate (evec, cst, snt, csf, snf)
           call rotate (xivec, cst, snt, csf, snf)
         endif
      endif

      return
      end

      subroutine rotate (vec, cst, snt, csf, snf)
      implicit double precision (a-h, o-z)
!     rotates vector to a new coordinate system
!     Euler angles: alpha=phi, beta=theta, gamma=0
      dimension vec(3), temp (3)

      temp(1) = vec(1)*cst*csf + vec(2)*cst*snf - vec(3)*snt
      temp(2) = -vec(1)*snf + vec(2)*csf
      temp(3) = vec(1)*csf*snt + vec(2)*snt*snf + vec(3)*cst
      do 10 i = 1,3
  10  vec(i) = temp(i)

      return
      end
