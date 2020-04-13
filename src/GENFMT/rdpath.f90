!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdpath.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rdpath (in, done, ipol)

      use dimsmod, only: legtot
	  use constants
      use pdata
      use str
      implicit double precision (a-h, o-z)

      logical done

      complex*16  alph, gamm
      dimension  alpha(0:legtot), gamma(legtot)
      character*512 slog
      external dist

      read(in,*,end=200)  ipath, nleg, deg
      if (nleg .gt. legtot)  then
         write(slog,'(a,2i6)')                                          &
     &         ' nleg .gt. legtot, nleg, legtot ', nleg, legtot
         call wlog(slog)
         call wlog(' ERROR')
         goto 200
      endif
!     skip label (x y z ipot rleg beta eta)
      read(in,*)
      do 20  ileg = 1, nleg
         read(in,*,end=999)  (rat(j,ileg),j=1,3), ipot(ileg),           &
     &                       potlbl(ipot(ileg))
!        convert to code units
         do 10  j = 1, 3
            rat(j,ileg) = rat(j,ileg)/bohr
   10    continue
         if (ipot(ileg) .gt. npot)  then
            write(slog,'(a,3i8)')                                       &
     &              ' ipot(ileg) too big, ipot, ileg, npot ',           &
     &               ipot(ileg), ileg, npot
            call wlog(' ERROR')
            goto 200
         endif
   20 continue
      nsc = nleg-1

!     We need the 'z' atom so we can use it below.  Put
!     it in rat(nleg+1).  No physical significance, just a handy
!     place to put it.
      if (ipol.gt.0) then
         rat(1,nleg+1) = rat(1,nleg)
         rat(2,nleg+1) = rat(2,nleg)
         rat(3,nleg+1) = rat(3,nleg) + 1.0
      endif

!     add rat(0) and ipot(0) (makes writing output easier)
      do 22 j = 1, 3
         rat(j,0) = rat(j,nleg)
   22 continue
      ipot(0) = ipot(nleg)

      nangle = nleg
      if (ipol.gt.0) then 
!        in polarization case we need one more rotation
         nangle = nleg + 1
      endif
      do 100  j = 1, nangle

!        for euler angles at point i, need th and ph (theta and phi)
!        from rat(i+1)-rat(i)  and  thp and php
!        (theta prime and phi prime) from rat(i)-rat(i-1)
!
!        Actually, we need cos(th), sin(th), cos(phi), sin(phi) and
!        also for angles prime.  Call these  ct,  st,  cp,  sp

!        i = (j)
!        ip1 = (j+1)
!        im1 = (j-1)
!        except for special cases...
         ifix = 0
         if (j .eq. nsc+1)  then
!           j+1 'z' atom, j central atom, j-1 last path atom
            i = 0
            ip1 = 1
            if (ipol.gt.0) then
               ip1 = nleg+1
            endif
            im1 = nsc

         elseif (j .eq. nsc+2)  then
!           j central atom, j+1 first path atom, j-1 'z' atom
            i = 0
            ip1 = 1
            im1 = nleg+1
            ifix = 1
         else
            i = j
            ip1 = j+1
            im1 = j-1
         endif

         x = rat(1,ip1) - rat(1,i)
         y = rat(2,ip1) - rat(2,i)
         z = rat(3,ip1) - rat(3,i)
         call trig (x, y, z, ctp, stp, cpp, spp)
         x = rat(1,i) - rat(1,im1)
         y = rat(2,i) - rat(2,im1)
         z = rat(3,i) - rat(3,im1)
         call trig (x, y, z, ct, st, cp, sp)

!        Handle special case, j=central atom, j+1 first
!        path atom, j-1 is 'z' atom.  Need minus sign
!        for location of 'z' atom to get signs right.
         if (ifix .eq. 1)  then
            x = 0
            y = 0
            z = 1.0
            call trig (x, y, z, ct, st, cp, sp)
            ifix = 0
         endif

!        cppp = cos (phi prime - phi)
!        sppp = sin (phi prime - phi)
         cppp = cp*cpp + sp*spp
         sppp = spp*cp - cpp*sp
         phi  = atan2(sp,cp)
         phip = atan2(spp,cpp)

!        alph = exp(i alpha)  in ref eqs 18
!        beta = cos (beta)         
!        gamm = exp(i gamma)
         alph = -(st*ctp - ct*stp*cppp - coni*stp*sppp)
         beta(j) = ct*ctp + st*stp*cppp
!        watch out for roundoff errors
         if (beta(j) .lt. -1) beta(j) = -1
         if (beta(j) .gt.  1) beta(j) =  1
         gamm = -(st*ctp*cppp - ct*stp + coni*st*sppp)
         call arg(alph,phip-phi,alpha(j))
         beta(j) = acos(beta(j))
         call arg(gamm,phi-phi,gamma(j))
!       Convert from the rotation of FRAME used before to the rotation 
!       of VECTORS used in ref.
         dumm = alpha(j)
         alpha(j) =  pi- gamma(j)
         gamma(j) =  pi- dumm

         if (j .le. nleg)  then
            ri(j) = dist (rat(1,i), rat(1,im1))
         endif
  100 continue

!     Make eta(i) = alpha(i-1) + gamma(i). 
!     We'll need alph(nangle)=alph(0)
      alpha(0) = alpha(nangle)
      do 150  j = 1, nleg
         eta(j) = alpha(j-1) + gamma(j)
  150 continue
      if (ipol.gt.0) then
         eta(0) = gamma(nleg+1)
         eta(nleg+1) = alpha(nleg)
      endif

!     eta and beta in radians at this point.
      done = .false.
      return

!     If no more data, tell genfmt we're done
  200 continue
      done = .true.
      return

!     If unexpected end of file, die
  999 continue
      call wlog(' Unexpected end of file')
      call par_stop('ERROR')
      end
      subroutine trig (x, y, z, ct, st, cp, sp)
      implicit double precision (a-h, o-z)
!     returns cos(theta), sin(theta), cos(phi), sin(ph) for (x,y,z)
!     convention - if x=y=0 and z>0, phi=0, cp=1, sp=0
!                  if x=y=0 and z<0, phi=180, cp=-1,sp=0
!                - if x=y=z=0, theta=0, ct=1, st=0
      parameter (eps = 1.0e-6)
      r = sqrt (x**2 + y**2 + z**2)
      rxy = sqrt (x**2 + y**2)
      if (r .lt. eps)  then
         ct = 1
         st = 0
      else
         ct = z/r
         st = rxy/r
      endif
      if (rxy .lt. eps)  then
         cp = 1
         if (ct .lt. 0) cp = -1
         sp = 0
      else
         cp = x / rxy
         sp = y / rxy
      endif
      return
      end
      subroutine arg(c,fi,th)
      implicit double precision (a-h, o-z)
      complex*16  c
      parameter (eps = 1.0e-6)
      x = dble(c)
      y = dimag(c)
      if (abs(x) .lt. eps) x = 0
      if (abs(y) .lt. eps) y = 0
      if (abs(x) .lt. eps  .and.  abs(y) .lt. eps) then
        th = fi
      else
        th = atan2(y,x)
      endif
      return
      end
