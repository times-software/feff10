!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: phmeshjas.f90,v $:
! $Revision: 1.5 $
! $Author: jorissen $
! $Date: 2010/11/30 19:41:54 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     make e mesh for phase
!     input:  iprint, ispec, edge, vi0, gamach, xkmax, xkstep
!     output: ne, ne1, em(ne), ik0 [grid point with k=0]
!             ne -  total number of points in array em
!             ne1 - number of points on horizontal grid 
!
!     JAS made a great equlizer version of the original.
!     constant energy step for all 
!

      subroutine phmeshjas (iprint, ispec, edge, emu, vi0, gamach, ecv, &
                       xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3)
	use constants
	use dimsmod, only: nex 
      implicit double precision (a-h, o-z)
      complex*16 em(nex), tempc
      double precision emax,emin,de
      integer net
      integer iecrit(9),ii
      integer iold,inew
      double precision eold,krange(9)
      integer ikcheck,ikk
!
!     Aleksi tries to make a constant delta e for all 

       external getxk
!
!     These needed special points in the original.
!     At this point JAS does not want to figure out 
!     why. Just check and use the old phmesh for these
!
       IF(.FALSE.) THEN
          PRINT*, 'nex', nex
          PRINT*, 'xkmax', xkmax
          PRINT*, 'ispec', ispec
       END IF
       krange(1)=0.d0
       krange(2)=0.5123d0*bohr
       krange(3)=1.0123d0*bohr
       krange(4)=1.5123d0*bohr
       krange(5)=2.0123d0*bohr
       krange(6)=3.0123d0*bohr
       krange(7)=4.0123d0*bohr
       krange(8)=5.0123d0*bohr
       krange(9)=6.0123d0*bohr
       ikcheck=1
       do ie=1,9
          if (abs(xkmax).gt.krange(ie)) then
             ikcheck=ie
          end if 
       end do
!
!     this is hom many points are needed after kmax
!
       ikcheck=9-ikcheck

       if (ispec.eq.2 .or. ispec.ge.3) then 
          write(6,*) "The constant energy step is not supposed to be used with DANES of XES."
          call phmesh (iprint, ispec, edge, emu, vi0, gamach, ecv, &
              xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3)
          return
       endif
       
       xloss = gamach/2 + vi0
       if (xloss.lt.0) xloss = 0
       xvert = max(xloss, 0.02/hart)
       xloss = xvert
       aa = 0.5d0
       ne3 = 0
       xim = xloss*aa
       if (vixan.gt.0.0001) xim = vixan

!     nemax - max number of points on horizontal axis

      nemax =0
      if (ispec.le.3)  then
!        make energy mesh for XANES with FMS
!        around fermi level step is regular in energy (xloss/2)
!        and regular in k at high energies

!        10 points below Fermi level
         nemax = 10
         dk = 2*xkstep
         n1 = int (xim/2/dk**2)
         n2 = int ( sqrt(n1*2*xim) / dk )
         
         if ( (dk*(n2+1))**2 .gt. (n1+1)*2*xim ) n1 = n1+1
         n1 = min (n1,nemax)
         do i = 1, n1
            em(nemax+1-i) = -xim*i + edge + coni*xloss
         end do
         nb = nemax-n1
         do i = 1, nb
            em(nb + 1 -i) = -(dk*(n2+i))**2/2 + edge + coni*xloss
         end do
         nmin = nemax
         ik0 = nemax+1
      endif
!
!
!
      
      iecrit(1) = ik0
      iecrit(2) = ik0 + 5
      iecrit(3) = ik0 + 10
      iecrit(4) = ik0 + 15
      iecrit(5) = ik0 + 20
      iecrit(6) = ik0 + 30
      iecrit(7) = ik0 + 34
      iecrit(8) = ik0 + 38
      iecrit(9) = ik0 + 40
      
      if (ispec .gt. 0 .and. ispec.le.3)  then
!        make energy mesh for XANES with FMS
!        around fermi level step is regular in energy (xloss/2)
!        and regular in k at high energies
!        90 points above Fermi level
         nemax = nex - nemax
         emax = xkmax*xkmax
         emax=xkmax*xkmax/2.0d0
         emin=0.0d0
         de = (emax-emin)/dble(nemax)
         n1=nemax
         do 30 i = 1, n1
  30         em(nmin+i) = de*(i-1)

         ne1 = nmin+n1
         do i = ik0, ne1
            em(i) = em(i) + edge + coni*xloss
         end do
      else
!        energy mesh for EXAFS or XANES without FMS
!        20 pts (0 le k le 1.9, delk=0.1 ang(-1) )
!        20 pts (2 le k le 5.8, delk=0.2 ang(-1) )
!         9 pts (6 le k le 10., delk=0.5 ang(-1) )
!        10 pts (11 le k le 20.0, delk=1.0 ang(-1) )
         ne = 0
         ik0=1
         if (ispec.lt.0) ne = 10
         nemax = nex-50-ne
         delk = bohr/10
         emax = xkmax*xkmax/2.0d0
         emin=0.0d0
         de = (emax-emin)/dble(nemax-9)
         n1=nemax
!
!     This is 'too' complicated method to make 
!     sure that the 'correct' k-values are
!     1. included into the grid
!     2. at points given by iecrit written to emesh.dat
! 
!
         iold=0
         i=1
         ii=2
         iecrit(1)=ik0
         do while (i.le.n1-ikcheck)
            inew=iold+1
            tempk=sqrt(2.0d0*de*dble(inew-1))
            ikk=-1
            if (ii.le.9) then 
               if (tempk.gt.krange(ii)) then
                  ikk=1
               end if
            end if
            if (ikk.gt.0) then
!     
!     next energy point is the k^2/2+edge
!     
               iecrit(ii)=ne+i
               em(ne+i)=krange(ii)*krange(ii)/2.0d0+edge+coni*xloss
               ii=ii+1
            else
               em(ne+i)=de*dble(inew-1)+edge+coni*xloss
               iold=inew
            end if
            i=i+1
         end do
!
!     It can happen 
!     1. that the highest points 
!        iecrit have not been found. 
!         So we have to replace the last points with 
!     2. Last point is less than kmax*kmax/2.0d0
!
!         inew=iold+1
!         tempk=sqrt(2.0d0*de*dble(inew-1))
         ii=9-ikcheck+1
         iold=n1-ikcheck
         do i=1,ikcheck
            iold=iold+1
            em(ne+iold)=krange(ii)*krange(ii)/2.0d0+edge+coni*xloss
            iecrit(ii)=ne+iold
            ii=ii+1
         end do


        
!
!
!         
         ne = ne+n1
         ne1 = ne  
      endif

      if (ispec.le.3)  then
!        make the vertical grid in energy plane
!        first point is at 0.005 ev, second at 0.01 ev and
!        exponential grid with step 0.4 after that up to 50 eV
         tempk = 0.005/hart
         em(ne1+1) = edge + coni*tempk
         tempk = tempk*2
         em(ne1+2) = edge + coni*tempk
!        chose delk that point edge+coni*xloss is in the middle of step
!        delk = 0.6 is ok for Cu K edge, but needs more testing
         delk = 0.4
         n1 = nint ( log(xloss/tempk)/delk - 0.5)
         if (n1.le.0) n1 = 1
         bb = exp(delk)
         aa = 2*xloss /(1+bb)
         aa = aa/bb**n1
         if (aa.le. tempk) aa = aa*bb
         if (aa.le.tempk .or. aa.ge. xloss) stop ' Bad mesh in phmesh'
!        delk = log (xloss/tempk) /(n1+0.5)
!        n1 = nint( log(1000/hart/tempk) / delk )
         n1 = nint( log(50/hart/aa) / delk )
         do i = 0, n1
            em(ne1+3+i) = edge +coni*aa*exp(delk*i)
         end do
         ne = ne1 + n1 + 3

!        for DANES need additional points
      endif



!
!     for constant delta eenergy mesh one should set 
!     up iecrit earlier 
!

!      IF (iprint .ge. 1)  THEN
!KJ 7-09 Aleksi used to have different format for emesh.dat.  However I don't see any use for it so I've resorted to the
!        format already used in phmesh2.  In particular, I've removed iecrit from the output.
         OPEN (unit=44, file='emesh.dat', status='unknown')
         WRITE(44,'(a,3(1x,f12.5))') '# edge, bohr, edge*hart ', edge, bohr, edge*hart
         WRITE(44,'(a,2(1x,i5))') '# ispec, ik0 ', ispec, ik0
         WRITE(44,*) '# ie, em(ie)*hart, xk(ie)'
         DO ie = 1, ne
           WRITE (44,'(i5, 3f20.5)') ie, dble(em(ie))*hart, getxk(dble(em(ie))-edge)/bohr
         END DO
         CLOSE(unit=44)
!      END IF



      return
      end
