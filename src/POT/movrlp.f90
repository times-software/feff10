!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: movrlp.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine movrlp ( nph, nat, iphat, rat, iatph, xnatph,          &
                     novr, iphovr, nnovr, rovr,                        &
                     imt, rmt, rnrm, ri, lnear,                        &
                     cmovp, ipiv, volint, inters)

!     Constructs overlap matrix based on geometry of overlapped
!     muffin-tin spheres. Uses LU decomposition for inversion of matrix

      use constants, only: pi
      use DimsMod, only: natx, nphx=>nphu, novrx, novp
      implicit none

      integer, intent(in) :: iphat(natx)
      real*8, intent(in) :: rat(3,natx)
      integer, intent(in) :: iatph(0:nphx)
      real*8, intent(in) :: xnatph(0:nphx)
      integer, intent(in) :: novr(0:nphx)
      integer, intent(in) :: iphovr(novrx,0:nphx)
      integer, intent(in) :: nnovr(novrx,0:nphx)
      real*8, intent(in) :: rovr(novrx,0:nphx)
      integer, intent(in) :: imt(0:nphx)
      real*8, intent(in) :: rmt(0:nphx)
      real*8, intent(in) :: rnrm(0:nphx)
      logical, intent(in) :: lnear(0:nphx)
      integer, intent(in) :: nat, nph, inters
      real*8, intent(out) :: ri(251)
      complex, intent(out) :: cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      real*8, intent(out) :: volint
      integer, intent(out) :: ipiv(novp*(nphx+1)+1)
!     local
      character*512 slog
!     work space for linear algebra
!      parameter (novp=40)
!      parameter (novp=100)
      real bmat(nphx+1,novp*(nphx+1))
      integer, external :: ii
      integer iat, ip1, i1, i2, ipot, irav, i,ncp, nlast, ntmp, iat0, ip2, &
        ix2, ind2, ix1, iph, istx, imin1, imin2, info
      real*8, external :: dist, calcvl
      real*8 xn, exphx, rnn, ri1, ri2, temp, xr, r1, r2, rav, aa

!     get ipot and irav from inters
      ipot = mod(inters,2)
      irav = (inters-ipot) / 2
      do i=1,251
         ri(i)=exp(-8.8d0+(i-1)*0.05d0)
      enddo
      exphx=exp(0.025d0)
!     initially cmovp is a unit matrix up to ncp
      ncp = novp*(nph+1)+1
	  do i1=1,ncp
	     do i2=1,ncp-1
	        cmovp(i1,i2)=0.d0
         enddo
	     cmovp(i1,i1)=1.d0
		 cmovp(i1,ncp)=0.01d0
	  enddo
      xn = 0.d0
      bmat(:,:) = 0.d0

      do 200 ip1=0,nph
        if (novr(ip1) .gt. 0 ) then
           nlast = novr(ip1)
        else
           iat0 = iatph(ip1)
           ntmp = 1
           nlast = nat
        endif
        if (irav.eq.1) then
          rav = (rmt(ip1) + rnrm(ip1)) / 2
        elseif (irav.eq.0) then
          rav =  rnrm(ip1)
        else
          rav=ri(imt(ip1)+1)
        endif
        if (lnear(ip1)) rav=ri(imt(ip1)+1)

        do 190 iat = 1,nlast
          if (novr(ip1) .gt. 0 ) then
             ntmp = nnovr(iat,ip1)
             ip2 = iphovr(iat,ip1)
             rnn = rovr(iat,ip1)
          else
            if (iat.eq.iat0) goto 190
            ip2 = iphat(iat)
            rnn = dist (rat(1,iat0), rat(1,iat))
          endif

!         correct for double counting volume and area
          if (rnn .lt. rmt(ip1)+rmt(ip2)) then
!            correct interstitial volume
             volint = volint + xnatph(ip1) * ntmp * (calcvl( rmt(ip1), rmt(ip2), rnn) + calcvl(rmt(ip1), rmt(ip2), rnn)) / 2.d0
          endif

!         using expression for vtot(jri) ,(jri=i1)
!         first fill  matrix bmat
          ix1 = ip1+1

          if (rav+rmt(ip2) .le. rnn) goto 100
          imin2 = ii( rnn-rav )
		  
          if (imt(ip2)-imin2 .ge. novp-1) then
             write(slog,132) ip1
  132        format(' FOLP for POTENTIAL type ',i3,' is too big.')
             call wlog (slog)
             write(slog,'(a)') ' Reduce overlap using FOLP and rerun'
             call wlog (slog)
			 write(*,*) 'ip1,ip2,imt(ip2),imin2,novp,rnn,rav,ii(rnn-rav)',ip1,ip2,imt(ip2),imin2,novp,rnn,rav,ii(rnn-rav)
             call par_stop('MOVRLP-1')
          endif
          imin2=imt(ip2)-novp+1

          do 80 i2 = imin2,imt(ip2)
             r1=ri(i2)/exphx
             r2=ri(i2)*exphx
             if (i2.eq.imt(ip2)) r2=rmt(ip2)
             if (i2.eq.imt(ip2))   r1=(r1+2*ri(imt(ip2))-rmt(ip2))/2.d0
             if (i2.eq.imt(ip2)-1) r2=(r2+2*ri(imt(ip2))-rmt(ip2))/2.d0
             if (r2+rav .lt. rnn) goto 80
             if (r1+rav .lt. rnn) then
!               use linear interpolation between cases xr=0, xr=1
                xr = (rnn-rav-r1)/ (r2-r1)
                r1 = rnn-rav   
                temp =  (r2**2 - r1**2) / (4*rnn*rav) * ntmp
                ind2=i2+1
                if (i2.eq.imt(ip2))  ind2=i2-1
                xr = xr * (r2-ri(i2)) / (ri(ind2)-ri(i2))
                ix2 = ip2*novp + i2 - imin2 + 1
                bmat (ix1,ix2) = bmat (ix1,ix2) + real(temp*(1-xr))
                ix2 = ip2*novp + ind2 - imin2 + 1
                bmat (ix1,ix2)=bmat (ix1,ix2) + real(temp*xr)
             else
                temp = (r2**2 - r1**2) / (4*rnn*rav   ) * ntmp
                ix2 = ip2*novp + i2 - imin2 + 1
                bmat (ix1,ix2) = bmat (ix1,ix2) + real( temp)
             endif
 80         continue !enddo !i2

!         using expression for vtot(i) ,(i<jri)
!         construct matrix  cmovp
 100      if (rmt(ip1)+rmt(ip2) .le. rnn) goto 190

          imin1=ii(rnn-rmt(ip2))
          imin2=ii(rnn-rmt(ip1))
          if (imt(ip1)-imin1.ge.novp-1 .or. imt(ip2)-imin2.ge.novp-1) then
                call wlog('Overlapping problem in potentials/movrlp.')
                call wlog('This sometimes happens in anisotropic systems where default overlapping of muffin-tin radii fails.')
                call wlog('You can try to use the FOLP card to solve this problem.  In particular, if you have H atoms, try using')
                call wlog('FOLP n 0.8')
                call wlog('where "n" is the potential index of H in ' // &
                          'feff.inp or pot.inp. Rerun FEFF. If you ' // &
                          'cannot solve the problem, ask the authors for help.')
                call wlog('Below is debugging output that you may send to the authors:')
                write(*,*) "overlap problem for ip1,ip2=",ip1,",",ip2,"."
                write(*,*) "rnn=",rnn,"rmt(ip1),imt(ip1)=",rmt(ip1),imt(ip1),"rmt(ip2),imt(ip2)=",rmt(ip2),imt(ip2),"."
                write(*,*) "novp=",novp,", crit1=",imt(ip1)-imin1,", crit2=",imt(ip2)-imin2,"."
                call par_stop('tell authors to INCREASE NOVP')
          endif
          imin1=imt(ip1)-novp+1
          imin2=imt(ip2)-novp+1

          do 180 i1 = imin1,imt(ip1)
            ri1=ri(i1)/exphx
            ri2=ri(i1)*exphx
            if (i1.eq.imt(ip1)) ri2=rmt(ip1)
            if (i1.eq.imt(ip1)) ri1=(ri1+2*ri(imt(ip1))-rmt(ip1))/2.d0
            if (i1.eq.imt(ip1)-1) ri2=(ri2+2*ri(imt(ip1))-rmt(ip1))/2.d0
            ix1 = i1-imin1+1  + ip1*novp
            do 170 i2 = imin2,imt(ip2)
              r1=ri(i2)/exphx
              r2=ri(i2)*exphx
              if (i2.eq.imt(ip2)) r2=rmt(ip2)
              if (i2.eq.imt(ip2))   r1=(r1+2*ri(imt(ip2))-rmt(ip2))/2.d0
              if (i2.eq.imt(ip2)-1) r2=(r2+2*ri(imt(ip2))-rmt(ip2))/2.d0
              if (r2+ri2.lt.rnn) goto 170

!             calculate volume of intersection
              temp = calcvl(ri2,r2,rnn) + calcvl(r2,ri2,rnn)
              if (ri1+r2.gt.rnn) temp = temp - calcvl(ri1,r2,rnn) - calcvl(r2,ri1,rnn)
              if (ri2+r1.gt.rnn) temp = temp - calcvl(ri2,r1,rnn) - calcvl(r1,ri2,rnn)
              if (ri1+r1.gt.rnn) temp = temp + calcvl(ri1,r1,rnn) + calcvl(r1,ri1,rnn)
!             volume of intersection (temp) should be divided by volume between spheres ri1 and ri2
              temp=temp / ( 4.d0/3.d0*pi * (ri2**3-ri1**3) ) * ntmp

              if (r1+ri2.lt.rnn) then
!               use linear interpolation between cases xr=0, xr=1
                xr = (rnn-ri(i1)-r1)/ (r2-r1)

                ind2=i2+1
                if (i2.eq.imt(ip2))  ind2=i2-1
                xr = xr * (r2-ri(i2)) / (ri(ind2)-ri(i2))
                ix2 = i2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2) +cmplx (temp*(1-xr))
                ix2 = ind2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2) +cmplx (temp*xr)
                r1=rnn-ri2
              else
                ix1 = i1-imin1+1 + ip1*novp
                ix2 = i2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2)  +cmplx (temp)
              endif
 170        continue
 180      continue

 190     continue
         xn = xn + xnatph(ip1)
  200 continue

!     using matrix bmat fill in the last row of matrix cmvovp
!     this is additional equation to find Vint.
!     switch to local equation from average over all atoms
      if (ipot .eq. 0) then
         do 260 iph=0, nph
!          xn may differ from nat, if atom list have more natx atoms
!          see rdinp.f
           aa = xnatph(iph)/xn
           do 250 ix1 = 1, ncp-1
  250      cmovp(ncp,ix1) = cmovp(ncp, ix1) + aa*bmat(iph+1,ix1)
  260    continue
      else  
         iph=0
         do 270 ix1 = 1, ncp-1
  270    cmovp(ncp,ix1) = cmovp(ncp, ix1) + bmat(iph+1,ix1)
      endif

! --- invert matrices by LU decomposition
!     call cgetrf from lapack.  this performs an LU decomposition on the matrix 
      istx=novp*(nphx+1) + 1
      call cgetrf( ncp, ncp, cmovp, istx, ipiv, info )
      if (info.ne.0) call wlog('    *** Error in cgetrf when computing cmovp')

!     have to check that the last was not permuted, otherwise the density calculation will be wrong
!     this is also why we put 0.01 in last column and not 1.0
      if (ipiv(ncp).ne.ncp) call par_stop('illegal permutation in ipiv ')

      return
      end
