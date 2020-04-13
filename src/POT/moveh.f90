!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: moveh.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine moveh (nat, iphat, iz, rath)
!    Increase length of bonds with hydrogen atoms
!    Move Hydrogens for potentials. Otherwise MT geometry is screwed up.
        use DimsMod, only: natx, nphx=>nphu

      implicit double precision (a-h, o-z)

!     input is everything; output is modified  atomic coordinates (rath) 
!     nat is number of atoms in cluster
      dimension iphat(natx),  iz(0:nphx)
      dimension rath(3,natx)

      do 970 iat = 1, nat
        if (iz(iphat(iat)) .eq. 1) then
!         find the nearest atom A, units for rat are bohr.
          rah = 100
          ia = 0
          do i = 1,nat
              rattmp = dist(rath(1,iat), rath(1,i) )
              if (rattmp.lt. rah .and. i.ne. iat) then
                 ia = i
                 rah = rattmp
              endif
          enddo
          if (iz(iphat(ia)).eq.1) goto 970

!         set max distance as function of rah ( set by calculations
!         for H2O and GeH4)
          ratmax = rah + 4.d0/rah**2 

!        find shortest AB bond (neither A or B are H)
          rab = 10
          ib = 0
          do i = 1,nat
              rattmp = dist(rath(1,ia), rath(1,i))
              if (i.ne.ia .and. iz(iphat(i)).ne.1 .and.                 &
     &            rab.gt.rattmp) then
                 rab = rattmp
                 ib = i
              endif
          enddo
          if (rab.lt.ratmax) ratmax = 0.95*rab + 0.05*rah
          if (rah .gt. ratmax) goto 970

!         increase rah to ratmax and check that A is still closest to H
          ratmin = rah
  960     do i = 1,3
           rath(i,iat)=rath(i,ia)+ratmax/ratmin*(rath(i,iat)-rath(i,ia))
          enddo
          rbh = 10
          ib = 0
          do i = 1,nat
              rattmp = dist(rath(1,iat), rath(1,i))
              if (i.ne.iat .and. rbh.gt.rattmp) then
                 rbh = rattmp
                 ib = i
              endif
          enddo

          if (ia.ne.ib) then
             rab = dist(rath(1,ia),rath(1,ib))
             rattmp = ratmax*rab**2/(ratmax**2+rab**2-rbh**2)
             ratmin = ratmax
             ratmax = 0.95*rattmp +0.05*rah
             goto 960
          endif
        endif
  970 continue

      return
      end
