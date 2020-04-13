!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: setedg.f90,v $:
! $Revision: 1.2 $
! $Author: jorissen $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setedg (a2, ihole)
!KJ Dec 2013: this routine converts between the "name" and "number" that identifies the edge we work with; e.g. "K" == "1".
!KJ I made the routine a 2-way affair: now, if a2 doesn't seem to contain a valid edge name to convert into a number, it will attempt to translate from number to name instead.
!KJ I remain suspicious of the "NO" option here.  I don't think that's how we normally execute a NOHOLE calculation ...
!KJ Will leave it for now.

      implicit none
      character*2,intent(inout) :: a2
      integer,intent(inout) :: ihole
      integer i
      character*2 edgelabel, edgenumber
      character*1 uppercase, lowercase
      dimension edgelabel(0:29), edgenumber(0:29)
      dimension uppercase(0:6),lowercase(0:6)
      logical label_to_number

      data edgelabel / 'NO', 'K ', 'L1', 'L2', 'L3',                       &
     &            'M1','M2','M3','M4','M5',                             &
     &            'N1','N2','N3','N4','N5','N6','N7',                   &
     &            'O1','O2','O3','O4','O5','O6','O7',                   &
     &            'P1','P2','P3','P4','P5','R1' /
      data edgenumber / '0', '1 ', '2', '3', '4',                           &
     &            '5','6','7','8','9',                                  &
     &            '10','11','12','13','14','15','16',                   &
     &            '17','18','19','20','21','22','23',                   &
     &            '24','25','26','27','28','29' /
     data uppercase /'K','L','M','N','O','P','R'/
     data lowercase /'k','l','m','n','o','p','r'/

      !Make uppercase
      do i=0,6
         if(a2(1:1).eq.lowercase(i)) a2(1:1)=uppercase(i)
      enddo
      if(a2(2:2).eq.'o') a2(2:2)='O'

      label_to_number=.false.
      do i = 0,29
         if (a2 .eq. edgelabel(i) .or. a2 .eq. edgenumber(i) ) then
            label_to_number=.true.
            ihole  = i
         endif
      enddo

      if(.not.label_to_number) then
         if(ihole.ge.0 .and. ihole.le.29) then
            a2=edgelabel(ihole)
         else
            call par_stop('Unknown EDGE in RDINP-setedg: check feff.inp for EDGE or HOLE errors.')
         endif
      endif


      !if (ihole  .lt. 0) call par_stop('unknown EDGE')

      return
      end
