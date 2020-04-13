
      subroutine stdnm (str)
!     Converts edge number to label if the edge was specified by number.
      implicit none
      integer i
      logical,external :: isedge
      character*2 str, edglbl, edglbp
      dimension edglbl(0:29), edglbp(0:29)

      data edglbl / 'NO', 'K ', 'L1', 'L2', 'L3',                       &
     &            'M1','M2','M3','M4','M5',                             &
     &            'N1','N2','N3','N4','N5','N6','N7',                   &
     &            'O1','O2','O3','O4','O5','O6','O7',                   &
     &            'P1','P2','P3','P4','P5','R1' /
      data edglbp / '0', '1 ', '2', '3', '4',                           &
     &            '5','6','7','8','9',                                  &
     &            '10','11','12','13','14','15','16',                   &
     &            '17','18','19','20','21','22','23',                   &
     &            '24','25','26','27','28','29' /

      if (isedge(str)) then
        do  i = 0,29
          if (str .eq. edglbp(i) ) str=edglbl(i) 
        enddo
      else
        str='  '
      endif
   
      return
      end
