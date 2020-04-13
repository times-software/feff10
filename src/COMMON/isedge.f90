
      logical function isedge (str)
!     checks str to see if it will be recognized as an edge
!     by rdinp.
      implicit none
      integer i
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

      isedge=.false.
      do  i = 0,29
        if (str .eq. edglbl(i) .or. str .eq. edglbp(i) ) isedge=.true.
      enddo
   
      return
      end

