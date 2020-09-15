
      logical function isedge (str)
!     checks str to see if it will be recognized as an edge
!     by rdinp.
      implicit none
! Josh Kas - Changed array dimensions from 30 to 40 (and others) for high Z elements
! according to Pavlo Baranov's changes.
      integer i
      character*2 str, edglbl, edglbp
      dimension edglbl(0:40), edglbp(0:40)

      data edglbl / 'NO', 'K ', 'L1', 'L2', 'L3',                       &
     &            'M1','M2','M3','M4','M5',                             &
     &            'N1','N2','N3','N4','N5','N6','N7',                   &
     &            'O1','O2','O3','O4','O5','O6','O7','O8','O9',         &
     &            'P1','P2','P3','P4','P5','P6','P7',                   &
     &            'R1','R2','R3','R4','R5',                             &
     &            'S1','S2','S3' /
      data edglbp / '0', '1 ', '2', '3', '4',                           &
     &            '5','6','7','8','9',                                  &
     &            '10','11','12','13','14','15','16',                   &
     &            '17','18','19','20','21','22','23',                   &
     &            '24','25','26','27','28','29','30',                   &
     &            '31','32','33','34','35','36','37',                   &
     &            '38','39','40' /
      isedge=.false.
      do  i = 0,40
        if (str .eq. edglbl(i) .or. str .eq. edglbp(i) ) isedge=.true.
      enddo
   
      return
      end

