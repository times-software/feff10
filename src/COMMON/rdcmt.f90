!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdcmt.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE rdcmt(iUnt,Cmt,ierr)
      INTEGER iUnt, i1, ierr
      CHARACTER(300) line
      CHARACTER(4) Cmt
      CHARACTER TmpCmt(4), ch
      LOGICAL CmtLin
      
      ierr = 0
      CmtLin = .true.
      DO i1 = 1, 4
         TmpCmt(i1) = Cmt(i1:i1)
      END DO
 5    CONTINUE
      READ(iUnt,*,END=10) ch
      DO i1 = 1, 4
         IF(ch.eq.TmpCmt(i1)) goto 5
      END DO
      
      BACKSPACE(iUnt)
      RETURN
 
 10   CONTINUE
      ierr = -1
      RETURN
      END
