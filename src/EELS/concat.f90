!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: concat.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine concat (str1,str2,strout,length)
      !places the concatenation of str1 (less trailing blanks) and str2
      !in strout and returns the length of the non-blank portion of strout
      !in length.
      implicit none
      character*(*) str1, str2, strout
      integer length, istrln
      external istrln
      
      length=istrln(str1) 
      strout=str1(1:length)//str2
      length=istrln(strout)
      return
      end

  





