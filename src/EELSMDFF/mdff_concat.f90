!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: mdff_concat.f90,v $:
! $Revision: 1.2 $
! $Author: jorissen $
! $Date: 2011/02/11 04:14:53 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mdff_concat (str1,str2,strout,length)
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

  





