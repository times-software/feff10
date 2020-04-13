      subroutine clearfp
!
!     This is a null version of a
!     subroutine to clear IEEE floating point exceptions
!     for inexact and underflow under SUN OS 4 f77
!     For most other systems, no action is needed.
!
!     character*1 out
!     ii = ieee_flags('clear','exception','underflow',out)
!     ii = ieee_flags('clear','execption','inexact',out)
      return
      end
