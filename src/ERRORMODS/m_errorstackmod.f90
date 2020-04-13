!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_errorstackmod.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE ErrorStackMod
      USE ErrorMod
      USE StringLLMod

      IMPLICIT NONE
      
      TYPE(StringLL),PRIVATE :: ErrorStack
      LOGICAL,PRIVATE :: PrintErrorStack
      INTEGER,PRIVATE :: NumSpaces

      CONTAINS

!     Messages will print as the stack is being filled/depleted
      SUBROUTINE InitErrorStack(PrintingIsOn)
      LOGICAL,INTENT(IN),OPTIONAL :: PrintingIsOn

!     Initialize string linked list
      CALL InitStringLL(ErrorStack)
      NumSpaces = 0
!     Set printing status
      IF(PRESENT(PrintingIsOn)) THEN
         PrintErrorStack = PrintingIsOn
      ELSE
         PrintErrorStack = .FALSE.
      END IF

      END SUBROUTINE InitErrorStack     

 
!     Fills the next error in the stack.
!     Prints the message if PrintErrorStack is true.
      SUBROUTINE NextError(Str)
      CHARACTER*(*) Str
      TYPE(StringNode),POINTER :: NextErrorNode

      ! Make the next node with Str as the message.
      CALL MakeStringNode(Str,NextErrorNode)

      ! Append NextErrorNode to the end of the stack.
      CALL AppendStringNode(ErrorStack,NextErrorNode)
      
      IF(PrintErrorStack) THEN
         
         PRINT '(A)', REPEAT(' ',NumSpaces) // '-->' // &
              & TRIM(ErrorStack%Tail%String)
         PRINT '(A)', REPEAT(' ',NumSpaces+2) // '|'
         NumSpaces = NumSpaces + 2
      END IF

      END SUBROUTINE NextError


      
!     Deletes the last error in the stack. Prints last error if PrintErrorStack 
!     is true.
      SUBROUTINE DeleteError()
      IF(PrintErrorStack) THEN
         PRINT '(A)', REPEAT(' ',NumSpaces-2) // '<--' // &
              & TRIM(ErrorStack%Tail%String)
         PRINT '(A)', REPEAT(' ',NumSpaces-4) // '|'
         NumSpaces = NumSpaces - 2
      END IF

      CALL DeleteStringNode( ErrorStack,StringLLLength(ErrorStack) )

      END SUBROUTINE DeleteError


!     Dump error stack to file error.stack
      SUBROUTINE DumpErrorStack()
      USE IOMod
      INTEGER :: iErrorUnit
      INTEGER iLen, i1, iSpace
      TYPE(StringNode),POINTER :: p
      LOGICAL FileIsOpen
      iSpace = 0

      ! Open file error.stack for writing if is not open.
      CALL OpenFl('error.stack', iErrorUnit, 'replace',FileIsOpen)
      ! Write errorstack to file, and to screen.
      iLen = StringLLLength(ErrorStack) ! iLen is number of messages
      p => ErrorStack%Head              ! Initialize pointer to first elements
      DO i1 = 1, iLen
         ! Print a formatting string '--> ' and the error message
         PRINT '(A)', REPEAT(' ',iSpace) // '--> ' // TRIM(p%String)
         WRITE(iErrorUnit,'(A)') REPEAT(' ',iSpace) // '--> ' // &
              & TRIM(p%String)
         
         ! increment the number of spaces
         iSpace = iSpace + 4 + LEN_TRIM(p%String)-2
         
         ! Print another formatting string (spaces followed by |)
         PRINT '(A)', REPEAT(' ',iSpace) // '|'
         WRITE(iErrorUnit,'(A)') REPEAT(' ',iSpace) // '|'

         ! Update pointer
         p => p%Next
      END DO
      
      END SUBROUTINE DumpErrorStack
         

      END MODULE ErrorStackMod
