!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_errormod.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE ErrorMod
        IMPLICIT NONE
        INTEGER :: MxErrLines
        PARAMETER( MxErrLines = 5 )
        LOGICAL,PRIVATE :: FirstErrorCall = .TRUE.

        TYPE ErrorType
           ! ErrorName - name of errror
           CHARACTER(30)  ErrorName
           CHARACTER(150) Message(MxErrLines)
           CHARACTER(10)  Action
           LOGICAL        ErrorOccured
           INTEGER        NErrLines
        END TYPE ErrorType

        CHARACTER(300),PRIVATE :: ErrorMessage
        LOGICAL,PRIVATE :: MessageIsSet = .FALSE.

        INTERFACE Error
           MODULE PROCEDURE ErrorSM
           MODULE PROCEDURE ErrorMM
        END INTERFACE
      CONTAINS
        
        ! Set the message to report when an error occurs.
        SUBROUTINE SetErrorMessage(Message)
          CHARACTER*(*) Message
          
          ErrorMessage = Message
          MessageIsSet = .TRUE.
          
          RETURN 
        END SUBROUTINE SetErrorMessage
        
        
        ! Returns true if error message is set, false otherwise.
        FUNCTION ErrorMessageIsSet()
          LOGICAL ErrorMessageIsSet
          ErrorMessageIsSet = MessageIsSet
          RETURN
        END FUNCTION ErrorMessageIsSet


        ! Print ErrorMessage if Message is not supplied,
        ! otherwise, Message is printed. Then STOP.
        SUBROUTINE ErrorSM(Message,StopProgram)
          CHARACTER*(*),INTENT(IN),OPTIONAL :: Message
          LOGICAL,INTENT(IN),OPTIONAL :: StopProgram
          ! If this is the first call to Error, print buffer.
          IF(FirstErrorCall) THEN
             PRINT*
             PRINT '(A)', '************************************************************************************'
             PRINT '(A)', '************************************************************************************'
             FirstErrorCall = .FALSE.
          END IF

          IF(PRESENT(Message)) THEN         
             PRINT*, TRIM(Message)
          ELSE
             IF(MessageIsSet) THEN
                PRINT*, TRIM(ErrorMessage)
             ELSE
                PRINT*, 'Unknown Error'
             END IF
          END IF

          IF(PRESENT(StopProgram)) THEN
             IF(StopProgram) THEN
                PRINT '(A)', '************************************************************************************'
                PRINT '(A)', '************************************************************************************'
                PRINT*
             END IF
          ELSE
             PRINT '(A)', '************************************************************************************'
             PRINT '(A)', '************************************************************************************'
             PRINT*
             STOP
          END IF

        END SUBROUTINE ErrorSM

        SUBROUTINE ErrorMM(Messg,StopProgram)
          CHARACTER*(*),INTENT(IN) :: Messg(:)
          INTEGER iMessg, NMessg
          LOGICAL,OPTIONAL,INTENT(IN) :: StopProgram
          
          NMessg = SIZE(Messg)
          DO iMessg = 1, NMessg - 1
             IF(PRESENT(StopProgram)) THEN
                CALL ErrorSM(Messg(iMessg), StopProgram = StopProgram)
             ELSE
                CALL ErrorSM(Messg(iMessg), StopProgram = .FALSE.)
             END IF
          END DO

          CALL ErrorSM(Messg(NMessg), StopProgram = .TRUE.)

        END SUBROUTINE ErrorMM
        
        SUBROUTINE CheckError(ErrType)
          TYPE(ErrorType),INTENT(IN) :: ErrType
          LOGICAL StopOnError
          CHARACTER(10) Action
          INTEGER i1

          ! If the error did not occur, return
          IF(.not.ErrType%ErrorOccured) RETURN

          ! Set the action. Stop, warn, or do nothing.
          Action = ADJUSTL(ErrType%Action)
          CALL Upper(Action)
          IF(Action(1:4).eq.'STOP') THEN              
             StopOnError = .TRUE.
             CALL Error('ERROR:', StopProgram = .FALSE.)
          ELSEIF(Action(1:4).eq.'WARN') THEN
             StopOnError = .FALSE.
             CALL Error('WARNING:', StopProgram = .FALSE.)
          ELSE
             RETURN
          END IF

          ! Print the errror messages
          DO i1 = 1, ErrType%NErrLines
             CALL Error(ErrType%Message(i1), StopProgram = .FALSE.)
          END DO
          
          ! If Action = 'STOP' Stop.
          IF(StopOnError) STOP
        END SUBROUTINE CheckError


        SUBROUTINE CheckErrors(ErrorTypes)
          TYPE(ErrorType) ErrorTypes(:)
          INTEGER i1
          
          DO i1 = 1, SIZE(ErrorTypes)
             CALL CheckError(ErrorTypes(i1))
          END DO
        END SUBROUTINE CheckErrors


        ! Check Allocation.
        SUBROUTINE CheckAllocation(iErr,Message)
          INTEGER,INTENT(IN) :: iErr
          CHARACTER*(*),INTENT(IN) :: Message

          IF(iErr.ne.0) CALL Error(Message)

          RETURN
        END SUBROUTINE CheckAllocation

      END MODULE ErrorMod
