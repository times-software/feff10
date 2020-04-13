!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_work.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/07/02 05:32:11 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        module work
!     some parameters that are needed in several subroutines
      real*8 thpart
        integer npos

        CONTAINS
          subroutine init_work(acoll,aconv,nqr,nqf,qmodus)
           implicit none
            real*8 acoll,aconv
            integer nqr,nqf
            character*1 qmodus
      IF ((acoll.GT.1.0D-6).OR.(aconv.GT.1.0D-6)) THEN
         ThPart = (acoll + aconv) / DBLE(2*nqr)
      ELSEIF((nqr+nqf).gt.2) then
	     write(11,*) ':INFO - You specified (almost?) zero collection and convergence angle.'
		 write(11,*) 'Your NQR and NQF values were reset to 1.'
         ThPart = dble(0)
         nqr = 1
         nqf = 1
      ENDIF
      if(qmodus.eq.'U'.or.qmodus.eq.'L') then
            npos=nqr*nqr*nqf
	  elseif(qmodus.eq.'1') then
            nqf=1   ! it is not used anyway
			npos=nqr
	  endif

      return
	  end subroutine init_work

	  end module work
