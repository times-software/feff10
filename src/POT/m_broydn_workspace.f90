module broydn_workspace
!KJ 12-2011
!KJ I really just wanted to use nscmt as dimensioning parameter.
!KJ But that meant the SAVE statement no longer worked.
!KJ So I just stuffed it in a module.

!     work space for Broyden algorithm
!      parameter (nbr=30)  !KJ 12-2011 A parameter by the same name exists in rdinp.f90 .  I'm not sure if they need to be identical.  I rather think this nbr may have to be
	  ! identical to the maximum value of nmix (whereas in rdinp.f90 nbr is used to limit nscmt).  Yucky programming.
      !parameter (nbr=100)  !KJ changed for debugging
	  !KJ I think a better solution is this:
	  !KJ 12-2011 use nscmt (from feff.inp) as dimensioning parameter; see use statement below
	  	  
      real*8,allocatable :: cmi(:,:), frho(:,:,:), urho(:,:,:)
      real*8,allocatable :: xnorm(:), wt(:), rhoold(:,:), ri05(:)
!KJ here's the sucker:      save cmi, frho, urho, xnorm, wt, rhoold, ri05

      contains
	    subroutine broydn_workspace_init
	       use potential_inp,only:nbr=>nscmt
	       use DimsMod,only: nphx=>nphu
		   allocate(cmi(nbr,nbr), frho(251,0:nphx,nbr), urho(251,0:nphx,nbr), xnorm(nbr), wt(251), rhoold(251,0:nphx), ri05(251))
		end subroutine broydn_workspace_init

end module broydn_workspace
