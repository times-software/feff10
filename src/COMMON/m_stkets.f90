!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_stkets.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module stkets
  !--------------------------------------------------------------------
  ! Module to contain state kets
  !
  ! Variables:
  !   lrstat:  Array of state kets at current energy
  !   istate:  ??? (just an index)
  ! Subroutines:
  !   init_stkets: Allocates memory for lrstat
  !   kill_stkets: Deallocates memory for lrstat
  !        getkts: Constructs state kets |iat,l,m> at this energy
  !--------------------------------------------------------------------
  implicit none
  integer ,allocatable :: lrstat(:,:)
  integer :: istate  

contains

    subroutine init_stkets(istatx)
    integer,intent(in) :: istatx
    allocate(lrstat(4,istatx))
    end subroutine init_stkets

    subroutine kill_stkets
    deallocate(lrstat)
    end subroutine kill_stkets

    subroutine getkts(nsp, nat, npot, iphx, lipotx, i0)
    !--------------------------------------------------------------------
    !  Construct state kets |iat,l,m> at this energy
    !--------------------------------------------------------------------
    !  input
    !    nsp:    ???
    !    nat:    number of atoms in cluster
    !    npot:   number of unique potentials
    !    iphx:   (nclusx) potential index of each atom in the cluster
    !    lipotx: (nphasx) maximum angular momentum to consider for each
    !            ipot
    !  output
    !    i0:     index shift for each potential representative
    !
    !  Passed in module
    !    istate: number of states 
    !    lrstat: (4, istatx) state kets |iat,l,m> 
    !--------------------------------------------------------------------
    use DimsMod, only: nclusx, nphasx, nphx=>nphu, lx, istatx

    integer, intent(in) :: nsp,nat,npot
    integer, intent(in) ::iphx(nclusx),lipotx(0:nphasx)
    integer, intent(out) :: i0(0:nphx)
    integer :: iat,l,m,isp,ip,lim


    istate = 0
    i0(:) = -1
    IAT_LOOP: do iat=1,nat
       ip = iphx(iat)
       ! i0(ip) - index for the ip-representative atom
       ! need for simple find of states for ip-representative.
       if (i0(ip).lt.0) i0(ip) = istate
       lim = min(lx, lipotx(ip))
       L_LOOP: do l=0,lim
          M_LOOP: do m = -l, l
             ISP_LOOP:  do isp = 1, nsp
                istate = istate + 1
                if (istate.gt.istatx) then
                   call wlog('Exceeded maximum number of LR states.  Stopping')
                   call par_stop('GETKTS-1')
                endif
                lrstat(1,istate) = iat
                lrstat(2,istate) = l
                lrstat(3,istate) = m
                lrstat(4,istate) = isp
             enddo ISP_LOOP
          enddo M_LOOP
       enddo L_LOOP
    enddo IAT_LOOP

    return
  end subroutine getkts

end module stkets
