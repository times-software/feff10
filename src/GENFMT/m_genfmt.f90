!module genfmt_workspace

! Replaces the old COMMON blocks, which don't allow proper dynamic allocation\
! Kevin Jorissen October 2014 for FEFF10

      module clmz
      use dimsmod, only: ltot,mtot,ntot,legtot
      implicit none
      complex*16 clmi(ltot+1,mtot+ntot+1,legtot)
      end module clmz

      module fmatrx
      use dimsmod, only: lamtot,legtot
      complex*16  fmati(lamtot,lamtot,legtot)
      end module fmatrx


      module lambda
      use dimsmod, only: lamtot
      implicit none
      integer mlam(lamtot), nlam(lamtot), lamx, laml0x, mmaxp1, nmax
!        mlam(lamtot), 	!mu for each lambda
!        nlam(lamtot),	!nu for each lambda
!        lamx, 		!max lambda in problem
!        laml0x, 	!max lambda for vectors involving absorbing atom
!        mmaxp1, nmax 	!max mu in problem + 1, max nu in problem
      end module lambda


      module nlm
      use dimsmod, only: ltot,mtot
      implicit none
      real*8 xnlm(ltot+1,mtot+1)
      end module nlm

!     Note that leg nleg is the leg ending at the central atom, so that
!     ipot(nleg) is central atom potential, rat(nleg) position of 
!     central atom.
!     Central atom has ipot=0
!     For later convience, rat(,0) and ipot(0) refer to the central
!     atom, and are the same as rat(,nleg), ipot(nleg).

      module str
      use dimsmod, only: nphx=>nphu
      implicit none
      character*80 text (5)
      character*6, allocatable ::  potlbl(:)
!     text from paths.dat, potential labels
      contains
         subroutine init_str
            implicit none
            allocate(potlbl(0:nphx))
         end subroutine init_str
      end module str

      module pdata
!     ph(nex,-ltot:ltot,0:nphx), !complex phase shifts ipot=0
!     eref(nex),		!complex energy reference
!     rat(3,0:legtot+1),	!position of each atom, code units(bohr)
!     em(nex),		!energy mesh
!     ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
!     deg, rnrmav, xmu, edge,	!(output only)
!     lmax(nex,0:nphx),	!max l with non-zero phase for each energy
!     ipot(0:legtot),	!potential for each atom in path
!     iz(0:nphx),	!atomic number (output only)
!     ltext (5),	!length of each string
!     nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
!     npot, ne,	!number of potentials, energy points
!     ik0,		!index of energy grid corresponding to k=0 (edge)
!     ipath, ihole, 	!index of current path  and hole (output only)
!     kinit, linit, ilinit,  ! initial state kappa and ang. mom.
!     lmaxp1,	!largest lmax in problem + 1
!     ntext 	!number of text  lines
      use dimsmod, only: nex, ltot, nphx=>nphu, nex, legtot
      implicit none
      complex*16, allocatable :: ph(:,:,:)
      complex*16 eref(nex), em(nex)
      integer nsc,nleg,npot,ne,ik0,ipath,ihole,linit,ilinit, &
        ipot(0:legtot), ltext (5), ntext, lmaxp1, kinit
      integer, allocatable :: lmax(:,:), iz(:)
      real*8  rat(3,0:legtot+1),                                               &
      ri(legtot), beta(legtot+1), eta(0:legtot+1),                     &
      deg, rnrmav, xmu, edge
      contains
         subroutine init_pdata
            implicit none
            allocate(ph(nex,-ltot:ltot,0:nphx),lmax(nex,0:nphx),iz(0:nphx))
         end subroutine init_pdata
      end module pdata

      module rotmat
      use dimsmod, only: ltot, mtot, legtot
      implicit none
      real*8 dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)
      end module rotmat

!end module genfmt_workspace
