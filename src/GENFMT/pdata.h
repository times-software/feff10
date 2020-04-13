!     Note that leg nleg is the leg ending at the central atom, so that
!     ipot(nleg) is central atom potential, rat(nleg) position of 
!     central atom.
!     Central atom has ipot=0
!     For later convience, rat(,0) and ipot(0) refer to the central
!     atom, and are the same as rat(,nleg), ipot(nleg).

!     text arrays 
      character*80 text 
      character*6  potlbl
!     text from paths.dat, potential labels
      common /str/ text (5), potlbl(0:nphx)

!     common /pdata/ ph(nex,-ltot:ltot,0:nphx), !complex phase shifts ipot=0
!    1 eref(nex),		!complex energy reference
!    1 rat(3,0:legtot+1),	!position of each atom, code units(bohr)
!    1 em(nex),		!energy mesh
!    1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
!    1 deg, rnrmav, xmu, edge,	!(output only)
!    1 lmax(nex,0:nphx),	!max l with non-zero phase for each energy
!    1 ipot(0:legtot),	!potential for each atom in path
!    1 iz(0:nphx),	!atomic number (output only)
!    1 ltext (5),	!length of each string
!    1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
!    1 npot, ne,	!number of potentials, energy points
!    1 ik0,		!index of energy grid corresponding to k=0 (edge)
!    1 ipath, ihole, 	!index of current path  and hole (output only)
!    1 kinit, linit, ilinit,  ! initial state kappa and ang. mom.
!    1 lmaxp1,	!largest lmax in problem + 1
!    1 ntext 	!number of text  lines

      complex*16 ph(nex,-ltot:ltot,0:nphx), eref(nex), em(nex)
      integer nsc,nleg,npot,ne,ik0,ipath,ihole,linit,ilinit,lmax(nex,0:nphx), &
        ipot(0:legtot), iz(0:nphx), ltext (5), ntext, lmaxp1
      real*8  rat(3,0:legtot+1),                                               &
     & ri(legtot), beta(legtot+1), eta(0:legtot+1),                     &
     & deg, rnrmav, xmu, edge


      common /pdata/ ph, eref, em,     &
     & rat, ri, beta, eta, deg, rnrmav, xmu, edge, lmax, ipot,        &
     & iz, ltext, nsc, nleg, npot, ne, ik0, ipath, ihole, kinit, linit, ilinit, lmaxp1, ntext 
