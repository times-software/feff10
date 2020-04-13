!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: vbh.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine vbh(rs,xmag,vxc)

      implicit double precision (a-h, o-z)

!   INPUT: density parameter rs, 2* fraction of given spin orientation.
!   OUTPUT: xc potential for given spin orientation.
!   Reference: Von Barth, Hedin, J.Phys.C, 5, 1629, (1972). eq.6.2
!   xmag is twice larger than 'x' in their paper
!   effect of tau was also found to be small. thus tau is not used

!     parameter (asm = 2.0**(-1.0/3.0) )
!     parameter (gamma = 4.0/3.0*asm/(1-asm) )
! APS parameter (gamma = 5.129762496709890 ) changed
      parameter (gamma = 5.129762802484097 )

      vxc = 0.0
      if (rs.gt.1000) goto 999
      epc = -0.0504 * flarge(rs/30)
      efc = -0.0254 * flarge(rs/75)
      xmup = -0.0504*log(1.0+30.0/rs)
!     xmuf = -0.0254*log(1.0+75.0/rs)
      vu = gamma*(efc - epc)
!     tau = xmuf-xmup-(efc-epc)*4.0/3.0
     
      alg = -1.22177412/rs + vu
      blg = xmup - vu
      vxc = alg*xmag**(1.0/3.0) + blg
!     vxc = alg*xmag**(1.0/3.0) + blg +tau*fsmall(xmag/2.0)

 999  continue
!     transform to code units (Hartrees) from Rydbergs
      vxc = vxc / 2.d0

      return
      end

      double precision function flarge(x)
      implicit double precision (a-h, o-z)
        flarge = (1+x**3)*log(1+1/x) + x/2 - x**2 - 1.0/3.0
      return
      end

!     double precision function fsmall(x)
!     implicit double precision (a-h, o-z)
!     parameter (a = 2.0**(-1.0/3.0) )
!       fsmall = ( x**(4/3) + (1.0-x)**(4/3) - a ) / (1.0-a)
!     return
!     end
