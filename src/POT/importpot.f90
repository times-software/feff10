

         subroutine importpot(vtot,vint,edens,rhoint,xmu)  
!         subroutine importpot(rnrmav, xmu, vint, rhoint, emu, erelax, ecv,rs,xf,        &
!                       edens, vclap, vtot, edenvl, vvalgs, qnrm, xnmues, inters, totvol)  

! For now, this is little more than a wrapper for rdpot.  However, more functionality may be added later.
! Kevin Jorissen 12/2010.
! Reads potentials from a pot.bin file and sets them up as starting point for a scf cycle.


      use dimsmod, only: nphx=>nphu, lx, nheadx
      implicit double precision (a-h, o-z)

      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(30,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
      dimension eorb(30), kappa(30)
      dimension iorb(-4:3,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)

      character*80 title(nheadx)


      call rdpot(ntitle, title, rnrmav, xmu, vint, rhoint, emu, s02, erelax, wp, ecv,rs,xf, qtotel,        &
                       imt, rmt, inrm, rnrm, folp, folpx, xnatph, dgc0, dpc0, dgc, dpc, adgc, adpc,               &
                       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,  eorb, kappa, iorb, &
                        qnrm, xnmues, nohole, ihole,inters, totvol, iafolp, xion, iunf, iz, jumprm)


      return
	  end
	  
