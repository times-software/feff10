*234567890
      program testedgec
      implicit double precision (a-h, o-z)
      complex*16 senergy,phshift
      complex*16 xintpm2
      complex*16 xmatel,pmm,graterc,cme
      double precision qkmin,qkmax,abr,rlr,xsing,error,qsing
c     maximum number of atoms for FMS. Reduce nclusx if you need
c     smaller executable.
      parameter (nclusx=175)
c     max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c     max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c     max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =10000)
c     max orbital momentum for FMS module.
      parameter (lx=3)
c     max number of unique potentials (potph)
      parameter (nphx = 7)
c     max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c     Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c     Number of energy points genfmt, etc.
      parameter (nex = 150)
c     Max number of distinct lambda's for genfmt
c     15 handles iord 2 and exact ss
      parameter (lamtot=15)
c     vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c     max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c     matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c     max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c     max number of header lines
      parameter (nheadx=30)
      character*80 text
      character*6  potlblf,potlbli
      dimension text(40),  potlblf(0:nphx),potlbli(0:nphx)
      dimension ltext(40)
      complex*16 phf(nex,-ltot:ltot,nspx,0:nphx), 
     2           ereff(nex,nspx), emf(nex)
      complex*16 phi(nex,-ltot:ltot,nspx,0:nphx), 
     2           erefi(nex,nspx), emi(nex)
      complex*16 deltaph(nex)
      complex*16 rkkf(nex,8,nspx),rkki(nex,8,nspx),
     2           rkk(nex,8,nspx)
      dimension lmaxf(0:nphx),lmaxi(0:nphx)
      dimension izf(0:nphx),izi(0:nphx)
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      parameter (xkb=3.1666643e-6)
      complex*16 pk
      double precision w,T,ef,qf,xmui,xmuf
      double precision xnormi(nex),xnormf(nex),xnorm(nex)
      complex*16 xsecti(nex),xsectf(nex),xsect(nex)
      character*10 filename

      call rdxspi (ne, ne1, ne3, nph, ihole, rnrmai,xmui,edgei,
     1               ik0i, emi, erefi, izi, potlbli, phi, rkki, 
     2               lmaxi, lmaxp1i)
      filename='xsecti.bin'
      call rdnorm (filename,emi,xnormi,xsecti)
      call rdxspf (ne, ne1, ne3, nph, ihole, rnrmavf,xmuf,edgef,
     1               ik0f, emf, ereff, izf, potlblf, phf, rkkf, 
     2               lmaxf, lmaxp1f)
      filename='xsectf.bin'
      call rdnorm (filename,emf,xnormf,xsectf)
*      write(6,*) ne,ne1,ne-ne1-ne3,ne3,nex-ne3
*      write(6,*) ((lmaxf(i,j),j=0,nphx),i=1,5)
*      write(6,*) (izf(i),i=0,nphx)
*      write(6,*) (potlblf(i),i=0,nphx)
      ef=dble(emf(11))
      qf=dble(sqrt(2.d0*(emf(11)-ereff(11,1))))
* angular momentum of absorption edge
      labs=1
      do i=1,ne1
        deltaph(i)=phf(i,labs+1,1,0)-phi(i,labs+1,1,0)
*        pk=sqrt(2.d0*(emf(i)-ereff(i,1)))
*        write(21,700) dble(pk)/qf,phi(i,labs+1,1,0)
      enddo
* Temperature in Kelvins
*      T=300.d0
      T=1.d0
      nsp=1
      call edgec(ne1,T,deltaph,emf,ereff,rkkf,rkki,rkk,
     2  xnormi,xnormf,xnorm,xsecti,xsectf,xsect,xmuf)
      call wrxsph (nsp, ne, ne1, ne3, nph, ihole, rnrmavf,xmuf,edgef,
     1              ik0f, emf, ereff, lmaxf, izf, potlblf, phf, rkk)
      call wrnorm(ne,emf,xnorm,xsect)
500   format(1x,5(e12.5,1x))
501   format(1x,5(f12.5,1x))
700   format(1x,7(f10.5,1x))
      stop
      end

      subroutine edgec(ne1,T,deltaph,emf,ereff,rkkf,rkki,rkk,
     2  xnormi,xnormf,xnorm,xsecti,xsectf,xsect,xmu)
      implicit double precision (a-h, o-z)
      complex*16 senergy,phshift
      complex*16 xmatel,pmm,graterc,cme
      double precision qkmin,qkmax,abr,rlr,xsing,error,qsing
c     max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c     Number of energy points genfmt, etc.
      parameter (nex = 150)
c     max number of unique potentials (potph)
      parameter (nphx = 7)
c     max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
      complex*16 deltaph(nex)
      complex*16 ereff(nex,nspx),emf(nex)
      complex*16 rkkf(nex,8,nspx),rkki(nex,8,nspx),
     2           rkk(nex,8,nspx)
      complex*16 cmatelf(nex,8,nspx),cmateli(nex,8,nspx),
     2           cmatel(nex,8,nspx)
      complex*16 ereff2(nex,nspx),emf2(nex)
      complex*16 rkki2(nex,8,nspx), coni
      parameter (coni=(0.d0,1.d0))
      double precision xnormi(nex),xnormf(nex),xnorm(nex)
      complex*16 xsecti(nex),xsectf(nex),xsect(nex)
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      parameter (xkb=3.1666643e-6)
      complex*16 pk
      double precision w,T,ef,qf,xmu
      common /econe/ pk,w,T2,ef,qf
      common /ectwo/ emf2,ereff2
      common /ecfour/ cmateli
      common /ecfive/ ipol
      external xmatel,graterc,xintpmr,xintpmi

      T2=T
      ifermi=1
      xmin=1.d30
      do i=1,nex
        emf2(i)=emf(i)
        if(xmin.gt.abs(dble(xmu-emf(i)))) then
           xmin = abs(dble(xmu-emf(i)))
           ifermi = i
        end if
        do j=1,nspx
          ereff2(i,j)=ereff(i,j)
          do k=1,8
            cmateli(i,k,j)=rkki(i,k,j)*dsqrt(xnormi(i))
            cmatelf(i,k,j)=rkkf(i,k,j)*dsqrt(xnormf(i))
          enddo
        enddo
      enddo
      print*, ifermi
      ef=xmu
      qf=dble(sqrt(2.d0*(xmur-ereff(ifermi,1))))
      qkmax=qf+dsqrt(10.d0*T*xkb)
      qkmin=0.d0
*      abr=1.d-4
*      rlr=1.d-4
      abr=3.d-5
      rlr=3.d-5
      nsing=0
      pi=dacos(-1.d0)

      iwrite=1
      do i=1,nex
*      do i=1,ne1
*      do i=iwrite,iwrite
        pk=sqrt(2.d0*(emf(i)-ereff(i,1)))
        w=dble(emf(i))
        qsing=dble(pk)
        rkpolavg=0.d0
        cmpolavg=0.d0
        do ipol=1,8

          if (i.eq.0.and.ipol.eq.1) then
            npts=400
            pmr=0.d0
            pmi=0.d0
            dq=(qkmax-qkmin)/npts
            do j=1,npts
              qk=qkmin+j*dq
              addr=xintpmr(qk)
              addi=xintpmi(qk)
              pmr=pmr+addr*dq
              pmi=pmi+addi*dq
              write(10,500) qk/qf,addr,pmr
              write(11,500) qk/qf,addi,pmi
            enddo
            close(10)
            close(11)
          endif

          if (qsing.lt.qkmax) then
            pmr=grater(xintpmr,qkmin,qsing,
     2                 abr,rlr,nsing,xsing,error,numcal,maxns)
            pmi=grater(xintpmi,qkmin,qsing,
     4                 abr,rlr,nsing,xsing,error,numcal,maxns)
            pmr=pmr+grater(xintpmr,qsing,qkmax,
     2                 abr,rlr,nsing,xsing,error,numcal,maxns)
            pmi=pmi+grater(xintpmi,qsing,qkmax,
     4                 abr,rlr,nsing,xsing,error,numcal,maxns)
          else
            pmr=grater(xintpmr,qkmin,qkmax,
     2                 abr,rlr,nsing,xsing,error,numcal,maxns)
            pmi=grater(xintpmi,qkmin,qkmax,
     4                 abr,rlr,nsing,xsing,error,numcal,maxns)
          endif
*
          pmm=cmatelf(i,ipol,1)-(pmr+coni*pmi)*sin(deltaph(i))
          cmatel(i,ipol,1)=pmm
          if (ipol.eq.1.and.i.le.100) then
            write(13,700) real(i),emf(i)-ereff(i,1)/eV,ereff(i,1)/eV,
     2                    emf(i)/eV
            write(14,500) real(i),dble(pmm),dimag(pmm),
     2                  dble(cmatelf(i,ipol,1)),
     3                  dimag(cmatelf(i,ipol,1))
            write(15,500) real(emf(i)-emf(11)), abs(pmm)**2, 
     2                  abs(cmatelf(i,ipol,1))**2,
     3                  abs(cmateli(i,ipol,1))**2
          endif
          rkpolavg=rkpolavg+abs(rkkf(i,ipol,1))
          cmpolavg=cmpolavg+abs(cmatel(i,ipol,1))
        enddo
        if (rkpolavg.ne.0.d0) then
          xnorm(i)=(cmpolavg/rkpolavg)**2
        else
          xnorm(i)=xnorm(i-1)
        endif
        do ipol=1,8
          rkk(i,ipol,1)=cmatel(i,ipol,1)/dsqrt(xnorm(i))
        enddo
        xsect(i)=xsectf(i)*dsqrt(xnorm(i)/xnormf(i))
        write(16,500) real(i),xnorm(i),xnormf(i)
        write(17,500) real(i),xsect(i),xsectf(i)
        write(18,500) real(i),rkk(i,1,1),rkkf(i,1,1)
      enddo

500   format(1x,5(e12.5,1x))
501   format(1x,5(f12.5,1x))
700   format(1x,7(f10.5,1x))
      return
      end

      double precision function xintpmr(qk)
      implicit double precision (a-h, o-z)
      double precision qk,ek
      complex*16 pk,senergy
      complex*16 xmatel, coni
      parameter (coni=(0.d0,1.d0))
      parameter (nex = 150)
      parameter (nspx=1)
      double precision w,T,ef,qf
      common /econe/ pk,w,T,ef,qf
      complex*16 cmateli(nex,8,nspx)
      common /ecfour/ cmateli
      common /ecfive/ ipol
      external fermi,senergy,xmatel
      ek=qk**2/2.d0+dble(senergy(qk))
      xintpmr=dble(fermi(ek,T,ef)*xmatel(qk,cmateli,ipol)
     2        /(qk-(dble(pk)+coni*qf*1.d-2)))
*     2        /((qk,0.d0)-pk-(0,qf*1.d-2)))
      return
      end

      double precision function xintpmi(qk)
      implicit double precision (a-h, o-z)
      double precision qk,ek
      complex*16 pk,senergy
      complex*16 xmatel, coni
      parameter (coni=(0.d0,1.d0))
      parameter (nex = 150)
      parameter (nspx=1)
      double precision w,T,ef,qf
      common /econe/ pk,w,T,ef,qf
      complex*16 cmateli(nex,8,nspx)
      common /ecfour/ cmateli
      common /ecfive/ ipol
      external fermi,senergy,xmatel
      ek=qk**2/2.d0+dble(senergy(qk))
      xintpmi=dimag(fermi(ek,T,ef)*xmatel(qk,cmateli,ipol)
     2        /(qk-(dble(pk)+coni*qf*1.d-2)))
*     2        /((qk,0.d0)-pk-(0,qf*1.d-2)))
      return
      end

      double precision function fermi(w,T,ef)
      implicit double precision (a-h, o-z)
      double precision w,T,ef,xkb
*     Boltzman constant in Hartree/K
      parameter (xkb=3.1666643e-6)
      fermi=1.d0/(exp((w-ef)/(xkb*T))+1.d0)
500   format(1x,5(e12.5,1x))
      return
      end

      complex*16 function senergy(qk)
      implicit double precision (a-h, o-z)
      double precision qk
      parameter (nex = 150)
      complex*16 eref(nex,1), em(nex)
      common /ectwo/ em,eref
      do i=1,99
        q1=sqrt(2.d0*(em(i)-eref(i,1)))
        q2=sqrt(2.d0*(em(i+1)-eref(i+1,1)))
        if(qk.gt.q1.and.qk.le.q2) then
          senergy=eref(i,1)+(eref(i+1,1)-eref(i,1))*(qk-q1)/(q2-q1)
          return
        endif
      enddo
      q1=sqrt(2.d0*(em(1)-eref(1,1)))
      if(qk.le.q1) then
        senergy = eref(1,1)
        return
      endif
      q1=sqrt(2.d0*(em(100)-eref(100,1)))
      if(qk.gt.q1) then
        senergy = eref(100,1)
        return
      endif
      return
      end
      
      complex*16 function xmatel(qk,cmatel,ipol)
      implicit double precision (a-h, o-z)
      double precision qk
      parameter (nex = 150)
      parameter (nspx=1)
      complex*16 cmatel(nex,8,nspx)
      complex*16 eref(nex,1), em(nex)
      common /ectwo/ em,eref
      do i=1,99
        q1=sqrt(2.d0*(em(i)-eref(i,1)))
        q2=sqrt(2.d0*(em(i+1)-eref(i+1,1)))
        if(qk.gt.q1.and.qk.le.q2) then
          xmatel=cmatel(i,ipol,1)+(cmatel(i+1,ipol,1)-cmatel(i,ipol,1))
     2           *(qk-q1)/(q2-q1)
          return
        endif
      enddo
      q1=sqrt(2.d0*(em(1)-eref(1,1)))
      if(qk.le.q1) then
        xmatel = cmatel(1,ipol,1)
        return
      endif
      q1=sqrt(2.d0*(em(100)-eref(100,1)))
      if(qk.gt.q1) then
        xmatel = cmatel(100,ipol,1)
        return
      endif
      return
      end

      subroutine rdnorm(filename,em,xnorm,xsect)
      implicit double precision (a-h, o-z)
      parameter (nex = 150)
      double precision xnorm(nex)
      complex*16 xsect(nex),em(nex),emt(nex),ediff, coni
      parameter (coni=(0.d0,1.d0))
      character*10 filename
      character*96 cline
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      parameter (xkb=3.1666643e-6)

      open (unit=31,status='old',file=filename)
      open (unit=32,status='unknown',file='xsectheader.dat')
 12   read(31,'(A)') cline
      last=80
* Trim off excess white space.
 13   if (cline(last:last).eq.' ') then
        last=last-1
        goto 13
      endif
      write(32,'(A)') cline(1:last)
      if (cline(6:14).ne.'---------') goto 12
      do i=1,3
        read(31,'(A)') cline
        write(32,'(A)') cline(1:last)
      enddo
      close(32)
      do i=1,nex
        xsect(i)=(0.d0,0.d0)
        xnorm(i)=1.d0
      enddo
      do i=1,nex
        read(31,*,end=30) er,ei,xnorm(i),xsr,xsi
        emt(i)=(er+coni*ei)
        xsect(i)=(xsr+coni*xsi)
        ediff=emt(i)*eV-em(i)
        if (real(ediff)/real(em(i)).gt.0.0001
     2  .or.imag(ediff)/imag(em(i)).gt.0.0001) then
          write(6,*) 'energy grids in xsect.bin, phase.bin do not match'
          write(6,*) emt(i),em(i)/eV
          stop
        endif
        ii=i
      enddo
 30   continue
      do i=ii+1,nex
        xnorm(i)=xnorm(i-1)
      enddo
      return
      end

      subroutine wrnorm(ne,em,xnorm,xsect)
      implicit double precision (a-h, o-z)
      parameter (nex = 150)
      double precision xnorm(nex)
      complex*16 xsect(nex),em(nex)
      character*96 cline
      parameter (aangstrom=1.d0/0.52917706d0,eV=1.d0/27.21160d0)
      parameter (xkb=3.1666643e-6)

      open (unit=31,status='old',file='xsectheader.dat')
      open (unit=32,status='unknown',file='xsect.bin')
 12   read(31,'(A)',end=20) cline
      last=80
* Trim off excess white space.
 13   if (cline(last:last).eq.' ') then
        last=last-1
        goto 13
      endif
      write(32,'(A)') cline(1:last)
      goto 12
 20   close(unit=31)
      do i=1,ne
        write(32,100) dble(em(i))/eV,dimag(em(i))/eV,xnorm(i),
     2                dble(xsect(i)),dimag(xsect(i)) 
      enddo
 100  format (e17.9, 4e13.5)
      close(unit=32)
      return
      end
       
      subroutine rdxspf (ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,
     1               ik0, em, eref, iz, potlbl, ph, rkk, lmax0, lmaxp1)
      implicit double precision (a-h, o-z)
c     reads file 'phasef.bin'
c  Energy grid information
c     em   - complex energy grid
c     eref - V_int + i*gamach/2 + self-energy correction
c     ne   - total number of points in complex energy grid
c     ne1  - number of points on main horizontal axis
c     ne2  - number of points on vertical vertical axis ne2=ne-ne1-ne3
c     ne3  - number of points on auxilary horizontal axis (need for f')
c     xmu  - Fermi energy
c     edge - x-ray frequency for final state at Fermi level
c     ik0  - grid point index at Fermi level
c  Potential type information
c     nph - number of potential types
c     iz  - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax - max orb momentum for each potential type
c     ihole - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph  - complex scattering phase shifts
c     rkk - complex multipole matrix elements

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
       parameter (nclusx=175)
c      max number of spins: 1 for spin average; 2 for spin-dep
       parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
       parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
       parameter (nattx =10000)
c      max orbital momentum for FMS module.
       parameter (lx=3)
c      max number of unique potentials (potph)
       parameter (nphx = 7)
c      max number of ang mom (arrays 1:ltot+1)
       parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
       parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
       parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
       parameter (lamtot=15)
c      vary mmax and nmax independently
       parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
       parameter (npatx = 8)
c      matches path finder, used in GENFMT
       parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
       parameter (novrx=8)
c      max number of header lines
       parameter (nheadx=30)
c= ../HEADERS/dim.h}

      character*6  potlbl
      dimension  potlbl(0:nphx)

      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex,8,nspx)
      dimension lmax0(0:nphx), lmax(nex,0:nphx)
      dimension iz(0:nphx)
c     kinit, linit, ilinit,  - initial state kappa and ang. mom.
c     lmaxp1  -largest lmax in problem + 1

c     phmin is min value to use for |phase shift|
      parameter (phmin = 1.0d-7)

c     Local staff
c     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)
      character*128 str
      parameter (nwordx=20)
      character*20 words(nwordx)

      open (unit=1, file='phasef.bin', status='old', iostat=ios)
      call chopen (ios, 'phasef.bin', 'rdxspf')

      read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
  10  format (8(1x,i4))

      call rdpadd(1, npadx, dum(1), 3)
      rnrmav = dum(1)
      xmu    = dum(2)
      edge   = dum(3)

      call rdpadx(1, npadx, em(1), ne)
c     call rdpadx(1, npadx, eref(1), ne)
      call rdpadx (1, npadx, temp(1), ne*nsp)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        eref (ie, isp) = temp(ii)
  60  continue

      do 80  iph = 0, nph
         read(1, 20)  lmax0(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)

         do 75 isp = 1,nsp
            ii = ne * (2*lmax0(iph)+1)
            call rdpadx (1, npadx, temp(1), ii )
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax0(iph), lmax0(iph)
               ii = ii+ 1
               ph(ie,ll,isp,iph) = temp(ii)
   70       continue
   75    continue
   80 continue

      call rdpadx (1, npadx, temp(1), ne*8*nsp)
      ii = 0
      do 90 isp = 1,nsp
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        rkk (ie, kdif, isp) = temp(ii)
  90  continue

      close (unit=1)

c     make additional data for output
      lmaxp1 = 0
      do 180  iph = 0, nph
      do 180  ie = 1, ne
c        Set lmax to include only non-zero phases
         do 160  il =  lmax0(iph), 0, -1
            lmax(ie,iph) = il
            if (abs(sin(ph(ie, il, 1, iph))) .gt. phmin .or.
     3          abs(sin(ph(ie, il,nsp,iph))) .gt. phmin)  goto 161
  160    continue
  161    continue
         if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
  180 continue

      return
      end

      subroutine rdxspi (ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,
     1               ik0, em, eref, iz, potlbl, ph, rkk, lmax0, lmaxp1)
      implicit double precision (a-h, o-z)
c     reads file 'phase.bin'
c  Energy grid information
c     em   - complex energy grid
c     eref - V_int + i*gamach/2 + self-energy correction
c     ne   - total number of points in complex energy grid
c     ne1  - number of points on main horizontal axis
c     ne2  - number of points on vertical vertical axis ne2=ne-ne1-ne3
c     ne3  - number of points on auxilary horizontal axis (need for f')
c     xmu  - Fermi energy
c     edge - x-ray frequency for final state at Fermi level
c     ik0  - grid point index at Fermi level
c  Potential type information
c     nph - number of potential types
c     iz  - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax - max orb momentum for each potential type
c     ihole - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph  - complex scattering phase shifts
c     rkk - complex multipole matrix elements

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
       parameter (nclusx=175)
c      max number of spins: 1 for spin average; 2 for spin-dep
       parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
       parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
       parameter (nattx =10000)
c      max orbital momentum for FMS module.
       parameter (lx=3)
c      max number of unique potentials (potph)
       parameter (nphx = 7)
c      max number of ang mom (arrays 1:ltot+1)
       parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
       parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
       parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
       parameter (lamtot=15)
c      vary mmax and nmax independently
       parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
       parameter (npatx = 8)
c      matches path finder, used in GENFMT
       parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
       parameter (novrx=8)
c      max number of header lines
       parameter (nheadx=30)
c= ../HEADERS/dim.h}

      character*6  potlbl
      dimension  potlbl(0:nphx)

      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex,8,nspx)
      dimension lmax0(0:nphx), lmax(nex,0:nphx)
      dimension iz(0:nphx)
c     kinit, linit, ilinit,  - initial state kappa and ang. mom.
c     lmaxp1  -largest lmax in problem + 1

c     phmin is min value to use for |phase shift|
      parameter (phmin = 1.0d-7)

c     Local staff
c     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)
      character*128 str
      parameter (nwordx=20)
      character*20 words(nwordx)

      open (unit=1, file='phasei.bin', status='old', iostat=ios)
      call chopen (ios, 'phasei.bin', 'rdxspi')

      read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
  10  format (8(1x,i4))

      call rdpadd(1, npadx, dum(1), 3)
      rnrmav = dum(1)
      xmu    = dum(2)
      edge   = dum(3)

      call rdpadx(1, npadx, em(1), ne)
c     call rdpadx(1, npadx, eref(1), ne)
      call rdpadx (1, npadx, temp(1), ne*nsp)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        eref (ie, isp) = temp(ii)
  60  continue

      do 80  iph = 0, nph
         read(1, 20)  lmax0(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)

         do 75 isp = 1,nsp
            ii = ne * (2*lmax0(iph)+1)
            call rdpadx (1, npadx, temp(1), ii )
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax0(iph), lmax0(iph)
               ii = ii+ 1
               ph(ie,ll,isp,iph) = temp(ii)
   70       continue
   75    continue
   80 continue

      call rdpadx (1, npadx, temp(1), ne*8*nsp)
      ii = 0
      do 90 isp = 1,nsp
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        rkk (ie, kdif, isp) = temp(ii)
  90  continue

      close (unit=1)

c     make additional data for output
      lmaxp1 = 0
      do 180  iph = 0, nph
      do 180  ie = 1, ne
c        Set lmax to include only non-zero phases
         do 160  il =  lmax0(iph), 0, -1
            lmax(ie,iph) = il
            if (abs(sin(ph(ie, il, 1, iph))) .gt. phmin .or.
     3          abs(sin(ph(ie, il,nsp,iph))) .gt. phmin)  goto 161
  160    continue
  161    continue
         if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
  180 continue

      return
      end

      subroutine wrxsph (nsp, ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,
     1                   ik0, em, eref, lmax, iz, potlbl, ph, rkk)
      implicit double precision (a-h, o-z)
c     writes down file 'phase.bin' to be read by rphbin
c  Energy grid information
c     em   - complex energy grid
c     eref - V_int + i*gamach/2 + self-energy correction
c     ne   - total number of points in complex energy grid
c     ne1  - number of points on main horizontal axis
c     ne2  - number of points on vertical vertical axis ne2=ne-ne1-ne3
c     ne3  - number of points on auxilary horizontal axis (need for f')
c     xmu  - Fermi energy
c     edge - x-ray frequency for final state at Fermi level
c     ik0  - grid point index at Fermi level
c  Potential type information
c     nph - number of potential types
c     iz  - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax - max orb momentum for each potential type
c     ihole - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph  - complex scattering phase shifts
c     rkk - complex multipole matrix elements

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
       parameter (nclusx=175)
c      max number of spins: 1 for spin average; 2 for spin-dep
       parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
       parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
       parameter (nattx =10000)
c      max orbital momentum for FMS module.
       parameter (lx=3)
c      max number of unique potentials (potph)
       parameter (nphx = 7)
c      max number of ang mom (arrays 1:ltot+1)
       parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
       parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
       parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
       parameter (lamtot=15)
c      vary mmax and nmax independently
       parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
       parameter (npatx = 8)
c      matches path finder, used in GENFMT
       parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
       parameter (novrx=8)
c      max number of header lines
       parameter (nheadx=30)
c= ../HEADERS/dim.h}

      character*6  potlbl
      dimension  potlbl(0:nphx)

      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex, 8, nspx)
      dimension lmax(0:nphx)
      dimension iz(0:nphx)

c     Local staff
c     npadx control padlib precision (see padlib package)
      parameter (npadx=8)
c     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)

      open (unit=1, file='phase.bin', status='unknown', iostat=ios)
      call chopen (ios, 'phase.bin', 'wrxsph')

      write(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
  10  format (8(1x,i4))

      dum(1) = rnrmav
      dum(2) = xmu
      dum(3) = edge
      call wrpadd(1, npadx, dum(1), 3)

      call wrpadx(1, npadx, em(1), ne)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        temp(ii) = eref (ie, isp)
  60  continue
      call wrpadx (1, npadx, temp(1), ii)

      do 80  iph = 0, nph
         write(1, 20) lmax(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)
         do 75  isp = 1, nsp
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax(iph), lmax(iph)
               ii = ii+ 1
               temp(ii) = ph(ie, ll, isp, iph)
   70       continue
            call wrpadx (1, npadx, temp(1), ii )
   75    continue
   80 continue

      ii = 0
      do 90 isp = 1, nsp
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        temp(ii) = rkk (ie, kdif, isp)
  90  continue
      call wrpadx (1, npadx, temp(1), ii)

      close (unit=1)

      return
      end

       subroutine wrpadx(iout,npack,array,npts)
c write complex*16 array as pad string
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       complex*16 array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = dimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end

       subroutine pad(xreal,npack,str)
c  convert dp number *xreal* to packed-ascii-data string *str*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer  iexp, itmp, isgn, i, npack, iok, j
       double precision xreal, xwork, xsave,onem, tenth
       parameter (onem  =  0.99999999997d0)
       parameter (tenth =  0.099999999994d0)
       character str*(*)
c
       str      = ' '
       xsave    = min(huge, max(-huge, xreal))
       isgn     = 1
       if (xsave.le.0) isgn = 0
c
       xwork    = dabs( xsave )
       iexp     = 0
       if ((xwork.lt.huge).and.(xwork.gt.tiny))  then
          iexp  =   1 + int(log(xwork) / tenlog  )
       else if (xwork.ge.huge) then
          iexp  = ihuge
          xwork = one
       else if (xwork.le.tiny)  then
          xwork = zero
       end if
c force xwork between ~0.1 and ~1
c note: this causes a loss of precision, but
c allows backward compatibility
       xwork    = xwork / (ten ** iexp)
 20    continue
       if (xwork.ge.one) then
          xwork = xwork * 0.100000000000000d0
          iexp  = iexp + 1
       else if (xwork.le.tenth) then
          xwork = xwork * ten
          iexp  = iexp - 1
       endif
       if (xwork.ge.one) go to 20

       itmp     = int ( ibas2 * xwork )
       str(1:1) = char(iexp  + ioff + ibas2 )
       str(2:2) = char( 2 * itmp + isgn + ioff)
       xwork    = xwork * ibas2 - itmp
       if (npack.gt.2) then
          do 100 i = 3, npack
             itmp     = int( base * xwork + 1.d-9)
             str(i:i) = char(itmp + ioff)
             xwork    = xwork * base - itmp
 100      continue
       end if
       if (xwork.ge.0.5d0) then
          i = itmp + ioff + 1
          if (i.le.126) then
             str(npack:npack)= char(i)
          else
             j = ichar(str(npack-1:npack-1))
             if (j.lt.126) then
                str(npack-1:npack-1) = char(j+1)
                str(npack:npack)     = char(37)
             endif
          endif
       endif
       return
       end


       subroutine wrpadd(iout,npack,array,npts)
c
c write a dp array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       double precision array(*), xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = array(i)
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end

       double precision function unpad(str,npack)
c
c  convert packed-ascii-data string *str* to dp number *unpad*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       double precision sum
       integer   iexp, itmp, isgn, i, npack
       character str*(*)
       unpad = zero
       if (npack.le.2) return
       iexp  =     (ichar(str(1:1)) - ioff   ) - ibas2
       isgn  = mod (ichar(str(2:2)) - ioff, 2) * 2 - 1
       itmp  =     (ichar(str(2:2)) - ioff   ) / 2
       sum   = dble(itmp/(base*base))
c       do 100 i = 3, npack
c          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
c 100   continue
       do 100 i = npack, 3, -1
          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
 100   continue
       unpad = 2 * isgn * ibase * sum * (ten ** iexp)
cc       print*, sum, iexp,unpad
       return
       end


       subroutine rdpadd(iou,npack,array,npts)
c read dparray from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       double precision    array(*), unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = tmp
             if (ipts.ge.npts) go to 50
 30       continue
          go to 10
 50    continue
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end

C SUBROUTINE UPPER (STRING)  Changes a-z to upper case.

      SUBROUTINE UPPER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 97)  .OR.  (IC .GT. 122))  GOTO 10
         STRING (I:I) = CHAR (IC - 32)
   10 CONTINUE

      RETURN
      END


      subroutine wlog (string)
      character*(*) string

c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process,
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c     This output routine is used to replace the PRINT statement
c     for output that "goes to the terminal", or to the log file.
c     If you use a window based system, you can modify this routine
c     to handle the running output elegantly.
c     Handle carriage control in the string you pass to wlog.
c
c     The log file is also written here, hard coded here.

c     The log file is unit 11.  The log file is opened in the
c     main program, program feff.

c     make sure not to write trailing blanks

   10 format (a)

c     Suppress output in sequential loops
      if (par_type .eq. 2) return

      il = istrln (string)
      if (il .eq. 0)  then
         print10
         if (par_type .ne. 3) write(11,10)
      else
         print10, string(1:il)
         if (par_type .ne. 3) write(11,10) string(1:il)
      endif
      return
      end
      subroutine lblank (string)
      character*(*) string
c     add a leading blank, useful for carriage control
      string = ' ' // string
      return
      end
      double precision function xx (j)
      implicit double precision (a-h, o-z)
c     x grid point at index j, x = log(r), r=exp(x)
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     xx = -8.8 + (j-1)*0.05
      xx = -c88 + (j-1)*delta
      return
      end


       subroutine rdpadx(iou,npack,array,npts)
c read complex*16 array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex*16  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue
          go to 10
 50    continue
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end

       integer function iread(lun,string)
c
c generalized internal read:
c    read a string the next line of an opened file
c    unit, returning the real length of string
c
c inputs:  
c   lun     opened file unit number
c outputs:
c   string  string read from file
c returns:
c   iread   useful length of string, as found from
c                  sending string to 'sclean' to
c                  remove non-printable characters
c                   and then istrln
c           or
c              -1   on 'end-of-file'
c              -2   on 'error'
c
c copyright (c) 1999  Matthew Newville
       implicit none
       character*(*) string
       integer    lun, istrln
       external   istrln
       string = ' '
 10    format(a)
       read(lun, 10, end = 40, err = 50) string
       call sclean(string)
       iread = istrln(string)
       return
 40    continue
       string = ' '
       iread = -1
       return
 50    continue
       string = ' '
       iread = -2
       return
       end

C FUNCTION ISTRLN (STRING)  Returns index of last non-blank
C                           character.  Returns zero if string is
C                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '   ')
C     there is a tab character here  ^

C  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

C  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
C SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*200  TMP
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '   ')
C     there is a tab character here  ^

      JLEN = ISTRLN (STRING)

C  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

C  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 20
   10 CONTINUE
   20 CONTINUE

C  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

C  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END

       subroutine sclean(str)
c
c  clean a string, especially for strings passed between
c  different file systems, or from C functions:
c
c   1. characters in the range char(0), or char(10)...char(15)
c      are interpreted as end-of-line characters, so that all
c      remaining characters are explicitly blanked.
c   2. all other characters below char(31) (including tab) are
c      replaced by a single blank
c
c  this is mostly useful when getting a string generated by a C
c  function and for handling dos/unix/max line-endings.
c
c copyright (c) 1999  Matthew Newville
       character*(*) str, blank*1
       parameter (blank = ' ')
       integer i,j,is
       do 20 i = 1, len(str)
          is = ichar(str(i:i))
          if ((is.eq.0) .or. ((is.ge.10) .and. (is.le.15))) then
             do 10 j= i, len(str)
                str(j:j) = blank
 10          continue
             return
          elseif (is.le.31)  then
             str(i:i)  = blank
          end if
 20    continue
       return
c end subroutine sclean
       end

      subroutine setkap(ihole, kinit, linit)
      implicit double precision (a-h, o-z)

c     Set initial state ang mom and quantum number kappa
c     ihole  initial state from ihole
c     1      K    1s      L=0 -> linit=0
c     2      LI   2s      L=0 -> linit=0
c     3      LII  2p 1/2  L=1 -> linit=1
c     4      LIII 2p 3/2  L=1 -> linit=1
c     5+     etc.
      if (ihole.le. 2 .or. ihole.eq. 5 .or. ihole.eq.10 .or.
     1    ihole.eq.17 .or. ihole.eq.24 .or. ihole.eq.27)  then
c        hole in s state
         linit = 0
         kinit = -1
      elseif (ihole.eq. 3 .or. ihole.eq. 6 .or. ihole.eq.11 .or.
     1        ihole.eq.18 .or. ihole.eq.25 .or. ihole.eq.30)  then
c        hole in p 1/2 state
         linit = 1
         kinit = 1
      elseif (ihole.eq. 4 .or. ihole.eq. 7 .or. ihole.eq.12 .or.
     1        ihole.eq.19 .or. ihole.eq.26)  then
c        hole in p 3/2 state
         linit = 1
         kinit = -2
      elseif (ihole.eq. 8 .or. ihole.eq.13 .or.
     1        ihole.eq.20 .or. ihole.eq.27)  then
c        hole in d 3/2 state
         linit = 2
         kinit = 2
      elseif (ihole.eq. 9 .or. ihole.eq.14 .or.
     1        ihole.eq.21 .or. ihole.eq.28)  then
c        hole in d 5/2 state
         linit = 2
         kinit = -3
      elseif (ihole.eq.15 .or. ihole.eq.22)  then
c        hole in  f 5/2 state
         linit = 3
         kinit = 3
      elseif (ihole.eq.16 .or. ihole.eq.23)  then
c        hole in  f 7/2 state
         linit = 3
         kinit = -4
      else
c        some unknown hole
         stop 'invalid hole number in setkap'
      endif

      return
      end


      subroutine chopen (ios, fname, mod)
c     Writes error msg and stops if error in ios flag from open
c     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod
      character*512 slog

c     open successful
      if (ios .le. 0)  return

c     error opening file, tell user and die.
*      i = istrln(fname)
*      j = istrln(mod)
*      write(slog,100)  fname(1:i), mod(1:j)
*      call wlog(slog)
      write(slog,100)  fname, mod

  100 format (' Error opening file, ', a,
     2        ' in module ', a)

*      call wlog(' Fatal error')
      stop 'CHOPEN'
      end


c**********************************************************************
c   This is Steve White's rewrite of Mike Teter's integration routine.
c   Modified by J. Rehr for complex integration.
c   The following is a listing of the arguments in the initial function
c   statement:
c      fn    -- routine requires external function statement in MAIN
c      xmin  -- lower limit
c      xmax  -- upper limit
c      abr   -- absolute tolerable error
c      rlr   -- relative tolerable error
c      nsing -- number of singularities or regions requiring
c                   special attention
c      xsing -- array of locations of singularities or endpoints
c                   of special regions
c      error -- output for routine error messages
c      numcal-- the number of times fn was called
c      maxns -- the maximum number of regions being considered simultaneously
c       function grater(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
c       fn declared double precision
c       double precision function grater(fn,xmin,xmax,abr,rlr,
c       fn declared complex*16
c      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1

       double precision function grater(fn,xmin,xmax,abr,rlr,
     1 nsing,xsing,error,numcal,maxns)

       implicit double precision (a-h,o-z)
       parameter (mx=1500)
       dimension xleft(mx),fval(3,mx),dx(3),wt(3)
       dimension wt9(9), xsing(20)
       external fn
        logical atsing
        save dx,wt,wt9
        data dx/0.1127016653792583  ,0.5  ,0.8872983346207417  /
        data wt/0.277777777777777778  ,0.4444444444444444444  ,
     1                               0.2777777777777777778  /
        data wt9/0.0616938806304841571  ,0.108384229110206161  ,
     1           0.0398463603260281088  ,0.175209035316976464  ,
     2           0.229732989232610220  ,0.175209035316976464  ,
     3           0.0398463603260281088  ,0.108384229110206161  ,
     4           0.0616938806304841571  /
c nstack is the number of different intervals into which the
c integration region is currently divided. The number of regions can
c grow if more accuracy is needed by dividing the right-most region
c into three regions. The number of regions shrinks when the integral
c over the right-most region is accurate enough, in which case that
c integral is added to the total (stored in grater) and the region
c is removed from consideration (and a new region is the right-most).
        nstack=nsing+1
        maxns = nstack
        error=0.
        grater=0.
c The array xleft stores the boundary points of the regions.
c The singular points just govern the initial placement of the regions.
        xleft(1)=xmin
        xleft(nsing+2)=xmax
        if(nsing.gt.0) then
          do 9 j=1,nsing
9           xleft(j+1)=xsing(j)
        endif
c For each region, calculate the function and store at three selected points.
        do 1 jj=1,nstack
          del=xleft(jj+1)-xleft(jj)
c         print*, 'fn call j= ,'
          do 1 j=1,3
c         print*, 'fn call in grater j= ',j
1           fval(j,jj)=fn(xleft(jj)+del*dx(j))
c         print*, 'output of fn call, fval(j,jj)',fval(j,jj)
        numcal = nstack * 3
6       continue
          if(nstack+3.ge.mx) then
            write(*,*) 'TOO MANY REGIONS'
            stop 0006
          endif
c Divide the rightmost region into three subregions.
          del=xleft(nstack+1)-xleft(nstack)
          xleft(nstack+3)=xleft(nstack+1)
          xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
          xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
c The three data points already found for the region become the
c middle data points (number 2 in first index of fval) for each region.
          fval(2,nstack+2)=fval(3,nstack)
          fval(2,nstack+1)=fval(2,nstack)
          fval(2,nstack)=fval(1,nstack)
c Now do the integral over the right-most region in two different ways-
c a three point integral (valu) over each of the three subregions
c and a more accurate nine-point integral (value) over whole region.
c valu is used only for the error estimate.
          icount=0
          value=0.
          valu=0.
          do 3 j=nstack,nstack+2
            del1=xleft(j+1)-xleft(j)
c         print*, 'fn call 2'
            fval(1,j)=fn(xleft(j)+dx(1)*del1)
            fval(3,j)=fn(xleft(j)+dx(3)*del1)
c         print*, 'fn call 2'
            numcal = numcal + 2
            do 5 k=1,3
              icount=icount+1
              value=value+wt9(icount)*fval(k,j)*del
5             valu=valu+fval(k,j)*wt(k)*del1
3         continue
          dif=abs(value-valu)
c If the following condition is true, add in this integral to the total,
c and reduce the number of regions under consideration.
          frac = del / (xmax - xmin)
          atsing = .false.
          if(frac .le. 1.0e-8) atsing = .true.
          if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.
     1       (atsing .and.
     2     (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
c The following commented out line is Teeter's old error criterion.
c          if(dif.le.abr.or.dif.le.rlr*abs(value))then
            grater=grater+value
            error=error+abs(dif)
            nstack=nstack-1
c If no more regions, we are done.
            if(nstack.le.0) return
          else
c If the integration is insufficiently accurate, make each of the
c three subregions of the right-most region into regions.
c On next pass the right-most of these is the new current region.
            nstack=nstack+2
            maxns = max(maxns,nstack)
          endif
        go to 6
        end

c**********************************************************************
c   This is Steve White's rewrite of Mike Teter's integration routine.
c   Modified by J. Rehr for complex integration.
c   The following is a listing of the arguments in the initial function
c   statement:
c      fn    -- routine requires external function statement in MAIN
c      xmin  -- lower limit
c      xmax  -- upper limit
c      abr   -- absolute tolerable error
c      rlr   -- relative tolerable error
c      nsing -- number of singularities or regions requiring
c                   special attention
c      xsing -- array of locations of singularities or endpoints
c                   of special regions
c      error -- output for routine error messages
c      numcal-- the number of times fn was called
c      maxns -- the maximum number of regions being considered simultaneously
c       function grater(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
c       fn declared double precision
c       double precision function grater(fn,xmin,xmax,abr,rlr,
c       fn declared complex*16
c      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1

       complex*16 function graterc(fn,xmin,xmax,abr,rlr,
     1 nsing,xsing,error,numcal,maxns)

       implicit double precision (a-h,o-z)
c      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1
       parameter (mx=1500)
       dimension xleft(mx),fval(3,mx),dx(3),wt(3)
       dimension wt9(9), xsing(20)
       external fn
        logical atsing
        save dx,wt,wt9
        data dx/0.1127016653792583  ,0.5  ,0.8872983346207417  /
        data wt/0.277777777777777778  ,0.4444444444444444444  ,
     1                               0.2777777777777777778  /
        data wt9/0.0616938806304841571  ,0.108384229110206161  ,
     1           0.0398463603260281088  ,0.175209035316976464  ,
     2           0.229732989232610220  ,0.175209035316976464  ,
     3           0.0398463603260281088  ,0.108384229110206161  ,
     4           0.0616938806304841571  /
c nstack is the number of different intervals into which the
c integration region is currently divided. The number of regions can
c grow if more accuracy is needed by dividing the right-most region
c into three regions. The number of regions shrinks when the integral
c over the right-most region is accurate enough, in which case that
c integral is added to the total (stored in graterc) and the region
c is removed from consideration (and a new region is the right-most).
        nstack=nsing+1
        maxns = nstack
        error=0.
        graterc=0.
c The array xleft stores the boundary points of the regions.
c The singular points just govern the initial placement of the regions.
        xleft(1)=xmin
        xleft(nsing+2)=xmax
        if(nsing.gt.0) then
          do 9 j=1,nsing
9           xleft(j+1)=xsing(j)
        endif
c For each region, calculate the function and store at three selected points.
        do 1 jj=1,nstack
          del=xleft(jj+1)-xleft(jj)
c         print*, 'fn call j= ,'
          do 1 j=1,3
c         print*, 'fn call in graterc j= ',j
1           fval(j,jj)=fn(xleft(jj)+del*dx(j))
c         print*, 'output of fn call, fval(j,jj)',fval(j,jj)
        numcal = nstack * 3
6       continue
          if(nstack+3.ge.mx) then
            write(*,*) 'TOO MANY REGIONS'
            stop 0006
          endif
c Divide the rightmost region into three subregions.
          del=xleft(nstack+1)-xleft(nstack)
          xleft(nstack+3)=xleft(nstack+1)
          xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
          xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
c The three data points already found for the region become the
c middle data points (number 2 in first index of fval) for each region.
          fval(2,nstack+2)=fval(3,nstack)
          fval(2,nstack+1)=fval(2,nstack)
          fval(2,nstack)=fval(1,nstack)
c Now do the integral over the right-most region in two different ways-
c a three point integral (valu) over each of the three subregions
c and a more accurate nine-point integral (value) over whole region.
c valu is used only for the error estimate.
          icount=0
          value=0.
          valu=0.
          do 3 j=nstack,nstack+2
            del1=xleft(j+1)-xleft(j)
c         print*, 'fn call 2'
            fval(1,j)=fn(xleft(j)+dx(1)*del1)
            fval(3,j)=fn(xleft(j)+dx(3)*del1)
c         print*, 'fn call 2'
            numcal = numcal + 2
            do 5 k=1,3
              icount=icount+1
              value=value+wt9(icount)*fval(k,j)*del
5             valu=valu+fval(k,j)*wt(k)*del1
3         continue
          dif=abs(value-valu)
c If the following condition is true, add in this integral to the total,
c and reduce the number of regions under consideration.
          frac = del / (xmax - xmin)
          atsing = .false.
          if(frac .le. 1.0e-8) atsing = .true.
          if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.
     1       (atsing .and.
     2     (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
c The following commented out line is Teeter's old error criterion.
c          if(dif.le.abr.or.dif.le.rlr*abs(value))then
            graterc=graterc+value
            error=error+abs(dif)
            nstack=nstack-1
c If no more regions, we are done.
            if(nstack.le.0) return
          else
c If the integration is insufficiently accurate, make each of the
c three subregions of the right-most region into regions.
c On next pass the right-most of these is the new current region.
            nstack=nstack+2
            maxns = max(maxns,nstack)
          endif
        go to 6
        end

