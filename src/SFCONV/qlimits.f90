!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: qlimits.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine qlimits(w,pk,omp,Ap,qh,nq,q1,q2,q3)
! finds the limiting q values of the inequalities
! omega(q)+(q-k)^2/2-omega < 0 
! omega(q)+(q+k)^2/2-omega > 0 and
! omega(q)-(q-k)^2/2+omega > 0
! for omega(q)^2 = omp^2+Ap q^2+q^4/4
! qh is the upper allowable value for q (usually set by omega(q)+E_f-omega < 0
! for the first two inequalities.
! input: energy w in atomic hartrees.
!        momentum pk in atomic units.
!        resonance frequency omp in hartrees.
!        dispersion parameter Ap.
!        upper limit qh.
! output: number of limiting q values (either 1 or 3).
!         limiting q values q1, q2,and q3
!         q1 and q2 bracket a region of allowed q values given by
!         the first two inequalities listed above,
!         q3 is the upper bound of the third inequality
      implicit none
      integer i,j,nar,nq
      double precision w,pk,omp,Ap,qh,q1,q2,q3,a,b,c,d
      double precision dev1,dev2,dev3,qa1,qa2,qa3,wdisp
      complex*16 aa1,aa2,aa3
      external wdisp
!      write(6,*) 'qlimits'
! find q for which omega(q)+(q-k)^2/2-omega = 0
      a=pk
      b=w+Ap-3.d0*pk**2/2.d0
      c=pk**3-2.d0*w*pk
      d=omp**2-w**2+w*pk**2-pk**4/4.d0
      call croots(a,b,c,d,aa1,aa2,aa3,nar)
! One of these roots is the solution of
! omega(q)-(q-k)^2/2+omega = 0.  This is q3.  The other two roots, if present,
! solve omega(q)+(q-k)^2/2-omega = 0
! It must hold that q>=0.  Fortunately, the symmetry of the inequalities
! means that solving omega(q)+(q+k)^2/2-omega = 0 gives roots which are the 
! negatives of the roots from the cubic we just solved.  It suffices to
! take the absolute value of the roots (if they are real) to give our
! limiting values of q.
      if (nar.eq.3) then
        qa1=dble(aa1)
        qa2=dble(aa2)
        qa3=dble(aa3)
        dev1=dabs(wdisp(qa1)+(qa1-pk)**2/2.d0-w)
        dev2=dabs(wdisp(qa2)+(qa2-pk)**2/2.d0-w)
        dev3=dabs(wdisp(qa3)+(qa3-pk)**2/2.d0-w)
        if (dev1.gt.dev2.and.dev1.gt.dev3) then
          q1=min(dabs(qa2),dabs(qa3))
          q2=max(dabs(qa2),dabs(qa3))
          q3=dabs(qa1)
        elseif(dev2.gt.dev3) then
          q1=min(dabs(qa1),dabs(qa3))
          q2=max(dabs(qa1),dabs(qa3))
          q3=dabs(qa2)
        else
          q1=min(dabs(qa1),dabs(qa2))
          q2=max(dabs(qa1),dabs(qa2))
          q3=dabs(qa3)
        endif
        q1=min(q1,qh)
        q2=min(q2,qh)
        nq=3
      else
! The equation omega(q)+(q-k)^2/2-omega = 0 has no real solutions.  The
! one real root of the cubic is q3. 
        q1=0.d0
        q2=0.d0
        qa1=dabs(dimag(aa1))
        qa2=dabs(dimag(aa2))
        qa3=dabs(dimag(aa3))
        if (qa1.lt.qa2.and.qa1.lt.qa3) then
          q3=dabs(dble(aa1))
        elseif (qa2.lt.qa3) then
          q3=dabs(dble(aa2))
        else
          q3=dabs(dble(aa3))
        endif
        nq=1
      endif
      return
      end
