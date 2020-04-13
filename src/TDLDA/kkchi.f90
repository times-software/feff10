!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: kkchi.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine kkchi ( emu, edge, refsh, kinitm, kfinm, matsize,      &
     &                  em, ne, ne1, chi0im, chi0r)
   
      use dimsmod, only: nex 
	  use constants     
      implicit double precision (a-h, o-z)

      parameter (maxsize = 78)
      dimension chi0im(nex,maxsize,maxsize)
      parameter (nexfk = 2000)
      dimension chi0fk(nexfk,maxsize,maxsize)
      dimension chi0r(nex,maxsize,maxsize)
      dimension kinitm(maxsize), kfinm(maxsize)
      dimension em(nex), emfk(nexfk)
      dimension refsh(maxsize)

!     E1 transitions, mult = 0

!     do 700 iemain = 1, ne
      do 700 iemain = 1, ne1
        w = em(iemain)   

!       make "fake" grid - a constant step energy grid around w

        eleft = em(1)
        eright = em(ne1) + 2.0/hart
        nefk = 2000
        step = (eright - eleft)/(nefk - 1)
        nel = int((w - eleft)/step)
!       w is between poits nel and nel+1
        delta = step/2.0 - ( (w - eleft) - step*nel )
        emfk(1) = eleft - delta         
        if (delta .gt. 0) emfk(1) = emfk(1) + step 
        do 11 ie = 2, nefk
          emfk(ie) = emfk(ie-1) + step
   11   continue       


        do 20 im = 1, matsize
        do 20 imp = 1, matsize

          chi0r(iemain,im,imp) = 0.0

          kinit = kinitm(im)
          kfin = kfinm(im)
          kdif = kfin - kinit
!         selection rules (buggy?)
!          ks = 1
!          if ( abs(kdif) .ne. ks .and. kfin .ne. - kinit ) goto 20    

!         find chi0im on the fake grid. Use simple linear interpolation

          do 16 ief = 1, nefk
           do 15 ie = 1, ne1
            if (emfk(ief) .eq. em(ie)) then
              chi0fk(ief,im,imp) = chi0im(ie,im,imp)
            else
              del1 = emfk(ief) - em(ie)
              del2 = emfk(ief) - em(ie+1)   
              if (del1*del2 .lt. 0.0 .or.ie.eq.ne1-1) then
                t1 = chi0im(ie,im,imp)*(em(ie+1) - emfk(ief))
                t2 = chi0im(ie+1,im,imp)*(emfk(ief) - em(ie))
                chi0fk(ief,im,imp) = (t1+t2)/(em(ie+1)-em(ie))
                goto 16
              endif
            endif
   15      continue
   16     continue

!         make integration
!         trapeozoid method 

!         here we use fake grid

          xint = 0.0
          do 920 ie = 1, (nefk-1)
            pint = 0.0
            e1 = emfk(ie)
            e2 = emfk(ie + 1)
!            w = wm(iw)
!c           integration starts at Ef (that is, em = edge)
            if (e1 .lt. (edge - refsh(im))) go to 920

!           add second pole contribution, if is.ne.0
            is = 1
!           is = 0
            nc = 0
            es = emu

            if ( e2 .gt. w .and. e1 .lt. w ) then
!             chi0(e) at E=w
              a1 = chi0fk(ie,im,imp) * (e2 - w)
!             if (is.ne.0) a1 = a1* 2 *(w +es)**2 /(2*es+w+e2)/(es+e2)
!             if (is.ne.0) a1 = a1* 2 /(2*es+w+e2)*(es+e2)
             if (is.ne.0) a1=a1*2/(w+es)**nc/(2*es+w+e2)*(es+e2)**(nc+1)
              a2 = chi0fk(ie+1,im,imp) * (w - e1)
!             if (is.ne.0) a2 = a2* 2 *(w +es)**2 /(2*es+w+e1)/(es+e1)
!             if (is.ne.0) a2 = a2* 2 /(2*es+w+e1)*(es+e1)
             if (is.ne.0) a2=a2*2/(w+es)**nc/(2*es+w+e1)*(es+e1)**(nc+1)
              xchiint = (a1 + a2) / (e2 - e1)
!             first term
              arg = (w - e1) / (e2 - w)
              pint = - xchiint * log(arg) 
              pint = pint + (chi0fk(ie+1,im,imp) - chi0fk(ie,im,imp))   
            else
              a1 = chi0fk(ie+1,im,imp) / (e2 - w)
!             if (is.ne.0) a1 = a1* 2 *(w +es)**2 /(2*es+w+e2)/(es+e2)
!             if (is.ne.0) a1 = a1* 2 /(2*es+w+e2)*(es+e2)
             if (is.ne.0) a1=a1*2/(w+es)**nc/(2*es+w+e2)*(es+e2)**(nc+1)
              a2 = chi0fk(ie,im,imp) / (e1 - w)       
!             if (is.ne.0) a2 = a2* 2 *(w +es)**2 /(2*es+w+e1)/(es+e1)
!             if (is.ne.0) a2 = a2* 2 /(2*es+w+e1)*(es+e1)
             if (is.ne.0) a2=a2*2/(w+es)**nc/(2*es+w+e1)*(es+e1)**(nc+1)
              pint = 0.5 * (e2 - e1) * (a2 + a1)
            endif
!           extra factor 1/pi
            pint =  pint / pi
 
            xint = xint + pint 
  920     continue  

! test constant shift
          shift = 0.0
          chi0r(iemain,im,imp) =  xint + shift

!         test - put real part of chi to 0
!         chi0r(iemain,im,imp) = 0

   20   continue
!     end of loop over im
   
  700 continue
!     end of energy loop            
       
      return 
      end
