      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc) !KJ ,func)
	use fitting, x=>fitx,y=>fity,n=>fitmeshsize,m=>fitorder  !KJ
      implicit none !KJ
      DOUBLE PRECISION ax,bx,cx,fa,fb,fc,GOLD,GLIMIT,TINY  !KJ ,func
!KJ      EXTERNAL func
      PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
!	Given a function func, and given distinct initial points ax and bx,
!	this routine searches in the downhill direction (defined by the function
!	as evaluated at the initial points) and returns new points ax, bx, cx
!	that bracket a minimum of the function.  Also returned are the function
!	values at the three points, fa, fb, and fc.
!	Parameters : GOLD is the default ration by which successive intervals are
!	magnified ; GLIMIT is the maximum magnification allowed for a parabolic-fit
!	step.
!!KJ 2/2007 : modified to use terp instead of func, which has been removed from argument list.

      DOUBLE PRECISION dum,fu,q,r,u,ulim
      call terp (x,y,n,m,ax,fa) !KJ fa=func(ax)
      call terp (x,y,n,m,bx,fb) !KJ fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      call terp (x,y,n,m,cx,fc) !KJ fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.d0)then
          call terp (x,y,n,m,u,fu) !KJ fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          call terp (x,y,n,m,u,fu) !KJ fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.d0)then
          call terp (x,y,n,m,u,fu) !KJ fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
             call terp (x,y,n,m,u,fu) !KJ fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
          call terp (x,y,n,m,u,fu) !KJ fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          call terp (x,y,n,m,u,fu) !KJ fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
!        write(6,*) 'nmbrak',ax,fa,bx,fb,cx,fc
        goto 1
      endif
      return
      END

