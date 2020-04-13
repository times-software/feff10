!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: yzkrdc.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine yzkrdc (i,k,flps,ps,qs,aps,aqs,p2, norb)
!       * calculate  function yk *
! yk = r * integral of f(s)*uk(r,s)
! uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)
! j=norb for photoelectron
! f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)
! f(s) is constructed by the calling programm  if i < or =0
! in the last case a function f (lies in the block dg) is supposedly
! tabulated untill point dr(j), and its' devlopment coefficients
! at the origin are in ag and the power in r of the first term is k+2

! the output functions yk and zk are in the blocks dp and dg.
! at the origin  yk = cte * r**(k+1) - developement limit,
! cte lies in ap(1) and development coefficients in ag.
!        this programm uses aprdec and yzktec
      use dimsmod, only: nrptx
      implicit double precision (a-h,o-z)

      complex*16 aprdec,p2, dyzk
!     complex*16 a1,a2,b1,b2,coni
!     complex*16 xck, temp, ck, phx
      parameter (coni=(0.d0,1.d0))
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10)
      common/dff/cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),      &
     &             fl(30), fix(30), ibgp
      complex*16 dg,ag,dp,ap,bidcom, chg(10)
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),            &
     &   bidcom(3*nrptx+30)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),            &
     &   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgi(10),bpi(10)
!#mn
       external aprdec
 
!     construction of the function f
      do  5 l= 1,ibgp
        bgi(l) = bg(l,i)
  5     bpi(l) = bp(l,i)
      id=min(nmax(i),np)
      ap(1)=fl(i)+flps
      do 11 l=1,id
 11   dg(l)=cg(l,i)*ps(l)+cp(l,i)*qs(l)
      do 12 l = id+1,idim
 12    dg(l) = 0.0d0
      do 21 l=1,ndor
 21   ag(l) = aprdec(aps,bgi,l) + aprdec(aqs,bpi,l)

      dyzk = 0
!     if (id .ge. nmax(norb)) then
!        id = nmax(norb)-1
!        ck0 = log(cg(id,i)/cg(id+1,i))  / (dr(id+1)-dr(id))
!        ck = sqrt(2*p2)
!        xck = ck/cl
!        xck = -xck/(1+sqrt(1+xck**2))
!        temp = -ps(id+1) / qs(id+1) *xck
!        xx = dble (temp)
!        yy = dimag(temp)
!        if (xx .ne. 0)  then
!            alph = (1 - xx**2 - yy**2)
!            alph = sqrt(alph**2 + 4*xx**2) - alph
!            alph = alph / (2 * xx)
!            alph = atan (alph)
!        else
!            alph = 0
!        endif
!        beta = (xx**2 + (yy+1)**2) / (xx**2 + (yy-1)**2)
!        beta = log(beta) / 4

!        phx = dcmplx (alph, beta)
!        a1 =   ps(id+1) / sin(phx)
!        a2 = - qs(id+1) / cos(phx)
!        xck=ck*dr(id+1)
!        phx = phx -xck
!        a1 = a1*cg(id+1,i)/2/coni
!        a2 = a2*cp(id+1,i)/2
!        b1=exp(coni*phx) * (a1 - a2)
!        b2=exp(-coni*phx) * (-a1 - a2)
!        xck = (ck0 - coni*ck)*dr(id+1)
!        n = k +1
!        dyzk = dyzk + b1*exp(-xck)/xck
!        dyzk = dyzk + b1*expint(n,xck)
!        xck = (ck0 + coni*ck)*dr(id+1)
!        dyzk = dyzk + b2*exp(-xck)/xck
!        dyzk = dyzk + b2*expint(n,xck)
!        dyzk = dyzk*dr(id+1)
!     endif

      call yzktec (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim, dyzk)
      return
      end

!     complex*16 function expint(n,x)
!     implicit double precision (a-h,o-z)
!     integer n, maxit
!     complex*16 x, b, c, d, h, del, fact, zero
!     parameter (zero=(0.d0,0.d0))
!     parameter (maxit=100, eps=1.d-7, fpmin=1.d-30, euler=.5772156649)

!     nm1 = n - 1
!     if (n.lt.0 .or. (dble(x).lt.0.d0 .and. dimag(x).eq.0.d0) .or.
!    1     (x.eq.zero .and. (n.eq.0.or.n.eq.1))) then
!        call par_stop('Bad arguments in expint')
!     elseif (n.eq.0) then
!        expint = exp(-x) / x
!     elseif (x.eq.0) then
!        expint = 1.d0 /nm1
!     elseif (dble(x).gt.1) then
!        b = x + n
!        c = 1/fpmin
!        d = 1/b
!        h = d
!        do 10 i=1,maxit
!           a = -i*(nm1+i)
!           b = b + 2
!           d = 1 / (a*d+b)
!           c = b + a/c
!           del = c*d
!           h = h*del
!           if (abs(del-1) .lt. eps) then
!              expint = h * exp(-x)
!              return
!           endif
! 10     continue
!        call par_stop(' continued fraction failed in expint')
!     else
!        if (nm1.ne.0) then
!           expint = 1/nm1
!        else
!           expint = -log(x) - euler
!        endif
!        fact = 1
!        do 30 i=1,maxit
!           fact = - fact *x / i
!           if (i.ne.nm1) then
!              del = - fact / (i-nm1)
!           else
!              psi = - euler
!              do 20 ii=1,nm1
!                 psi = psi + 1.d0 / ii
! 20           continue
!              del = fact*(-log(x)+psi)
!           endif
!           expint = expint + del
!           if (abs(del).lt.abs(expint)*eps) return
! 30     continue
!        call par_stop('series failed in expint')
!     endif
!     return
!     end
