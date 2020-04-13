

 subroutine apply_cifsymop(symname,x,xn)
! Given coordinates x(1:3) and a symmetry operation in .cif notation (symname),
! generate the coordinates xn by applying symname to x.
! The resulting numbers are 0 <= xn < 1

! The syntax for a symmetry operation is:
! 'f1, f2, f3'
! where f is a sum of several terms of the form
! "x", "-z" ...
! "1/2", "2/3" ...
! "0.5", "0", "0.6667" ...
!
! Some examples of valid syntax:
! 'x, y, z'
! '-x+y, y, 0.5+z'
! '0.66667-x+y, 0.33333+y, 0.83333+z'

   integer irec
   character symname*(*)
   real*8    x(3),xn(3)
   integer n,i,j,lp
   character*10 part(3)
   real*8 a,b
   logical, external :: try_read

   n=len_trim(symname)
   do jj=1,n
        if(symname(jj:jj).eq.'X')symname(jj:jj)='x'
        if(symname(jj:jj).eq.'Y')symname(jj:jj)='y'
        if(symname(jj:jj).eq.'Z')symname(jj:jj)='z'
   enddo

   i=index(symname,',',.false.)
   j=index(symname,',',.true.)
   part(1)=adjustl(symname(1:i-1))
   part(2)=adjustl(symname(i+1:j-1))
   part(3)=adjustl(symname(j+1:n))


   do i=1,3
      ! i loops over the new coordinates x, y, z

      xn(i)=0.0d0  
      lp=0
      n=len_trim(part(i))

      ! First find chunks of the form "x", "-y" etc.
      j=index(part(i),'-x')
      if (j.ne.0) then
         xn(i)=xn(i)-x(1)
         lp=lp+2 
      else
         j=index(part(i),'x')
         if (j.ne.0) then
            xn(i)=xn(i)+x(1)
            lp=lp+1 
            if (j.ne.1) lp=lp+1
         endif
      endif
      j=index(part(i),'-y')
      if (j.ne.0) then
         xn(i)=xn(i)-x(2) 
         lp=lp+2 
      else
         j=index(part(i),'y')
         if (j.ne.0) then
            xn(i)=xn(i)+x(2)
            lp=lp+1
            if (j.ne.1) lp=lp+1
         endif
      endif
      j=index(part(i),'-z')
      if (j.ne.0) then
         xn(i)=xn(i)-x(3) 
         lp=lp+2 
      else
         j=index(part(i),'z')
         if (j.ne.0) then
            xn(i)=xn(i)+x(3)
            lp=lp+1
            if (j.ne.1) lp=lp+1
         endif
      endif
      ! Next find chunks of the form "1/3", "-2/3" etc.
      j=index(part(i),'/')
      if (j.ne.0) then
         if (j.eq.1.or.j.eq.n) call stop_syntax_error(symname)
         read(part(i)(j-1:j-1),*) a
         read(part(i)(j+1:j+1),*) b
         if (j.eq.2) then
            xn(i)=xn(i)+a/b                  
            lp=lp+3
         else 
            if (part(i)(j-2:j-2).eq.'-') then
               xn(i)=xn(i)-a/b                  
               lp=lp+4
            else
               xn(i)=xn(i)+a/b                  
               lp=lp+4
            endif
         endif
      endif
      ! Finally, find chunks of the form "0.5", "0.6667" etc.
      j=index(part(i),'.')
      if(j.ne.0) then
         ! We have a decimal number somewhere.  Let's try to read it using format statements.
         if(.not.try_read(part(i),a)) then
            call stop_syntax_error(symname)
         endif
         !Correct "a" as a fraction - cif's are often less accurate than theoretical symmetry routines
         do j=-24,24
            if (dabs(a -dble(j)/dble(12)).lt.0.0002) a =dble(j)/dble(12)  ! so 0.6667 becomes 0.66666666666...
         enddo
         do j=-16,16
            if (dabs(a -dble(j)/dble(8)).lt.0.0002) a =dble(j)/dble(8)  ! so 0.6667 becomes 0.66666666666...
         enddo
         !write(*,*) 'HA: ',a, ' --- ', part(i)

         xn(i)=xn(i)+a
      else
         ! This sanity check only works if we didn't get a decimal number.
         if (lp.ne.len_trim(part(i))) call stop_syntax_error(symname)
      endif

   enddo

   do i=1,3
      if (xn(i).lt.0.0d0) xn(i)=xn(i)+1.0d0
      if (xn(i).gt.1.0d0) xn(i)=xn(i)-1.0d0
   enddo
   
   do i=1,3
      if (xn(i).lt.0.0d0) xn(i)=xn(i)+1.0d0
      if (xn(i).gt.1.0d0) xn(i)=xn(i)-1.0d0
   enddo
   
   do i=1,3
      if (xn(i).lt.0.0d0) xn(i)=xn(i)+1.0d0
      if (xn(i).gt.1.0d0) xn(i)=xn(i)-1.0d0
   enddo
   
   return
 end subroutine apply_cifsymop
 


 logical function try_read(s,a)
    ! Try to extract a number "a" from a string "s" by applying a series of format statements
    implicit none
    real*8, intent(out) :: a
    character*(*), intent(in) :: s
    character*10 s1,s2 ! forgettable dummies
    character*100 ss ! copy of "s" modified for easy reading
    integer i,j

    ! Copy array into a new array for reading
    ! Add a space before "-" and around "+"
    i=1
    do j=1,len_trim(s)
       if(s(j:j).eq.'-') then
          ss(i:i+1)=' -'
          i=i+2
       elseif(s(j:j).eq.'+') then
          ss(i:i+1)=' + '
          i=i+3
       else
          ss(i:i)=s(j:j)
          i=i+1
       endif
    enddo

    a=0.d0
!write(*,*) 'try_read: ',ss
    read(ss,*, err=901, end=901) a, s1
    goto 904
    901 continue   ! if this one fails, try the next:
    read(ss,*,err=902, end=902) s1, a
    goto 904
    902 continue  ! if this one fails, try the next:
    read(ss,*,err=903, end=903) s1, a, s2
    goto 904

    904 continue  ! succesful read:
    try_read = .true.
    return

    903 continue  ! all reads were unsuccessful :
    try_read=.false.
    return
 end function try_read

 
 subroutine stop_syntax_error(symname)
   implicit none
   character*(*) symname
      call wlog('wrong syntax in _symmetry_equiv_pos_as_xyz: '//symname)
      stop ' choked on symmetry operations in cif file'
 end subroutine stop_syntax_error



 subroutine gen_equiv(nsym,xn,lattyp,indequiv,mult)

   integer nsym
   real*8  xn(3,nsym)
   character lattyp*(*)
   integer   indequiv(nsym),mult
   integer i,j,k
   logical thesame    

   mult=0
   do i=1,nsym
      do j=1,mult
         k=indequiv(j)
         if (thesame(lattyp,xn(1:3,k),xn(1:3,i))) goto 1
      enddo
      mult=mult+1
      indequiv(mult)=i
1     continue
   enddo

 end subroutine gen_equiv







 logical function thesame(lattyp,x1,x2)

    character lattyp*(*)
    real*8    x1(3),x2(3)
    real*8 tv(3,1000),dx(3),dxl,small,small2
    integer ind,i,j,k


    small=1.0d-5
    small2=small/2.0d0

    ind=0
    do i=-1,1
       do j=-1,1
          do k=-1,1
             ind=ind+1
             tv(1,ind)=i
             tv(2,ind)=j
             tv(3,ind)=k
          enddo 
        enddo
     enddo

     if (lattyp(1:1).eq.'B' .or. lattyp(1:1).eq.'I') then

        do i=-1,1,2
           do j=-1,1,2
              do k=-1,1,2
                 ind=ind+1
                 tv(1,ind)=0.5d0*i
                 tv(2,ind)=0.5d0*j
                 tv(3,ind)=0.5d0*k
              enddo
           enddo
        enddo

     else if (lattyp(1:1).eq.'F') then

        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.0d0
           enddo
        enddo
        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.0d0
              tv(3,ind)=0.5d0*j
           enddo
        enddo
        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.0d0
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.5d0*j
           enddo
        enddo

     else if (lattyp(1:3).eq.'CXY') then
          
        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.0d0
           enddo
        enddo

     else if (lattyp(1:3).eq.'CXZ') then

        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.0d0
              tv(3,ind)=0.5d0*j
           enddo
        enddo

     else if (lattyp(1:3).eq.'CYZ') then

        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.0d0
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.5d0*j
           enddo
        enddo

     endif

     thesame=.false.
     do i=1,ind
        dx(1:3)=x1(1:3)-x2(1:3)-tv(1:3,i)
        dx=dx+10.0d0
        dx=mod(dx+small2,1.0d0)-small2
        if (abs(dx(1)-1.0d0).lt.small2) dx(1)=0.0d0
        if (abs(dx(2)-1.0d0).lt.small2) dx(2)=0.0d0
        if (abs(dx(3)-1.0d0).lt.small2) dx(3)=0.0d0

        dxl=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
        if (dxl.lt.small) then
           thesame=.true.
           return
        endif
     enddo

  end function thesame
