module planchon_2001

implicit none
!private
!public :: remove_pits,next_cell
integer :: npos = 8
integer :: pos,k,l,positions(8,2)
real :: undef = -9999
real :: wmax = 1000000
real :: eps
integer :: r0(8),c0(8),dr(8),dc(8),fr(8),fc(8)

contains

subroutine remove_pits(z,w,res,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: z(nx,ny),res
 real,intent(out) :: w(nx,ny)
 integer :: b(nx,ny)
 real :: minslope = 0.0001
 
 !Define minimum elevation change
 eps = minslope*res

 !Initialize others 
 r0 = (/1,ny,1,ny,1,ny,1,ny/)
 c0 = (/1,nx,nx,1,nx,1,1,nx/)
 dr = (/0,0,1,-1,0,0,1,-1/)
 dc = (/1,-1,0,0,-1,1,0,0/)
 fr = (/1,-1,-ny,ny,1,-1,-ny,ny/)!careful
 fc = (/-nx,nx,-1,1,nx,-nx,1,-1/)!careful

 !Construct positions array
 npos = 8
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !Stage 1
 call stage1(z,w,b,nx,ny)

 !Stage 2 (section 1)
 call stage2section1(z,w,b,nx,ny)

 !Stage 2 (section 2)
 call stage2section2(z,w,b,nx,ny)
 
end subroutine

subroutine next_cell(i,r,c,nx,ny,ncflag)

 implicit none
 integer :: r,c,i,nx,ny
 logical :: ncflag
 ncflag = .True.

 r = r + dr(i)
 c = c + dc(i)
 if ((r .lt. 1) .or. (c .lt. 1) .or. (r .gt. ny) .or. (c .gt. nx))then
  r = r + fr(i)
  c = c + fc(i)
  if ((r .lt. 1) .or. (c .lt. 1) .or. (r .gt. ny) .or. (c .gt. nx))then
   ncflag = .False.
  endif
 endif

end subroutine

subroutine stage2section2(z,w,b,nx,ny)

 implicit none
 integer :: nx,ny,x,y,it,scan,xnew,ynew
 real ::w(nx,ny),z(nx,ny)
 integer :: b(nx,ny),r,c,i
 logical :: something_done,ncflag

 do it=1,1000
  do scan=1,8
   r = r0(scan)
   c = c0(scan)
   something_done = .False.
   ncflag = .True.
   do while (ncflag .eqv. .True.)
    if ((z(c,r) .ne. undef) .and. (w(c,r) .gt. z(c,r)))then
     do i=1,npos
      xnew = c+positions(i,1)
      ynew = r+positions(i,2)
      !in grid?
      if ((xnew .le. 0) .or. (xnew .gt. nx) .or. (ynew .le. 0) .or. &
        (ynew .gt. ny) .or. (w(xnew,ynew) .eq. undef))cycle
      !operation 1
      if(z(c,r) .ge. (w(xnew,ynew) + eps))then
       w(c,r) = z(c,r)
       something_done = .True.
       call dry_upward_cell(z,w,c,r,nx,ny)
       exit
      endif
      !operation 2
      if(w(c,r) .gt. (w(xnew,ynew) + eps))then
       w(c,r) = w(xnew,ynew) + eps
       something_done = .True.
      endif
     enddo
    endif
    call next_cell(scan,r,c,nx,ny,ncflag)
   enddo
   if (something_done .eqv. .False.)exit
  enddo
  if (something_done .eqv. .False.)exit
 enddo

end subroutine

subroutine stage2section1(z,w,b,nx,ny)

 implicit none
 integer :: nx,ny,x,y
 real ::w(nx,ny),z(nx,ny)
 integer :: b(nx,ny)

 do x = 1,nx
  do y = 1,ny
   if (b(x,y) .eq. 1) then
    call dry_upward_cell(z,w,x,y,nx,ny)
   endif
  enddo
 enddo


end subroutine

recursive subroutine dry_upward_cell(z,w,x,y,nx,ny)

 implicit none
 integer :: i,x,y,xnew,ynew,nx,ny
 real :: w(nx,ny),z(nx,ny)
 
 do i=1,npos
  xnew = x + positions(i,1)
  ynew = y + positions(i,2)
  if ((xnew .le. 0) .or. (xnew .gt. nx) .or. (ynew .le. 0) .or. &
        (ynew .gt. ny) .or. (w(xnew,ynew) .eq. undef))cycle
  !if (w(xnew,ynew) .ne. wmax)cycle
  if ((z(xnew,ynew) .ge. (w(x,y) + eps)) .and. (w(xnew,ynew) .eq. wmax)) then
   w(xnew,ynew) = z(xnew,ynew)
   call dry_upward_cell(z,w,xnew,ynew,nx,ny)
  endif
 enddo

end subroutine

subroutine stage1(z,w,b,nx,ny)

 implicit none
 integer :: nx,ny,x,y,i,xnew,ynew,b(nx,ny)
 real :: w(nx,ny),z(nx,ny)
 w = z
 b = 0

 do x = 1,nx
  do y = 1,ny
   if (w(x,y) .eq. undef)cycle
   do i = 1,npos
    xnew = x + positions(i,1)
    ynew = y + positions(i,2)
    if ((xnew .le. 0) .or. (xnew .gt. nx) .or. (ynew .le. 0) .or. &
        (ynew .gt. ny) .or. (w(xnew,ynew) .eq. undef))then
     b(x,y) = 1
     w(x,y) = z(x,y)
     exit
    else
     w(x,y) = wmax
    endif
   enddo
  enddo
 enddo

end subroutine   

end module
