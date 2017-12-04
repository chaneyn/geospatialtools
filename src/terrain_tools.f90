subroutine remove_pits_planchon(dem,res,demns,nx,ny)

 use planchon_2001, only: remove_pits
 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: res
 real,intent(in),dimension(nx,ny) :: dem
 real,intent(out),dimension(nx,ny) :: demns

 call remove_pits(dem,demns,res,nx,ny)

end subroutine

subroutine calculate_slope_and_aspect(dem,dx,dy,slope,aspect,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in),dimension(nx,ny) :: dem,dx,dy
 real,intent(out),dimension(nx,ny) :: slope,aspect
 integer :: imin,imax,jmin,jmax,i,j
 real :: undef,dzdx,dzdy
 undef = -9999.0
 slope = dem

 do i=1,nx
  do j=1,ny
   if (slope(i,j) .eq. undef) cycle
   imin = i-1
   imax = i+1
   jmin = j-1
   jmax = j+1
   if (imin .le. 0)imin = 1
   if (imax .gt. nx)imax = nx
   if (jmin .le. 0)jmin = 1
   if (jmax .gt. ny)jmax = ny
   !dz/dy
   dzdy = ((dem(imin,jmin) + 2*dem(i,jmin) + dem(imax,jmin)) &
          - (dem(imin,jmax) + 2*dem(i,jmax) + dem(imax,jmax))) &
          /((dy(imin,jmin) + 2*dy(i,jmin) + dy(imax,jmin)) &
          + (dy(imin,jmax) + 2*dy(i,jmax) + dy(imax,jmax)))
   !dz/dx
   dzdx = ((dem(imin,jmin) + 2*dem(imin,j) + dem(imin,jmax)) &
          - (dem(imax,jmin) + 2*dem(imax,j) + dem(imax,jmax))) &
          /((dx(imin,jmin) + 2*dx(imin,j) + dx(imin,jmax)) &
          + (dx(imax,jmin) + 2*dx(imax,j) + dx(imax,jmax))) 
   !slope
   slope(i,j) = (dzdx**2 + dzdy**2)**0.5
   !aspect
   aspect(i,j) = atan2(dzdy,-dzdx)
   !aspect(i,j) = atan2(dxdy,dzdy)
  enddo
 enddo
 
 end subroutine

subroutine remove_pits(dem,demns,res,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in),dimension(nx,ny) :: dem
 real,intent(out),dimension(nx,ny) :: demns
 real,intent(in) :: res
 integer :: positions(8,2)
 integer :: i,j,k,l,pos,count,p,maxp
 integer :: imin,jmin,imax,jmax,imin0,imax0,jmin0,jmax0
 imin = 1
 imax = nx
 jmin = 1
 jmax = ny
 maxp = 10000

 !Copy the dem
 demns = dem

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !Iterate through all cells and find ones that cannot flow anywhere
 !get flow direction map
 do p=1,100000
  imin0 = imin
  imax0 = imax
  jmin0 = jmin
  jmax0 = jmax
  imin = 99999
  imax = -99999
  jmin = 99999
  jmax = -99999
  count = 0
  do i=imin0,imax0
   do j=jmin0,jmax0
    call check_remove_pit(i,j,demns,positions,nx,ny,8,count,imax,imin,jmin,jmax,res)
   enddo
  enddo
  !print*,p,count
  if (count .eq. 0)then
   !zoom out to check 
   imin = 99999
   imax = -99999
   jmin = 99999
   jmax = -99999
   count = 0
   do i=1,nx
    do j=1,ny
     call check_remove_pit(i,j,demns,positions,nx,ny,8,count,imax,imin,jmin,jmax,res)
    enddo
   enddo
   if (count .eq. 0)exit
  endif
  if (p .eq. maxp)then
   print*,"pit removal did not converge"
   exit
  endif
  
 enddo

end subroutine 

recursive subroutine check_remove_pit(i,j,demns,positions,nx,ny,npos,count,&
                                      imax,imin,jmin,jmax,res)
 
 implicit none
 integer,intent(in) :: i,j,npos,nx,ny
 integer,intent(inout) :: count,imin,imax,jmin,jmax
 real,intent(inout),dimension(nx,ny) :: demns
 integer,intent(in),dimension(npos,2) :: positions
 real,intent(in) :: res
 integer :: inew,jnew,pos,k,l,tmp(1)
 real :: slopes(8),length
 real :: minslope = 0.01

 if (demns(i,j) .eq. -9999)return
 if ((i .eq. 1) .or. (i .eq. nx) .or. (j .eq. 1) .or. (j .eq. ny))return
 slopes = -9999.0
 do pos=1,8
  k = positions(pos,1)
  l = positions(pos,2)
  if ((i+k .lt. 1) .or. (j+l .lt. 1) .or. (i+k .gt. nx) .or. &
      (j+l .gt. ny)) cycle !skip due to on boundary
  if ((k + l .eq. -2) .or. (k + l .eq. 2) .or. (k + l .eq. 0))then
      length = 1.41421356237*res
  else
      length = res
  endif
  slopes(pos) = (demns(i,j) - demns(i+k,j+l))/length
 enddo
 if (maxval(slopes) .le. 0)then
  if (i .gt. imax)imax = i
  if (j .gt. jmax)jmax = j
  if (i .lt. imin)imin = i
  if (j .lt. jmin)jmin = j
  tmp = maxloc(slopes)
  inew = i+positions(tmp(1),1)
  jnew = j+positions(tmp(1),2)
  demns(i,j) = demns(inew,jnew)+minslope*res
  !demns(inew,jnew) = demns(i,j) - minslope*res
  count = count + 1
  call check_remove_pit(inew,jnew,demns,positions,nx,ny,npos,count,imax,imin,jmin,jmax,res)
 endif
    
end subroutine
    
subroutine calculate_d8_acc(dem,mask,res,area,fdir,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in),dimension(nx,ny) :: dem,mask
 real,intent(in) :: res
 real,intent(out),dimension(nx,ny) :: area
 integer,intent(out),dimension(nx,ny,2) :: fdir
 integer,allocatable,dimension(:,:) :: positions
 real,allocatable,dimension(:) :: slopes
 real :: length,undef,demns(nx,ny)
 integer :: catchment(nx,ny)
 integer :: i,j,k,l,pos,tmp(1),npos=8
 allocate(positions(npos,2),slopes(npos))
 undef = -9999.0
 demns = dem

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !get flow direction map
 do i=1,nx
  do j=1,ny
   slopes = -9999.0
   do pos=1,npos
    k = positions(pos,1)
    l = positions(pos,2)
    if ((i+k .lt. 1) .or. (j+l .lt. 1) .or. (i+k .gt. nx) .or. &
        (j+l .gt. ny)) then
      cycle !skip due to on boundary
    endif
    if ((k + l .eq. -2) .or. (k + l .eq. 2) .or. (k + l .eq. 0))then
        length = 1.41421356237*res
    else 
       length = res
    endif
    slopes(pos) = (demns(i,j) - demns(i+k,j+l))/length
   enddo
   if (maxval(slopes) .gt. 0) then
    tmp = maxloc(slopes)
    fdir(i,j,1) = i+positions(tmp(1),1)
    fdir(i,j,2) = j+positions(tmp(1),2)
   else if(minval(slopes) .eq. -9999.0) then
    tmp = minloc(slopes)
    fdir(i,j,1) = i+positions(tmp(1),1)
    fdir(i,j,2) = j+positions(tmp(1),2)
   else
    fdir(i,j,:) = int(undef)
   endif
  enddo
 enddo

 !get the cell count
 catchment(:,:) = 0
 do i=1,nx
  do j=1,ny
   call neighbr_check_d8(i,j,demns,catchment,fdir,positions,nx,ny,npos)
  enddo
 enddo

 !Calculate accumulation area
 area = res**2*catchment

 !Where the mask is 0 set area to undefined
 where (mask .eq. 0)
  area = undef
 endwhere
 where (fdir(:,:,1) .eq. -9999)
  area = undef
 endwhere

end subroutine

subroutine calculate_d8_acc_wipoints(dem,mask,ipoints,res,area,fdir,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in),dimension(nx,ny) :: dem,mask
 integer,intent(in),dimension(nx,ny) :: ipoints
 real,intent(in) :: res
 real,intent(out),dimension(nx,ny) :: area
 integer,intent(out),dimension(nx,ny,2) :: fdir
 integer,allocatable,dimension(:,:) :: positions
 real,allocatable,dimension(:) :: slopes
 real :: length,undef,demns(nx,ny)
 integer :: catchment(nx,ny),pc
 integer :: i,j,k,l,pos,tmp(1),npos=8
 allocate(positions(npos,2),slopes(npos))
 undef = -9999.0
 demns = dem

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !get flow direction map
 do i=1,nx
  do j=1,ny
   slopes = -9999.0
   do pos=1,npos
    k = positions(pos,1)
    l = positions(pos,2)
    if ((i+k .lt. 1) .or. (j+l .lt. 1) .or. (i+k .gt. nx) .or. &
        (j+l .gt. ny)) then
      cycle !skip due to on boundary
    endif
    if ((k + l .eq. -2) .or. (k + l .eq. 2) .or. (k + l .eq. 0))then
        length = 1.41421356237*res
    else 
       length = res
    endif
    slopes(pos) = (demns(i,j) - demns(i+k,j+l))/length
   enddo
   if (maxval(slopes) .gt. 0) then
    tmp = maxloc(slopes)
    fdir(i,j,1) = i+positions(tmp(1),1)
    fdir(i,j,2) = j+positions(tmp(1),2)
   else if(minval(slopes) .eq. -9999.0) then
    tmp = minloc(slopes)
    fdir(i,j,1) = i+positions(tmp(1),1)
    fdir(i,j,2) = j+positions(tmp(1),2)
   else
    fdir(i,j,:) = int(undef)
   endif
  enddo
 enddo

 !get the cell count
 catchment(:,:) = 0
 do i=1,nx
  do j=1,ny
   if (ipoints(i,j) .eq. undef)cycle
   pc = 0
   call neighbr_check_d8_wipoints(i,j,demns,catchment,fdir,positions,nx,ny,npos,pc)
  enddo
 enddo

 !Calculate accumulation area
 area = res**2*catchment

 !Where the mask is 0 set area to undefined
 where (mask .eq. 0)
  area = undef
 endwhere
 where (fdir(:,:,1) .eq. -9999)
  area = undef
 endwhere

end subroutine

recursive subroutine neighbr_check_d8_wipoints(i,j,dem,catchment,fdir,positions,nx,ny,npos,pc)
 
 implicit none
 integer,intent(in) :: i,j,npos,nx,ny
 real,intent(in),dimension(nx,ny) :: dem
 integer,intent(inout),dimension(nx,ny) :: catchment
 integer,intent(in),dimension(nx,ny,2) :: fdir
 integer,intent(in),dimension(npos,2) :: positions
 integer,intent(inout) :: pc
 integer :: ipos,inew,jnew

 if (catchment(i,j) .le. 0)then
  catchment(i,j) = catchment(i,j) + pc + 1
  pc = catchment(i,j)
  do ipos=1,npos
   inew = i+positions(ipos,1)
   jnew = j+positions(ipos,2)
   if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny))cycle
   if (dem(i,j) .gt. dem(inew,jnew))then
    if ((fdir(i,j,1) .eq. inew) .and. (fdir(i,j,2) .eq. jnew))then
     call neighbr_check_d8_wipoints(inew,jnew,dem,catchment,fdir,positions,nx,ny,npos,pc)
     !catchment(inew,jnew) = catchment(inew,jnew) + catchment(i,j)
    endif
   endif
  enddo
 endif
    
end subroutine

recursive subroutine neighbr_check_d8(i,j,dem,catchment,fdir,positions,nx,ny,npos)
 
 implicit none
 integer,intent(in) :: i,j,npos,nx,ny
 real,intent(in),dimension(nx,ny) :: dem
 integer,intent(inout),dimension(nx,ny) :: catchment
 integer,intent(in),dimension(nx,ny,2) :: fdir
 integer,intent(in),dimension(npos,2) :: positions
 integer :: ipos,inew,jnew

 if (catchment(i,j) .le. 0)then
  catchment(i,j) = 1
  do ipos=1,npos
   inew = i+positions(ipos,1)
   jnew = j+positions(ipos,2)
   if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny))cycle
   if (dem(inew,jnew) .gt. dem(i,j))then
    if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
     call neighbr_check_d8(inew,jnew,dem,catchment,fdir,positions,nx,ny,npos)
     catchment(i,j) = catchment(i,j) + catchment(inew,jnew)
    endif
   endif
  enddo
 endif
    
end subroutine

subroutine calculate_d8_acc_neighbors(dem,res,variable,area,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: res
 real,intent(in),dimension(nx,ny) :: dem,variable
 real,intent(out),dimension(nx,ny) :: area
 integer,dimension(nx,ny,2) :: fdir
 integer,allocatable,dimension(:,:) :: positions
 real,allocatable,dimension(:) :: slopes
 real :: length,catchment(nx,ny)
 integer :: i,j,k,l,pos,tmp(1),npos=8
 integer :: ipos,inew,jnew
 allocate(positions(npos,2),slopes(npos))

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !get flow direction map
 do i=1,nx
  do j=1,ny
   slopes = 0.0
   do pos=1,npos
    k = positions(pos,1)
    l = positions(pos,2)
    if ((i+k .lt. 1) .or. (j+l .lt. 1) .or. (i+k .gt. nx) .or. &
        (j+l .gt. ny)) cycle !skip due to on boundary
    if ((k + l .eq. -2) .or. (k + l .eq. 2) .or. (k + l .eq. 0))then
        length = 1.41421356237*res
    else 
       length = res
    endif
    slopes(pos) = (dem(i,j) - dem(i+k,j+l))/length
   enddo
   if (maxval(slopes) .gt. 0) then
    tmp = maxloc(slopes)
    fdir(i,j,1) = i+positions(tmp(1),1)
    fdir(i,j,2) = j+positions(tmp(1),2)
   else
    fdir(i,j,:) = -9999
   endif
  enddo
 enddo

 !get the cell count
 catchment(:,:) = 0
 do i=1,nx
  do j=1,ny
   do ipos=1,npos
    inew = i+positions(ipos,1)
    jnew = j+positions(ipos,2)
    if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny))cycle
    if (dem(inew,jnew) .gt. dem(i,j))then
     if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
      catchment(i,j) = catchment(i,j) + variable(i,j)
     endif
    endif
   enddo
  enddo
 enddo

 !Calculate accumulation area
 area = catchment

end subroutine

subroutine calculate_mfd_acc(dem,res,p,area,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: res,p
 real,intent(in),dimension(nx,ny) :: dem
 real,intent(out),dimension(nx,ny) :: area
 integer,allocatable,dimension(:,:) :: positions
 real,allocatable,dimension(:) :: slopes
 integer :: i,j,k,l,pos,npos=8
 real :: catchment(nx,ny)
 allocate(positions(npos,2),slopes(npos))

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !get the cell count
 catchment(:,:) = 0
 do i=1,nx
  do j=1,ny
   call neighbr_check_mfd(i,j,dem,catchment,positions,nx,ny,npos,p,res)
  enddo
 enddo

 !Calculate accumulation area
 area = res**2*catchment

end subroutine

recursive subroutine neighbr_check_mfd(i,j,dem,catchment,positions,nx,ny,npos,p,res)
 
 implicit none
 integer,intent(in) :: i,j,npos,nx,ny
 real,intent(in) :: p,res
 real,intent(in),dimension(nx,ny) :: dem
 real,intent(inout),dimension(nx,ny) :: catchment
 integer,intent(in),dimension(npos,2) :: positions
 integer :: ipos,inew,jnew
 real :: fract

 if (catchment(i,j) .le. 0)then
  catchment(i,j) = 1
  do ipos=1,npos
   inew = i+positions(ipos,1)
   jnew = j+positions(ipos,2)
   if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny)) cycle
   if (dem(inew,jnew) .gt. dem(i,j))then
    call neighbr_check_mfd(inew,jnew,dem,catchment,positions,nx,ny,npos,p,res)
    call fract_flow_mfd(i,j,inew,jnew,dem,fract,p,positions,nx,ny,npos,res)
    catchment(i,j) = catchment(i,j) + fract*catchment(inew,jnew)
   endif
  enddo
 endif
    
end subroutine

subroutine fract_flow_mfd(iorg,jorg,i,j,dem,fract,p,positions,nx,ny,npos,res)

 implicit none
 integer,intent(in) :: iorg,jorg,i,j,nx,ny,npos
 real,intent(in) :: res,p
 real,intent(out) :: fract
 real,intent(in),dimension(nx,ny) :: dem
 integer,intent(in),dimension(npos,2) :: positions
 integer :: ipos,inew,jnew,k,l
 real :: angle_sum,slope,length,slopes(npos)
 slope = 0.0
 slopes = 0.0
 !Calculate all the slopes for the surrounding cells
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny))cycle
  if (dem(i,j) .gt. dem(inew,jnew))then
   k = inew - i
   l = jnew - j
   if ((k + l .eq. -2) .or. (k + l .eq. 2) .or. (k + l .eq. 0))then
    length = 1.41421356237*res
   else
    length = res
   endif
   slopes(ipos) = (dem(i,j) - dem(inew,jnew))/length
  endif
 enddo
 !Calculate sum of angles
 angle_sum = sum(slopes**p)
 !Calculate fraction
 k = iorg - i
 l = jorg - j
 if ((k + l .eq. -2) .or. (k + l .eq. 2) .or. (k + l .eq. 0))then
  length = 1.41421356237*res
 else
  length = res
 endif
 if (angle_sum .eq. 0.0)then
  fract = 0.0
 else
  slope = (dem(i,j) - dem(iorg,jorg))/length
  fract = slope**p/angle_sum
 endif

end subroutine

subroutine calculate_depth2channel_mfd(channels,mask,p,dem,res,depth2channel,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny,p
 integer,intent(in) :: channels(nx,ny)
 integer,intent(inout) :: mask(nx,ny)
 real,intent(in) :: dem(nx,ny),res
 real,intent(out) :: depth2channel(nx,ny)
 real :: channeldepth(nx,ny),cd
 integer :: i,j,positions(8,2),k,pos,l
 real :: undef = -9999.0

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !Initialize the channel depth array to -9999.0
 channeldepth = dem
 depth2channel = undef

 !If the channeldepth is below 0 then set to 0
 where(channeldepth .lt. 0)
  channeldepth = 0.0
 endwhere

 !Mask out all the elevation elements that are not channels
 where ((mask .le. 0) .or. (channels .le. 0))
  channeldepth = undef
 endwhere

 !Iterate cell by cell
 do i=1,nx
  do j=1,ny
   !Only work on this cell if the basin id is unknown and the mask is positive
   if ((channeldepth(i,j) .eq. undef) .and. (mask(i,j) .ge. 1)) then
    call determine_channel_depth_mfd(i,j,channeldepth,cd,mask,nx,ny,p,positions,dem,res)
   endif
  enddo
 enddo

 !Calculate the depth2channel by subtracting the channel depth
 depth2channel = dem - channeldepth

 !Set any channel depths below 0 to 0 (This is a hack)
 where (depth2channel .lt. 0) 
  depth2channel = undef
 endwhere

 !Set all values that are outside of the mask to undef
 where (mask .le. 0)
  depth2channel = undef
 endwhere

end subroutine

recursive subroutine determine_channel_depth_mfd(i,j,channeldepth,cd,mask,nx,ny,p,positions,dem,res)

 implicit none
 integer,intent(in) :: i,j,nx,ny,p
 integer,intent(inout) :: mask(nx,ny)
 real,intent(inout) :: cd,channeldepth(nx,ny)
 real,intent(in) :: res
 real,intent(in),dimension(nx,ny) :: dem
 integer,intent(in) :: positions(8,2)
 integer :: ipos,inew,jnew,k,l,npos
 real :: angle_sum,slope,length,slopes(8),fract
 if (mask(i,j) .eq. 0)return
 npos = 8
 slopes = 0.0
 !Calculate all the slopes for the surrounding cells
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny))cycle
  if (mask(inew,jnew) .eq. 0)cycle
  if (dem(i,j) .gt. dem(inew,jnew))then
   k = inew - i
   l = jnew - j
   if ((k + l .eq. -2) .or. (k + l .eq. 2) .or. (k + l .eq. 0))then
    length = 1.41421356237*res
   else
    length = res
   endif
   slopes(ipos) = (dem(i,j) - dem(inew,jnew))/length
  endif
 enddo
 !Calculate sum of angles
 angle_sum = sum(slopes**p)
 !Iterate through all the positions
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  slope = slopes(ipos)
  if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny))cycle
  if (slope .eq. 0.0)cycle
  if (mask(inew,jnew) .eq. 0)cycle
  fract = slope**p/angle_sum
  !Figure out if downhill has a value if not then recurse. If it does then
  if (channeldepth(inew,jnew) .ge. 0)then
   cd = channeldepth(inew,jnew)
   if (channeldepth(i,j) .eq. -9999)channeldepth(i,j) = 0.0
   channeldepth(i,j) = channeldepth(i,j) + fract*cd
  else
   call determine_channel_depth_mfd(inew,jnew,channeldepth,cd,mask,nx,ny,p,positions,dem,res)
   if (channeldepth(i,j) .eq. -9999)channeldepth(i,j) = 0.0
   channeldepth(i,j) = channeldepth(i,j) + fract*cd
  endif
 enddo
 !Return the cd
 cd = channeldepth(i,j)
  
end subroutine

subroutine calculate_channels(area_in,threshold,basin_threshold,fdir,channels,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: threshold,basin_threshold
 real,intent(in),dimension(nx,ny) :: area_in
 integer,intent(in),dimension(nx,ny,2) :: fdir
 integer,intent(out),dimension(nx,ny) :: channels
 real,dimension(nx,ny) :: area
 integer,dimension(nx,ny) :: mask
 integer,dimension(2) :: placement
 integer,dimension(:,:),allocatable :: positions
 integer :: i,j,pos,cid,k,l,npos
 logical :: bool
 real :: undef
 undef = -9999.0
 npos = 8
 allocate(positions(npos,2))
 !Copy the area array
 area = area_in
 !Initialize channels to undef
 channels = int(undef)

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !Define the channels mask
 where (area .gt. threshold)
  mask = 1
 elsewhere
  mask = 0
 endwhere

 !Differentiate the channels by segments
 cid = 1
 bool = .False.
 do while (bool .eqv. .False.)
  
  !Determine if there are still are cells
  if (maxval(mask) .eq. 0) bool = .True.

  !Maskout the area
  where (mask .eq. 0) 
   area = 0
  endwhere

  !Find the highest accumulation area
  placement = maxloc(area)
  i = placement(1)
  j = placement(2)
  !Set the channel id
  if ((mask(i,j) .eq. 1) .and. (area(i,j) .ge. basin_threshold))then
   channels(i,j) = cid
  endif
  mask(i,j) = 0

  !Go upstream
  call channels_upstream(i,j,fdir,channels,positions,nx,ny,cid,npos,&
                         mask,basin_threshold,area)

 enddo

end subroutine

subroutine gap_fill_hrus(hrus_in,channels,hrus_out,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(in) :: channels(nx,ny)
 integer,intent(in) :: hrus_in(nx,ny)
 integer,intent(out) :: hrus_out(nx,ny)
 integer :: i,j,ii,jj,imin,imax,jmin,jmax
 integer :: hru_count(8),hru_id(8),hid,hru,iloc(1)

 !Copy the original hrus map
 hrus_out = hrus_in

 !Iterate through each point (cid = 999999 is shoreline)
 do i = 1,nx
  do j = 1,ny
   if ((channels(i,j) .gt. 0) .and. (channels(i,j) .lt. 999999))then
    imin = i-1
    imax = i+1
    jmin = j-1
    jmax = j+1
    if (i .eq. 1)imin = 1
    if (i .eq. nx)imax = nx
    if (j .eq. 1)jmin = 1
    if (j .eq. ny)jmax = ny
    hru_count = 0
    hru_id = -9999
    hid = 1
    do ii = imin,imax
     do jj = jmin,jmax
      hru = hrus_in(ii,jj)
      if (hru .eq. -9999)cycle
      if (any(hru_id .eq. hru)) then
       where (hru_id .eq. hru)
        hru_count = hru_count + 1
       endwhere
      else
       hru_id(hid) = hru
       hru_count(hid) = hru_count(hid) + 1
       hid = hid + 1
      endif
     enddo
    enddo
    !Compute the most frequent
    iloc = maxloc(hru_count)
    !Set it to the out array
    hrus_out(i,j) = hru_id(iloc(1))
   endif
  enddo
 enddo

end subroutine gap_fill_hrus
 

subroutine calculate_channels_wocean(area_in,threshold,basin_threshold,fdir,mask,channels,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: threshold,basin_threshold
 real,intent(in),dimension(nx,ny) :: area_in,mask
 integer,intent(in),dimension(nx,ny,2) :: fdir
 integer,intent(out),dimension(nx,ny) :: channels
 real,dimension(nx,ny) :: area
 integer,dimension(nx,ny) :: cmask
 integer,dimension(2) :: placement
 integer,dimension(:,:),allocatable :: positions
 integer :: i,j,pos,cid,k,l,npos,imin,imax,jmin,jmax
 logical :: bool
 real :: undef
 undef = -9999.0
 npos = 8
 allocate(positions(npos,2))
 !Copy the area array
 area = area_in
 !Initialize channels 
 channels = 0

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !Define the channels mask
 where (area .gt. threshold)
  cmask = 1
 elsewhere
  cmask = 0
 endwhere

 !Differentiate the channels by segments
 cid = 1
 bool = .False.
 do while (bool .eqv. .False.)
  
  !Determine if there are still are cells
  if (maxval(cmask) .eq. 0) bool = .True.

  !Maskout the area
  where (cmask .eq. 0) 
   area = 0
  endwhere

  !Find the highest accumulation area
  placement = maxloc(area)
  i = placement(1)
  j = placement(2)
  !Set the channel id
  if ((cmask(i,j) .eq. 1) .and. (area(i,j) .ge. basin_threshold))then
   channels(i,j) = cid
  endif
  cmask(i,j) = 0

  !Go upstream
  call channels_upstream(i,j,fdir,channels,positions,nx,ny,cid,npos,&
                         cmask,basin_threshold,area)

 enddo

 !Set the ocean/land, lake/land, and glacier/land boundaries as "channels"
 cid = 999999!max(maxval(channels),1)
 do i = 1,nx
  do j = 1,ny
   !Determine if this point is not "land"
   if (mask(i,j) == 0.0) then
    imin = i-1
    imax = i+1
    jmin = j-1
    jmax = j+1
    !Determine if any of the surrounding points are land
    if (i .eq. 1)imin = 1
    if (i .eq. nx)imax = nx
    if (j .eq. 1)jmin = 1
    if (j .eq. ny)jmax = ny
    if (maxval(mask(imin:imax,jmin:jmax)) .gt. 0)channels(i,j) = cid
    cid = cid + 1
   endif
  enddo
 enddo

 !Where the mask is 0 set area to undefined
 where (((mask .eq. 0) .and. (channels .eq. 0)))! .or. (channels .eq. 0))
  channels = int(undef)
 endwhere

end subroutine

recursive subroutine channels_upstream(i,j,fdir,channels,positions,nx,ny,cid,npos,&
                             mask,basin_threshold,area)

 implicit none
 integer,intent(in) :: npos,i,j,nx,ny
 integer,intent(in) :: positions(npos,2),fdir(nx,ny,2)
 real,intent(in) :: basin_threshold,area(nx,ny)
 integer,intent(inout) :: cid,channels(nx,ny),mask(nx,ny)
 integer :: inew,jnew,count,ipos,cid_org
 !Memorize the channel id
 cid_org = cid

 !Determine how many cells flow into this cell 
 count = 0
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt.ny))cycle
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   if (mask(inew,jnew) .eq. 1) then
    !Make sure we are above the threshold
    if (area(inew,jnew) .ge. basin_threshold) then
     count = count + 1
    endif
   endif
  endif
 enddo
 !Decide the path to take
 !1.Only one upstream cell
 if (count .le. 1)then 
  do ipos=1,npos
   inew = i+positions(ipos,1)
   jnew = j+positions(ipos,2)
   if ((inew.lt.1).or.(jnew.lt.1).or.(inew.gt.nx).or.(jnew.gt.ny)) cycle
   if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
    if (mask(inew,jnew) .eq. 1) then
     mask(inew,jnew) = 0
     channels(inew,jnew) = channels(i,j)
     call channels_upstream(inew,jnew,fdir,channels,positions,nx,ny,&
                             cid,npos,mask,basin_threshold,area)
    endif
   endif
  enddo
 !2.More than one upstream cell
 elseif (count .gt. 1)then
  do ipos=1,npos
   inew = i+positions(ipos,1)
   jnew = j+positions(ipos,2)
   if ((inew.lt.1).or.(jnew.lt.1).or.(inew.gt.nx).or.(jnew.gt.ny)) cycle
   if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
    if (mask(inew,jnew) .eq. 1) then
     if (area(inew,jnew) .ge. basin_threshold)then
      cid = cid + 1
      mask(inew,jnew) = 0
      channels(inew,jnew) = cid
      call channels_upstream(inew,jnew,fdir,channels,positions,nx,ny,&
                             cid,npos,mask,basin_threshold,area)
     else
      mask(inew,jnew) = 0
      channels(inew,jnew) = cid_org
      call channels_upstream(inew,jnew,fdir,channels,positions,nx,ny,&
                             cid_org,npos,mask,basin_threshold,area)
     endif
    endif
   endif
  enddo
 endif

end subroutine

subroutine delineate_basins(channels,basins,mask,fdir,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(in) :: channels(nx,ny),fdir(nx,ny,2)
 integer,intent(in) :: mask(nx,ny)
 integer,intent(out) :: basins(nx,ny)
 integer :: i,j,basin_id,cmask(nx,ny)
 cmask = mask

 !Initialize the basin delineation to the channel network (everythin else 0)
 basins = channels

 !Iterate cell by cell
 do i=1,nx
  do j=1,ny
   !Only work on this cell if the basin id is unknown and the mask is positive
   if ((basins(i,j) .eq. 0) .and. (cmask(i,j) .ge. 1)) then
    !Find the id
    call determine_basin_id(i,j,basins,basin_id,fdir,cmask,nx,ny)
   endif
  enddo
 enddo
 !Go downstream to find the correct id
 !When it reaches a value then stop recursing
 !Clean up the "lost ones"

end subroutine

recursive subroutine determine_basin_id(i,j,basins,basin_id,fdir,mask,nx,ny)

 implicit none
 integer,intent(in) :: i,j,nx,ny,fdir(nx,ny,2)
 integer,intent(inout) :: basin_id,basins(nx,ny),mask(nx,ny)
 integer :: inew,jnew
 basin_id = 0
 !Determine which is way down
 inew = fdir(i,j,1)
 jnew = fdir(i,j,2)
 if ((inew.lt.1).or.(jnew.lt.1).or.(inew.gt.nx).or.(jnew.gt.ny))return
 if (mask(i,j) .eq. 0)return
 !Figure out if downhill has a value if not then recurse. If it does then
 if (basins(inew,jnew) .gt. 0)then
  basin_id = basins(inew,jnew)
  basins(i,j) = basin_id
 else
  call determine_basin_id(inew,jnew,basins,basin_id,fdir,mask,nx,ny)
  basins(i,j) = basin_id
 endif

end subroutine

subroutine delineate_hillslopes(channels,area_in,fdir,mask,hillslopes,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(in) :: channels(nx,ny),fdir(nx,ny,2),mask(nx,ny)
 real,intent(in) :: area_in(nx,ny)
 integer,intent(out) :: hillslopes(nx,ny)
 real :: maxarea,area(nx,ny)
 integer :: i,j,placement(2),hillslope_id,count
 area = area_in

 !Mask out all the elements that are not channels
 hillslopes = 0
 where ((mask .le. 0) .or. (channels .gt. 0))
  hillslopes = int(-9999.0)
 endwhere

 !Initialize the hillslope_id
 hillslope_id = 1

 !Iterate through all the catchments starting at the outlet
 count = 0
 do while (maxarea .ne. -9999.0)
  maxarea = maxval(area)
  if (maxarea .eq. -9999.0)exit
  placement = maxloc(area)
  i = placement(1)
  j = placement(2)
  !Go inside
  !print*,count,maxarea,placement,channels(i,j),area(i,j),fdir(i,j,:),mask(i,j)
  call delineate_hillslopes_catchment(channels,area,fdir,mask,hillslopes,nx,ny,i,j,hillslope_id)
  !Set the area value to undefined
  area(i,j) = -9999
  count = count + 1
  !if (count .eq. 100)exit
 enddo

 !Cleanup hillslopes
 call cleanup_hillslopes(hillslopes,nx,ny)

end subroutine

subroutine delineate_hillslopes_catchment(channels,area,fdir,mask,hillslopes,nx,ny,i,j,hillslope_id)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(in) :: channels(nx,ny),fdir(nx,ny,2),mask(nx,ny)
 real,intent(inout) :: area(nx,ny)
 integer,intent(inout) :: hillslope_id,i,j
 integer,intent(out) :: hillslopes(nx,ny)
 !integer,dimension(2) :: placement
 integer,dimension(:,:),allocatable :: positions
 integer :: ipos,pos,k,l,npos,inew,jnew,iold,jold,ipos_old
 integer :: cid
 npos = 8
 ipos_old = -9999
 allocate(positions(npos,2))

 !Construct positions array
 pos = 0
 do k=-1,1
  do l=-1,1
   if ((k == 0) .and. (l == 0)) cycle
   pos = pos + 1
   positions(pos,1) = k
   positions(pos,2) = l
  enddo
 enddo

 !Calculate the position of the maximum area
 !placement = maxloc(area)
 !i = placement(1)
 !j = placement(2)

 !Initialize the hillslope_id
 !hillslope_id = 1

 !Set the channel
 cid = channels(i,j)

 !Figure out the origin position
 iold = fdir(i,j,1)
 jold = fdir(i,j,2)
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((inew .eq. iold) .and. (jnew .eq. jold))then
   ipos_old = ipos
  endif
 enddo

 !Define the positions (clock-wise orientation)
 positions(:,:) = -9999
 positions(1,1:2) = (/-1,0/)
 positions(2,1:2) = (/-1,-1/)
 positions(3,1:2) = (/0,-1/)
 positions(4,1:2) = (/1,-1/)
 positions(5,1:2) = (/1,0/)
 positions(6,1:2) = (/1,1/)
 positions(7,1:2) = (/0,1/)
 positions(8,1:2) = (/-1,1/)

 !Start from the downstream cell
 do ipos=ipos_old,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((inew.lt.1).or.(jnew.lt.1).or.(inew.gt.nx).or.(jnew.gt.ny))cycle
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !If it is a channel then move upstream
   if (channels(inew,jnew) .gt. 0)then
    call move_upstream(inew,jnew,hillslope_id,hillslopes,fdir,&
                       channels,nx,ny,positions,i,j,cid,mask,area)
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask,area)
   endif
  endif
 enddo

 !Continue to close the loop
 do ipos=1,ipos_old-1
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !If it is a channel then move upstream
   if (channels(inew,jnew) .gt. 0)then
    call move_upstream(inew,jnew,hillslope_id,hillslopes,fdir,&
                       channels,nx,ny,positions,i,j,cid,mask,area)
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask,area)
   endif
  endif
 enddo
 
end subroutine

recursive subroutine move_upstream(i,j,hillslope_id,hillslopes,fdir,&
                                   channels,nx,ny,positions,&
                                   iold,jold,cid,mask,area)

 implicit none
 integer,intent(in) :: i,j,nx,ny,positions(8,2)
 integer,intent(in) :: channels(nx,ny),fdir(nx,ny,2),mask(nx,ny)
 real,intent(inout) :: area(nx,ny)
 integer,intent(inout) :: hillslopes(nx,ny),hillslope_id,cid
 integer :: inew,jnew,ipos,npos=8,iold,jold,ipos_old,channel_count

 ipos_old = -9999
 !Initialize channel count
 channel_count = 0
 area(iold,jold) = -9999
 
 !Figure out the origin position
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((inew .eq. iold) .and. (jnew .eq. jold))then
   ipos_old = ipos
  endif
 enddo

 !Determine if it is a node
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !If it is a channel then move upstream
   if (channels(inew,jnew) .gt. 0)then
    !Count how many channels there are that flow into here
    channel_count = channel_count + 1
   endif
  endif
 enddo

 !Start from the downstream cell
 do ipos=ipos_old,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !If it is a channel then move upstream
   if (channels(inew,jnew) .gt. 0)then
    !If we are on a different channel link then update
    if (channel_count .gt. 1) then
     cid = channels(i,j)
     hillslope_id = hillslope_id + 1
    endif
    call move_upstream(inew,jnew,hillslope_id,hillslopes,fdir,&
                       channels,nx,ny,positions,i,j,cid,mask,area)
    !If we are on a different channel link then update
    if (channel_count .gt. 1)then
     cid = channels(i,j)
     hillslope_id = hillslope_id + 1
    endif
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask,area)
   endif
  endif
 enddo

 !Continue to close the loop
 do ipos=1,ipos_old-1
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !If it is a channel then move upstream
   if (channels(inew,jnew) .gt. 0)then
    !If we are on a different channel link then update
    if (channel_count .gt. 1) then
     cid = channels(i,j)
     hillslope_id = hillslope_id + 1
    endif
    call move_upstream(inew,jnew,hillslope_id,hillslopes,fdir,&
                       channels,nx,ny,positions,i,j,cid,mask,area)
    !If we are on a different channel link then update
    if (channel_count .gt. 1) then
     cid = channels(i,j)
     hillslope_id = hillslope_id + 1
    endif
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask,area)
   endif
  endif
 enddo

 !If we are at a terminal point then create own hillslope
 if (channel_count == 0)then
  hillslope_id = hillslope_id + 1
  do ipos=1,npos
   inew = i+positions(ipos,1)
   jnew = j+positions(ipos,2)
   if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask,area)
   endif
  enddo
  !Update the hillslope id for the other side slope
  hillslope_id = hillslope_id + 1
 endif
 
 !Set channel area to undefined
 area(i,j) = -9999.0

end subroutine

recursive subroutine define_hillslope_id(i,j,hillslope_id,hillslopes,&
                                         fdir,nx,ny,positions,mask,area)

 implicit none
 integer,intent(in) :: i,j,nx,ny,fdir(nx,ny,2),positions(8,2),mask(nx,ny)
 integer,intent(inout) :: hillslope_id,hillslopes(nx,ny)
 real,intent(inout) :: area(nx,ny)
 integer :: inew,jnew,ipos,npos=8
 !Determine if the point is in the mask
 if (mask(i,j) .eq. 0)return
 !Define the id
 hillslopes(i,j) = hillslope_id
 area(i,j) = -9999
 !Determine the cells that flow into this cell
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !Recurse to place hillslope id 
   call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,fdir&
                            ,nx,ny,positions,mask,area)
  endif
 enddo

end subroutine

subroutine assign_properties_to_hillslopes(hillslopes,hp_slope,hp_area,slope,area,nx,ny,nh)

 implicit none
 integer,intent(in) :: nx,ny,nh
 integer,intent(in) :: hillslopes(nx,ny)
 real,intent(in) :: hp_slope(nh),hp_area(nh)
 real,intent(out) :: slope(nx,ny)
 real,intent(out) :: area(nx,ny)
 integer :: i,j,undef
 undef = 0

 !Initialize the new array
 slope(:,:) = -9999.0
 area(:,:) = -9999.0

 !Go through and set the ids
 do i=1,nx
  do j=1,ny
   if (hillslopes(i,j) .ne. -9999)then
    slope(i,j) = hp_slope(hillslopes(i,j))
    area(i,j) = hp_area(hillslopes(i,j))
   endif
  enddo
 enddo

end subroutine

subroutine cleanup_hillslopes(hillslopes,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(inout) :: hillslopes(nx,ny)
 integer,allocatable,dimension(:) :: hids,hcounts,mapping
 integer :: hid,i,j

 !Create array with hillslopes ids and count
 allocate(hids(maxval(hillslopes)))
 allocate(hcounts(size(hids)))
 allocate(mapping(size(hids)))
 hcounts = 0
 mapping = 0
 do i=1,size(hids)
  hids(i) = i
 enddo

 !Count the hillslopes 
 do i=1,nx
  do j=1,ny
   if (hillslopes(i,j) .gt. 0)then
    hid = hillslopes(i,j)
    hcounts(hid) = hcounts(hid) + 1
   endif
  enddo
 enddo

 !Assign mapping
 hid = 0
 do i=1,size(hids)
  if (hcounts(i) .gt. 0)then
   mapping(i) = hid
   hid = hid + 1
  endif
 enddo

 !Assign new hillslope ids
 do i=1,nx
  do j=1,ny
   if (hillslopes(i,j) .gt. 0)then
    hillslopes(i,j) = mapping(hillslopes(i,j))
   endif
  enddo
 enddo

end subroutine

subroutine calculate_depth2channel(channels,mask,fdir,dem,depth2channel,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(in) :: channels(nx,ny),fdir(nx,ny,2)
 integer,intent(inout) :: mask(nx,ny)
 real,intent(in) :: dem(nx,ny)
 real,intent(out) :: depth2channel(nx,ny)
 real :: channeldepth(nx,ny),cd
 integer :: i,j
 real :: undef = -9999.0

 !Initialize the channel depth array to -9999.0
 channeldepth = dem
 depth2channel = undef

 !If the channeldepth is below 0 then set to 0
 where(channeldepth .lt. 0)
  channeldepth = 0.0
 endwhere

 !Mask out all the elevation elements that are not channels
 where ((mask .le. 0) .or. (channels .le. 0))
  channeldepth = undef
 endwhere

 !Iterate cell by cell
 cd = 0.0
 do i=1,nx
  do j=1,ny
   !Only work on this cell if the basin id is unknown and the mask is positive
   if ((channeldepth(i,j) .eq. undef) .and. (mask(i,j) .ge. 1)) then
    call determine_channel_depth(i,j,channeldepth,cd,fdir,mask,nx,ny)
   endif
  enddo
 enddo

 !Calculate the depth2channel by subtracting the channel depth
 depth2channel = dem - channeldepth

 !Set any channel depths below 0 to 0 (This is a hack)
 where (depth2channel .lt. 0) 
  depth2channel = undef
 endwhere

 !Set all values that are outside of the mask to undef
 where (mask .le. 0)
  depth2channel = undef
 endwhere

end subroutine

recursive subroutine determine_channel_depth(i,j,channeldepth,cd,fdir,mask,nx,ny)

 implicit none
 integer,intent(in) :: i,j,nx,ny,fdir(nx,ny,2)
 integer,intent(inout) :: mask(nx,ny)
 real,intent(inout) :: cd,channeldepth(nx,ny)
 integer :: inew,jnew
 !Determine which way is down
 inew = fdir(i,j,1)
 jnew = fdir(i,j,2)
 if ((inew.lt.1).or.(jnew.lt.1).or.(inew.gt.nx).or.(jnew.gt.ny))return
 if (mask(i,j) .eq. 0)return
 !Figure out if downhill has a value if not then recurse. If it does then
 if (channeldepth(inew,jnew) .gt. 0)then
  cd = channeldepth(inew,jnew)
  channeldepth(i,j) = cd
 else
  call determine_channel_depth(inew,jnew,channeldepth,cd,fdir,mask,nx,ny)
  channeldepth(i,j) = cd
 endif

end subroutine

subroutine assign_clusters_to_hillslopes(hillslopes_org,clusters,hillslopes_new,nx,ny,nh)

 implicit none
 integer,intent(in) :: nx,ny,nh
 integer,intent(in) :: hillslopes_org(nx,ny),clusters(nh)
 integer,intent(out) :: hillslopes_new(nx,ny)
 integer :: i,j,undef
 undef = -9999

 !Initialize the new array
 hillslopes_new = undef!hillslopes_org

 !Go through and set the ids
 do i=1,nx
  do j=1,ny
   if (hillslopes_org(i,j) .ne. undef)then
    hillslopes_new(i,j) = clusters(hillslopes_org(i,j) + 1)
   endif
  enddo
 enddo

end subroutine
