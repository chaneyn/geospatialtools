subroutine calculate_d8_acc(dem,res,area,fdir,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: res
 real,intent(in),dimension(nx,ny) :: dem
 real,intent(out),dimension(nx,ny) :: area
 integer,intent(out),dimension(nx,ny,2) :: fdir
 integer,allocatable,dimension(:,:) :: positions
 real,allocatable,dimension(:) :: slopes
 real :: length
 integer :: i,j,k,l,pos,tmp(1),npos=8,catchment(nx,ny)
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
   call neighbr_check_d8(i,j,dem,catchment,fdir,positions,nx,ny,npos)
  enddo
 enddo

 !Calculate accumulation area
 area = res**2*catchment

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
   if ((inew .lt. 1) .or. (jnew .lt. 1) .or. (inew .gt. nx) .or. (jnew .gt. ny)) cycle
   if (dem(inew,jnew) .gt. dem(i,j))then
    if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
     call neighbr_check_d8(inew,jnew,dem,catchment,fdir,positions,nx,ny,npos)
     catchment(i,j) = catchment(i,j) + catchment(inew,jnew)
    endif
   endif
  enddo
 endif
    
end subroutine

subroutine calculate_mfd_acc(dem,res,area,nx,ny,p)

 implicit none
 integer,intent(in) :: nx,ny
 real,intent(in) :: res,p
 real,intent(in),dimension(nx,ny) :: dem
 real,intent(out),dimension(nx,ny) :: area
 integer,allocatable,dimension(:,:) :: positions
 real,allocatable,dimension(:) :: slopes
 real :: length
 integer :: i,j,k,l,pos,tmp(1),npos=8
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
 real,intent(in) :: p,res
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
 angle_sum = sum(slopes)
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
  fract = slope/angle_sum
 endif

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
 npos = 8
 allocate(positions(npos,2))
 !Copy the area array
 area = area_in

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

subroutine channels_upstream(i,j,fdir,channels,positions,nx,ny,cid,npos,&
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
 integer,intent(in) :: channels(nx,ny),mask(nx,ny),fdir(nx,ny,2)
 integer,intent(out) :: basins(nx,ny)
 integer :: i,j,basin_id

 !Initialize the basin delineation to the channel network (everythin else 0)
 basins = channels

 !Iterate cell by cell
 do i=1,nx
  do j=1,ny
   !Only work on this cell if the basin id is unknown and the mask is positive
   if ((basins(i,j) .eq. 0) .and. (mask(i,j) .ge. 1)) then
    !Find the id
    call determine_basin_id(i,j,basins,basin_id,fdir,mask,nx,ny)
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
