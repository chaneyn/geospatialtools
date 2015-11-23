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

subroutine delineate_hillslopes(channels,area,fdir,mask,hillslopes,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(in) :: channels(nx,ny),fdir(nx,ny,2),mask(nx,ny)
 real,intent(in) :: area(nx,ny)
 integer,intent(out) :: hillslopes(nx,ny)
 integer,dimension(2) :: placement
 integer,dimension(:,:),allocatable :: positions
 integer :: i,j,ipos,pos,k,l,npos,inew,jnew,iold,jold,ipos_old
 integer :: hillslope_id,cid
 npos = 8
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
 placement = maxloc(area)
 i = placement(1)
 j = placement(2)

 !Initialize the hillslope_id
 hillslope_id = 1

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
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !If it is a channel then move upstream
   if (channels(inew,jnew) .gt. 0)then
    call move_upstream(inew,jnew,hillslope_id,hillslopes,fdir,&
                       channels,nx,ny,positions,i,j,cid,mask)
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask)
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
                       channels,nx,ny,positions,i,j,cid,mask)
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask)
   endif
  endif
 enddo
 
 !Cleanup hillslopes
 call cleanup_hillslopes(hillslopes,nx,ny)

end subroutine

recursive subroutine move_upstream(i,j,hillslope_id,hillslopes,fdir,&
                                   channels,nx,ny,positions,&
                                   iold,jold,cid,mask)

 implicit none
 integer,intent(in) :: i,j,nx,ny,positions(8,2)
 integer,intent(in) :: channels(nx,ny),fdir(nx,ny,2),mask(nx,ny)
 integer,intent(inout) :: hillslopes(nx,ny),hillslope_id,cid
 integer :: inew,jnew,ipos,npos=8,iold,jold,ipos_old,channel_count

 !Initialize channel count
 channel_count = 0
 
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
                       channels,nx,ny,positions,i,j,cid,mask)
    !If we are on a different channel link then update
    if (channel_count .gt. 1)then
     cid = channels(i,j)
     hillslope_id = hillslope_id + 1
    endif
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask)
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
                       channels,nx,ny,positions,i,j,cid,mask)
    !If we are on a different channel link then update
    if (channel_count .gt. 1) then
     cid = channels(i,j)
     hillslope_id = hillslope_id + 1
    endif
   !If it not a channel then recurse to define the id
   else
    !Recurse to place hillslope id 
    call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,&
                             fdir,nx,ny,positions,mask)
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
                             fdir,nx,ny,positions,mask)
   endif
  enddo
  !Update the hillslope id for the other side slope
  hillslope_id = hillslope_id + 1
 endif

end subroutine

recursive subroutine define_hillslope_id(i,j,hillslope_id,hillslopes,&
                                         fdir,nx,ny,positions,mask)

 implicit none
 integer,intent(in) :: i,j,nx,ny,fdir(nx,ny,2),positions(8,2),mask(nx,ny)
 integer,intent(inout) :: hillslope_id,hillslopes(nx,ny)
 integer :: inew,jnew,ipos,npos=8
 !Determine if the point is in the mask
 if (mask(i,j) .eq. 0)return
 !Define the id
 hillslopes(i,j) = hillslope_id
 !Determine the cells that flow into this cell
 do ipos=1,npos
  inew = i+positions(ipos,1)
  jnew = j+positions(ipos,2)
  if ((fdir(inew,jnew,1) .eq. i) .and. (fdir(inew,jnew,2) .eq. j))then
   !Recurse to place hillslope id 
   call define_hillslope_id(inew,jnew,hillslope_id,hillslopes,fdir&
                            ,nx,ny,positions,mask)
  endif
 enddo

end subroutine

subroutine calculate_hillslope_properties(hillslopes,dem,basins,res,nh,&
           latitude,longitude,depth2channel,&
           hillslopes_elevation,hillslopes_area,hillslopes_basin,&
           hillslopes_latitude,hillslopes_longitude,hillslopes_range,&
           hillslopes_id,hillslopes_depth2channel,&
           nx,ny)

 implicit none
 integer,intent(in) :: nx,ny,hillslopes(nx,ny),nh,basins(nx,ny)
 real,intent(in) :: dem(nx,ny),res,latitude(nx,ny),longitude(nx,ny)
 real,intent(in) :: depth2channel(nx,ny)
 integer,intent(out) :: hillslopes_basin(nh),hillslopes_id(nh)
 real,intent(out) :: hillslopes_elevation(nh),hillslopes_area(nh)
 real,intent(out) :: hillslopes_latitude(nh),hillslopes_longitude(nh)
 real,intent(out) :: hillslopes_range(nh),hillslopes_depth2channel(nh)
 real :: hillslopes_maxelev(nh),hillslopes_minelev(nh)
 integer :: hillslopes_count(nh)
 integer :: i,j,ih

 !Hillslope elevation,latitude,longitude
 hillslopes_minelev = 9.99e+08
 hillslopes_maxelev = -9.99e+08
 hillslopes_elevation = 0.0
 hillslopes_latitude = 0.0
 hillslopes_longitude = 0.0
 hillslopes_count = 0
 do i=1,nx
  do j=1,ny
   ih = hillslopes(i,j)
   if (ih .gt. 0)then
    hillslopes_elevation(ih) = hillslopes_elevation(ih) + dem(i,j)
    hillslopes_depth2channel(ih) = hillslopes_depth2channel(ih) + depth2channel(i,j)
    if (dem(i,j) .gt. hillslopes_maxelev(ih))hillslopes_maxelev(ih) = dem(i,j)
    if (dem(i,j) .lt. hillslopes_minelev(ih))hillslopes_minelev(ih) = dem(i,j)
    hillslopes_latitude(ih) = hillslopes_latitude(ih) + latitude(i,j)
    hillslopes_longitude(ih) = hillslopes_longitude(ih) + longitude(i,j)
    hillslopes_count(ih) = hillslopes_count(ih) + 1
   endif
  enddo
 enddo
 hillslopes_depth2channel = hillslopes_depth2channel/hillslopes_count
 hillslopes_elevation = hillslopes_elevation/hillslopes_count
 hillslopes_latitude = hillslopes_latitude/hillslopes_count
 hillslopes_longitude = hillslopes_longitude/hillslopes_count
 hillslopes_range = hillslopes_maxelev - hillslopes_minelev

 !Hillslope area
 hillslopes_area = res**2*hillslopes_count

 !Corresponding subbasin
 hillslopes_basin = 0
 do i=1,nx
  do j=1,ny
   ih = hillslopes(i,j)
   if (ih .gt. 0)then
    hillslopes_basin(ih) = basins(i,j)
   endif
  enddo
 enddo

 !Define the index
 hillslopes_id = 0
 do ih=1,nh
  hillslopes_id(ih) = ih
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
   hid = hid + 1
   mapping(i) = hid
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
 integer,intent(in) :: channels(nx,ny),mask(nx,ny),fdir(nx,ny,2)
 real,intent(in) :: dem(nx,ny)
 real,intent(out) :: depth2channel(nx,ny)
 real :: channeldepth(nx,ny),cd
 integer :: i,j
 real :: undef = -9.99e+08

 !Initialize the channel depth array to -9999.0
 channeldepth = dem
 depth2channel = undef

 !Mask out all the elevation elements that are not channels
 where ((mask .le. 0) .or. (channels .le. 0))
  channeldepth = undef
 endwhere

 !Iterate cell by cell
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

subroutine calculate_hillslopesd8(channels,mask,fdir,hillslopes,nx,ny)

 implicit none
 integer,intent(in) :: nx,ny
 integer,intent(in) :: channels(nx,ny),mask(nx,ny),fdir(nx,ny,2)
 integer,intent(out) :: hillslopes(nx,ny)
 integer :: i,j,ih,maxih
 real :: undef = -9999

 !Initialize the hillslopes array
 hillslopes = channels

 !Mask out all the elements that are not channels
 where ((mask .le. 0) .or. (channels .le. 0))
  hillslopes = undef
 endwhere

 !Iterate cell by cell
 maxih = 0
 do i=1,nx
  do j=1,ny
   !Only work on this cell if the basin id is unknown and the mask is positive
   if ((hillslopes(i,j) .eq. undef) .and. (mask(i,j) .ge. 1)) then
    ih = maxih
    call determine_hillslopesd8(i,j,hillslopes,ih,fdir,mask,channels,nx,ny)
    if (ih .gt. maxih) maxih = ih
   endif
  enddo
 enddo

 !Maskout the channels
 where (channels .gt. 0)
  hillslopes = undef
 endwhere

 !Clean up the hillslopes
 call cleanup_hillslopes(hillslopes,nx,ny)
 
end subroutine

recursive subroutine determine_hillslopesd8(i,j,hillslopes,ih,fdir,mask,&
                                            channels,nx,ny)

 implicit none
 integer,intent(in) :: i,j,nx,ny,fdir(nx,ny,2),channels(nx,ny)
 integer,intent(inout) :: mask(nx,ny),ih,hillslopes(nx,ny)
 integer :: inew,jnew
 !Determine which way is down
 inew = fdir(i,j,1)
 jnew = fdir(i,j,2)
 if ((inew.lt.1).or.(jnew.lt.1).or.(inew.gt.nx).or.(jnew.gt.ny))return
 if (mask(i,j) .eq. 0)return
 !Figure out if downhill has a value if not then recurse. If it does then
 if (hillslopes(inew,jnew) .gt. 0)then
  !If it is a channel then increase the id
  if (channels(inew,jnew) .gt. 0)then
   ih = ih + 1
  else 
   ih = hillslopes(inew,jnew)
  endif
  !ih = hillslopes(inew,jnew)
  hillslopes(i,j) = ih
 else
  call determine_hillslopesd8(inew,jnew,hillslopes,ih,fdir,mask,channels,nx,ny)
  hillslopes(i,j) = ih
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
 hillslopes_new = hillslopes_org

 !Go through and set the ids
 do i=1,nx
  do j=1,ny
   if (hillslopes_org(i,j) .ne. undef)then
    hillslopes_new(i,j) = clusters(hillslopes_org(i,j))
   endif
  enddo
 enddo

end subroutine

subroutine calculate_hru_properties(hillslopes,tiles,channels,nhru,res,nhillslope,&
           hrus,dem,slope,&
           hru_bwidth,hru_twidth,hru_length,hru_position,&
           hru_hid,hru_tid,hru_id,hru_area,hru_dem,hru_slope,&
           nx,ny,nc)

 implicit none
 integer,intent(in) :: nx,ny,nhru,nc,nhillslope(nc)
 integer,intent(in) :: hillslopes(nx,ny),tiles(nx,ny),channels(nx,ny),hrus(nx,ny)
 real,intent(in) :: dem(nx,ny),slope(nx,ny)
 real,intent(in) :: res
 real,intent(out) :: hru_bwidth(nhru),hru_twidth(nhru),hru_length(nhru)
 real,intent(out) ::hru_area(nhru),hru_dem(nhru),hru_position(nhru),hru_slope(nhru)
 integer,intent(out) :: hru_hid(nhru),hru_tid(nhru),hru_id(nhru)
 integer :: i,j,ntile,hru,pos,k,l,npos,inew,jnew,cluster,tid,hid,hillslope,ipos
 integer,allocatable,dimension(:,:) :: positions
 real,allocatable,dimension(:,:) :: tile_bwidth,tile_twidth,tile_length,tile_area
 real,allocatable,dimension(:,:) :: tile_dem,tile_slope,tile_position
 real :: dp(nhru),up(nhru),tile_pos
 npos = 8
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

 !Define some parameters
 ntile = maxval(tiles)
 
 !Allocate and initialize the tile arrays
 allocate(tile_bwidth(ntile,nc),tile_twidth(ntile,nc),tile_length(ntile,nc),&!
          tile_area(ntile,nc),tile_dem(ntile,nc),tile_slope(ntile,nc),&!
          tile_position(ntile,nc))
 tile_bwidth = 0.0
 tile_twidth = 0.0
 tile_length = 0.0
 tile_area = 0.0
 tile_dem = 0.0
 tile_slope = 0.0
 tile_position = 0.0

 !Initialize the properties for each tile
 hru_id = 0
 hru_hid = 0
 hru_tid = 0 
 
 !Iterate through each tile to compute the properties
 do i=1,nx
  do j=1,ny
   tid = tiles(i,j)
   hid = hillslopes(i,j)
   hru = hrus(i,j)
   if ((hid .gt. 0) .and. (tid .gt. 0))then
    !hillslope id
    hru_hid(hru) = hillslopes(i,j)
    !tile id
    hru_tid(hru) = tiles(i,j)
    !hru
    hru_id(hru) = hru
    !tile area
    tile_area(tid,hid) = tile_area(tid,hid) + 1
    !tile dem
    tile_dem(tid,hid) = tile_dem(tid,hid) + dem(i,j)
    !tile slope
    tile_slope(tid,hid) = tile_slope(tid,hid) + slope(i,j)
    !Add to downstream perimeter
    do pos=1,npos
     inew = i + positions(pos,1)
     jnew = j + positions(pos,2)
     if (((tid .gt. tiles(inew,jnew)) &
        .and. (hillslopes(inew,jnew) .eq. hid)) &
        .or. (channels(inew,jnew) .gt. 0)) then
       tile_bwidth(tid,hid) = tile_bwidth(tid,hid) + 1.0
       exit
     endif
    enddo
    !Add to upstream perimeter
    do pos=1,npos
     inew = i + positions(pos,1)
     jnew = j + positions(pos,2)
     if ((tid .lt. tiles(inew,jnew)) &
        .and. (hillslopes(inew,jnew) .eq. hid))then
       tile_twidth(tid,hid) = tile_twidth(tid,hid) + 1.0
       exit
     endif
    enddo
   endif
  enddo
 enddo

 do i = 1,ntile
  do j = 1,nc
   !Calculate the average hillslope properties
   tile_bwidth(i,j) = res*tile_bwidth(i,j)/nhillslope(j)
   tile_twidth(i,j) = res*tile_twidth(i,j)/nhillslope(j)
   tile_dem(i,j) = tile_dem(i,j)/tile_area(i,j)
   tile_slope(i,j) = tile_slope(i,j)/tile_area(i,j)
   tile_area(i,j) = res**2*tile_area(i,j)/nhillslope(j)
   !Where missing set it to the previous
   if (tile_twidth(i,j) .eq. 0.0) tile_twidth(i,j) = tile_twidth(i-1,j)
   if (tile_bwidth(i,j) .eq. 0.0) tile_bwidth(i,j) = tile_bwidth(i-1,j)
   !Figure out what to do with the undef lower and top width
   tile_length(i,j) = 2*tile_area(i,j)/(tile_bwidth(i,j) + tile_twidth(i,j))
  enddo
 enddo

 !Create the positions array
 do hid = 1,nc
  do tid = 1,ntile
   !Define the position
   if (tid .eq. 1) then
    tile_position(tid,hid) = tile_length(tid,hid)/2
   else
    tile_position(tid,hid) = tile_position(tid-1,hid) + &
    tile_length(tid,hid)/2 + tile_length(tid-1,hid)/2
   endif
  enddo
 enddo

 !Assign the info to each hru
 do i = 1,nhru
  tid = hru_tid(i)
  hid = hru_hid(i)
  hru_bwidth(i) = tile_bwidth(tid,hid)
  hru_twidth(i) = tile_twidth(tid,hid)
  hru_length(i) = tile_length(tid,hid)
  hru_area(i) = tile_area(tid,hid)
  hru_dem(i) = tile_dem(tid,hid)
  hru_slope(i) = tile_slope(tid,hid)
  hru_position(i) = tile_position(tid,hid)
 enddo

end subroutine
