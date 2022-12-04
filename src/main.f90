Program Main
  use, intrinsic :: iso_fortran_env, only: ou => error_unit
  use :: la_precision, only: wp => dp
  use utils
  use CSV
  use SortUnique,  only : unique
  use source_linker
  
  implicit none
  real(kind=wp), allocatable :: table(:,:), remap_table(:,:), plantlist(:,:)
  character(len=:), allocatable :: columns(:), filename, remap_columns(:), pl_columns(:)
  character(len=1024) value
  integer :: ra, dec, id, jd, rra, rdec, i, j, k, expnum, mag, length, status
  integer :: pl_id, pl_rate, pl_angle, pl_mag, rate_limit
  real(kind=wp) :: rlow, rhigh, med_ra, med_dec
  integer, allocatable :: object(:), plant(:), ids(:), links(:), link_length(:), idx(:)
  integer, allocatable :: primary(:), secondary(:), trial(:), expnums(:)
  real(kind=wp), allocatable :: rates(:), seps(:), rlink_length(:)
  real(kind=wp) r, angle, ra0, dec0, jd0, variance, plant_res(3)
  logical, allocatable :: cond(:)
  real(kind=wp), parameter :: &
       Pi = 3.141592653589793238d0, &! Pi
       drad = Pi/180.0d0             ! Degree to radian convertion: Pi/180

  ! get the name of the file to work on from the commandline
  CALL GET_COMMAND_ARGUMENT(1, value, length, status)
  if(status/=0)then
     write(ou,*) "Bad CL value ",value
     stop
  end if
  filename = value(1:length)
  write (ou,*)  "# Reading in ", filename
  call read_csv_file(filename, table, columns)
  write (ou,*) '# length of input table ',size(table,1),'x',size(table,2)
  call reshape_channels(table, columns, remap_table, remap_columns)
  table = remap_table
  columns = remap_columns
  write (ou,*) '# length of work table ',size(table,1),'x',size(table,2)
    ! get the name of the file to work on from the commandline

  CALL GET_COMMAND_ARGUMENT(2, value, length, status)
  if(status/=0)then
     write(ou,*) "Bad CL value ",value
     stop
  end if
  filename = value(1:length)
  write (ou,*) "# Reading in ", filename
  call read_csv_file(filename, plantlist, pl_columns)
  write (ou,*) '# length of input table ',size(plantlist,1),'x',size(plantlist,2)

  pl_id = column(pl_columns, 'id')
  pl_rate = column(pl_columns, 'rate')
  pl_angle = column(pl_columns, 'angle')
  pl_mag = column(pl_columns, 'mag')
  
  id = column(columns, 'kbo_id')
  mag = column(columns, 'pred_mag')
  ra = column(columns, 'real_ra')
  dec = column(columns, 'real_dec')
  jd = column(columns, 'time')
  expnum = column(columns, 'expnum')
  
  write(ou,*) '# Coverting inputs to radians'
  table(:, ra) = table(:, ra)*drad
  table(:, dec) = table(:, dec)*drad

  
  links = [(0, i=1,size(table,1))] ! will old value of object that source links to
  idx = [(i, i=1,size(table,1))] ! a convenience vector of indice values

  write(ou,*) '# merging observations of the same object on the same exposure'
  ! link sources that are the same object, per chip (separation between two measures less than 3")
  do j=1,size(table,1)
     ! when links for index is non zero then this source is part of an object
     if (links(j).ne.0) then
        cycle
     end if
     ! in the input catalog we can detect the same source on a single chip many times
     ! as the all expnum1-expnum2 combos are detected independently.
     ! find all the sources that are other measurs of this source in this expousre
     expnums = pack([idx],table(:,expnum).eq.table(j,expnum)) ! expnums is an index into all sources for this exposure
     seps = separation(table(j,ra),table(j,dec), table(expnums,ra), table(expnums,dec)) 
     cond = seps.lt.drad*3d0/3600d0.and.links(expnums).eq.0 ! considered close enough to be the same object
     expnums = pack([expnums],cond) ! select from index those sources that meet the condition
     if (size(expnums).lt.2) then
        ! no need to 'merge' with other measures.
        cycle
     endif
     ids = argsort(table(expnums,ra)) ! order the value to find the median.
     med_ra = table(expnums(ids(size(ids)/2)), ra)
     if (count(abs(table(expnums,ra)-med_ra).gt.5d0/3600d0/drad).gt.0) then ! scatter too large
        write(ou,*) j,table(expnums,ra)-med_ra
        stop
     endif
     table(expnums,ra) = med_ra
     ids = argsort(table(expnums,dec))
     med_dec = table(expnums(ids(size(ids)/2)), dec)
     table(expnums,dec) = med_dec
     links(expnums(2:)) = -1 ! a -ve here indicates that we can now ignore this source measure
  end do

  write(ou,*) "# Working with ",COUNT(links(:).eq.0)," source measurements"
  write(ou,*) '# Will now do search through the input table to find linkable sources'
  write(ou,'(A)') '# First Pass'
  write (ou,'(9A8)') "KBO_ID","N_OBS","RATE","ANGLE","STDEV","NCOR","RRATE","RANGLE","RMAG"

  idx = pack([idx],links.eq.0)
  do j=1,SIZE(idx)
     if (links(j).ne.0) then ! skip this source as it must linked to an object already
        cycle
     endif
     do rate_limit=1,10,2 ! search for possilbe sources that would be of this object
        object = pack([idx], links(idx).eq.0.and.table(idx,expnum).ne.table(index(j),expnum))
        rates = rate(table(idx(j), ra), table(idx(j), dec), table(idx(j), jd), &
             &table(object, ra), table(object, dec), table(object, jd))
        
        ! pick out sources that might link into an object
        rlow = (rate_limit * drad / 150d0)
        rhigh = (rate_limit + 4) * drad / 150d0
        cond = rates.lt.rhigh.and.rates.gt.rlow
        if (count(cond).le.3) then
           cycle
        endif
        ! look at sources within 'rate' of this source
        object = pack([object], cond)
        !write(ou,*) size(object)
        !write(ou,'(I0,2F16.10,F16.8)') (object(i), table(object(i),ra), table(object(i),dec), table(object(i),jd), i=1,size(object))
        !stop
        call skymotion_lsq(table(object, ra), table(object, dec), table(object, jd), &
             &r, angle, ra0, dec0, jd0, variance)

        ! accept all object measurements that are within 1 arcseconds of the
        ! lsq ra/dec rate/angle fits
        cond = (table(object, ra) - (ra0 + r * cos(angle) * (table(object, jd) - jd0)).lt.drad * (1. / 3600.0))
        cond = (cond.and.(table(object, dec) - (dec0 + r * sin(angle) * (table(object, jd) - jd0)).lt.drad * (1. / 3600.0)))

        if (count(cond).lt.3) then
           cycle
        endif

        object = pack([object], cond)

        call skymotion_lsq(table(object, ra), table(object, dec), table(object, jd), &
             &r, angle, ra0, dec0, jd0, variance)

        if (sqrt(variance) * 3600d0 / drad.gt.1.00) then
           ! fit not good enough
           cycle
        end if

        ! accept all object measurements that are within 1 arcseconds of the
        ! lsq ra/dec rate/angle fits
        cond = abs(table(object, ra) - (ra0 + r * cos(angle) * (table(object, jd) - jd0))).lt.drad * (1. / 3600.0)
        cond = (cond.and.abs(table(object, dec) - (dec0 + r * sin(angle) * (table(object, jd) - jd0))).lt.drad * (1. / 3600.0))
        object = pack([object], cond)
        if (count(cond).lt.3) then
           cycle
        endif

        links(object) = table(idx(j), j)
        write (ou, '(I8,2(I8,F8.2,F8.2,F8.2))') &
             &floor(table(object(1), id)), size(object), r * 150 / drad, angle / drad, sqrt(variance) * 3600d0 / drad, &
             &size(expnums), 0d0, 0d0, 0d0
     end do
  end do
  stop
  
  write(ou,'(A)') '# Second Pass'
  write (ou,'(9A8)') "KBO_ID","N_OBS","RATE","ANGLE","STDEV","NCOR","RRATE","RANGLE","RMAG"
  
  idx = pack([(i,i=1,size(links))], links.gt.0)
  do j=1,size(idx)
     primary = pack([(k,k=1,size(links)],links.eq.table(idx(j),id))
     call skymotion_lsq(table(object, ra), table(object, dec), table(object, jd), &
          &r, angle, ra0, dec0, jd0, variance)

     do i=size(idx),j+1,-1
        secondary = pack([(k,k=1,size(links)],links.eq.table(idx(i),id))
        cond = abs(table(secondary,ra) - (ra0 + r*cos(angle)*(table(secondary,jd)-jd0))).lt.drad*(1./3600.0)
        cond = cond.and.abs(table(secondary,dec) - (dec0 + r*sin(angle)*(table(secondary,jd)-jd0))).lt.drad*(1./3600.0)
        ! Assign objects that meet 'cond' to primary.
        links(pack([secondary],cond),2)=ids(j)
     enddo
  enddo
  
  ids = pack([links(:,2)],links(:,2).gt.0)
  if (size(ids).gt.1) then
     ids = unique(ids)
     ids = ids(:size(ids)-1)
  endif

  do j=1,size(ids)
     object = pack([links(:,1)],links(:,2).eq.ids(j))
     if (size(object).lt.3) then
        cycle
     endif
     if (table(object(1),id).gt.0) then
        plant = pack([(i,i=1,size(plantlist,1))],plantlist(:,pl_id).eq.table(object(1),id))
        plant_res = [plantlist(plant,pl_rate), plantlist(plant,pl_angle), plantlist(plant,pl_mag)]
     else
         plant_res = [0d0, 0d0, 0d0]
     end if
     call skymotion_lsq(table(object,ra), table(object,dec), table(object,jd),&
          &r, angle, ra0, dec0, jd0, variance)

     expnums = floor(table(object,id))
     if (size(expnums).gt.1) then
        expnums = unique(expnums)
        if (size(expnums).gt.1) then
            expnums = expnums(:size(expnums) - 1)
        end if
     endif
     write (ou, '(I8,2(I8,F8.2,F8.2,F8.2))') &
          &floor(table(object(1),id)), size(object), r*150/drad, angle/drad, sqrt(variance)*3600d0/drad,&
          &SIZE(expnums),&
             plant_res(1), plant_res(1), plant_res(3)

      do i=1,size(object)
          write(ou,*) '---> ',table(object(i),id), table(object(i), ra), table(object(i),dec), table(object(i),expnum)
      end do
  end do
        
end program main
