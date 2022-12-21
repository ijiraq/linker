Program Main
  use, intrinsic :: iso_fortran_env, only: ou => error_unit
  use :: la_precision, only: wp => dp
  use utils
  use CSV
  use SortUnique,  only : unique
  use source_linker
  
  implicit none
  real(kind=wp), allocatable :: table(:,:), remap_table(:,:), plantlist(:,:)
  character(len=:), allocatable :: columns(:), filename, remap_columns(:), pl_columns(:), out_filename
  character(len=1024) value
  integer :: ra, dec, id, jd, rra, rdec, i, j, k, expnum, mag, length, status, likelihood
  integer :: pl_id, pl_rate, pl_angle, pl_mag, rate_limit, obj_num, obj_num_pp, unit
  integer :: expnum0, expnum1, nexp, kbo_id
  real(kind=wp) :: rlow, rhigh, med_ra, med_dec
  integer, allocatable :: object(:), plant(:), ids(:), links(:), link_length(:), idx(:)
  integer, allocatable :: primary(:), secondary(:), trial(:), expnums(:), kboids(:)
  real(kind=wp), allocatable :: rates(:), seps(:), rlink_length(:), rand_expnums(:)
  real(kind=wp) r, angle, ra0, dec0, jd0, variance, plant_res(3)
  logical, allocatable :: cond(:), cond0(:), cond1(:)
  real(kind=wp), parameter :: &
       Pi = 3.141592653589793238d0, &! Pi
       drad = Pi/180.0d0,           &! Degree to radian convertion: Pi/180
       match = drad*3.0d0/3600.0d0   ! this is a matched source

  ! get the name of the file to work on from the commandline
  CALL GET_COMMAND_ARGUMENT(1, value, length, status)
  if(status/=0)then
     write(ou,*) "Bad CL value ",value
     stop
  end if
  filename = value(1:length)
  write (ou,*)  "# Reading in ", filename
  write (ou,*)  "# Writing results to ", out_filename
  call read_csv_file(filename, table, columns)
  write (ou,*) '# length of input table ',size(table,1),'x',size(table,2)

  CALL GET_COMMAND_ARGUMENT(3, value, length, status)
  read(value,*) nexp
  out_filename = filename(:len_trim(filename))//"_"//value(:len_trim(value))//".linked"

  expnums = [(i,i=1,44,44/nexp)]
  expnums = expnums(1:nexp)
  ! only use detections from exposures 1, 22, 43
  expnum0 = column(columns, 'expnum0')
  expnum1 = column(columns, 'expnum1')

  allocate(cond0(size(table,1)), cond1(size(table,1)))
  cond0 = .false.
  cond1 = .false.
  do i=1,size(expnums)
     cond0 = cond0.or.floor(table(:,expnum0)).eq.expnums(i)
     cond1 = cond1.or.floor(table(:,expnum1)).eq.expnums(i)
  end do
  cond = cond0.and.cond1
  object = pack([(i,i=1,size(table,1))],cond)
  table = table(object,:)
  write(ou,*) "# only consider 4 exposures ", size(table,1)

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
  ra = column(columns, 'ra')
  dec = column(columns, 'dec')
  jd = column(columns, 'time')
  expnum = column(columns, 'expnum')
  likelihood = column(columns, 'p')
  obj_num = column(columns, 'original_num')
  obj_num_pp = column(columns, 'original_num_among_positively_predicted')
  
  write(ou,*) '# Coverting inputs to radians'
  table(:, ra) = table(:, ra)*drad
  table(:, dec) = table(:, dec)*drad

  ! only accept objects with probabilty from detection model of 0.99 or higher
  object = pack([(i,i=1,size(table,1))], table(:,likelihood).gt.0.99)
  table = table(object,:)
  write(ou,*) '# only keep p>0.99 sources ',size(table,1)
  
  table = source_merge(table, ra, dec, expnum)
  idx = [(i,i=1,size(table,1))]
  links = 0*idx
  rlink_length = 0d0*idx
  
  do j=1,SIZE(idx)
     if (links(j).ne.0) then ! skip this source as it must linked to an object already
        cycle
     endif
     do rate_limit=1,10,2 ! search for possilbe sources that would be of this object
        object = pack([idx], links.eq.0)
        if (size(object).lt.3) then
           cycle
        end if
        rates = rate(table(idx(j), ra), table(idx(j), dec), table(idx(j), jd), &
             &table(object, ra), table(object, dec), table(object, jd))

        
        ! pick out sources that might link into an object require at least 1 hour arc
        rlow = (rate_limit * drad / 150d0)
        rhigh = (rate_limit + 10) * drad / 150d0
        cond = rates.lt.rhigh.and.rates.gt.rlow
        ! and always include the 'root' source measure.
        cond(findloc(object,j))=.true.
        object = pack([object], cond)
        expnums = table(object,expnum)
        expnums = unique(expnums)
        if (size(expnums).lt.3) then
           cycle
        endif
        ! look at sources within 'rate' of this source

        !write(ou,*) size(object)
        !write(ou,'(I0,2F16.10,F16.8)') (object(i), table(object(i),ra), table(object(i),dec), table(object(i),jd), i=1,size(object))
        !stop
        call skymotion_lsq(table(object, ra), table(object, dec), table(object, jd), &
             &r, angle, ra0, dec0, jd0, variance)

        ! accept all object measurements that are within 1 arcseconds of the
        ! lsq ra/dec rate/angle fits

        cond = (abs(table(object, ra) - (ra0 + r * cos(angle) * (table(object, jd) - jd0))).lt.match)
        cond = (cond.and.(abs(table(object, dec) - (dec0 + r * sin(angle) * (table(object, jd) - jd0))).lt.match))

        object = pack([object], cond)
        expnums = table(object,expnum)
        expnums = unique(expnums)
        if (size(expnums).lt.3) then
           cycle
        endif

        call skymotion_lsq(table(object, ra), table(object, dec), table(object, jd), &
             &r, angle, ra0, dec0, jd0, variance)

        if (sqrt(variance) * 3600d0 / drad.gt.1.00) then
           ! fit not good enough
           cycle
        end if

        ! accept all object measurements that are within 1 arcseconds of fit and aren't linked to an object
        object = pack([idx],links.eq.0)
        cond = abs(table(object, ra) - (ra0 + r * cos(angle) * (table(object, jd) - jd0))).lt.match
        cond = cond.and.(abs(table(object, dec) - (dec0 + r * sin(angle) * (table(object, jd) - jd0))).lt.match)
        object = pack([object], cond)
        
        expnums = table(object,expnum)
        expnums = unique(expnums)
        if (size(expnums).lt.3) then
           cycle
        endif

        links(object) = object(1)
        ids = table(object,id)
        ids = unique(ids)
        rlink_length(object(1)) = size(expnums)
        !write (ou, '(I8,2(I8,F8.2,F8.2,F8.2))') &
        !     &floor(table(object(1), id)), size(expnums), r * 150 / drad, angle / drad, sqrt(variance) * 3600d0 / drad, &
        !     &size(object), 0d0, 0d0, 0d0
        exit
     end do
  end do
  
  
  ids = pack([(i,i=1,size(links))],links.eq.0)
  links(ids) = ids
  ids = argsort(rlink_length)
  cond = (rlink_length.gt.0)
  
  write(ou,'(A)') '# Second Pass'
  write(ou, *) '# lookig at ', count(cond), ' link groups'
  write (ou,'(12A12)') "ID", "OBJNUM", "OBJNUMPP", "KBO_ID","N_OBS","RATE","ANGLE","STDEV","NCOR","RRATE","RANGLE","RMAG"

  ! loop over all the linked objects, starting with the ones with the longest link_length
  ! and check to seeif of ther sources are part of this object, starting with the shortest
  ! link_length
  do j=1,size(ids)
     if (rlink_length(ids(j)).lt.3) then
        cycle
     endif
     ! get the source id values for all the sources in object ids(j)
     primary = pack([idx],links.eq.ids(j))
     expnums = table(primary,expnum)
     expnums = unique(expnums)
     if (size(expnums).lt.3) then
        links(primary) = 0
        cycle
     endif

     call skymotion_lsq(table(primary, ra), table(primary, dec), table(primary, jd), &
          &r, angle, ra0, dec0, jd0, variance)

     do i=size(ids),j+1,-1
        if (rlink_length(ids(i)).lt.1) then
           cycle
        endif
        ! get the source id values for all the sources in object ids(i)
        secondary = pack([idx],links.eq.ids(i))
        cond = abs(table(secondary,ra) - (ra0 + r*cos(angle)*(table(secondary,jd)-jd0))).lt.match
        cond = cond.and.(abs(table(secondary,dec) - (dec0 + r*sin(angle)*(table(secondary,jd)-jd0))).lt.match)
        ! Assign sources that were part of object ids(i) to object ids(j) if cond is true
        links(pack([secondary], cond))=ids(j)
        rlink_length(ids(j)) = rlink_length(ids(j)) + COUNT(cond)
        rlink_length(ids(i)) = rlink_length(ids(i)) - COUNT(cond)
     enddo
  enddo

  ! write newly linked object information to output file
  open(newunit=unit,file=out_filename)

  ids = argsort(rlink_length)
  ! write (ou,'(A)') "# looking for link to artificial source and refining fit"
  do j=1,size(ids)
     object = pack([idx],links.eq.ids(j))
     if (size(object).lt.3) then
        cycle
     endif
     if (table(object(1),id).gt.0) then
        plant = pack([(i,i=1,size(plantlist,1))],plantlist(:,pl_id).eq.table(object(1),id))
        plant_res = [plantlist(plant,pl_rate), plantlist(plant,pl_angle), plantlist(plant,pl_mag)]
     else
         plant_res = [0d0, 0d0, 0d0]
     end if
     expnums = table(object,expnum)
     expnums = unique(expnums)
     if (size(expnums).lt.3) then
        links(object) = 0
        cycle
     endif
     call skymotion_lsq(table(object,ra), table(object,dec), table(object,jd),&
          &r, angle, ra0, dec0, jd0, variance)
     if (variance*3600d0/drad.gt.5) then
        cycle
     endif
     kboids = floor(table(object,id))
     if (size(kboids).gt.1) then
        kboids = unique(kboids)
        if (size(kboids).gt.1) then
            kboids = expnums(:size(kboids) - 1)
        end if
     endif
     
     expnums = pack([object],table(object,id).gt.0)
     if (size(expnums).gt.0) then
        kbo_id = floor(table(expnums(1),id))
     else
        kbo_id = -1
     end if
     expnums = table(object, expnum)
     expnums = unique(expnums)
     write(ou, '(4I10,2(I10,F12.2,F12.2,F12.2))') &
          object(1), floor(table(object(1),obj_num)), floor(table(object(1),obj_num_pp)),&
          kbo_id,&
          size(expnums), r*150/drad, 180+angle/drad, sqrt(variance)*3600d0/drad,&
          SIZE(kboids), plant_res(1), plant_res(1), plant_res(3)
     expnums = argsort(table(object, expnum))
     do i=1,size(expnums)
        k=object(expnums(i))
        write(unit, '(2F12.5,F16.8,I8,2F8.2,5I8)') &
             table(k, ra)/drad, table(k, dec)/drad,&
             table(k, jd), floor(table(k, expnum)), table(k, likelihood), &
             table(k, mag), floor(table(k,id)), &
             object(1), floor(table(k, obj_num)), floor(table(k, obj_num_pp)), floor(table(k, size(table,2)))
     end do
     !do i=1,size(object)
     !    write(ou,*) '---> ',table(object(i),id), table(object(i), ra), table(object(i),dec), table(object(i),expnum)
     !end do
  end do
  close(unit)
end program main
