Module CSV
    use :: utils, only: replace

contains

  integer Function column(cols, name)
    ! return the index in cols where name appears
    implicit none
    character(*), intent(in) :: cols(:), name
    do column=1,size(cols)
       if (trim(cols(column)).eq.trim(name)) then
          return
       end if
    end do
    column = -1
    return
  end Function column


  subroutine read_csv_file(filename, table, columns)
    ! read a CSV file into table with column headers in columns
    
    implicit none
    character(*), intent(in) :: filename
    real(kind=8), allocatable, intent(out) :: table(:,:)
    character(len=:), allocatable, intent(out) :: columns(:)
    real(kind=8), allocatable :: row(:,:)
    character(1024) line
    integer infile, max_name_size, name_size, i, j, ncol, tlen

    open(newunit=infile, file=filename, status='old')
    read(infile, '(A)') line
    ! CSV is 'comma separated' so count commas to get number of columns
    ! also get the maximum size of a column name
    ncol = 0
    max_name_size = 0
    name_size = 0
    j = len_trim(line)
    do i=1,j
       name_size = name_size + 1
       if ( line(i:i) .eq. ',' ) then
          if ( max_name_size.lt.name_size ) max_name_size=name_size
          name_size = 0
          ncol = ncol + 1
       end if
    end do
    if ( line(j:j) .ne. ',' ) then
       if ( max_name_size.lt.name_size) max_name_size=name_size
       ncol = ncol + 1
    end if

    allocate(character(max_name_size) :: columns(ncol))

    read(line, *) columns

    ! now deetermine the size of the file to allocate the table array
    tlen = 0
    do
       read(infile, *, end=100)
       tlen = tlen + 1
    end do

100 continue
    allocate(row(ncol,tlen))
    rewind(infile)
    ! pop off the header line
    read(infile, *) line
    read(infile, *, end=200) row(:,:)
200 continue
    close(infile)

    !send it back as a table that is row,column arranged.
    table = transpose(row)
    return

  end subroutine read_csv_file

  subroutine reshape_channels(in, columns, out, out_columns)
    !> Make a table where the '0' and '1' columns are merged to a common column
    implicit none
    real(kind=8), intent(inout) :: in(:,:)
    character(len=*), intent(inout) :: columns(:)
    real(kind=8), allocatable :: out(:,:)
    character(len=:), allocatable :: out_columns(:)
    integer in_col_idx, out_col_idx, remap_count

    ! e.g. if in is 10x4 and columns is ['ra1', 'dec1', 'ra0', dec0'] then
    ! out will be 20x2 with columns ['ra0', 'dec0']
    ! count number of columns that have a '1' in their name, these will be mapped
    ! to new rows in the output table
    remap_count=0
    do in_col_idx=1,size(columns)
       if (INDEX(columns(in_col_idx),'1').gt.0) remap_count = remap_count+1
    enddo
    ! make the out array twice as long as in but remove the columns that match '1'

    allocate(out(size(in,1)*2, size(columns)-remap_count))
    allocate(character(len=len(columns(1))) :: out_columns(size(columns)-remap_count))

    out_col_idx=0
    do in_col_idx=1,size(columns)
       if (INDEX(columns(in_col_idx),'0').gt.0.or.INDEX(columns(in_col_idx),'1').eq.0) then
          out_col_idx=out_col_idx+1
          ! this is a '0' labelled column or a generic column
          ! make a new table with the '1' labelled column appended to the '0' ones
          out(:size(in,1),out_col_idx) = in(:,in_col_idx)
          out(size(in,1)+1:,out_col_idx) = in(:,column(columns,replace(columns(in_col_idx),'0','1')))
          out_columns(out_col_idx) = replace(columns(in_col_idx),'0','')
       end if
    enddo
    return
    
  end subroutine reshape_channels

end Module CSV

