module source_linker
  use :: la_precision, only: wp => dp
  use, intrinsic :: iso_fortran_env, only: ou => error_unit
  use :: f95_lapack, only: la_gels
  implicit none
  private

  public :: rate, separation, skymotion_lsq, source_merge

contains

  function angle(a1, a2)
    ! compute the angular difference between two degree angles
    implicit none
    real(kind=wp), intent(in) :: a1, a2(:)
    real(kind=wp), allocatable :: angle(:)
    real(kind=wp), allocatable ::  da(:)

    da = (a1-a2)
    angle = DATAN2(dsin(da), dcos(da))
    return 
  end function angle

  function separation(ra1, dec1, ra2, dec2)
    implicit none
    real(kind=wp), allocatable :: separation(:), sep(:)
    real(kind=wp), intent(in) :: ra1, dec1, ra2(:), dec2(:)
    integer, allocatable :: idx(:)
    integer :: i
    sep = dsin(dec1)*dsin(dec2) + dcos(dec1)*dcos(dec2)*dcos(ra2-ra1)
    idx = pack([(i,i=1,size(sep))], sep.gt.1d0)
    sep(idx) = 1d0
    idx = pack([(i,i=1,size(sep))], sep.lt.-1d0)
    sep(idx) = -1d0
    separation = dacos(sep)
    return
  end function separation

  function rate(ra1, dec1, jd1, ra2, dec2, jd2)
    implicit none
    real(kind=wp), allocatable :: rate(:)
    real(kind=wp), intent(in) :: ra1, dec1, jd1
    real(kind=wp), intent(in) :: ra2(:), dec2(:), jd2(:)
    real(kind=wp), allocatable :: dt(:)

    dt = jd2
    WHERE ( jd2 < jd1 ) dt = jd1 - dt
    WHERE ( jd2 > jd1 ) dt = dt - jd1

    rate = separation(ra1, dec1, ra2, dec2)/dt
    return
  end function rate

  subroutine skymotion_lsq(ra, dec, jd, rate, angle, ra0, dec0, jd0, variance)
    !> use LAPACK least Squares analysis to determine the most likely rate/angle
    !> given the source measurements associated to the object
    implicit none
    real(kind=wp), intent(in) :: ra(:), dec(:), jd(:)
    real(kind=wp), intent(out) :: rate, angle, ra0, dec0, jd0, variance

    integer :: N, NRHS  ! number of values in optimization, number of parameters per value
    real(kind=wp), allocatable :: a(:, :)
    real(kind=wp), allocatable :: l(:), dt(:)
    integer                   :: i

    !
    ! We are determining the parameters to these two coupled equations
    ! ra  = ra_0  + rate*cos(angle)*(jd-jd0)
    ! dec = dec_0 + rate*sin(angle)*(jd-jd0)
    ! we have freedom to set jd0 to an arbitrary value, so use the middle
    ! of the sequence
    !
    ! Using LAPACK so need to express this as matrix
    !
    !> Matrix solution for rates of motion looks like this
    ! / dt 0 1 0 \   / rate*cos(angle) \  / ra_1 \
    ! | dt 0 1 0 |   | rate*sin(angle) |  | ra_2 |
    ! | .        | X | ra_0            |  |  .   |
    ! | .        |   | dec_0           | =|  .   |
    ! | 0 dt 0 1 |   |                 |  | de_1 |
    ! | 0 dt 0 1 |   |                 |  | de_2 |
    ! | .        |   |                 |  |.     |
    ! \ .        /   \                 /  \.     /
    ! based on example from: https://cyber.dabamos.de/programming/modernfortran/lapack.html

    if ( abs(maxval(jd)-minval(jd)) .lt. 1d0/24d0 .or. size(ra) .lt. 3 ) then
       rate = 0d0
       angle= 0d0
       ra0  = 0d0
       dec0 = 0d0
       jd0  = 0d0
       variance = 1d6
       return
    end if

    ! Allocate arrays.
    N = size(ra)+size(dec)
    NRHS = 4
    allocate (a(N, NRHS), l(N), dt(N) )

    ! set the time as relative to the middle observation
    ! Control points in source (t)
    jd0 = jd(size(jd)/2)
    dt = jd - jd0

    ! Fill coefficient matrix A.
    do i=1,size(ra)
       a(i,:) = [ dt(i), 0d0, 1d0, 0d0 ]
       a(i+size(ra),:) = [ 0d0, dt(i), 0d0, 1d0] 
    end do

    ! Control points in target system.
    ! contatinate the ra and dec lists
    l = [ra, dec]
    ! Get transformation parameters `rate`, `angle`
    call solve(a, l, rate, angle, ra0, dec0, variance)

    return
  end subroutine skymotion_lsq

  subroutine solve(a, l, rate, angle, ra0, dec0, variance)
    implicit none
    !! Solves AX = L, then calculates rotation angle and scale factor.
    real(kind=wp), intent(inout) :: a(:, :)  ! Input:  Coefficient matrix A.
    real(kind=wp), intent(inout) :: l(:)     ! Input:  Observations.
    real(kind=wp), intent(out) :: rate, angle, ra0, dec0, variance
    real(kind=wp), allocatable :: x(:), v(:)     ! Output: Solution [ rate*cos(a), rate*sin(a), ra0, dec0 ]
    real(kind=wp), allocatable   :: wa(:, :) ! Work array of `a`.
    real(kind=wp), allocatable   :: wl(:)    ! Work array of `l`.
    
    allocate (wa(size(a, 1), size(a, 2)), wl(size(l)))
    allocate (x(size(a,2)))
    wa = a
    wl = l
    
    ! Compute minimum-norm least squares solution to AX = L using LAPACK95.
    ! Work array `wl` contains both the solution and the residuals afterwards.
    ! The solution is stored in `wl(1:size(a, 2))`, and the sum of squares
    ! of the residuals in `wl(size(a, 2) + 1:size(a, 1))`.
    call la_gels(wa, wl)
    
    x = wl(1:size(a, 2))  ! Solution [ rate*cos(angle), rate*sin(angle), ra0, dec0 ].

    ! Calculate residuals of observations.
    v = matmul(a, x) - l
    
    ! Calculate adjustment's reference variance, given the
    ! N observations and NRHS degrees of freedom.
    variance = dot_product(v, v) / (size(a, 1) - size(a, 2))
    angle = atan2(x(2), x(1)) ! Calculate angle
    rate = x(1) / cos(angle)  ! rate
    ra0 = x(3)
    dec0 = x(4)

    return
  end subroutine solve

  function source_merge(table, ra, dec, expnum)
    use utils
    ! merge sources with same ra,dec,expnum
    real(kind=wp), allocatable :: source_merge(:,:)
    real(kind=wp), intent(inout) :: table(:,:)
    integer, intent(in) :: ra, dec, expnum ! column numbers
    real(kind=wp), allocatable :: seps(:)
    logical, allocatable :: cond(:)
    integer, allocatable :: idx(:), links(:), expnums(:), ids(:)
    real(kind=wp) med_ra, med_dec
    integer :: i, j, k
    
    real(kind=wp), parameter :: &
         Pi = 3.141592653589793238d0, &! Pi
         drad = Pi/180.0d0,           &! Degree to radian convertion: Pi/180
         match = 3.5d0*drad/3600d0
    
    links = [(0, i=1,size(table,1))] ! will old value of object that source links to
    idx = [(i, i=1,size(table,1))] ! a convenience vector of indice values

    write(ou,*) '# merging observations of the nearby souces on the same exposure'
    ! link sources that are the same object, per chip (separation between two measures less than 3")
    do j=1,size(table,1)
       ! when links for index is non zero then this source is part of an object
       if (links(j).ne.0) then
          cycle
       end if
       ! in the input catalog we can detect the same source on a single chip many times
       ! as the all expnum1-expnum2 combos are detected independently.
       ! find all the sources that are other measurs of this source in this expousre
       expnums = pack([idx],floor(table(:,expnum)).eq.floor(table(j,expnum))) ! expnums is an index into all sources for this exposure
       seps = separation(table(j,ra),table(j,dec), table(expnums,ra), table(expnums,dec))

       cond = seps.lt.match.and.links(expnums).eq.0 ! considered close enough to be the same object
       expnums = pack([expnums],cond) ! select from index those sources that meet the condition
       if (size(expnums).lt.2) then
          ! no need to 'merge' with other measures.
          cycle
       endif
       if (size(expnums).gt.4) then
          ids = argsort(table(expnums,ra)) ! order the value to find the median.
          med_ra = table(expnums(ids(size(ids)/2)), ra)
          ids = argsort(table(expnums,dec))
          med_dec = table(expnums(ids(size(ids)/2)), dec)
       else
          med_ra = SUM(table(expnums,ra))/SIZE(expnums)
          med_dec = SUM(table(expnums,dec))/SIZE(expnums)
       end if
       if (count(abs(table(expnums,ra)-med_ra).gt.match).gt.(1+floor(0.32*size(expnums)))) then ! scatter too large
          write(ou,*) j,size(expnums),match,(table(expnums,ra)-med_ra)*3600d0/drad
          cycle
       endif
       if (count(abs(table(expnums,dec)-med_dec).gt.match).gt.(1+floor(0.32*size(expnums)))) then ! scatter too large
          write(ou,*) j,size(expnums),match,(table(expnums,dec)-med_dec)*3600d0/drad
          cycle
       endif
       table(expnums,ra) = med_ra
       table(expnums,dec) = med_dec
       links(expnums(2:)) = -1 ! a -ve here indicates that we can now ignore this source measure
    end do
    
    write(ou,*) "# Working with ",COUNT(links(:).eq.0)," source measurements"
    !write(ou,*) '# Will now do search through the input table to find linkable sources'
    !write(ou,'(A)') '# First Pass'
    !write (ou,'(9A8)') "KBO_ID","N_OBS","RATE","ANGLE","STDEV","NCOR","RRATE","RANGLE","RMAG"
    
    idx = pack([idx],links.eq.0)
    source_merge = table(idx, :)
    return 
  end function source_merge

end module source_linker
