module source_linker
  use :: la_precision, only: wp => dp
  use :: f95_lapack, only: la_gels
  implicit none
  private

  public :: rate, separation, skymotion_lsq

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
    real(kind=wp), allocatable :: separation(:)
    real(kind=wp), intent(in) :: ra1, dec1, ra2(:), dec2(:)
    separation = dacos(dsin(dec1)*dsin(dec2) + dcos(dec1)*dcos(dec2)*dcos((ra2-ra1)))
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


end module source_linker
