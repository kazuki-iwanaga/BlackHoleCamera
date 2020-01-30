program black_hole_camera

  use parameters
  use calculation,    only : deg2rad, calc_energy, calc_angular_momentum, &
                             calc_carter_constant, &
                             calc_sigma, calc_delta
  use transformation, only : screen2bhcartesian, bhcartesian2boyerlindquist, &
                             calc_direction
  use initialization, only : initialize_4position, initialize_4momentum, &
                             initialize_ray
  use integration,    only : null_geodesic_RK4, update_ray

  implicit none
  double precision :: r_screen, theta_screen, phi_screen
  double precision :: i_pixel, j_pixel
  double precision :: x_bh, y_bh, z_bh
  double precision :: r_by, theta_by, phi_by
  double precision :: dot_r_by, dot_theta_by, dot_phi_by
  double precision :: E, L, Q
  double precision :: four_position_by(4), four_momentum_by(4)
  double precision :: ray(6), ray_update(6)
  character :: linebuf*128
  integer :: i, j
  double precision :: minimum, maximum
  double precision :: gass_range
  double precision :: intensity
  double precision :: x, y, z, tmp
  double precision :: dr, dtheta, dphi, ds
  double precision :: delta, sigma

  open(97, file="output.csv", status="replace")

  ! ...
  ! Set screen Position
  ! <-
  r_screen     = DISTANCE
  theta_screen = PI*0.5d0 - deg2rad(ELEVATION)
  ! theta_screen = PI*0.5d0 - deg2rad(30.0d0)
  phi_screen   = deg2rad(AZIMUTH)
  ! ->

  ! ...
  ! Select pixels
  ! <-
  do i = 1, 100
  ! do i = 1, 51
  do j = 1, 100
  ! i_pixel = X_PIXEL
  ! i_pixel = -15.0d0 + 0.6d0*dble(i-1)
  i_pixel = -30.0d0 + 0.6d0*dble(i-1)
  ! i_pixel = -8.75d0 + 0.5d0*dble(i-1)
  ! j_pixel = Y_PIXEL
  j_pixel = -30.0d0 + 0.6d0*dble(j-1)
  ! ->
  ! print *, i_pixel, j_pixel, r_screen, theta_screen, phi_screen

  minimum = sqrt(MIN_R*MIN_R - BH_A*BH_A)
  maximum = sqrt(MAX_R*MAX_R - BH_A*BH_A)
  ! gass_range = sqrt(4.0d0*4.0d0 - BH_A*BH_A)
  gass_range = 10.0d0
  intensity = 0.0d0

  ! ...
  ! Transform coordinates
  ! <-
  call screen2bhcartesian(i_pixel, j_pixel, &
                          r_screen, theta_screen, phi_screen, &
                          x_bh, y_bh, z_bh)

  call bhcartesian2boyerlindquist(x_bh, y_bh, z_bh, &
                                  r_by, theta_by, phi_by)
  ! ->
  ! print *, x_bh, y_bh, z_bh
  ! print *, r_by, theta_by, phi_by
  ! print *, r_by*cos(phi_by), r_by*sin(phi_by)

  ! ...
  ! Calculate light directions, conserved quantities(E, L, Q)
  ! <-
  call calc_direction(r_by, theta_by, phi_by, &
                      theta_screen, phi_screen, &
                      dot_r_by, dot_theta_by, dot_phi_by)

  E = calc_energy(r_by, theta_by, dot_r_by, dot_theta_by, dot_phi_by)
  L = calc_angular_momentum(r_by, theta_by, &
                            dot_r_by, dot_theta_by, dot_phi_by)
  Q = calc_carter_constant(r_by, theta_by, dot_r_by, dot_theta_by, dot_phi_by)
  ! ->
  ! print *, dot_r_by, dot_theta_by, dot_phi_by
  ! print *, E, L, Q

  ! ...
  ! Set 4-position and 4-momentum
  ! <-
  call initialize_4position(r_by, theta_by, phi_by, four_position_by)
  call initialize_4momentum(r_by, theta_by, dot_r_by, dot_theta_by, &
                            E, L, four_momentum_by)
  call initialize_ray(four_position_by, four_momentum_by, ray)
  ! ->
  ! print *, four_position_by
  ! print *, four_momentum_by
  ! print *, ray

  ! open(97, file="output.csv", status="replace")
  do
    if (ray(2) .le. MIN_R) then
      intensity = 0.0d0
      exit
    elseif (ray(2) .gt. MAX_R) then
      exit
    end if

    ! ...
    ! Integrate geodesic equations
    ! <-
    call null_geodesic_RK4(E, L, Q, ray, ray_update)
    call update_ray(ray_update, ray)
    ! ->

    ! if (ray(4) .gt. 2.0d0*PI) then
    !   ray(4) = ray(4) - 2.0d0*PI
    ! elseif ( ray(4) .lt. -2.0d0*PI ) then
    !   ray(4) = ray(4) + 2.0d0*PI
    ! end if

    tmp = ray(2)
    x   = tmp*sin(ray(3))*cos(ray(4))
    y   = tmp*sin(ray(3))*sin(ray(4))
    z   = ray(2)*cos(ray(3))
    tmp = sqrt(x*x + y*y)

    if (ray(2) .lt. 60.0d0) then
      ! if (ray(2) .lt. gass_range) then
      if ((tmp .lt. 40.0d0) .and. (abs(z) .lt. 0.5d0)) then

        ! write(97, '(f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6)') &
        !   ray
        ! write (linebuf, &
        !       '(f15.5,",",f15.5,",",f15.5)') &
        !       ray(2), ray(3), ray(4)
        ! ! write (linebuf, &
        ! !       '(i4,",",i4,",",f15.5)') &
        ! !       i_pixel, j_pixel, intensity
        ! call del_spaces(linebuf)
        ! write (97, '(a)') trim(linebuf)

        ! sigma = calc_sigma(ray(2), ray(3))
        ! delta = calc_delta(ray(2))

        ! dr     = ray_update(2) - ray(2)
        ! dtheta = ray_update(3) - ray(3)
        ! dphi   = ray_update(4) - ray(4)

        ! tmp = dr*dr*sigma/delta + dtheta*dtheta*sigma + dphi*dphi &
        !       *sin(ray(3))*sin(ray(3))*(ray(2)*ray(2) + BH_A*BH_A &
        !       + 2.0d0*BH_A*BH_A*BH_M*ray(2)*sin(ray(3))*sin(ray(3))/sigma)
        ! ds = sqrt(abs(tmp))

        ! intensity = intensity + ds
          intensity = intensity + 1.0d0

      end if
    end if

    ! print *, ray
  end do

  write (linebuf, &
        '(i4,",",i4,",",f15.5)') &
        i, j, intensity
  call del_spaces(linebuf)
  write (97, '(a)') trim(linebuf)

  print *, i, j
  end do ! j
  end do ! i

  close(97)

contains
  subroutine del_spaces(s)
    character (*), intent (inout) :: s
    character (len=len(s)) tmp
    integer i, j
    j = 1
    do i = 1, len(s)
      if (s(i:i)==' ') cycle
      tmp(j:j) = s(i:i)
      j = j + 1
    end do
    s = tmp(1:j-1)
  end subroutine del_spaces

end program black_hole_camera
