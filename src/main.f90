program main

  use parameters
  use calculation, only : deg2rad
  use initialization
  use integration

  implicit none
  real(8) :: i_pixel, j_pixel
  real(8) :: r_screen, theta_screen, phi_screen
  real(8) :: t, r, theta, phi
  real(8) :: p_t, p_r, p_theta, p_phi
  real(8) :: E, L, Q

  real(8) :: tmp
  real(8) :: t_upd, r_upd, theta_upd, phi_upd, p_r_upd, p_theta_upd

  ! loop variants
  i_pixel = SCREEN_X
  j_pixel = SCREEN_Y

  r_screen = DISTANCE
  call deg2rad(ELEVATION, tmp)
  theta_screen = PI * 0.5d0 - tmp
  call deg2rad(AZIMUTH, phi_screen)

  call initialize(i_pixel, j_pixel, &
                  r_screen, theta_screen, phi_screen, &
                  t, r, theta, phi, &
                  p_t, p_r, p_theta, p_phi, &
                  E, L, Q)

  open(97, file="output2.csv", status="replace")
  do
    ! if ((r .gt. MAX_R) .or. (r .le. MIN_R)) then
    !   exit
    ! end if

    write(97, '(f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6)') &
          t, r, theta, phi, p_r, p_theta
    print *, t, r, theta, phi, p_r, p_theta, p_phi

    call null_geodesic_RK4(E, L, Q, t, r, theta, phi, p_r, p_theta, &
          t_upd, r_upd, theta_upd, phi_upd, p_r_upd, p_theta_upd)

    call update(t, r, theta, phi, p_r, p_theta, &
                t_upd, r_upd, theta_upd, phi_upd, p_r_upd, p_theta_upd)

    if ((r .gt. MAX_R) .or. (r .le. MIN_R)) then
      exit
    end if
  end do
  close(97)

end program main
