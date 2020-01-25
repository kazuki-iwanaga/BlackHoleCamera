program main

  use integration
  use calculation
  use parameters
  use transformation

  implicit none
  real(8) :: i, j
  real(8) :: x, y, z
  real(8) :: t, r, theta, phi, p_t, p_r, p_theta, p_phi
  real(8) :: E, L, Q
  real(8) :: t_upd, r_upd, theta_upd, phi_upd, p_r_upd, p_theta_upd

  ! loop variants
  i = SCREEN_X
  j = SCREEN_Y

  call screen2cartesian(i, j, x, y, z)
  call cartesian2boyerlindquist(x, y, z, r, theta, phi)
  call initialize(t, r, theta, phi, p_t, p_r, p_theta, p_phi, E, L, Q)

  open(97, file="output.csv", status="replace")
  do
    if ((r .gt. MAX_R) .or. (r .le. MIN_R)) then
      exit
    end if

    write(97, '(f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6)') &
          t, r, theta, phi, p_r, p_theta, p_phi
    ! print *, t, r, theta, phi, p_r, p_theta, p_phi

    call null_geodesic_RK4(E, L, Q, t, r, theta, phi, p_r, p_theta, &
          t_upd, r_upd, theta_upd, phi_upd, p_r_upd, p_theta_upd)

    call update(t, r, theta, phi, p_r, p_theta, &
                t_upd, r_upd, theta_upd, phi_upd, p_r_upd, p_theta_upd)
  end do
  close(97)

end program main
