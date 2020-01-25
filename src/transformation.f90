module transformation

  implicit none
  public :: screen2cartesian, cartesian2boyerlindquist, initialize, update
  private

contains

!===============================================================================
  subroutine screen2cartesian(i, j, x, y, z)

    use parameters
    use calculation, only : deg2rad

    implicit none
    real(8), intent(in)  :: i, j
    real(8), intent(out) :: x, y, z

    real(8) :: theta, phi

    call deg2rad(ELEVATION, theta)
    theta = PI * 0.5d0 - theta
    call deg2rad(AZIMUTH, phi)

    x = -j * cos(theta) * cos(phi) - i * sin(phi) &
        + DISTANCE * sin(theta) * cos(phi)
    y = -j * cos(theta) * sin(phi) + i * cos(phi) &
        + DISTANCE * sin(theta) * sin(phi)
    z = DISTANCE * cos(theta) + j * sin(theta)

    return
  end subroutine screen2cartesian

!===============================================================================
  subroutine cartesian2boyerlindquist(x, y, z, r, theta, phi)

    use parameters

    implicit none
    real(8), intent(in)  :: x, y, z
    real(8), intent(out) :: r, theta, phi

    real(8) :: tmp

    tmp = x * x + y * y + z * z - BH_A * BH_A

    r     = sqrt((tmp + sqrt(tmp * tmp + 4.0d0 * BH_A * BH_A * z * z)) * 0.5d0)
    theta = acos(z / r)
    phi   = atan2(y, x)

    return
  end subroutine cartesian2boyerlindquist

!===============================================================================
  subroutine initialize(t, r, theta, phi, p_t, p_r, p_theta, p_phi, E, L, Q)

    use parameters
    use calculation, only : deg2rad, calc_Sigma, calc_Delta, calc_ELQ

    implicit none
    real(8), intent(in)  :: r, theta, phi
    real(8), intent(out) :: t, p_t, p_r, p_theta, p_phi, E, L, Q

    real(8) :: theta_screen, phi_screen
    real(8) :: tmp_r, tmp_phi
    real(8) :: sigma, delta
    real(8) :: dot_r, dot_theta, dot_phi

    call deg2rad(ELEVATION, theta_screen)
    theta_screen = PI * 0.5d0 - theta_screen
    call deg2rad(AZIMUTH, phi_screen)

    tmp_r   = sqrt(r * r + BH_A * BH_A)
    tmp_phi = phi - phi_screen

    call calc_Sigma(r, theta, sigma)
    call calc_Delta(r, delta)

    dot_r = -(r * tmp_r * sin(theta) * sin(theta_screen) * cos(tmp_phi) &
              + tmp_r * tmp_r * cos(theta) * cos(theta_screen)) / sigma
    dot_theta = (r * sin(theta) * cos(theta_screen) &
                - tmp_r * cos(theta) * sin(theta_screen) * cos(tmp_phi)) / sigma
    dot_phi = sin(theta_screen) * sin(tmp_phi) / (tmp_r * sin(theta))

    call calc_ELQ(r, theta, dot_r, dot_theta, dot_phi, E, L, Q)

    t       = 0.0d0
    p_t     = -E
    p_r     = dot_r * sigma / delta
    p_theta = dot_theta * sigma
    p_phi   = L

    return
  end subroutine initialize

!===============================================================================
  subroutine update(t, r, theta, phi, p_r, p_theta, &
                    t_upd, r_upd, theta_upd, phi_upd, p_r_upd, p_theta_upd)

    implicit none
    real(8), intent(in)  :: t_upd, r_upd, theta_upd, phi_upd, &
                            p_r_upd, p_theta_upd
    real(8), intent(out) :: t, r, theta, phi, p_r, p_theta

    t       = t_upd
    r       = r_upd
    theta   = theta_upd
    phi     = phi_upd
    p_r     = p_r_upd
    p_theta = p_theta_upd

    return
  end subroutine update

end module transformation
