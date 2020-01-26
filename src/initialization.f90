module initialization

  implicit none
  public  :: initialize
  private :: screen2bhcartesian, bhcartesian2boyerlindquist, &
             calc_direction

contains

!===============================================================================
  subroutine initialize(i_pixel, j_pixel, &
                        r_screen, theta_screen, phi_screen, &
                        t_by, r_by, theta_by, phi_by, &
                        p_t_by, p_r_by, p_theta_by, p_phi_by, &
                        E, L, Q)

    implicit none
    real(8), intent(in)  :: i_pixel, j_pixel
    real(8), intent(in)  :: r_screen, theta_screen, phi_screen
    real(8), intent(out) :: t_by, r_by, theta_by, phi_by
    real(8), intent(out) :: p_t_by, p_r_by, p_theta_by, p_phi_by
    real(8), intent(out) :: E, L, Q

    real(8) :: x_bh, y_bh, z_bh

    call screen2bhcartesian(i_pixel, j_pixel, &
                            r_screen, theta_screen, phi_screen, &
                            x_bh, y_bh, z_bh)

    call bhcartesian2boyerlindquist(x_bh, y_bh, z_bh, &
                                    r_by, theta_by, phi_by)

    call calc_direction(r_by, theta_by, phi_by, &
                        theta_screen, phi_screen, &
                        p_r_by, p_theta_by, &
                        E, L, Q)

    t_by     = 0.0d0
    p_t_by   = -E
    p_phi_by = L

    return
  end subroutine initialize

!===============================================================================
  subroutine screen2bhcartesian(i_pixel, j_pixel, &
                                r_screen, theta_screen, phi_screen, &
                                x_bh, y_bh, z_bh)

    use parameters

    implicit none
    real(8), intent(in)  :: i_pixel, j_pixel
    real(8), intent(in)  :: r_screen, theta_screen, phi_screen
    real(8), intent(out) :: x_bh, y_bh, z_bh

    real(8) :: tmp

    tmp = sin(theta_screen)*sqrt(r_screen*r_screen + BH_A*BH_A) &
          - j_pixel*cos(theta_screen)

    x_bh = tmp*cos(phi_screen) - i_pixel*sin(phi_screen)
    y_bh = tmp*sin(phi_screen) + i_pixel*cos(phi_screen)
    z_bh = r_screen*cos(theta_screen) + j_pixel*sin(theta_screen)

    return
  end subroutine screen2bhcartesian

!===============================================================================
  subroutine bhcartesian2boyerlindquist(x_bh, y_bh, z_bh, &
                                        r_by, theta_by, phi_by)

    use parameters

    implicit none
    real(8), intent(in)  :: x_bh, y_bh, z_bh
    real(8), intent(out) :: r_by, theta_by, phi_by

    real(8) :: tmp

    tmp = x_bh*x_bh + y_bh*y_bh + z_bh*z_bh - BH_A*BH_A

    r_by     = sqrt(0.5d0*(tmp + sqrt(tmp*tmp + 4.0d0*BH_A*BH_A*z_bh*z_bh)))
    theta_by = acos(z_bh/r_by)
    phi_by   = atan2(y_bh, x_bh)

    return
  end subroutine bhcartesian2boyerlindquist

!===============================================================================
  subroutine calc_direction(r_by, theta_by, phi_by, &
                            theta_screen, phi_screen, &
                            p_r_by, p_theta_by, &
                            E, L, Q)

    use parameters
    use calculation, only : calc_Sigma, calc_Delta, calc_ELQ

    implicit none
    real(8), intent(in)  :: r_by, theta_by, phi_by
    real(8), intent(in)  :: theta_screen, phi_screen
    real(8), intent(out) :: p_r_by, p_theta_by
    real(8), intent(out) :: E, L, Q

    real(8) :: tmp0, tmp1, sigma, delta
    real(8) :: dot_r_by, dot_theta_by, dot_phi_by

    tmp0 = sqrt(r_by*r_by + BH_A*BH_A)
    tmp1 = phi_by - phi_screen

    call calc_Sigma(r_by, theta_by, sigma)
    call calc_Delta(r_by, delta)

    dot_r_by     = -(r_by*tmp0*sin(theta_by)*sin(theta_screen)*cos(tmp1) &
                    + tmp0*tmp0*cos(theta_by)*cos(theta_screen)) / sigma
    dot_theta_by = (r_by*sin(theta_by)*cos(theta_screen) &
                    - tmp0*cos(theta_by)*sin(theta_screen)*cos(tmp1)) / sigma
    dot_phi_by   = sin(theta_screen)*sin(tmp1)/(tmp1*sin(theta_by))

    call calc_ELQ(r_by, theta_by, dot_r_by, dot_theta_by, dot_phi_by, &
                  E, L, Q)

    p_r_by     = dot_r_by*sigma/delta
    p_theta_by = dot_theta_by*sigma

    return
  end subroutine calc_direction

end module initialization
