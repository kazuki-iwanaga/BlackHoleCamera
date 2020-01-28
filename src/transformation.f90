module transformation

  implicit none
  public :: screen2bhcartesian, bhcartesian2boyerlindquist, calc_direction
  private

contains

!===============================================================================
  subroutine screen2bhcartesian(i_pixel, j_pixel, &
                                r_screen, theta_screen, phi_screen, &
                                x_bh, y_bh, z_bh)

    use parameters

    implicit none
    double precision, intent(in)  :: i_pixel, j_pixel
    double precision, intent(in)  :: r_screen, theta_screen, phi_screen
    double precision, intent(out) :: x_bh, y_bh, z_bh

    double precision :: tmp

    tmp = sqrt(r_screen*r_screen + BH_A*BH_A)*sin(theta_screen) &
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
    double precision, intent(in)  :: x_bh, y_bh, z_bh
    double precision, intent(out) :: r_by, theta_by, phi_by

    double precision :: tmp

    tmp = x_bh*x_bh + y_bh*y_bh + z_bh*z_bh - BH_A*BH_A

    r_by     = sqrt(0.5d0*(tmp + sqrt(tmp*tmp + 4.0d0*BH_A*BH_A*z_bh*z_bh)))
    theta_by = acos(z_bh/r_by)
    phi_by   = atan2(y_bh, x_bh)

    return
  end subroutine bhcartesian2boyerlindquist

!===============================================================================
  subroutine calc_direction(r_by, theta_by, phi_by, &
                            theta_screen, phi_screen, &
                            dot_r_by, dot_theta_by, dot_phi_by)

    use parameters
    use calculation, only : calc_sigma

    implicit none
    double precision, intent(in)  :: r_by, theta_by, phi_by
    double precision, intent(in)  :: theta_screen, phi_screen
    double precision, intent(out) :: dot_r_by, dot_theta_by, dot_phi_by

    double precision :: tmp_R, tmp_PHI

    tmp_R = sqrt(r_by*r_by + BH_A*BH_A)
    tmp_PHI = phi_by - phi_screen

    dot_r_by     = -(r_by*tmp_R*sin(theta_by)*sin(theta_screen)*cos(tmp_PHI) &
                    + tmp_R*tmp_R*cos(theta_by)*cos(theta_screen)) &
                    /calc_sigma(r_by, theta_by)
    dot_theta_by = (r_by*sin(theta_by)*cos(theta_screen) &
                    - tmp_R*cos(theta_by)*sin(theta_screen)*cos(tmp_PHI) &
                    /calc_sigma(r_by, theta_by))
    dot_phi_by   = sin(theta_screen)*sin(tmp_PHI)/(tmp_R*sin(theta_by))

    return
  end subroutine calc_direction

end module transformation
