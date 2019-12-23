module integration

  implicit none
  public  :: null_geodesic_RK4
  private :: null_geodesic_equation

contains

!===============================================================================
  subroutine null_geodesic_RK4(E, L, Q, &
                               t_0, r_0, theta_0, phi_0, p_r_0, p_theta_0, &
                               t_1, r_1, theta_1, phi_1, p_r_1, p_theta_1)

    use parameters
    use calculation, only : calc_Sigma, calc_Delta

    implicit none
    real(8), intent(in) :: E, L, Q
    real(8), intent(in) :: t_0, r_0, theta_0, phi_0, p_r_0, p_theta_0
    real(8), intent(out) :: t_1, r_1, theta_1, phi_1, p_r_1, p_theta_1
    real(8) :: Sigma, Delta, kappa
    real(8), dimension(6) :: k0, k1, k2, k3, RK_inc
    integer :: i

    ! ...
    ! for convenience
    ! <-
    call calc_Sigma(r_0, theta_0, Sigma)
    call calc_Delta(r_0, Delta)
    kappa = Q + L * L + BH_A * BH_A * E * E
    ! ->

    ! ...
    ! 4th order Runge-Kutta method
    ! <-
    call null_geodesic_equation(E, L, Sigma, Delta, kappa, &
                                r_0, &
                                theta_0, &
                                p_r_0, &
                                k0)

    call null_geodesic_equation(E, L, Sigma, Delta, kappa, &
                                r_0     + 0.5d0 * D_LAMBDA * k0(2), &
                                theta_0 + 0.5d0 * D_LAMBDA * k0(3), &
                                p_r_0   + 0.5d0 * D_LAMBDA * k0(5), &
                                k1)

    call null_geodesic_equation(E, L, Sigma, Delta, kappa, &
                                r_0     + 0.5d0 * D_LAMBDA * k1(2), &
                                theta_0 + 0.5d0 * D_LAMBDA * k1(3), &
                                p_r_0   + 0.5d0 * D_LAMBDA * k1(5), &
                                k2)

    call null_geodesic_equation(E, L, Sigma, Delta, kappa, &
                                r_0     + D_LAMBDA * k2(2), &
                                theta_0 + D_LAMBDA * k2(3), &
                                p_r_0   + D_LAMBDA * k2(5), &
                                k3)

    ! initialization
    do i = 1, 6
      RK_inc(i) = 0.0d0
    end do

    ! calculate RK4 gradients
    do i = 1, 6
      RK_inc(i) = RK_inc(i) &
                  + (k0(i) + 2.0d0 * k1(i) + 2.0d0 * k2(i) + k3(i)) / 6.0d0
    end do

    ! update variables
    t_1       = t_0       + D_LAMBDA * RK_inc(1)
    r_1       = r_0       + D_LAMBDA * RK_inc(2)
    theta_1   = theta_0   + D_LAMBDA * RK_inc(3)
    phi_1     = phi_0     + D_LAMBDA * RK_inc(4)
    p_r_1     = p_r_0     + D_LAMBDA * RK_inc(5)
    p_theta_1 = p_theta_0 + D_LAMBDA * RK_inc(6)
    ! ->

    return
  end subroutine null_geodesic_RK4

!===============================================================================
  subroutine null_geodesic_equation(E, L, Sigma, Delta, kappa, &
                                    r, theta, p_r, k)

    use parameters

    implicit none
    real(8), intent(in) :: E, L, Sigma, Delta, kappa
    real(8), intent(in) :: r, theta, p_r
    real(8), dimension(6), intent(out) :: k
    real(8) :: sintheta, Sigma_inverse, Sigma_Delta_inverse

    ! ...
    ! for convenience
    ! <-
    sintheta = sin(theta)
    Sigma_inverse = 1.0d0 / Sigma
    Sigma_Delta_inverse = 1.0d0 / (Sigma * Delta)
    ! ->

    ! ...
    ! Null Geodesic Equations
    ! <-
    ! dot_p_r
    k(5) = (2.0d0 * r * (r * r + BH_A * BH_A) * E * E &
            - 2.0d0 * BH_A * E * L &
            - kappa * (r - 1.0d0)) * Sigma_Delta_inverse &
            - 2.0d0 * p_r * p_r * (r - 1.0d0) * Sigma_inverse
    ! dot_p_theta
    k(6) = (L * L / (sintheta * sintheta * sintheta * sintheta) &
            - BH_A * BH_A * E * E) * sintheta * cos(theta) * Sigma_inverse
    ! dot_t
    k(1) = E + 2.0d0 * (r * (r * r + BH_A * BH_A) * E - BH_A * r * L) &
            * Sigma_Delta_inverse
    ! dot_r
    k(2) = k(5) * Delta * Sigma_inverse
    ! dot_theta
    k(3) = k(6) * Sigma_inverse
    ! dot_phi
    k(4) = (2.0d0 * BH_A * r * E + (Sigma - 2.0d0 * r) * L &
            / (sintheta * sintheta)) * Sigma_Delta_inverse
    ! ->

    return
  end subroutine null_geodesic_equation

end module integration
