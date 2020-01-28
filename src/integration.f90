module integration

  implicit none
  public  :: null_geodesic_RK4, update_ray
  private :: null_geodesic_equations, add_ray

contains

!===============================================================================
  subroutine null_geodesic_RK4(E, L, Q, ray, ray_update)

    use parameters

    implicit none
    double precision, intent(in)  :: E, L, Q
    double precision, intent(in)  :: ray(6)
    double precision, intent(out) :: ray_update(6)

    double precision :: tmp1(6), tmp2(6), tmp3(6), tmp4(6)
    double precision :: k1(6), k2(6), k3(6), k4(6), RK_inc(6)
    integer :: i

    ! Initialization
    do i = 1, 6
      RK_inc(i) = 0.0d0
    end do

    ! ...
    ! 4th Order Runge-Kutta Method
    ! <-

    ! Calculate RK4 gradients
    call update_ray(ray, tmp1)
    call null_geodesic_equations(E, L, Q, tmp1, k1)

    call add_ray(tmp1, k1, 0.5d0, tmp2)
    call null_geodesic_equations(E, L, Q, tmp2, k2)

    call add_ray(tmp1, k2, 0.5d0, tmp3)
    call null_geodesic_equations(E, L, Q, tmp3, k3)

    call add_ray(tmp1, k3, 1.0d0, tmp4)
    call null_geodesic_equations(E, L, Q, tmp4, k4)

    ! Calculate RK4 increments
    do i = 1, 6
      RK_inc(i) = RK_inc(i) &
                  + (k1(i) + 2.0d0*k2(i) + 2.0d0*k3(i) + k4(i))/6.0d0
    end do
    ! print *, RK_inc

    ! Update variables
    do i = 1, 6
      ray_update(i) = ray(i) + RK_STEP*RK_inc(i)
    end do

    ! ->

    return
  end subroutine null_geodesic_RK4

!===============================================================================
  subroutine null_geodesic_equations(E, L, Q, ray, k)

    use parameters
    use calculation, only : calc_sigma, calc_delta, calc_energy, &
                            calc_angular_momentum, calc_carter_constant

    implicit none
    double precision, intent(in)  :: E, L, Q
    double precision, intent(in)  :: ray(6)
    double precision, intent(out) :: k(6)

    double precision :: t, r, theta, phi, p_r, p_theta
    double precision :: sigma, delta, sigma_inverse, delta_inverse
    double precision :: kappa, sintheta

    ! tmp
    t       = ray(1)
    r       = ray(2)
    theta   = ray(3)
    phi     = ray(4)
    p_r     = ray(5)
    p_theta = ray(6)

    sigma         = calc_sigma(r, theta)
    sigma_inverse = 1.0d0/sigma
    delta         = calc_delta(r)
    delta_inverse = 1.0d0/delta

    kappa    = Q + L*L + BH_A*BH_A*E*E
    sintheta = sin(theta)

    ! ...
    ! Null Geodesic Equations
    ! <-
    k(1) = E + sigma_inverse*delta_inverse &
            *(2.0d0*r*(r*r + BH_A*BH_A)*E*E - 2.0d0*BH_A*r*L)
    k(2) = p_r*delta*sigma_inverse
    k(3) = p_theta*sigma_inverse
    k(4) = sigma_inverse*delta_inverse*(2.0d0*BH_A*r*E &
            + (sigma - 2.0d0*r)*L/(sintheta*sintheta))
    k(5) = sigma_inverse*delta_inverse*(-kappa*(r - 1.0d0) &
            + 2.0d0*r*(r*r + BH_A*BH_A)*E*E - 2.0d0*BH_A*E*L) &
            - 2.0d0*p_r*p_r*(r - 1.0d0)*sigma_inverse
    k(6) = sintheta*cos(theta)*sigma_inverse &
            *(L*L/(sintheta*sintheta*sintheta*sintheta) &
            - BH_A*BH_A*E*E)
    ! ->

    ! print *, ray
    ! print *, k

    return
  end subroutine null_geodesic_equations

!===============================================================================
  subroutine update_ray(ray_update, ray)

    implicit none
    double precision, intent(in)  :: ray_update(6)
    double precision, intent(out) :: ray(6)

    integer :: i

    do i = 1, 6
      ray(i) = ray_update(i)
    end do

    return
  end subroutine update_ray

!===============================================================================
  subroutine add_ray(ray, ray_add, weight, ray_added)

    use parameters

    implicit none
    double precision, intent(in)  :: ray(6), ray_add(6)
    double precision, intent(in)  :: weight
    double precision, intent(out) :: ray_added(6)

    integer :: i

    do i = 1, 6
      ray_added(i) = ray(i) + weight*RK_STEP*ray_add(i)
    end do

    return
  end subroutine add_ray

end module integration
