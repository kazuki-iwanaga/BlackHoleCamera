module initialization

  implicit none
  public :: initialize_4position, initialize_4momentum, initialize_ray
  private

contains

!===============================================================================
  subroutine initialize_4position(r_by, theta_by, phi_by, &
                                  four_position_by)
    implicit none
    double precision, intent(in)  :: r_by, theta_by, phi_by
    double precision, intent(out) :: four_position_by(4)

    four_position_by(1) = 0.0d0
    four_position_by(2) = r_by
    four_position_by(3) = theta_by
    four_position_by(4) = phi_by

    return
  end subroutine initialize_4position

!===============================================================================
  subroutine initialize_4momentum(r_by, theta_by, dot_r_by, dot_theta_by, &
                                  energy, angular_momentum, &
                                  four_momentum_by)

    use parameters
    use calculation, only : calc_sigma, calc_delta

    implicit none
    double precision, intent(in)  :: r_by, theta_by, dot_r_by, dot_theta_by
    double precision, intent(in)  :: energy, angular_momentum
    double precision, intent(out) :: four_momentum_by(4)

    double precision :: sigma, delta

    sigma = calc_sigma(r_by, theta_by)
    delta = calc_delta(r_by)

    four_momentum_by(1) = -energy
    four_momentum_by(2) = dot_r_by*sigma/delta
    four_momentum_by(3) = dot_theta_by*sigma
    four_momentum_by(4) = angular_momentum

    return
  end subroutine initialize_4momentum

!===============================================================================
  subroutine initialize_ray(four_position_by, four_momentum_by, ray)

    implicit none
    double precision, intent(in)  :: four_position_by(4), four_momentum_by(4)
    double precision, intent(out) :: ray(6)

    integer :: i

    do i = 1, 4
      ray(i) = four_position_by(i)
    end do

    ray(5) = four_momentum_by(2)
    ray(6) = four_momentum_by(3)

    return
  end subroutine initialize_ray

end module initialization
