module calculation

  implicit none
  public :: deg2rad, calc_sigma, calc_delta, &
            calc_energy, calc_angular_momentum, calc_carter_constant
  private

contains

!===============================================================================
  function deg2rad(deg)

    use parameters

    implicit none
    double precision :: deg
    double precision :: deg2rad

    deg2rad = deg * PI / 180.0d0

  end function deg2rad

!===============================================================================
  function calc_sigma(r_by, theta_by)

    use parameters

    implicit none
    double precision :: r_by, theta_by
    double precision :: calc_sigma

    calc_sigma = r_by*r_by + BH_A*BH_A*cos(theta_by)*cos(theta_by)

  end function calc_sigma

!===============================================================================
  function calc_delta(r_by)

    use parameters

    implicit none
    double precision :: r_by
    double precision :: calc_delta

    calc_delta = r_by*r_by - 2.0d0*BH_M*r_by + BH_A*BH_A

  end function calc_delta

!===============================================================================
  function calc_energy(r_by, theta_by, &
                       dot_r_by, dot_theta_by, dot_phi_by)

    use parameters

    implicit none
    double precision :: r_by, theta_by
    double precision :: dot_r_by, dot_theta_by, dot_phi_by
    double precision :: calc_energy

    double precision :: sigma, delta

    sigma = calc_sigma(r_by, theta_by)
    delta = calc_delta(r_by)

    calc_energy = sqrt((sigma - 2.0d0*r_by)*(dot_r_by*dot_r_by/delta &
                    + dot_theta_by*dot_theta_by) &
                    + delta*dot_phi_by*dot_phi_by*sin(theta_by)*sin(theta_by))

  end function calc_energy

!===============================================================================
  function calc_angular_momentum(r_by, theta_by, &
                                 dot_r_by, dot_theta_by, dot_phi_by)

    use parameters

    implicit none
    double precision :: r_by, theta_by
    double precision :: dot_r_by, dot_theta_by, dot_phi_by
    double precision :: calc_angular_momentum

    double precision :: sigma, delta

    sigma = calc_sigma(r_by, theta_by)
    delta = calc_delta(r_by)

    calc_angular_momentum = sin(theta_by)*sin(theta_by) &
                              *(sigma*delta*dot_phi_by - 2.0d0*BH_A*r_by &
                              *calc_energy(r_by, theta_by, dot_r_by, &
                                           dot_theta_by, dot_phi_by)) &
                              /(sigma - 2.0d0*r_by)

  end function calc_angular_momentum

!===============================================================================
  function calc_carter_constant(r_by, theta_by, &
                                dot_r_by, dot_theta_by, dot_phi_by)

    use parameters

    implicit none
    double precision :: r_by, theta_by
    double precision :: dot_r_by, dot_theta_by, dot_phi_by
    double precision :: calc_carter_constant

    double precision :: sigma, delta, E, L

    sigma = calc_sigma(r_by, theta_by)
    delta = calc_delta(r_by)

    E =  calc_energy(r_by, theta_by, &
                     dot_r_by, dot_theta_by, dot_phi_by)
    L =  calc_angular_momentum(r_by, theta_by, &
                               dot_r_by, dot_theta_by, dot_phi_by)

    calc_carter_constant = sigma*sigma*dot_theta_by*dot_theta_by &
                            + cos(theta_by)*cos(theta_by) &
                            *(L*L/sin(theta_by)/sin(theta_by) &
                            - BH_A*BH_A*E*E)

  end function calc_carter_constant

end module calculation
