module calculation

  implicit none
  public :: calc_ELQ, calc_Energy, calc_Angular_Momentum, &
            calc_Carter_Constant, calc_Sigma, calc_Delta
  private

contains

!===============================================================================
  subroutine calc_ELQ(r, theta, dot_r, dot_theta, dot_phi, E, L, Q)

    implicit none
    real(8), intent(in)  :: r, theta, dot_r, dot_theta, dot_phi
    real(8), intent(out) :: E, L, Q
    real(8) :: Sigma, Delta

    call calc_Sigma(r, theta, Sigma)
    call calc_Delta(r, Delta)

    call calc_Energy(r, theta, dot_r, dot_theta, dot_phi, Sigma, Delta, E)
    call calc_Angular_Momentum(r, theta, dot_phi, Sigma, Delta, E, L)
    call calc_Carter_Constant(theta, dot_theta, Sigma, E, L, Q)

    return
  end subroutine calc_ELQ

!===============================================================================
  subroutine calc_Energy(r, theta, dot_r, dot_theta, dot_phi, Sigma, Delta, E)

    implicit none
    real(8), intent(in)  :: r, theta, dot_r, dot_theta, dot_phi
    real(8), intent(in)  :: Sigma, Delta
    real(8), intent(out) :: E
    real(8) :: E_sq

    E_sq = (Sigma - 2.0d0 * r) * (Sigma * dot_r * dot_r &
            + Sigma * Delta * dot_theta * dot_theta) / (Sigma * Delta) &
            + Delta * dot_phi * dot_phi * sin(theta) * sin(theta)
    E = sqrt(E_sq)

    return
  end subroutine calc_Energy

!===============================================================================
  subroutine calc_Angular_Momentum(r, theta, dot_phi, Sigma, Delta, E, L)

    use parameters

    implicit none
    real(8), intent(in)  :: r, theta, dot_phi
    real(8), intent(in)  :: Sigma, Delta, E
    real(8), intent(out) :: L

    L = (Sigma * Delta * dot_phi - 2.0d0 * BH_A * r * E) &
          * sin(theta)*sin(theta) / (Sigma - 2.0d0 * r)

    return
  end subroutine calc_Angular_Momentum

!===============================================================================
  subroutine calc_Carter_Constant(theta, dot_theta, Sigma, E, L, Q)

    use parameters

    implicit none
    real(8), intent(in)  :: theta, dot_theta
    real(8), intent(in)  :: Sigma, E, L
    real(8), intent(out) :: Q

    Q = Sigma * Sigma * dot_theta * dot_theta + cos(theta) * cos(theta) &
        * (L * L / (sin(theta) * sin(theta)) - BH_A * BH_A * E * E)

    return
  end subroutine calc_Carter_Constant

!===============================================================================
  subroutine calc_Sigma(r, theta, Sigma)

    use parameters

    implicit none
    real(8), intent(in)  :: r, theta
    real(8), intent(out) :: Sigma

    Sigma = r * r + BH_A * BH_A * cos(theta) * cos(theta)

    return
  end subroutine calc_Sigma

!===============================================================================
  subroutine calc_Delta(r, Delta)

    use parameters

    implicit none
    real(8), intent(in)  :: r
    real(8), intent(out) :: Delta

    Delta = r * r - 2.0d0 * r * BH_M + BH_A * BH_A

    return
  end subroutine calc_Delta

end module calculation
