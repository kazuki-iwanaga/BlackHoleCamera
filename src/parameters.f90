module parameters

  implicit none

  ! ...
  ! Constants
  ! <-
  real(8), parameter :: PI = 4.0d0 * atan(1.0d0)
  ! ->

  ! ...
  ! Simulation Parameters
  ! <-
  real(8), parameter :: EPS = 1.0d-10

  real(8), parameter :: BH_M = 1.0d0
  real(8), parameter :: BH_R = 1.0d0
  real(8), parameter :: BH_A = 0.998d0

  real(8), parameter :: DISTANCE  = 100.0d0
  real(8), parameter :: ELEVATION = 0.0d0
  real(8), parameter :: AZIMUTH   = 0.0d0

  real(8), parameter :: SCREEN_X = 5.0d0
  real(8), parameter :: SCREEN_Y = 0.0d0

  real(8), parameter :: D_LAMBDA = 1.0d-2
  real(8), parameter :: MAX_R    = 110.0d0
  real(8), parameter :: MIN_R    = 1.0d0
  ! ->

end module parameters
