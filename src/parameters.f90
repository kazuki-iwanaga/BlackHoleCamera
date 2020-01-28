module parameters

  implicit none

  ! ...
  ! Constants
  ! <-
  double precision, parameter :: PI = 4.0d0*atan(1.0d0)
  ! ->

  ! ...
  ! Simulation Parameters
  ! <-
  double precision, parameter :: DISTANCE  = 100.0d0
  double precision, parameter :: ELEVATION = 0.0d0
  double precision, parameter :: AZIMUTH   = 0.0d0

  double precision, parameter :: X_PIXEL = 3.6d0
  double precision, parameter :: Y_PIXEL = 0.0d0

  double precision, parameter :: BH_M = 1.0d0
  double precision, parameter :: BH_R = 1.0d0
  ! double precision, parameter :: BH_A = 0.998d0
  ! double precision, parameter :: BH_A = 0.5d0
  double precision, parameter :: BH_A = 0.0d0

  double precision, parameter :: RK_STEP = 1.0d-3
  double precision, parameter :: MAX_R   = 150.0d0
  double precision, parameter :: MIN_R   = 1.0d0
  ! ->

end module parameters
