program hello

  use integration
  use calculation
  use parameters

  implicit none
  real(8) :: tmp

  call calc_Sigma(0.0d0, 0.0d0, tmp)
  print *, tmp

  call calc_Delta(0.0d0, tmp)
  print *, tmp

end program hello
