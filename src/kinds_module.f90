module kinds_module
  ! Central place for numeric precision. Using iso_fortran_env is the
  ! modern, portable way to request a specific real kind instead of
  ! relying on compiler-dependent "real*8" or "double precision".
  use iso_fortran_env, only: dp => real64, ip => int32
  implicit none
  public :: dp, ip
end module kinds_module
