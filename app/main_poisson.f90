program main_poisson
  use kinds_module, only: dp, ip
  use grid_module,  only: grid_t, make_grid_1d, make_grid_2d
  use poisson_module, only: solve_poisson_1d, solve_poisson_2d
  use problem_definitions_module, only: poisson1d_rhs_manufactured, poisson1d_bc_manufactured, &
                                         poisson2d_rhs_manufactured, poisson2d_bc_manufactured, &
                                         poisson2d_rhs_custom, poisson2d_bc_zero
  use io_module, only: write_field_1d, write_field_2d
  implicit none

  type(grid_t) :: g1, g2
  real(dp), allocatable :: u1(:), u2(:,:)
  integer(ip) :: iters

  print *, '=== 1D Poisson: -u'''' = pi^2 sin(pi x), manufactured solution ==='
  g1 = make_grid_1d(nx=101_ip, xmin=0.0_dp, xmax=1.0_dp)
  call solve_poisson_1d(g1, poisson1d_rhs_manufactured, poisson1d_bc_manufactured, u1)
  call write_field_1d('results/poisson1d_manufactured.dat', g1, u1)
  print *, 'max error vs sin(pi x):', maxval(abs(u1 - sin_pi(g1%x)))
  print *, 'written to results/poisson1d_manufactured.dat'
  print *

  print *, '=== 2D Poisson: manufactured solution sin(pi x) sin(pi y) ==='
  g2 = make_grid_2d(nx=81_ip, ny=81_ip, xmin=0.0_dp, xmax=1.0_dp, ymin=0.0_dp, ymax=1.0_dp)
  call solve_poisson_2d(g2, poisson2d_rhs_manufactured, poisson2d_bc_manufactured, u2, &
                         tol=1.0e-10_dp, maxiter=5000_ip, iters_used=iters)
  call write_field_2d('results/poisson2d_manufactured.dat', g2, u2)
  print *, 'CG iterations:', iters
  print *, 'written to results/poisson2d_manufactured.dat'
  print *

  print *, '=== 2D Poisson: custom source (two Gaussian point charges) ==='
  call solve_poisson_2d(g2, poisson2d_rhs_custom, poisson2d_bc_zero, u2, &
                         tol=1.0e-10_dp, maxiter=5000_ip, iters_used=iters)
  call write_field_2d('results/poisson2d_custom_source.dat', g2, u2)
  print *, 'CG iterations:', iters
  print *, 'written to results/poisson2d_custom_source.dat'

contains

  ! Small local helper, only used for the printed error diagnostic
  ! above (not part of the reusable library).
  pure function sin_pi(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: s(size(x))
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    s = sin(pi * x)
  end function sin_pi

end program main_poisson
