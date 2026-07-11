program check_poisson
  ! Validates the Poisson solver against a manufactured solution:
  !   u(x,y) = sin(pi x) sin(pi y)
  ! This is THE standard way to check a FEM/FVM/FDM PDE code: pick a
  ! known analytic solution, derive the source term f that produces
  ! it, solve numerically, and confirm (a) the error is small and
  ! (b) it shrinks at the expected rate as the grid is refined
  ! (2nd-order central differences => error should behave like h^2,
  ! i.e. halving h should quarter the error).
  use kinds_module, only: dp, ip
  use grid_module, only: grid_t, make_grid_2d
  use poisson_module, only: solve_poisson_2d
  use problem_definitions_module, only: poisson2d_rhs_manufactured, poisson2d_bc_manufactured, pi
  implicit none

  integer(ip), parameter :: sizes(3) = [21_ip, 41_ip, 81_ip]
  real(dp) :: errors(3)
  integer(ip) :: k, i, j, iters
  type(grid_t) :: g
  real(dp), allocatable :: u(:,:), uexact(:,:)
  logical :: passed

  do k = 1, 3
    g = make_grid_2d(sizes(k), sizes(k), 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call solve_poisson_2d(g, poisson2d_rhs_manufactured, poisson2d_bc_manufactured, u, &
                           tol=1.0e-12_dp, maxiter=20000_ip, iters_used=iters)
    allocate(uexact(g%nx, g%ny))
    do j = 1, g%ny
      do i = 1, g%nx
        uexact(i,j) = sin(pi * g%x(i)) * sin(pi * g%y(j))
      end do
    end do
    errors(k) = maxval(abs(u - uexact))
    print '(A,I4,A,ES12.4,A,I6)', 'n=', sizes(k), '  max error=', errors(k), '  CG iters=', iters
    deallocate(u, uexact)
  end do

  print *
  print '(A,F6.3)', 'observed order (n=21->41): ', log(errors(1)/errors(2)) / log(2.0_dp)
  print '(A,F6.3)', 'observed order (n=41->81): ', log(errors(2)/errors(3)) / log(2.0_dp)

  ! Central differences are 2nd order: expect the ratio to be close
  ! to 2.0. Allow a loose tolerance since this is a coarse grid test.
  passed = (log(errors(1)/errors(2)) / log(2.0_dp) > 1.7_dp) .and. &
           (log(errors(2)/errors(3)) / log(2.0_dp) > 1.7_dp)

  if (passed) then
    print *, 'PASS: convergence order is consistent with 2nd-order finite differences'
  else
    print *, 'FAIL: convergence order lower than expected'
    stop 1
  end if

end program check_poisson
