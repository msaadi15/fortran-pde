module problem_definitions_module
  ! This is the ONE file you edit to change "the function" being
  ! solved for. Every function here matches one of the abstract
  ! interfaces declared in poisson_module / transport_module, so you
  ! can pass any of them straight into the solvers, or add your own
  ! following the same pattern (just keep the same signature).
  use kinds_module, only: dp, ip
  implicit none
  private
  public :: pi
  public :: poisson1d_rhs_manufactured, poisson1d_bc_manufactured
  public :: poisson2d_rhs_manufactured, poisson2d_bc_manufactured
  public :: poisson2d_rhs_custom,       poisson2d_bc_zero
  public :: transport1d_init_gaussian,  transport1d_bc_zero
  public :: transport2d_init_gaussian,  transport2d_bc_zero

  real(dp), parameter :: pi = 3.14159265358979323846_dp

contains

  ! -----------------------------------------------------------------
  ! 1D Poisson, manufactured solution:  u(x) = sin(pi x)
  ! => -u'' = pi^2 sin(pi x) = f(x)
  ! Used in test/check_poisson.f90 to verify convergence order.
  ! -----------------------------------------------------------------
  pure function poisson1d_rhs_manufactured(x) result(f)
    real(dp), intent(in) :: x
    real(dp) :: f
    f = pi * pi * sin(pi * x)
  end function poisson1d_rhs_manufactured

  pure function poisson1d_bc_manufactured(x) result(g)
    real(dp), intent(in) :: x
    real(dp) :: g
    g = sin(pi * x)   ! = 0 at x=0 and x=1, but written generally
  end function poisson1d_bc_manufactured

  ! -----------------------------------------------------------------
  ! 2D Poisson, manufactured solution: u(x,y) = sin(pi x) sin(pi y)
  ! => -Laplacian(u) = 2 pi^2 sin(pi x) sin(pi y) = f(x,y)
  ! -----------------------------------------------------------------
  pure function poisson2d_rhs_manufactured(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f
    f = 2.0_dp * pi * pi * sin(pi * x) * sin(pi * y)
  end function poisson2d_rhs_manufactured

  pure function poisson2d_bc_manufactured(x, y) result(g)
    real(dp), intent(in) :: x, y
    real(dp) :: g
    g = sin(pi * x) * sin(pi * y)   ! = 0 on the unit-square boundary
  end function poisson2d_bc_manufactured

  ! -----------------------------------------------------------------
  ! 2D Poisson, "ANY function" example: two Gaussian point-sources.
  ! This is the one to copy/edit if you want to try your own f(x,y) --
  ! there is no restriction on what f can be (smooth, localized,
  ! oscillatory, etc.), the solver does not know or care.
  ! -----------------------------------------------------------------
  pure function poisson2d_rhs_custom(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f
    real(dp), parameter :: x1 = 0.3_dp, y1 = 0.3_dp
    real(dp), parameter :: x2 = 0.7_dp, y2 = 0.6_dp
    real(dp), parameter :: sigma = 0.06_dp
    f =  1.0_dp * exp(-((x-x1)**2 + (y-y1)**2) / (2.0_dp*sigma**2)) &
       - 1.0_dp * exp(-((x-x2)**2 + (y-y2)**2) / (2.0_dp*sigma**2))
  end function poisson2d_rhs_custom

  pure function poisson2d_bc_zero(x, y) result(g)
    real(dp), intent(in) :: x, y
    real(dp) :: g
    g = 0.0_dp
  end function poisson2d_bc_zero

  ! -----------------------------------------------------------------
  ! Transport equation initial conditions and BCs: a Gaussian "blob"
  ! that gets advected and diffused, with u=0 held on the boundary.
  ! -----------------------------------------------------------------
  pure function transport1d_init_gaussian(x) result(u0)
    real(dp), intent(in) :: x
    real(dp) :: u0
    real(dp), parameter :: x0 = 0.2_dp, sigma = 0.05_dp
    u0 = exp(-(x - x0)**2 / (2.0_dp * sigma**2))
  end function transport1d_init_gaussian

  pure function transport1d_bc_zero(x, t) result(g)
    real(dp), intent(in) :: x, t
    real(dp) :: g
    g = 0.0_dp
  end function transport1d_bc_zero

  pure function transport2d_init_gaussian(x, y) result(u0)
    real(dp), intent(in) :: x, y
    real(dp) :: u0
    real(dp), parameter :: x0 = 0.2_dp, y0 = 0.5_dp, sigma = 0.06_dp
    u0 = exp(-((x - x0)**2 + (y - y0)**2) / (2.0_dp * sigma**2))
  end function transport2d_init_gaussian

  pure function transport2d_bc_zero(x, y, t) result(g)
    real(dp), intent(in) :: x, y, t
    real(dp) :: g
    g = 0.0_dp
  end function transport2d_bc_zero

end module problem_definitions_module
