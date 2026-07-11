module grid_module
  ! A simple structured (Cartesian) grid, used by both the Poisson
  ! solver and the transport solver. We keep it dimension-agnostic by
  ! providing both 1D and 2D grid constructors that fill the same
  ! derived type; nz/z are simply unused (size 0) in the 1D case.
  use kinds_module, only: dp, ip
  implicit none
  private
  public :: grid_t, make_grid_1d, make_grid_2d

  type :: grid_t
    integer(ip) :: ndim = 0
    integer(ip) :: nx = 0, ny = 0
    real(dp)    :: xmin = 0.0_dp, xmax = 1.0_dp
    real(dp)    :: ymin = 0.0_dp, ymax = 1.0_dp
    real(dp)    :: dx = 0.0_dp, dy = 0.0_dp
    real(dp), allocatable :: x(:)   ! node coordinates, size nx
    real(dp), allocatable :: y(:)   ! node coordinates, size ny (2D only)
  end type grid_t

contains

  function make_grid_1d(nx, xmin, xmax) result(g)
    integer(ip), intent(in) :: nx
    real(dp), intent(in)    :: xmin, xmax
    type(grid_t) :: g
    integer(ip) :: i
    g%ndim = 1
    g%nx = nx
    g%xmin = xmin; g%xmax = xmax
    g%dx = (xmax - xmin) / real(nx - 1, dp)
    allocate(g%x(nx))
    do i = 1, nx
      g%x(i) = xmin + real(i - 1, dp) * g%dx
    end do
  end function make_grid_1d

  function make_grid_2d(nx, ny, xmin, xmax, ymin, ymax) result(g)
    integer(ip), intent(in) :: nx, ny
    real(dp), intent(in)    :: xmin, xmax, ymin, ymax
    type(grid_t) :: g
    integer(ip) :: i, j
    g%ndim = 2
    g%nx = nx; g%ny = ny
    g%xmin = xmin; g%xmax = xmax
    g%ymin = ymin; g%ymax = ymax
    g%dx = (xmax - xmin) / real(nx - 1, dp)
    g%dy = (ymax - ymin) / real(ny - 1, dp)
    allocate(g%x(nx), g%y(ny))
    do i = 1, nx
      g%x(i) = xmin + real(i - 1, dp) * g%dx
    end do
    do j = 1, ny
      g%y(j) = ymin + real(j - 1, dp) * g%dy
    end do
  end function make_grid_2d

end module grid_module
