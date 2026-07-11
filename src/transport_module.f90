module transport_module
  ! Solves the time-dependent advection-diffusion (transport) equation
  !
  !     du/dt + v . grad(u) = D * Laplacian(u)
  !
  ! in 1D:  u_t + vx*u_x           = D*u_xx
  ! in 2D:  u_t + vx*u_x + vy*u_y  = D*(u_xx + u_yy)
  !
  ! Discretization:
  !   time      -> explicit forward Euler
  !   advection -> first-order upwind (stable, diffusive but robust --
  !                the standard first solver to reach for)
  !   diffusion -> second-order central differences
  !
  ! Stability (explicit scheme, must be respected by the caller /
  ! chosen dt):
  !   CFL (advection):     |v| * dt / dx        <= 1
  !   Diffusion number:    2 * D * dt / dx^2    <= 1   (1D)
  !                        D*dt*(1/dx^2+1/dy^2) <= 0.5 (2D)
  !
  ! Boundary conditions: Dirichlet, supplied via the same bc_func
  ! interface as the Poisson module, evaluated at the current time.
  use kinds_module, only: dp, ip
  use grid_module,  only: grid_t
  implicit none
  private
  public :: init_func_1d, bc_func_1d_t, init_func_2d, bc_func_2d_t
  public :: step_transport_1d, step_transport_2d, max_stable_dt_1d, max_stable_dt_2d

  abstract interface
    pure function init_func_1d(x) result(u0)
      import :: dp
      real(dp), intent(in) :: x
      real(dp) :: u0
    end function init_func_1d

    pure function bc_func_1d_t(x, t) result(g)
      import :: dp
      real(dp), intent(in) :: x, t
      real(dp) :: g
    end function bc_func_1d_t

    pure function init_func_2d(x, y) result(u0)
      import :: dp
      real(dp), intent(in) :: x, y
      real(dp) :: u0
    end function init_func_2d

    pure function bc_func_2d_t(x, y, t) result(g)
      import :: dp
      real(dp), intent(in) :: x, y, t
      real(dp) :: g
    end function bc_func_2d_t
  end interface

contains

  ! Largest forward-Euler timestep satisfying both the CFL and
  ! diffusion-number stability limits, with a safety factor.
  pure function max_stable_dt_1d(g, vx, D, safety) result(dt)
    type(grid_t), intent(in) :: g
    real(dp), intent(in) :: vx, D
    real(dp), intent(in), optional :: safety
    real(dp) :: dt, dt_adv, dt_diff, s
    s = 0.9_dp
    if (present(safety)) s = safety
    dt_adv  = huge(1.0_dp)
    dt_diff = huge(1.0_dp)
    if (abs(vx) > 1.0e-14_dp) dt_adv = g%dx / abs(vx)
    if (D > 1.0e-14_dp)       dt_diff = 0.5_dp * g%dx**2 / D
    dt = s * min(dt_adv, dt_diff)
  end function max_stable_dt_1d

  pure function max_stable_dt_2d(g, vx, vy, D, safety) result(dt)
    type(grid_t), intent(in) :: g
    real(dp), intent(in) :: vx, vy, D
    real(dp), intent(in), optional :: safety
    real(dp) :: dt, dt_adv, dt_diff, s
    s = 0.9_dp
    if (present(safety)) s = safety
    dt_adv  = huge(1.0_dp)
    dt_diff = huge(1.0_dp)
    if (abs(vx)/g%dx + abs(vy)/g%dy > 1.0e-14_dp) &
      dt_adv = 1.0_dp / (abs(vx)/g%dx + abs(vy)/g%dy)
    if (D > 1.0e-14_dp) &
      dt_diff = 0.5_dp / (D * (1.0_dp/g%dx**2 + 1.0_dp/g%dy**2))
    dt = s * min(dt_adv, dt_diff)
  end function max_stable_dt_2d

  ! Advances u by one timestep dt, in place, for the 1D equation.
  subroutine step_transport_1d(g, u, vx, D, dt, t_new, bc)
    type(grid_t), intent(in)    :: g
    real(dp), intent(inout)     :: u(:)
    real(dp), intent(in)        :: vx, D, dt, t_new
    procedure(bc_func_1d_t)     :: bc

    integer(ip) :: i, n
    real(dp) :: idx, idx2, adv, diff
    real(dp), allocatable :: unew(:)

    n = size(u)
    idx  = 1.0_dp / g%dx
    idx2 = 1.0_dp / (g%dx * g%dx)
    allocate(unew(n))

    do i = 2, n - 1
      ! First-order upwind advection: pick the difference on the side
      ! the flow is coming FROM, which is what keeps this scheme
      ! stable for pure advection.
      if (vx >= 0.0_dp) then
        adv = vx * (u(i) - u(i-1)) * idx
      else
        adv = vx * (u(i+1) - u(i)) * idx
      end if
      diff = D * (u(i+1) - 2.0_dp*u(i) + u(i-1)) * idx2
      unew(i) = u(i) + dt * (diff - adv)
    end do

    unew(1) = bc(g%x(1), t_new)
    unew(n) = bc(g%x(n), t_new)
    u = unew
  end subroutine step_transport_1d

  ! Advances u by one timestep dt, in place, for the 2D equation.
  subroutine step_transport_2d(g, u, vx, vy, D, dt, t_new, bc)
    type(grid_t), intent(in)    :: g
    real(dp), intent(inout)     :: u(:,:)
    real(dp), intent(in)        :: vx, vy, D, dt, t_new
    procedure(bc_func_2d_t)     :: bc

    integer(ip) :: i, j, nx, ny
    real(dp) :: idx, idy, idx2, idy2, advx, advy, diff
    real(dp), allocatable :: unew(:,:)

    nx = size(u,1); ny = size(u,2)
    idx  = 1.0_dp / g%dx;         idy  = 1.0_dp / g%dy
    idx2 = 1.0_dp / (g%dx*g%dx);  idy2 = 1.0_dp / (g%dy*g%dy)
    allocate(unew(nx,ny))

    do j = 2, ny - 1
      do i = 2, nx - 1
        if (vx >= 0.0_dp) then
          advx = vx * (u(i,j) - u(i-1,j)) * idx
        else
          advx = vx * (u(i+1,j) - u(i,j)) * idx
        end if
        if (vy >= 0.0_dp) then
          advy = vy * (u(i,j) - u(i,j-1)) * idy
        else
          advy = vy * (u(i,j+1) - u(i,j)) * idy
        end if
        diff = D * ( (u(i+1,j) - 2.0_dp*u(i,j) + u(i-1,j)) * idx2 &
                   + (u(i,j+1) - 2.0_dp*u(i,j) + u(i,j-1)) * idy2 )
        unew(i,j) = u(i,j) + dt * (diff - advx - advy)
      end do
    end do

    do i = 1, nx
      unew(i,1)  = bc(g%x(i), g%y(1),  t_new)
      unew(i,ny) = bc(g%x(i), g%y(ny), t_new)
    end do
    do j = 1, ny
      unew(1,j)  = bc(g%x(1),  g%y(j), t_new)
      unew(nx,j) = bc(g%x(nx), g%y(j), t_new)
    end do
    u = unew
  end subroutine step_transport_2d

end module transport_module
