module poisson_module
  ! Solves the Poisson / Laplace equation
  !
  !     - Laplacian(u) = f(x)      in the domain
  !                 u  = g(x)      on the boundary  (Dirichlet BC)
  !
  ! in 1D:  -u''(x) = f(x)
  ! in 2D:  -(u_xx + u_yy) = f(x,y)
  !
  ! f is passed in as a procedure argument, so this module works for
  ! ANY right-hand-side function you supply (that is what "for any
  ! function" means here) -- see app/main_poisson.f90 for examples.
  !
  ! Discretization: standard second-order central finite differences
  ! (3-point stencil in 1D, 5-point stencil in 2D).
  !
  ! Solvers:
  !   1D -> direct tridiagonal solve (LAPACK dgtsv)
  !   2D -> matrix-free Conjugate Gradient (the matrix is symmetric
  !         positive definite, so CG is the natural choice, and
  !         "matrix-free" means we never build/store the huge sparse
  !         matrix -- we only ever apply the 5-point stencil).
  use kinds_module, only: dp, ip
  use grid_module,  only: grid_t
  implicit none
  private
  public :: rhs_func_1d, bc_func_1d, rhs_func_2d, bc_func_2d
  public :: solve_poisson_1d, solve_poisson_2d

  abstract interface
    pure function rhs_func_1d(x) result(f)
      import :: dp
      real(dp), intent(in) :: x
      real(dp) :: f
    end function rhs_func_1d

    pure function bc_func_1d(x) result(g)
      import :: dp
      real(dp), intent(in) :: x
      real(dp) :: g
    end function bc_func_1d

    pure function rhs_func_2d(x, y) result(f)
      import :: dp
      real(dp), intent(in) :: x, y
      real(dp) :: f
    end function rhs_func_2d

    pure function bc_func_2d(x, y) result(g)
      import :: dp
      real(dp), intent(in) :: x, y
      real(dp) :: g
    end function bc_func_2d
  end interface

contains

  ! ---------------------------------------------------------------
  ! 1D solve:  -u'' = f(x),  u(xmin) = g(xmin),  u(xmax) = g(xmax)
  ! ---------------------------------------------------------------
  subroutine solve_poisson_1d(g, f, bc, u)
    type(grid_t), intent(in)          :: g
    procedure(rhs_func_1d)            :: f
    procedure(bc_func_1d)             :: bc
    real(dp), allocatable, intent(out) :: u(:)

    integer(ip) :: n, i, ninterior, info
    real(dp) :: h2
    real(dp), allocatable :: dl(:), d(:), du(:), rhs(:)

    n = g%nx
    h2 = g%dx * g%dx
    allocate(u(n))

    ! Boundary values are known exactly.
    u(1) = bc(g%x(1))
    u(n) = bc(g%x(n))

    ninterior = n - 2
    if (ninterior <= 0) return

    ! Tridiagonal system for interior points:
    !   (-u(i-1) + 2 u(i) - u(i+1)) / h^2 = f(x_i)
    allocate(dl(ninterior), d(ninterior), du(ninterior), rhs(ninterior))
    d  = 2.0_dp / h2
    dl = -1.0_dp / h2
    du = -1.0_dp / h2
    do i = 1, ninterior
      rhs(i) = f(g%x(i + 1))
    end do
    ! Fold in known boundary contributions.
    rhs(1)         = rhs(1)         + u(1) / h2
    rhs(ninterior) = rhs(ninterior) + u(n) / h2

    ! LAPACK dgtsv: solves a general tridiagonal system A x = rhs in
    ! place (rhs is overwritten with the solution x).
    call dgtsv(ninterior, 1, dl(2:ninterior), d, du(1:ninterior-1), rhs, ninterior, info)
    if (info /= 0) then
      print *, 'poisson_module: dgtsv failed, info = ', info
      stop 1
    end if

    u(2:n-1) = rhs
  end subroutine solve_poisson_1d

  ! ---------------------------------------------------------------
  ! 2D solve:  -(u_xx + u_yy) = f(x,y),  u = g(x,y) on the boundary
  ! Matrix-free Conjugate Gradient.
  ! ---------------------------------------------------------------
  subroutine solve_poisson_2d(g, f, bc, u, tol, maxiter, iters_used)
    type(grid_t), intent(in)           :: g
    procedure(rhs_func_2d)             :: f
    procedure(bc_func_2d)              :: bc
    real(dp), allocatable, intent(out) :: u(:,:)
    real(dp), intent(in), optional     :: tol
    integer(ip), intent(in), optional  :: maxiter
    integer(ip), intent(out), optional :: iters_used

    integer(ip) :: nx, ny, i, j, it, itmax
    real(dp) :: tolerance, rs_old, rs_new, alpha, beta, pAp
    real(dp), allocatable :: rhs(:,:), r(:,:), p(:,:), Ap(:,:)

    nx = g%nx; ny = g%ny
    tolerance = 1.0e-10_dp
    if (present(tol)) tolerance = tol
    itmax = 5000
    if (present(maxiter)) itmax = maxiter

    allocate(u(nx,ny), rhs(nx,ny), r(nx,ny), p(nx,ny), Ap(nx,ny))
    u = 0.0_dp

    ! Impose Dirichlet boundary values exactly; they stay fixed
    ! throughout the CG iteration (see apply_laplacian_2d, which
    ! only touches interior points).
    do i = 1, nx
      u(i, 1)  = bc(g%x(i), g%y(1))
      u(i, ny) = bc(g%x(i), g%y(ny))
    end do
    do j = 1, ny
      u(1, j)  = bc(g%x(1),  g%y(j))
      u(nx, j) = bc(g%x(nx), g%y(j))
    end do

    do i = 2, nx - 1
      do j = 2, ny - 1
        rhs(i,j) = f(g%x(i), g%y(j))
      end do
    end do

    ! r0 = rhs - A u0   (A = discrete -Laplacian on interior points)
    call apply_laplacian_2d(g, u, Ap)
    r = 0.0_dp
    do i = 2, nx - 1
      do j = 2, ny - 1
        r(i,j) = rhs(i,j) - Ap(i,j)
      end do
    end do
    p = r
    rs_old = interior_dot(r, r, nx, ny)

    do it = 1, itmax
      if (sqrt(rs_old) < tolerance) exit
      call apply_laplacian_2d(g, p, Ap)
      pAp = interior_dot(p, Ap, nx, ny)
      alpha = rs_old / pAp
      do i = 2, nx - 1
        do j = 2, ny - 1
          u(i,j) = u(i,j) + alpha * p(i,j)
          r(i,j) = r(i,j) - alpha * Ap(i,j)
        end do
      end do
      rs_new = interior_dot(r, r, nx, ny)
      beta = rs_new / rs_old
      do i = 2, nx - 1
        do j = 2, ny - 1
          p(i,j) = r(i,j) + beta * p(i,j)
        end do
      end do
      rs_old = rs_new
    end do

    if (present(iters_used)) iters_used = it
  end subroutine solve_poisson_2d

  ! Applies the discrete operator A = -Laplacian (5-point stencil) to
  ! v, writing the result into interior points of Av. Boundary points
  ! of v are assumed fixed (Dirichlet) and are not touched/used as
  ! unknowns -- this is what makes the CG solver "matrix-free": we
  ! never form or store A, we only know how to multiply by it.
  subroutine apply_laplacian_2d(g, v, Av)
    type(grid_t), intent(in) :: g
    real(dp), intent(in)  :: v(:,:)
    real(dp), intent(out) :: Av(:,:)
    integer(ip) :: i, j, nx, ny
    real(dp) :: idx2, idy2

    nx = g%nx; ny = g%ny
    idx2 = 1.0_dp / (g%dx * g%dx)
    idy2 = 1.0_dp / (g%dy * g%dy)
    Av = 0.0_dp
    do j = 2, ny - 1
      do i = 2, nx - 1
        Av(i,j) = (2.0_dp * v(i,j) - v(i-1,j) - v(i+1,j)) * idx2 &
                + (2.0_dp * v(i,j) - v(i,j-1) - v(i,j+1)) * idy2
      end do
    end do
  end subroutine apply_laplacian_2d

  pure function interior_dot(a, b, nx, ny) result(s)
    real(dp), intent(in) :: a(:,:), b(:,:)
    integer(ip), intent(in) :: nx, ny
    real(dp) :: s
    integer(ip) :: i, j
    s = 0.0_dp
    do j = 2, ny - 1
      do i = 2, nx - 1
        s = s + a(i,j) * b(i,j)
      end do
    end do
  end function interior_dot

end module poisson_module
