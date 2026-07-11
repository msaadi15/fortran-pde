program main_transport
  use kinds_module, only: dp, ip
  use grid_module,  only: grid_t, make_grid_1d, make_grid_2d
  use transport_module, only: step_transport_1d, step_transport_2d, &
                               max_stable_dt_1d, max_stable_dt_2d
  use problem_definitions_module, only: transport1d_init_gaussian, transport1d_bc_zero, &
                                         transport2d_init_gaussian, transport2d_bc_zero
  use io_module, only: write_field_1d, write_field_2d
  implicit none

  type(grid_t) :: g1, g2
  real(dp), allocatable :: u1(:), u2(:,:)
  real(dp) :: vx, vy, D, dt, t, t_final
  integer(ip) :: i, j, nstep, nsnap, snap_every
  character(len=64) :: fname

  ! ------------------- 1D advection-diffusion -------------------
  print *, '=== 1D transport: du/dt + vx du/dx = D d2u/dx2 ==='
  g1 = make_grid_1d(nx=201_ip, xmin=0.0_dp, xmax=1.0_dp)
  vx = 0.5_dp
  D  = 0.001_dp
  t_final = 1.0_dp

  allocate(u1(g1%nx))
  do i = 1, g1%nx
    u1(i) = transport1d_init_gaussian(g1%x(i))
  end do

  dt = max_stable_dt_1d(g1, vx, D)
  nstep = ceiling(t_final / dt)
  dt = t_final / real(nstep, dp)   ! adjust so we land exactly on t_final
  snap_every = max(1_ip, nstep / 20_ip)

  print *, 'dt =', dt, ' nstep =', nstep
  call write_field_1d('results/transport1d_t0000.dat', g1, u1)

  t = 0.0_dp
  nsnap = 0
  do i = 1, nstep
    t = t + dt
    call step_transport_1d(g1, u1, vx, D, dt, t, transport1d_bc_zero)
    if (mod(i, snap_every) == 0 .or. i == nstep) then
      nsnap = nsnap + 1
      write(fname, '(A,I4.4,A)') 'results/transport1d_t', nsnap, '.dat'
      call write_field_1d(trim(fname), g1, u1)
    end if
  end do
  print *, 'wrote', nsnap + 1, '1D snapshots to results/transport1d_t*.dat'
  print *

  ! ------------------- 2D advection-diffusion -------------------
  print *, '=== 2D transport: du/dt + v.grad(u) = D*Laplacian(u) ==='
  g2 = make_grid_2d(nx=101_ip, ny=101_ip, xmin=0.0_dp, xmax=1.0_dp, ymin=0.0_dp, ymax=1.0_dp)
  vx = 0.4_dp; vy = 0.2_dp
  D  = 0.0015_dp
  t_final = 1.0_dp

  allocate(u2(g2%nx, g2%ny))
  do j = 1, g2%ny
    do i = 1, g2%nx
      u2(i,j) = transport2d_init_gaussian(g2%x(i), g2%y(j))
    end do
  end do

  dt = max_stable_dt_2d(g2, vx, vy, D)
  nstep = ceiling(t_final / dt)
  dt = t_final / real(nstep, dp)
  snap_every = max(1_ip, nstep / 20_ip)

  print *, 'dt =', dt, ' nstep =', nstep
  call write_field_2d('results/transport2d_t0000.dat', g2, u2)

  t = 0.0_dp
  nsnap = 0
  do i = 1, nstep
    t = t + dt
    call step_transport_2d(g2, u2, vx, vy, D, dt, t, transport2d_bc_zero)
    if (mod(i, snap_every) == 0 .or. i == nstep) then
      nsnap = nsnap + 1
      write(fname, '(A,I4.4,A)') 'results/transport2d_t', nsnap, '.dat'
      call write_field_2d(trim(fname), g2, u2)
    end if
  end do
  print *, 'wrote', nsnap + 1, '2D snapshots to results/transport2d_t*.dat'

end program main_transport
