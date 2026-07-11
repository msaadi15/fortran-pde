module io_module
  ! Minimal, dependency-free output: plain whitespace-separated text
  ! files. Easy to load from Python (numpy.loadtxt / pandas), gnuplot,
  ! or a Manim script for animation, without needing any Fortran-side
  ! plotting library.
  use kinds_module, only: dp, ip
  use grid_module,  only: grid_t
  implicit none
  private
  public :: write_field_1d, write_field_2d

contains

  ! Writes "x  u(x)" pairs, one per line.
  subroutine write_field_1d(filename, g, u)
    character(len=*), intent(in) :: filename
    type(grid_t), intent(in) :: g
    real(dp), intent(in) :: u(:)
    integer(ip) :: i, unit_id

    open(newunit=unit_id, file=filename, status='replace', action='write')
    do i = 1, g%nx
      write(unit_id, '(ES16.8, 1X, ES16.8)') g%x(i), u(i)
    end do
    close(unit_id)
  end subroutine write_field_1d

  ! Writes "x  y  u(x,y)" triples, with a blank line between each row
  ! of constant y -- this is the format gnuplot's `splot` and
  ! numpy-based loaders both expect for a structured grid.
  subroutine write_field_2d(filename, g, u)
    character(len=*), intent(in) :: filename
    type(grid_t), intent(in) :: g
    real(dp), intent(in) :: u(:,:)
    integer(ip) :: i, j, unit_id

    open(newunit=unit_id, file=filename, status='replace', action='write')
    do j = 1, g%ny
      do i = 1, g%nx
        write(unit_id, '(ES16.8, 1X, ES16.8, 1X, ES16.8)') g%x(i), g%y(j), u(i,j)
      end do
      write(unit_id, *)   ! blank line separates scanlines
    end do
    close(unit_id)
  end subroutine write_field_2d

end module io_module
