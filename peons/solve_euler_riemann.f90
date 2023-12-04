! Copyright (c) 2013, 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2020 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
!
! Permission to use, copy, modify, and distribute this software for any
! purpose with or without fee is hereby granted, provided that the above
! copyright notice and this permission notice appear in all copies.
!
! THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
! WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
! ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
! WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
! ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
! OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
! **************************************************************************** !

program solve_euler_riemann
  use env_module, only: rk
  use aotus_module, only: flu_state, open_config_file, close_config, &
    &                     aot_get_val

  use atl_exact_riemann_euler_module, only: atl_ere_solState1D_type, &
    &                                   atl_ere_init, &
    &                                   atl_ere_sample, atl_ere_dump_solState

  implicit none

  character(len=256) :: filename
  type(flu_state) :: conf
  real(kind=rk) :: isen_coef
  real(kind=rk) :: prim_left(3)
  real(kind=rk) :: prim_right(3)
  real(kind=rk) :: tolerance
  real(kind=rk) :: origin_offset
  real(kind=rk) :: t_q
  real(kind=rk), allocatable :: t(:)
  real(kind=rk), allocatable :: x(:)
  real(kind=rk), allocatable :: rho(:)
  real(kind=rk), allocatable :: v(:)
  real(kind=rk), allocatable :: p(:)
  integer :: nSlices
  integer :: nPoints
  integer :: max_iter
  integer :: max_data_length
  integer :: err
  integer :: errv(3)
  integer, allocatable :: err_vec(:)
  integer :: i,j
  logical :: createFiles
  character(len=64) :: FilePrefix
  character(len=12) :: timestamp
  character(len=128) :: fname

  type(atl_ere_solState1D_type) :: sol

  write(*,*) 'Welcome to the Riemann solver for euler equations.'

  call get_command_argument(1,filename)
  if (trim(filename) .eq. '') then
    ! Warn about usage of default configuration file name.
    write(*,*) 'No Lua configuration file specified, using'
    write(*,*) 'riemann.lua as default input file.'
    filename = 'riemann.lua'
  end if

  call open_config_file(L = conf, filename = trim(filename))

  ! Description of the Riemann problem:
  call aot_get_val(L = conf, key = 'isen_coef', default=1.4_rk, &
    &              val = isen_coef, errCode = err)
  call aot_get_val(L = conf, key = 'left', &
    &              val = prim_left, errCode = errv)
  call aot_get_val(L = conf, key = 'right', &
    &              val = prim_right, errCode = errv)
  call aot_get_val(L = conf, key = 'max_iter', default = 100000, &
    &              val = max_iter, errCode = err)
  call aot_get_val(L = conf, key = 'tolerance', &
    &              default = 8*epsilon(tolerance), &
    &              val = tolerance, errCode = err)

  call aot_get_val(L = conf, key = 'origin_offset', &
    &              default = 0.0_rk, &
    &              val = origin_offset, errCode = err)

  call aot_get_val(L = conf, key = 'max_data_length', default=1000000000, &
    &              val = max_data_length, errCode = err)

  call aot_get_val(L = conf, key = 'create_files', default = .false., &
    &              val = createFiles, errCode = err)

  call aot_get_val(L = conf, key = 'file_prefix', default = 'exact', &
    &              val = FilePrefix, errCode = err)

  call atl_ere_init( me = sol, prim_left = prim_left, prim_right = prim_right, &
    &                gam = isen_coef, tolerance = tolerance,                   &
    &                nMaxIter = max_iter)

  call atl_ere_dump_solState(sol)

  call aot_get_val(L = conf, key = 't', maxlength=max_data_length, &
    &              val = t, errCode = err_vec)
  deallocate(err_vec)
  call aot_get_val(L = conf, key = 'x', maxlength=max_data_length, &
    &              val = x, errCode = err_vec)
  deallocate(err_vec)
  call close_config(conf)

  nPoints = size(x)
  nSlices = size(t)

  allocate(rho(nPoints))
  allocate(v(nPoints))
  allocate(p(nPoints))

  do i=1,nSlices
    if (t(i)>tiny(tolerance)) then
      t_q = 1.0_rk/t(i)
      call atl_ere_sample(me = sol, s = (x-origin_offset)*t_q, rho=rho, v=v, p=p)
    else
      where ((x-origin_offset)>0.0_rk)
        rho = prim_right(1)
        v = prim_right(2)
        p = prim_right(3)
      elsewhere
        rho = prim_left(1)
        v = prim_left(2)
        p = prim_left(3)
      end where
    end if

    if (createFiles) then
      write(timestamp, '(en12.3)') t(i)
      fname = trim(adjustl(filePrefix)) // '_t' &
        &     // trim(adjustl(timestamp)) // '.ascii'
      open(unit=22, file=trim(fname), recl=128)
      write(22,*) '# t=', t(i)
      write(22,*) '# x, rho, v, p'
      do j=1,nPoints
        write(22,*) x(j), rho(j), v(j), p(j)
      end do
      close(22)
    else
      write(*,*) '# t=', t(i)
      write(*,*) '# x, rho, v, p'
      do j=1,nPoints
        write(*,*) x(j), rho(j), v(j), p(j)
      end do
      write(*,*)
    end if

  end do

end program solve_euler_riemann
