! Copyright (c) 2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2018 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> This program acts as a unit test to check the proper functionality of
!! the atl_initialize_state_der subroutine in atl_compute_module.
program atl_state_der_init_test

  use env_module,               only: rk, eps, fin_env
  use tem_general_module,       only: tem_start
  use tem_logging_module,       only: logUnit
  use atl_kerneldata_module,    only: atl_kerneldata_type
  use atl_solver_param_module,  only: atl_solver_param_type
  use atl_compute_module,       only: atl_initialize_state_der

  implicit none

  type(atl_kerneldata_type) :: kerneldata
  integer :: x_l = 2, x_u = 7, y_l = 3, y_u = 9, z_l = 1, z_u = 15
  integer :: x, y, z
  real(kind=rk) :: state_sum
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codename = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = params%general      )

  allocate(kerneldata%state_der(x_l:x_u, y_l:y_u, z_l:z_u))

  kerneldata%state_der = 1.0_rk

  state_sum = 0.0_rk
  do x = x_l, x_u
    do y = y_l, y_u
      do z = z_l, z_u
        state_sum = state_sum + kerneldata%state_der(x, y, z)
      end do
    end do
  end do
!  write(logUnit(1),*)                                                      &
!    & 'Before initializing the array, the sum of all array elements is: ', &
!    & state_sum

  call atl_initialize_state_der(kerneldata%state_der)

  state_sum = 0.0_rk
  do x = x_l, x_u
    do y = y_l, y_u
      do z = z_l, z_u
        state_sum = state_sum + kerneldata%state_der(x, y, z)
      end do
    end do
  end do
!  write(logUnit(1),*)                                                     &
!    & 'After initializing the array, the sum of all array elements is: ', &
!    & state_sum

  if (state_sum < eps) then
    write(logUnit(1),*) 'PASSED'
  end if

  call fin_env()

end program atl_state_der_init_test
