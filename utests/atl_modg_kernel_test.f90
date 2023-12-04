! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
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

! This is the unit test for the flux function 
program atl_modg_kernel_test
  use env_module,               only: fin_env
  use tem_logging_module,       only: logUnit
  use tem_general_module,       only: tem_start

  use atl_solver_param_module,  only: atl_solver_param_type 
  use atl_modg_kernel_module,   only: atl_modg_kernel_utests

implicit none

  type(atl_solver_param_type) :: params
  logical :: passed

  ! Init the Treelm environment, needed to init the log Unit
  ! Init the Treelm environment, needed to init the log Unit
  call tem_start( codeName = 'MODG_Kernel_test', &
    &             version  = 'utest',            &
    &             general  = params%general      )

  call atl_modg_kernel_utests(passed)

  if (passed) write(logUnit(1),*) 'PASSED'

  call fin_env()


end program atl_modg_kernel_test
