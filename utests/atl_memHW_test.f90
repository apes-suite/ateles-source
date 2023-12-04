! Copyright (c) 2015 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
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

program atl_memHW_test

  use env_module,  only: rk, my_status_int
  use tem_aux_module, only: tem_global_vmhwm
  use tem_logging_module, only: logUnit
  use tem_general_module, only: tem_general_type, tem_start

  implicit none
 
  integer :: mem
  real(kind=rk) :: mem_global(3)
  type(tem_general_type) :: tgen

  call tem_start(codename = 'memHWM_test', &
    &            version  = 'v1',          &
    &            general  = tgen           )

  mem = my_status_int('VmHWM:')
  write(*,'(a,i0,a)') 'MemHWM: ', mem, ' kByte'

  mem_global = tem_global_vmhwm()
  write(logUnit(1),'(a,f0.3,a)') 'global MemHWM min: ', mem_global(1), ' MByte'
  write(logUnit(1),'(a,f0.3,a)') 'global MemHWM max: ', mem_global(2), ' MByte'
  write(logUnit(1),'(a,f0.3,a)') 'global MemHWM sum: ', mem_global(3), ' MByte'
  write(*,*) 'PASSED'

end program atl_memHW_test

