! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
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

!> Unit test to check the functionality of the
!! Exact Riemann solver for Euler equations.
program ere_test
  use env_module, only: rk
  use atl_exact_riemann_euler_module, only: atl_ere_solState1D_type, &
    &                                       atl_ere_init, &
    &                                       atl_ere_sample, &
    &                                       atl_ere_eval_onEdge, &
    &                                       atl_ere_dump_solState

  implicit none

  type(atl_ere_solState1D_type) :: sol, solM

  real(kind=rk) :: simple_state(3)
  real(kind=rk) :: probe_vel(3)
  real(kind=rk) :: rho(3), rhoM(3)
  real(kind=rk) :: v(3), vM(3)
  real(kind=rk) :: p(3), pM(3)
  real(kind=rk) :: state_onEdge(5)


  ! Check if a constant state is maintained by the exact riemann solver:
  simple_state = [ 1.0_rk, 0.25_rk, 2.0_rk ]

  call atl_ere_init( me         = sol,          &
    &                prim_left  = simple_state, &
    &                prim_right = simple_state, &
    &                gam        = 1.4_rk        )

  probe_vel = [ -0.5_rk, 0.0_rk, 0.5_rk ]

  call atl_ere_sample( me  = sol,       &
    &                  s   = probe_vel, &
    &                  rho = rho,       &
    &                  v   = v,         &
    &                  p   = p          )

  if ( any(rho < (simple_state(1)-epsilon(rho))) &
    &  .or. any(rho > simple_state(1)+epsilon(rho))) then
    write(*,*) 'Density is not maintained for a constant state'
    write(*,*) 'is       : ', rho
    write(*,*) 'should be: ', [1.0_rk, 1.0_rk, 1.0_rk]
    write(*,*)
    write(*,*) '.....................................'
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( any(v < (simple_state(2)-epsilon(v))) &
    &  .or. any(v > simple_state(2)+epsilon(v))) then
    write(*,*) 'Velocity is not maintained for a constant state:'
    write(*,*) 'is       : ', v
    write(*,*) 'should be: ', [0.25_rk, 0.25_rk, 0.25_rk]
    write(*,*)
    write(*,*) '.....................................'
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( any(p < (simple_state(3)-epsilon(p))) &
    &  .or. any(p > simple_state(3)+epsilon(p))) then
    write(*,*) 'Pressure is not maintained for a constant state'
    write(*,*) 'is       : ', p
    write(*,*) 'should be: ', [2.0_rk, 2.0_rk, 2.0_rk]
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if

  call atl_ere_eval_onEdge( me       = sol,             &
    &                      rho_left  = simple_state(1), &
    &                      vn_left   = simple_state(2), &
    &                      v1_left   = 0.0_rk,          &
    &                      v2_left   = 0.0_rk,          &
    &                      p_left    = simple_state(3), &
    &                      rho_right = simple_state(1), &
    &                      vn_right  = simple_state(2), &
    &                      v1_right  = 0.0_rk,          &
    &                      v2_right  = 0.0_rk,          &
    &                      p_right   = simple_state(3), &
    &                      rho       = state_onEdge(1), &
    &                      vn        = state_onEdge(2), &
    &                      v1        = state_onEdge(3), &
    &                      v2        = state_onEdge(4), &
    &                      p         = state_onEdge(5)  )
  if ( (state_onEdge(1) < (simple_state(1)-epsilon(rho(2)))) &
    &  .or. (state_onEdge(1) > simple_state(1)+epsilon(rho(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct density'
    write(*,*) 'is       : ', state_onEdge(1)
    write(*,*) 'should be: ', simple_state(1)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (state_onEdge(2) < (simple_state(2)-epsilon(v(2)))) &
    &  .or. (state_onEdge(2) > simple_state(2)+epsilon(v(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct velocity'
    write(*,*) 'is       : ', state_onEdge(2)
    write(*,*) 'should be: ', simple_state(2)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (state_onEdge(5) < (simple_state(3)-epsilon(p(2)))) &
    &  .or. (state_onEdge(5) > simple_state(3)+epsilon(p(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct pressure'
    write(*,*) 'is       : ', state_onEdge(5)
    write(*,*) 'should be: ', simple_state(3)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if

  ! Check for symmetry

  ! shock to the right
  call atl_ere_init( me         = sol,                      &
    &                prim_left  = [1.0_rk, 0.0_rk, 2.0_rk], &
    &                prim_right = [1.0_rk, 0.0_rk, 1.0_rk], &
    &                gam        = 1.4_rk                    )

  call atl_ere_sample( me  = sol,       &
    &                  s   = probe_vel, &
    &                  rho = rho,       &
    &                  v   = v,         &
    &                  p   = p          )

  call atl_ere_eval_onEdge( me        = sol,             &
    &                       rho_left  = 1.0_rk,          &
    &                       vn_left   = 0.0_rk,          &
    &                       v1_left   = 0.0_rk,          &
    &                       v2_left   = 0.0_rk,          &
    &                       p_left    = 2.0_rk,          &
    &                       rho_right = 1.0_rk,          &
    &                       vn_right  = 0.0_rk,          &
    &                       v1_right  = 0.0_rk,          &
    &                       v2_right  = 0.0_rk,          &
    &                       p_right   = 1.0_rk,          &
    &                       rho       = state_onEdge(1), &
    &                       vn        = state_onEdge(2), &
    &                       v1        = state_onEdge(3), &
    &                       v2        = state_onEdge(4), &
    &                       p         = state_onEdge(5)  )
  if ( (state_onEdge(1) < (rho(2)-epsilon(rho(2)))) &
    &  .or. (state_onEdge(1) > rho(2)+epsilon(rho(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct density'
    write(*,*) 'is       : ', state_onEdge(1)
    write(*,*) 'should be: ', rho(2)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (state_onEdge(2) < (v(2)-epsilon(v(2)))) &
    &  .or. (state_onEdge(2) > v(2)+epsilon(v(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct velocity'
    write(*,*) 'is       : ', state_onEdge(2)
    write(*,*) 'should be: ', v(2)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (state_onEdge(5) < (p(2)-epsilon(p(2)))) &
    &  .or. (state_onEdge(5) > p(2)+epsilon(p(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct pressure'
    write(*,*) 'is       : ', state_onEdge(5)
    write(*,*) 'should be: ', p(2)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if

  ! shock to the left
  call atl_ere_init( me         = solM,                     &
    &                prim_left  = [1.0_rk, 0.0_rk, 1.0_rk], &
    &                prim_right = [1.0_rk, 0.0_rk, 2.0_rk], &
    &                gam        = 1.4_rk                    )

  call atl_ere_sample( me  = solM,      &
    &                  s   = probe_vel, &
    &                  rho = rhoM,      &
    &                  v   = vM,        &
    &                  p   = pM         )

  call atl_ere_eval_onEdge( me        = solM,            &
    &                       rho_left  = 1.0_rk,          &
    &                       vn_left   = 0.0_rk,          &
    &                       v1_left   = 0.0_rk,          &
    &                       v2_left   = 0.0_rk,          &
    &                       p_left    = 1.0_rk,          &
    &                       rho_right = 1.0_rk,          &
    &                       vn_right  = 0.0_rk,          &
    &                       v1_right  = 0.0_rk,          &
    &                       v2_right  = 0.0_rk,          &
    &                       p_right   = 2.0_rk,          &
    &                       rho       = state_onEdge(1), &
    &                       vn        = state_onEdge(2), &
    &                       v1        = state_onEdge(3), &
    &                       v2        = state_onEdge(4), &
    &                       p         = state_onEdge(5)  )
  if ( (state_onEdge(1) < (rhoM(2)-epsilon(rhoM(2)))) &
    &  .or. (state_onEdge(1) > rhoM(2)+epsilon(rhoM(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct density'
    write(*,*) 'is       : ', state_onEdge(1)
    write(*,*) 'should be: ', rhoM(2)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (state_onEdge(2) < (vM(2)-epsilon(vM(2)))) &
    &  .or. (state_onEdge(2) > vM(2)+epsilon(vM(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct velocity'
    write(*,*) 'is       : ', state_onEdge(2)
    write(*,*) 'should be: ', vM(2)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (state_onEdge(5) < (pM(2)-epsilon(pM(2)))) &
    &  .or. (state_onEdge(5) > pM(2)+epsilon(pM(2))) ) then
    write(*,*) 'eval_onEdge did not return the correct pressure'
    write(*,*) 'is       : ', state_onEdge(5)
    write(*,*) 'should be: ', pM(2)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if

  if ( (rho(1) < (rhoM(3)-epsilon(rho))) &
    &  .or. (rho(1) > rhoM(3)+epsilon(rho))) then
    write(*,*) 'Density is not the same with mirrored characteristics'
    write(*,*) 'is       : ', rhoM(3)
    write(*,*) 'should be: ', rho(1)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'MIRRORED:'
    call atl_ere_dump_solState(solM)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (v(1) < (-vM(3)-epsilon(v))) &
    &  .or. (v(1) > -vM(3)+epsilon(v))) then
    write(*,*) 'Velocity is not the same with mirrored characteristics'
    write(*,*) 'is       : ', vM(3)
    write(*,*) 'should be: ', -v(1)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'MIRRORED:'
    call atl_ere_dump_solState(solM)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if
  if ( (p(1) < (pM(3)-epsilon(v))) &
    &  .or. (p(1) > pM(3)+epsilon(v))) then
    write(*,*) 'Pressure is not the same with mirrored characteristics'
    write(*,*) 'is       : ', pM(3)
    write(*,*) 'should be: ', p(1)
    write(*,*) '.....................................'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'MIRRORED:'
    call atl_ere_dump_solState(solM)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if

  ! Check for critical state
  call atl_ere_init( me         = sol,                       &
    &                prim_left  = [1.4_rk, -6.0_rk, 1.0_rk], &
    &                prim_right = [1.4_rk, 6.0_rk, 1.0_rk],  &
    &                gam        = 1.4_rk                     )
  if (.not. sol%critical) then
    write(*,*) 'Did not detect a critical state!'
    write(*,*)
    call atl_ere_dump_solState(sol)
    write(*,*)
    write(*,*) 'FAILED'
    STOP
  end if

  write(*,*) 'PASSED'

end program ere_test
