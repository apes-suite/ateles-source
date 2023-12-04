! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> Checking the viscous numerical fluxes in Navier-Stokes
!! with a simple one-dimensional state.
!! 2D and 3D should yield the same results for this, and this utest
!! just checks for that.
program atl_visc2v3_test
  use env_module, only: rk
  use atl_viscNumFlux_Nvrstk_2d_module, only: atl_viscNavierStokes_2d
  use atl_viscNumFlux_Nvrstk_module, only: atl_viscNavierStokes

  implicit none

  integer :: i

  logical :: success

  real(kind=rk) :: stateL_3D(5)
  real(kind=rk) :: stateR_3D(5)
  real(kind=rk) :: gradL_3D(5,3)
  real(kind=rk) :: gradR_3D(5,3)
  real(kind=rk) :: flux_3D(5)

  real(kind=rk) :: stateL_2D(4)
  real(kind=rk) :: stateR_2D(4)
  real(kind=rk) :: gradL_2D(4,2)
  real(kind=rk) :: gradR_2D(4,2)
  real(kind=rk) :: flux_2D(4)
  real(kind=rk) :: noise(3)

  success = .true.

  write(*,*) 'Constant state with 0 velocity, no gradient.'
  stateL_3D = [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  stateR_3D = [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  gradL_3D = 0.0_rk
  gradR_3D = 0.0_rk

  stateL_2D = [1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  stateR_2D = [1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  gradL_2D = 0.0_rk
  gradR_2D = 0.0_rk

  call atl_viscNavierStokes(     &
    &    left       = stateL_3D, &
    &    left_grad  = gradL_3D,  &
    &    right      = stateR_3D, &
    &    right_grad = gradR_3D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_3D    )

  write(*,*) flux_3D

  call atl_viscNavierStokes_2D(  &
    &    left       = stateL_2D, &
    &    left_grad  = gradL_2D,  &
    &    right      = stateR_2D, &
    &    right_grad = gradR_2D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_2D    )

  write(*,*) flux_2D

  success = success                                                &
    &     .and. (maxval(abs(flux_3D([1,2,5]) - flux_2D([1,2,4])))  &
    &            < epsilon(flux_2D(1))                             )

  write(*,*) 'Constant state with 1 velocity, no gradient.'
  stateL_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  stateR_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  gradL_3D = 0.0_rk
  gradR_3D = 0.0_rk

  stateL_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  stateR_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  gradL_2D = 0.0_rk
  gradR_2D = 0.0_rk

  call atl_viscNavierStokes(     &
    &    left       = stateL_3D, &
    &    left_grad  = gradL_3D,  &
    &    right      = stateR_3D, &
    &    right_grad = gradR_3D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_3D    )

  write(*,*) flux_3D

  call atl_viscNavierStokes_2D(  &
    &    left       = stateL_2D, &
    &    left_grad  = gradL_2D,  &
    &    right      = stateR_2D, &
    &    right_grad = gradR_2D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_2D    )

  write(*,*) flux_2D

  success = success                                                &
    &     .and. (maxval(abs(flux_3D([1,2,5]) - flux_2D([1,2,4])))  &
    &            < epsilon(flux_2D(1))                             )

  write(*,*) 'Constant state with 1 velocity, 0.5 gradient in X for density.'
  stateL_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  stateR_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  gradL_3D = 0.0_rk
  gradR_3D = 0.0_rk
  gradL_3D(1,1) = 0.5_rk
  gradR_3D(1,1) = 0.5_rk

  stateL_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  stateR_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  gradL_2D = 0.0_rk
  gradR_2D = 0.0_rk
  gradL_2D(1,1) = 0.5_rk
  gradR_2D(1,1) = 0.5_rk

  call atl_viscNavierStokes(     &
    &    left       = stateL_3D, &
    &    left_grad  = gradL_3D,  &
    &    right      = stateR_3D, &
    &    right_grad = gradR_3D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_3D    )

  write(*,*) flux_3D

  call atl_viscNavierStokes_2D(  &
    &    left       = stateL_2D, &
    &    left_grad  = gradL_2D,  &
    &    right      = stateR_2D, &
    &    right_grad = gradR_2D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_2D    )

  write(*,*) flux_2D

  success = success                                                &
    &     .and. (maxval(abs(flux_3D([1,2,5]) - flux_2D([1,2,4])))  &
    &            < epsilon(flux_2D(1))                             )

  write(*,*) 'Constant state with 1 velocity, 0.5 gradient in X for velocity.'
  stateL_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  stateR_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  gradL_3D = 0.0_rk
  gradR_3D = 0.0_rk
  gradL_3D(2,1) = 0.5_rk
  gradR_3D(2,1) = 0.5_rk

  stateL_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  stateR_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  gradL_2D = 0.0_rk
  gradR_2D = 0.0_rk
  gradL_2D(2,1) = 0.5_rk
  gradR_2D(2,1) = 0.5_rk

  call atl_viscNavierStokes(     &
    &    left       = stateL_3D, &
    &    left_grad  = gradL_3D,  &
    &    right      = stateR_3D, &
    &    right_grad = gradR_3D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_3D    )

  write(*,*) flux_3D

  call atl_viscNavierStokes_2D(  &
    &    left       = stateL_2D, &
    &    left_grad  = gradL_2D,  &
    &    right      = stateR_2D, &
    &    right_grad = gradR_2D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_2D    )

  write(*,*) flux_2D

  success = success                                                &
    &     .and. (maxval(abs(flux_3D([1,2,5]) - flux_2D([1,2,4])))  &
    &            < epsilon(flux_2D(1))                             )

  write(*,*) 'Constant state with 1 velocity, 0.5 gradient in X for energy.'
  stateL_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  stateR_3D = [1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk]
  gradL_3D = 0.0_rk
  gradR_3D = 0.0_rk
  gradL_3D(5,1) = 0.5_rk
  gradR_3D(5,1) = 0.5_rk

  stateL_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  stateR_2D = [1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk]
  gradL_2D = 0.0_rk
  gradR_2D = 0.0_rk
  gradL_2D(4,1) = 0.5_rk
  gradR_2D(4,1) = 0.5_rk

  call atl_viscNavierStokes(     &
    &    left       = stateL_3D, &
    &    left_grad  = gradL_3D,  &
    &    right      = stateR_3D, &
    &    right_grad = gradR_3D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_3D    )

  write(*,*) flux_3D

  call atl_viscNavierStokes_2D(  &
    &    left       = stateL_2D, &
    &    left_grad  = gradL_2D,  &
    &    right      = stateR_2D, &
    &    right_grad = gradR_2D,  &
    &    mu         = 1.0_rk,    &
    &    lambda     = 1.0_rk,    &
    &    thermCond  = 1.0_rk,    &
    &    heatCap    = 1.0_rk,    &
    &    penaltyIP  = 4.0_rk,    &
    &    flux       = flux_2D    )

  write(*,*) flux_2D

  success = success                                                &
    &     .and. (maxval(abs(flux_3D([1,2,5]) - flux_2D([1,2,4])))  &
    &            < epsilon(flux_2D(1))                             )

  write(*,*) 'Randomized.'

  stateL_3D = 0.0_rk
  stateR_3D = 0.0_rk
  gradL_3D = 0.0_rk
  gradR_3D = 0.0_rk

  do i=1,10
    write(*,*) '  random run ', i
    call random_number(noise)
    stateL_3D(1) = 1.0_rk + 0.5_rk * (noise(1) - 0.5_rk)
    stateL_3D(2) = noise(2)
    stateL_3D(5) = 1.0_rk + 0.5_rk * (noise(3) - 0.5_rk)
    stateL_2D(1) = 1.0_rk + 0.5_rk * (noise(1) - 0.5_rk)
    stateL_2D(2) = noise(2)
    stateL_2D(4) = 1.0_rk + 0.5_rk * (noise(3) - 0.5_rk)

    call random_number(noise)
    stateR_3D(1) = 1.0_rk + 0.5_rk * (noise(1) - 0.5_rk)
    stateR_3D(2) = noise(2)
    stateR_3D(5) = 1.0_rk + 0.5_rk * (noise(3) - 0.5_rk)
    stateR_2D(1) = 1.0_rk + 0.5_rk * (noise(1) - 0.5_rk)
    stateR_2D(2) = noise(2)
    stateR_2D(4) = 1.0_rk + 0.5_rk * (noise(3) - 0.5_rk)

    call random_number(noise)
    gradL_3D([1,2,5],1) = noise
    gradL_2D([1,2,4],1) = noise
    call random_number(noise)
    gradR_3D([1,2,5],1) = noise
    gradR_2D([1,2,4],1) = noise

    call atl_viscNavierStokes(     &
      &    left       = stateL_3D, &
      &    left_grad  = gradL_3D,  &
      &    right      = stateR_3D, &
      &    right_grad = gradR_3D,  &
      &    mu         = 1.0_rk,    &
      &    lambda     = 1.0_rk,    &
      &    thermCond  = 1.0_rk,    &
      &    heatCap    = 1.0_rk,    &
      &    penaltyIP  = 4.0_rk,    &
      &    flux       = flux_3D    )

    write(*,*) flux_3D

    call atl_viscNavierStokes_2D(  &
      &    left       = stateL_2D, &
      &    left_grad  = gradL_2D,  &
      &    right      = stateR_2D, &
      &    right_grad = gradR_2D,  &
      &    mu         = 1.0_rk,    &
      &    lambda     = 1.0_rk,    &
      &    thermCond  = 1.0_rk,    &
      &    heatCap    = 1.0_rk,    &
      &    penaltyIP  = 4.0_rk,    &
      &    flux       = flux_2D    )

    write(*,*) flux_2D

    success = success                                                &
      &     .and. (maxval(abs(flux_3D([1,2,5]) - flux_2D([1,2,4])))  &
      &            < epsilon(flux_2D(1))                             )
  end do

  if (success) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program atl_visc2v3_test
