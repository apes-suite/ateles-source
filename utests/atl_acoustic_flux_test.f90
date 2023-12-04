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

! This is the unit test for the flux function 
program atl_acoustic_flux_test
  use env_module,                   only: rk, fin_env
  use tem_logging_module,           only: logUnit
  use atl_eqn_acoustic_module,      only: atl_acoustic_type

  use atl_solver_param_module,      only: atl_solver_param_type 
  use tem_general_module,           only: tem_start
  use atl_acoustic_physFlux_module, only: atl_acoustic_physFlux
  use atl_acoustic_numflux_module,     only: atl_acoustic_numFlux_oneDir

implicit none

  ! physical flux check
  real(kind=rk) :: F_res,G_res,H_res
  real(kind=rk) :: F_ref(4),G_ref(4),H_ref(4)
  real(kind=rk) :: F(4),G(4),H(4), state(4)
  ! numerical flux check
  real(kind=rk) :: F_num_res,G_num_res,H_num_res
  real(kind=rk) :: F_num_ref(4),G_num_ref(4),H_num_ref(4)
  real(kind=rk) :: F_num(4),G_num(4),H_num(4), state_left(4), state_right(4)

  type(atl_solver_param_type) :: params
  type(atl_acoustic_type) :: acoustic

  !init the acoustic type
  acoustic%density_0 = 1.225     
  allocate(acoustic%velocity_0(3))
  acoustic%velocity_0(1) = 0.0 
  acoustic%velocity_0(2) = 0.0
  acoustic%velocity_0(3) = 0.0
  acoustic%pressure_0 = 100000.0    
  acoustic%speedOfSound=sqrt(acoustic%pressure_0 / acoustic%density_0) 


  ! Init the Treelm environment, needed to init the log Unit
  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = params%general      )

  ! init res to 0
  F_res = 0.0_rk
  G_res = 0.0_rk
  H_res = 0.0_rk

  ! init the state vector
  state(1) = 1.0_rk
  state(2) = 1.0_rk
  state(3) = 1.0_rk
  state(4) = 1.0_rk

  ! physcial flux function
  call check_physical_flux(state,acoustic,F,G,H)

  ! calcualte the reference physical flux
  F_ref(1) = acoustic%density_0*state(2)
  F_ref(2) = acoustic%speedOfSound**2/acoustic%density_0 * state(1)
  F_ref(3) = state(3)
  F_ref(4) = state(4)

  G_ref(1) = acoustic%density_0*state(3)
  G_ref(2) = state(2)
  G_ref(3) = acoustic%speedOfSound**2/acoustic%density_0 *state(1)
  G_ref(4) = state(4)

  H_ref(1) = acoustic%density_0*state(4)
  H_ref(2) = state(2)
  H_ref(3) = state(3)
  H_ref(4)= acoustic%speedOfSound**2/acoustic%density_0 * state(1)

  ! calcualte the residuals
  F_res = maxval(abs(F_ref-F) ) 
  G_res = maxval(abs(G_ref-G) ) 
  H_res = maxval(abs(H_ref-H) ) 

  ! write the results
  write(logUnit(1),*) ' Check of the physicalflux'
  write(logUnit(1),*) 'state= ', state
  write(logUnit(1),*) 'F= ', F
  write(logUnit(1),*) 'G= ', G
  write(logUnit(1),*) 'H= ', H
  write(logUnit(1),*) 'F_ref = ', F_ref
  write(logUnit(1),*) 'G_ref = ', G_ref
  write(logUnit(1),*) 'H_ref = ', H_ref
  write(logUnit(1),*) 'F_res = ', F_res
  write(logUnit(1),*) 'G_res = ', G_res
  write(logUnit(1),*) 'H_res = ', H_res
  write(logUnit(10),*) '--------------------------- DONE'
  
  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  if( (F_res.lt.1e-08) .and. (G_res.lt.1e-08) .and. (H_res.lt.1e-08) ) then
    write(logUnit(1),*) 'PASSED'
  end if
 
  !!! numerical flux function

  ! init res to 0
  F_num_res = 0.0_rk
  G_num_res = 0.0_rk
  H_num_res = 0.0_rk

  ! init the state vector
  state_left(1) = 1.0_rk
  state_left(2) = 1.0_rk
  state_left(3) = 1.0_rk
  state_left(4) = 1.0_rk
  state_right(1) = 1.0_rk
  state_right(2) = 1.0_rk
  state_right(3) = 1.0_rk
  state_right(4) = 1.0_rk
  
  ! numerical flux function
  call check_numerical_flux(state_left, state_right,acoustic,F_num,G_num,H_num)

  ! calculate the reference numerical flux
  F_num_ref(1) = 1.225*0.5*2
  F_num_ref(2) = acoustic%speedOfSound**2/acoustic%density_0 * 0.5 * 2
  F_num_ref(3) = state_left(3)
  F_num_ref(4) = state_left(4)

  G_num_ref(1) = 1.225*0.5*2
  G_num_ref(2) = state_left(2)
  G_num_ref(3) = acoustic%speedOfSound**2/acoustic%density_0 * 0.5 * 2
  G_num_ref(4) = state_left(4)

  H_num_ref(1) =1.225*0.5*2
  H_num_ref(2) = state_left(2)
  H_num_ref(3) = state_left(3)
  H_num_ref(4)= acoustic%speedOfSound**2/acoustic%density_0 * 0.5 *2

  ! calculate the residual
  F_num_res = maxval(abs(F_num_ref-F_num) ) 
  G_num_res = maxval(abs(G_num_ref-G_num) ) 
  H_num_res = maxval(abs(H_num_ref-H_num) ) 

  ! write the results
  write(logUnit(1),*) ' Check of the numerical flux'
  write(logUnit(1),*) 'state_left= ', state_left
  write(logUnit(1),*) 'state_right= ', state_right
  write(logUnit(1),*) 'F_num= ', F_num
  write(logUnit(1),*) 'G_num= ', G_num
  write(logUnit(1),*) 'H_num= ', H_num
  write(logUnit(1),*) 'F_num_ref = ', F_num_ref
  write(logUnit(1),*) 'G_num_ref = ', G_num_ref
  write(logUnit(1),*) 'H_num_ref = ', H_num_ref
  write(logUnit(1),*) 'F_num_res = ', F_num_res
  write(logUnit(1),*) 'G_num_res = ', G_num_res
  write(logUnit(1),*) 'H_num_res = ', H_num_res
  write(logUnit(10),*) '--------------------------- DONE'

  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  if( (F_num_res.lt.1e-08) .and. (G_num_res.lt.1e-08) &
    &                      .and. (H_num_res.lt.1e-08) ) then
    write(logUnit(1),*) 'PASSED'
  end if
  
  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  call fin_env()

contains

  subroutine check_physical_flux(state, acoustic, F,G,H)
    !-------------------------------------------------------------------------
    !> State to compute the fluxes (rho, u, v, w)
    real(kind=rk), intent(in) :: state(4)
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    !> The resulting flux in x direction
    real(kind=rk) :: F(4)
    !> The resulting flux in y direction
    real(kind=rk) :: G(4)
    !> The resulting flux in z direction
    real(kind=rk) :: H(4)
    !--------------------------------------------------------------------------
    integer :: rotX(4), rotY(4), rotZ(4)
    !--------------------------------------------------------------------------

    ! get the rotation for the physical flux calculation in y and z direction.
    rotX = [ 1, 2, 3, 4]
    rotY = [ 1, 3, 4, 2]
    rotZ = [ 1, 4, 2, 3] 

     ! Calculate the physical flux point by point within this cell - xdirection
     F(rotX) = atl_acoustic_physFlux(&
       &       state = state(rotX),  &
       &       acoustic = acoustic,  &
       &       idir = 1              )

     ! Calculate the physical flux point by point within this cell - ydirection
     G(rotY) = atl_acoustic_physFlux( & 
       &       state = state(rotY),   &
       &       acoustic = acoustic,   &
       &       idir = 2               )

     ! Calculate the physical flux point by point within this cell - zdirection
     H(rotZ) = atl_acoustic_physFlux( &
       &       state = state(rotZ),   &
       &       acoustic = acoustic,   &
       &       idir = 3               )

  end subroutine check_physical_flux

  subroutine check_numerical_flux(state_left,state_right,acoustic,F,G,H)
    !--------------------------------------------------------------------------
    !> State to compute the fluxes (rho, u, v, w)
    real(kind=rk), intent(in) :: state_left(4),state_right(4)
    !> Datatype for acoustic equation include all background data
    type(atl_acoustic_type), intent(in) :: acoustic
    !> The resulting flux in x direction
    real(kind=rk) :: F(4)
    !> The resulting flux in y direction
    real(kind=rk) :: G(4)
    !> The resulting flux in z direction
    real(kind=rk) :: H(4)
    !--------------------------------------------------------------------------
    integer :: rotX(4), rotY(4), rotZ(4)
    !--------------------------------------------------------------------------

    ! get the rotation for the physical flux calculation in y and z direction.
    rotX = [ 1, 2, 3, 4]
    rotY = [ 1, 3, 4, 2]
    rotZ = [ 1, 4, 2, 3] 
    
    F(rotX) = atl_acoustic_numFlux_oneDir(left = state_left(rotX), &
      &                      right = state_right(rotX),            &
      &                      acoustic = acoustic,                  &
      &                      idir = 1                              ) 
    
    G(rotY) = atl_acoustic_numFlux_oneDir(left = state_left(rotY), &
      &                      right = state_right(rotY),            &
      &                      acoustic = acoustic,                  &
      &                      idir = 2                              ) 

    H(rotZ) = atl_acoustic_numFlux_oneDir(left = state_left(rotZ), &
      &                      right = state_right(rotZ),            &
      &                      acoustic = acoustic,                  &
      &                      idir = 3                              ) 

  end subroutine check_numerical_flux

end program atl_acoustic_flux_test
