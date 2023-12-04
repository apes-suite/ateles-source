! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> Unit test to check functionallity of primal
!! to co-volume projections.
!! \author{Jens Zudrop}
program atl_covolumeToPrimal_1d_test
  use env_module,               only: rk, fin_env
  use tem_general_module,       only: tem_start
  use tem_logging_module,       only: logUnit
  use atl_solver_param_module,  only: atl_solver_param_type
  use atl_stabilize_module,     only: atl_covolume_to_primal_projection_1d
  use atl_covolume_module,      only: atl_covolume_type
  use atl_scheme_module,        only: atl_scheme_type, atl_modg_1d_scheme_prp
  use atl_space_basis,          only: atl_init_spacebasis

  implicit none

  real(kind=rk) :: res, newRes
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codename = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = params%general      )

  ! Check the transformation for different number of polynomial coefficients.
  res = 0.0_rk
  write(logUnit(10),*) 'CHECKING 1D COVOLUME->PRIMAL PROJECTION FOR'
  call check(power = 2, res = newRes)
  if(newRes.gt.res) then
    res = newRes
  end if
  write(logUnit(10),*) '... DONE!'

  if(res.lt.1e-08_rk) then
    write(logUnit(1),*) 'PASSED'
  end if

  call fin_env()

contains

  subroutine check(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    real(kind=rk), allocatable :: left(:,:), right(:,:), state(:,:), &
                                & primal_x(:,:), ref_primal_x(:,:)
    integer :: maxPolyDegree
    type(atl_covolume_type) :: filter
    integer :: level
    type(atl_scheme_type), allocatable :: scheme_list(:)

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree = 2**power - 1  ! maxPolyDegree+1 has to be a power of 2

    level = 4
    allocate(scheme_list(level:level))
    scheme_list(level)%scheme = atl_modg_1d_scheme_prp
    scheme_list(level)%modg_1d%maxPolyDegree = maxPolyDegree
    call atl_init_spacebasis(scheme_list = scheme_list, minlevel = level, maxlevel = level)

    ! Create the Legendre expansion coefficients
    allocate(left((maxPolyDegree+1),1))
    allocate(right((maxPolyDegree+1),1))
    allocate(state((maxPolyDegree+1),1))
    allocate(ref_primal_x((maxPolyDegree+1),1))
    allocate(primal_x((maxPolyDegree+1),1))
    left(:,:) = 1.0_rk
    right(:,:) = 0.5_rk
    state(:,:) = 0.34_rk

    ! Set up the filter
    filter%alpha = 36.0_rk
    filter%order = 14.0_rk
    filter%beta = 1.0_rk

    ! Co-volume projection in x direction
    primal_x = atl_covolume_to_primal_projection_1d( &
      & left       = left,                           &
      & right      = right,                          &
      & filter     = filter,                         &
      & scheme     = scheme_list(level),             &
      & maxPolyDeg = maxPolyDegree,                  &
      & state      = state,                          &
      & nScalars   = 1                               )

    !do degCo = 1, maxPolyDegree+1
    !    pos_x_cov = posOfModgCoeffQTens(degCo, 1, 1, maxPolyDegree)
    !    write(*,*) 'ref_primal_x(',pos_x_cov,',1)=', primal_x(pos_x_cov,1),'_rk'
    !end do

    ! Define the reference result
    ref_primal_x(1 ,1)=  0.84374999999995415_rk
    ref_primal_x(2 ,1)= -0.37499999999997785_rk
    ref_primal_x(3 ,1)= -0.89062499999995426_rk
    ref_primal_x(4 ,1)=  0.74999999999995148_rk

    ! Find the maximal deviation from the reference
    res = maxval(abs(primal_x - ref_primal_x))

  end subroutine


end program atl_covolumeToPrimal_1d_test
