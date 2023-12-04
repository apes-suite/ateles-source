! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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

?? include "ply_dof_module.inc"
!> Unit test to check functionallity of primal
!! to co-volume projections.
!! \author{Jens Zudrop}
program atl_covolumeToPrimal_test
  use env_module,               only: rk, fin_env
  use tem_general_module,       only: tem_start
  use tem_logging_module,       only: logUnit
  use atl_solver_param_module,  only: atl_solver_param_type
  use atl_stabilize_module,     only: atl_covolume_to_primal_projection
  use atl_covolume_module,      only: atl_covolume_type
  use atl_scheme_module,        only: atl_scheme_type, atl_modg_scheme_prp
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
  write(logUnit(10),*) 'CHECKING COVOLUME->PRIMAL PROJECTION FOR'
  call check(power = 2, res = newRes)
  if(newRes.gt.res) then
    res = newRes
  end if
  write(logUnit(10),*) '... DONE!'

  if(res.lt.1e-08_rk) then
    write(logUnit(1),*) 'PASSED'
  else
    write (logUnit(1),*) 'res = ', res
    write (*,*) 'res = ', res
  end if


  call fin_env()

contains

  subroutine check(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    real(kind=rk), allocatable :: left(:,:), right(:,:), state(:,:), &
                                & primal_x(:,:), primal_y(:,:), primal_z(:,:), &
                                & res_x_y(:,:), res_x_z(:,:), ref_primal_x(:,:)
    integer :: maxPolyDegree
    type(atl_covolume_type) :: filter
    integer :: level, pos_x_cov, pos_y_cov, pos_z_cov, deg1, deg2, degCo
    type(atl_scheme_type), allocatable :: scheme_list(:)

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree = 2**power - 1  ! maxPolyDegree+1 has to be a power of 2

    level = 4
    allocate(scheme_list(level:level))
    scheme_list(level)%scheme = atl_modg_scheme_prp
    scheme_list(level)%modg%maxPolyDegree = maxPolyDegree
    call atl_init_spacebasis(scheme_list = scheme_list, minlevel = level, maxlevel = level)

    ! Create the Legendre expansion coefficients
    allocate(left((maxPolyDegree+1)**3,1))
    allocate(right((maxPolyDegree+1)**3,1))
    allocate(state((maxPolyDegree+1)**3,1))
    allocate(ref_primal_x((maxPolyDegree+1)**3,1))
    allocate(primal_x((maxPolyDegree+1)**3,1))
    allocate(primal_y((maxPolyDegree+1)**3,1))
    allocate(primal_z((maxPolyDegree+1)**3,1))
    left(:,:) = 1.0_rk
    right(:,:) = 0.5_rk
    state(:,:) = 0.34_rk

    ! Set up the filter
    filter%alpha = 36.0_rk
    filter%order = 14.0_rk
    filter%beta = 1.0_rk

    ! Co-volume projection in x direction
    primal_x = atl_covolume_to_primal_projection( &
      & left       = left,                        &
      & right      = right,                       &
      & filter     = filter,                      &
      & scheme     = scheme_list(level),          &
      & maxPolyDeg = maxPolyDegree,               &
      & dir        = 1,                           &
      & state      = state,                       &
      & nScalars   = 1                            )

    ! Co-volume projection in y direction
    primal_y = atl_covolume_to_primal_projection( &
      & left       = left,                        &
      & right      = right,                       &
      & filter     = filter,                      &
      & scheme     = scheme_list(level),          &
      & maxPolyDeg = maxPolyDegree,               &
      & dir        = 2,                           &
      & state      = state,                       &
      & nScalars   = 1                            )

    ! Co-volume projection in z direction
    primal_z = atl_covolume_to_primal_projection( &
      & left       = left,                        &
      & right      = right,                       &
      & filter     = filter,                      &
      & scheme     = scheme_list(level),          &
      & maxPolyDeg = maxPolyDegree,               &
      & dir        = 3,                           &
      & state      = state,                       &
      & nScalars   = 1                            )

    ! Check rotational invariance
    allocate(res_x_y((maxPolyDegree+1)**3,1))
    allocate(res_x_z((maxPolyDegree+1)**3,1))
    do deg1 = 1, maxPolyDegree+1
      do deg2 = 1, maxPolyDegree+1
        do degCo = 1, maxPolyDegree+1
?? copy :: posOfModgCoeffQTens(degCo, deg1, deg2, maxPolyDegree, pos_x_cov)
?? copy :: posOfModgCoeffQTens(deg1, degCo, deg2, maxPolyDegree, pos_y_cov)
?? copy :: posOfModgCoeffQTens(deg1, deg2, degCo, maxPolyDegree, pos_z_cov)
          res_x_y(pos_x_cov,:) = primal_x(pos_x_cov,:) - primal_y(pos_y_cov,:)
          res_x_z(pos_x_cov,:) = primal_x(pos_x_cov,:) - primal_z(pos_z_cov,:)
        end do
      end do
    end do

    !do degCo = 1, maxPolyDegree+1
    !  do deg1 = 1, maxPolyDegree+1
    !    do deg2 = 1, maxPolyDegree+1
!?? copy :: posOfModgCoeffQTens(degCo, deg1, deg2, maxPolyDegree, pos_x_cov)
    !      write(*,*) 'ref_primal_x(',pos_x_cov,',1)=', primal_x(pos_x_cov,1),'_rk'
    !    end do
    !  end do
    !end do

    ! Define the reference result
    ref_primal_x( 1 ,1)=  0.84374999999995415_rk
    ref_primal_x(17 ,1)=  0.84374999999995415_rk
    ref_primal_x(33 ,1)=  0.84374999999995415_rk
    ref_primal_x(49 ,1)=  0.84374999999995415_rk
    ref_primal_x( 5 ,1)=  0.84374999999995415_rk
    ref_primal_x(21 ,1)=  0.84374999999995415_rk
    ref_primal_x(37 ,1)=  0.84374999999995415_rk
    ref_primal_x(53 ,1)=  0.84374999999995415_rk
    ref_primal_x( 9 ,1)=  0.84374999999995415_rk
    ref_primal_x(25 ,1)=  0.84374999999995415_rk
    ref_primal_x(41 ,1)=  0.84374999999995415_rk
    ref_primal_x(57 ,1)=  0.84374999999995415_rk
    ref_primal_x(13 ,1)=  0.84374999999995415_rk
    ref_primal_x(29 ,1)=  0.84374999999995415_rk
    ref_primal_x(45 ,1)=  0.84374999999995415_rk
    ref_primal_x(61 ,1)=  0.84374999999995415_rk
    ref_primal_x( 2 ,1)= -0.37499999999997785_rk
    ref_primal_x(18 ,1)= -0.37499999999997785_rk
    ref_primal_x(34 ,1)= -0.37499999999997785_rk
    ref_primal_x(50 ,1)= -0.37499999999997785_rk
    ref_primal_x( 6 ,1)= -0.37499999999997785_rk
    ref_primal_x(22 ,1)= -0.37499999999997785_rk
    ref_primal_x(38 ,1)= -0.37499999999997785_rk
    ref_primal_x(54 ,1)= -0.37499999999997785_rk
    ref_primal_x(10 ,1)= -0.37499999999997785_rk
    ref_primal_x(26 ,1)= -0.37499999999997785_rk
    ref_primal_x(42 ,1)= -0.37499999999997785_rk
    ref_primal_x(58 ,1)= -0.37499999999997785_rk
    ref_primal_x(14 ,1)= -0.37499999999997785_rk
    ref_primal_x(30 ,1)= -0.37499999999997785_rk
    ref_primal_x(46 ,1)= -0.37499999999997785_rk
    ref_primal_x(62 ,1)= -0.37499999999997785_rk
    ref_primal_x( 3 ,1)= -0.89062499999995426_rk
    ref_primal_x(19 ,1)= -0.89062499999995426_rk
    ref_primal_x(35 ,1)= -0.89062499999995426_rk
    ref_primal_x(51 ,1)= -0.89062499999995426_rk
    ref_primal_x( 7 ,1)= -0.89062499999995426_rk
    ref_primal_x(23 ,1)= -0.89062499999995426_rk
    ref_primal_x(39 ,1)= -0.89062499999995426_rk
    ref_primal_x(55 ,1)= -0.89062499999995426_rk
    ref_primal_x(11 ,1)= -0.89062499999995426_rk
    ref_primal_x(27 ,1)= -0.89062499999995426_rk
    ref_primal_x(43 ,1)= -0.89062499999995426_rk
    ref_primal_x(59 ,1)= -0.89062499999995426_rk
    ref_primal_x(15 ,1)= -0.89062499999995426_rk
    ref_primal_x(31 ,1)= -0.89062499999995426_rk
    ref_primal_x(47 ,1)= -0.89062499999995426_rk
    ref_primal_x(63 ,1)= -0.89062499999995426_rk
    ref_primal_x( 4 ,1)=  0.74999999999995148_rk
    ref_primal_x(20 ,1)=  0.74999999999995148_rk
    ref_primal_x(36 ,1)=  0.74999999999995148_rk
    ref_primal_x(52 ,1)=  0.74999999999995148_rk
    ref_primal_x( 8 ,1)=  0.74999999999995148_rk
    ref_primal_x(24 ,1)=  0.74999999999995148_rk
    ref_primal_x(40 ,1)=  0.74999999999995148_rk
    ref_primal_x(56 ,1)=  0.74999999999995148_rk
    ref_primal_x(12 ,1)=  0.74999999999995148_rk
    ref_primal_x(28 ,1)=  0.74999999999995148_rk
    ref_primal_x(44 ,1)=  0.74999999999995148_rk
    ref_primal_x(60 ,1)=  0.74999999999995148_rk
    ref_primal_x(16 ,1)=  0.74999999999995148_rk
    ref_primal_x(32 ,1)=  0.74999999999995148_rk
    ref_primal_x(48 ,1)=  0.74999999999995148_rk
    ref_primal_x(64 ,1)=  0.74999999999995148_rk


    ! Find the maximal deviation from the reference
    res = max(maxval(abs(res_x_y)), maxval(abs(res_x_z)))
    res = max(res, maxval(abs(primal_x - ref_primal_x)))

  end subroutine


end program atl_covolumeToPrimal_test
