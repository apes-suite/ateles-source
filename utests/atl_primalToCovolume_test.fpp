! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014-2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
program atl_primalToCovolume_test
  use env_module,               only: rk, fin_env
  use tem_general_module,       only: tem_start
  use tem_logging_module,       only: logUnit
  use atl_solver_param_module,  only: atl_solver_param_type
  use atl_stabilize_module,     only: atl_primal_to_covolume_projection
  use atl_covolume_module,      only: atl_covolume_type
  use atl_scheme_module,        only: atl_scheme_type, atl_modg_scheme_prp
  use atl_space_basis,          only: atl_init_spacebasis

  implicit none

  real(kind=rk) :: res, newRes
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codename = 'Ateles unit test', &
    &            general  = params%general      )

  ! Check the transformation for different number of polynomial coefficients.
  res = 0.0_rk
  write(logUnit(10),*) 'CHECKING PRIMAL->COVOLUME PROJECTION FOR '
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
    real(kind=rk), allocatable :: left(:,:), right(:,:), &
                                & covol_x(:,:), covol_y(:,:), covol_z(:,:), &
                                & res_x_y(:,:), res_x_z(:,:), ref_covol_x(:,:)
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
    allocate(ref_covol_x((maxPolyDegree+1)**3,1))
    allocate(covol_x((maxPolyDegree+1)**3,1))
    allocate(covol_y((maxPolyDegree+1)**3,1))
    allocate(covol_z((maxPolyDegree+1)**3,1))
    left(:,:) = 1.0_rk
    right(:,:) = 0.5_rk

    ! Set up the filter
    filter%alpha = 36.0_rk
    filter%order = 14.0_rk
    filter%beta = 1.0_rk

    ! Co-volume projection in x direction
    call atl_primal_to_covolume_projection( &
      & left       = left,                  &
      & right      = right,                 &
      & filter     = filter,                &
      & scheme     = scheme_list(level),    &
      & order      = filter%order,          &
      & maxPolyDeg = maxPolyDegree,         &
      & dir        = 1,                     &
      & covolume   = covol_x                )

    ! Co-volume projection in y direction
    call atl_primal_to_covolume_projection( &
      & left       = left,                  &
      & right      = right,                 &
      & filter     = filter,                &
      & scheme     = scheme_list(level),    &
      & order      = filter%order,          &
      & maxPolyDeg = maxPolyDegree,         &
      & dir        = 2,                     &
      & covolume   = covol_y                )

    ! Co-volume projection in z direction
    call atl_primal_to_covolume_projection( &
      & left       = left,                  &
      & right      = right,                 &
      & filter     = filter,                &
      & scheme     = scheme_list(level),    &
      & order      = filter%order,          &
      & maxPolyDeg = maxPolyDegree,         &
      & dir        = 3,                     &
      & covolume   = covol_z                )

    ! Check rotational invariance
    allocate(res_x_y((maxPolyDegree+1)**3,1))
    allocate(res_x_z((maxPolyDegree+1)**3,1))
    do deg1 = 1, maxPolyDegree+1
      do deg2 = 1, maxPolyDegree+1
        do degCo = 1, maxPolyDegree+1
?? copy :: posOfModgCoeffQTens(degCo, deg1, deg2, maxPolyDegree, pos_x_cov)
?? copy :: posOfModgCoeffQTens(deg1, degCo, deg2, maxPolyDegree, pos_y_cov)
?? copy :: posOfModgCoeffQTens(deg1, deg2, degCo, maxPolyDegree, pos_z_cov)
          res_x_y(pos_x_cov,:) = covol_x(pos_x_cov,:) - covol_y(pos_y_cov,:)
          res_x_z(pos_x_cov,:) = covol_x(pos_x_cov,:) - covol_z(pos_z_cov,:)
        end do
      end do
    end do

    !do degCo = 1, maxPolyDegree+1
    !  do deg1 = 1, maxPolyDegree+1
    !    do deg2 = 1, maxPolyDegree+1
!?? copy :: posOfModgCoeffQTens(degCo, deg1, deg2, maxPolyDegree, pos_x_cov)
    !      write(*,*) 'ref_covol_x(',pos_x_cov,',1)=', covol_x(pos_x_cov,1),'_rk'
    !    end do
    !  end do
    !end do

    ! Define the reference values
    ref_covol_x( 1 ,1)=  0.84374999999995415_rk
    ref_covol_x(17 ,1)=  0.84374364936592572_rk
    ref_covol_x(33 ,1)=  0.74586058258167542_rk
    ref_covol_x(49 ,1)=   1.9570973880179054E-016_rk
    ref_covol_x( 5 ,1)=  0.84374364936592572_rk
    ref_covol_x(21 ,1)=  0.84293750723320582_rk
    ref_covol_x(37 ,1)=  0.46863881756189496_rk
    ref_covol_x(53 ,1)=   1.7303991520946168E-033_rk
    ref_covol_x( 9 ,1)=  0.74586058258167542_rk
    ref_covol_x(25 ,1)=  0.46863881756189496_rk
    ref_covol_x(41 ,1)=   1.1776850828459816E-007_rk
    ref_covol_x(57 ,1)=   6.5078556689605132E-206_rk
    ref_covol_x(13 ,1)=   1.9570973880179054E-016_rk
    ref_covol_x(29 ,1)=   1.7303991520946168E-033_rk
    ref_covol_x(45 ,1)=   6.5078556689605132E-206_rk
    ref_covol_x(61 ,1)=   0.0000000000000000_rk
    ref_covol_x( 2 ,1)= -0.37499717749596523_rk
    ref_covol_x(18 ,1)= -0.37463889210364526_rk
    ref_covol_x(34 ,1)= -0.20828391891639678_rk
    ref_covol_x(50 ,1)=  -7.6906628981982597E-034_rk
    ref_covol_x( 6 ,1)= -0.37463889210364526_rk
    ref_covol_x(22 ,1)= -0.36887768813377553_rk
    ref_covol_x(38 ,1)=  -4.5601022576706930E-002_rk
    ref_covol_x(54 ,1)=  -7.4857475991950530E-065_rk
    ref_covol_x(10 ,1)= -0.20828391891639678_rk
    ref_covol_x(26 ,1)=  -4.5601022576706930E-002_rk
    ref_covol_x(42 ,1)=  -8.6982106134128721E-017_rk
    ref_covol_x(58 ,1)=  -0.0000000000000000_rk
    ref_covol_x(14 ,1)=  -7.6906628981982597E-034_rk
    ref_covol_x(30 ,1)=  -7.4857475991950530E-065_rk
    ref_covol_x(46 ,1)=  -0.0000000000000000_rk
    ref_covol_x(62 ,1)=  -0.0000000000000000_rk
    ref_covol_x( 3 ,1)= -0.78729728161399304_rk
    ref_covol_x(19 ,1)= -0.49467430742644619_rk
    ref_covol_x(35 ,1)=  -1.2431120318929845E-007_rk
    ref_covol_x(51 ,1)=  -6.8694032061250073E-206_rk
    ref_covol_x( 7 ,1)= -0.49467430742644619_rk
    ref_covol_x(23 ,1)= -0.10830242861967980_rk
    ref_covol_x(39 ,1)=  -2.0658250206855732E-016_rk
    ref_covol_x(55 ,1)=  -0.0000000000000000_rk
    ref_covol_x(11 ,1)=  -1.2431120318929845E-007_rk
    ref_covol_x(27 ,1)=  -2.0658250206855732E-016_rk
    ref_covol_x(43 ,1)=  -6.6441442759770857E-118_rk
    ref_covol_x(59 ,1)=  -0.0000000000000000_rk
    ref_covol_x(15 ,1)=  -6.8694032061250073E-206_rk
    ref_covol_x(31 ,1)=  -0.0000000000000000_rk
    ref_covol_x(47 ,1)=  -0.0000000000000000_rk
    ref_covol_x(63 ,1)=  -0.0000000000000000_rk
    ref_covol_x( 4 ,1)=   1.7396421226825648E-016_rk
    ref_covol_x(20 ,1)=   1.5381325796396434E-033_rk
    ref_covol_x(36 ,1)=   5.7847605946315072E-206_rk
    ref_covol_x(52 ,1)=   0.0000000000000000_rk
    ref_covol_x( 8 ,1)=   1.5381325796396434E-033_rk
    ref_covol_x(24 ,1)=   1.4971495198390022E-064_rk
    ref_covol_x(40 ,1)=   0.0000000000000000_rk
    ref_covol_x(56 ,1)=   0.0000000000000000_rk
    ref_covol_x(12 ,1)=   5.7847605946315072E-206_rk
    ref_covol_x(28 ,1)=   0.0000000000000000_rk
    ref_covol_x(44 ,1)=   0.0000000000000000_rk
    ref_covol_x(60 ,1)=   0.0000000000000000_rk
    ref_covol_x(16 ,1)=   0.0000000000000000_rk
    ref_covol_x(32 ,1)=   0.0000000000000000_rk
    ref_covol_x(48 ,1)=   0.0000000000000000_rk
    ref_covol_x(64 ,1)=   0.0000000000000000_rk

    ! Find the maximal deviation from the reference
    res = max(maxval(abs(res_x_y)), maxval(abs(res_x_z)))
    res = max(res, maxval(abs(covol_x - ref_covol_x)))

  end subroutine


end program atl_primalToCovolume_test
