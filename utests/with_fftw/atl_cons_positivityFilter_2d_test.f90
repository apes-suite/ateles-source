! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014, 2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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

!> Unit test to check functionallity of the conservative 2d positivity limiter
!! \author{Jens Zudrop}
program atl_cons_positivityFilter_2d_test
  use env_module,                         only: rk, fin_env

  use tem_logging_module,                 only: logUnit
  use tem_aux_module,                     only: tem_abort
  use tem_general_module,                 only: tem_start
  use tem_element_module,                 only: eT_fluid

  use ply_dof_module,                     only: Q_space
  use ply_dynArray_project_module,        only: ply_prj_init_define, &
    &                                           ply_prj_init_type
  use ply_poly_project_module,            only: ply_poly_project_fillbody, &
    &                                           ply_poly_project_m2n,      &
    &                                           ply_poly_project_n2m,      &
    &                                           ply_poly_project_type
  use ply_oversample_module,              only: ply_convert2oversample,  &
    &                                           ply_convertFromoversample
  use ply_prj_header_module,              only: ply_prj_header_type
  use ply_fpt_header_module,              only: ply_fpt_header_define

  use atl_cube_elem_module,               only: atl_cube_elem_type
  use atl_kerneldata_module,              only: atl_statedata_type
  use atl_cons_positivity_preserv_module, only: atl_cons_positivity_preserv_type
  use atl_stabilize_module,               only: atl_cons_positivity_preserv_2d
  use atl_eqn_euler_module,               only: atl_euler_type
  use atl_modg_2d_scheme_module,          only: atl_modg_2d_scheme_type
  use atl_solver_param_module,            only: atl_solver_param_type

  implicit none

  real(kind=rk), parameter :: tol = 1.0e-9_rk
  real(kind=rk), parameter :: lbfact = 1.0_rk - epsilon(tol)
  integer :: iPower
  real(kind=rk) :: res, newRes, nonCons, newCons
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = params%general      )

  res = 0.0_rk
  nonCons = 0.0_rk
  write(logunit(1),*) 'Checking with a tolerance of', tol
  do iPower = 1,5
    call atl_check_pntPosLimiter(iPower, newRes, newCons)
    write(logunit(1),*) 'minimal pressure or density:', newres
    write(logunit(1),*) 'maximal diff in integral mean:', newCons
    if(newRes.gt.res) then
      res = newRes
    end if
    if(newCons.gt.nonCons) then
      nonCons = newCons
    end if
  end do

  ! The smallest value should be just above the tolerance,
  ! but always larger than the tolerance.
  ! Conservativity should be preserved.
  if ( (res < 10*tol) .and. (res >= tol*lbfact) &
    &  .and. (nonCons < 1.e-12)                 ) then
    write(logUnit(1),*) 'PASSED'
  else
    write(logUnit(1),*) 'Res: ', res
    write(logUnit(1),*) 'Non-Cons: ', nonCons
    write(logUnit(1),*) 'Failed ...'
  end if

  call fin_env()

contains

  subroutine atl_check_pntPosLimiter(power, res, nonCons)
    !--------------------------------------------------------------------------
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    real(kind=rk), intent(out) :: nonCons
    !--------------------------------------------------------------------------
    integer :: nDim, nVars
    type(ply_prj_init_type) :: prj_init
    integer :: maxPolyDegree, iVar, basisType, iPnt
    type(atl_statedata_type) :: state
    type(atl_cube_elem_type) :: mesh
    type(atl_cons_positivity_preserv_type) :: filter
    type(atl_modg_2d_scheme_type) :: modg
    type(atl_Euler_type) :: euler
    integer :: nSeeds
    integer, allocatable :: rand_seed(:)
    real(kind=rk) :: rfac
    real(kind=rk) :: press
    real(kind=rk), allocatable :: pointVal(:,:), pressure(:), tmpState(:,:)
    real(kind=rk) :: minDensity, minPressure
    real(kind=rk) :: meanVal(4)
    type(ply_poly_project_type) :: poly_proj
    type(ply_prj_header_type) :: poly_proj_header
    !--------------------------------------------------------------------------
    maxPolyDegree =  2**power-1
    nVars = 4
    nDim = 2

    ! setup the filter
    filter%eps = 1.e-09_rk

    ! setup the equation info
    euler%isen_coef = 1.4_rk

    ! set up the projection scheme
    basisType = Q_space
    poly_proj_header%kind = 'fpt'
    call ply_fpt_header_define( me = poly_proj_header%fpt_header, &
      &                         lobattoPoints = .true.            )

    ! define my poly projection type
    call ply_prj_init_define( me            = prj_init,         &
      &                       header        = poly_proj_header, &
      &                       maxPolyDegree = maxPolydegree,    &
      &                       basisType     = basistype         )
    ! fill the projection body according to the header
    call ply_poly_project_fillbody( me         = poly_proj, &
      &                             proj_init  = prj_init,  &
      &                             scheme_dim = nDim       )
    ! setup the info for modg
    modg%maxPolyDegree = maxPolyDegree

    ! Setup the mesh info
    mesh%descriptor%nElems = 1
    mesh%descriptor%elem%nElems(eT_fluid) = 1

    ! transform to point values
    allocate( pressure( poly_proj%body_2d%nquadPoints ) )
    allocate( tmpState( poly_proj%body_2d%nquadPoints, nVars) )
    allocate( pointVal( poly_proj%body_2d%nquadPoints, nVars) )

    ! setup the state
    allocate( state%state(1, (maxPolyDegree+1)**nDim, nVars) )
    call random_seed(size=nSeeds)
    allocate(rand_seed(nSeeds))
    rand_seed = 0
    rand_seed(1) = 8345
    call random_seed(put=rand_seed)

    ! Ensure three points which violate the tolerance:
    ! Density < eps
    pointval(1,1) = filter%eps*0.999_rk !/(maxPolyDegree+1)
    do iVar = 2, nVars-1
      call random_number(rfac)
      pointval(1, iVar) = (rfac-0.5_rk)*pointval(1,1)
    end do
    press = filter%eps*1.01_rk
    pointval(1,4) = press/(euler%isen_coef-1._rk) &
      &           + 0.5_rk/pointval(1,1) &
      &             * (pointval(1,2)**2 + pointval(1,3)**2)

    ! Pressure < eps
    call random_number(rfac)
    pointval(2, 1) = filter%eps+0.01_rk*rfac
    do iVar = 2, nVars-1
      call random_number(rfac)
      pointval(2, iVar) = pointVal(1,iVar)*(1.0_rk + (rfac-0.5_rk)*0.01_rk)
    end do
    press = filter%eps*0.999_rk
    pointval(2,4) = press/(euler%isen_coef-1._rk) &
      &           + 0.5_rk/pointval(2,1) &
      &             * (pointval(2,2)**2 + pointval(2,3)**2)

    ! Density and pressure < eps
    pointval(3,1) = filter%eps*0.999_rk !/(maxPolyDegree+1)
    do iVar = 2, nVars-1
      call random_number(rfac)
      pointval(3, iVar) = pointVal(1,iVar)*(1.0_rk + (rfac-0.5_rk)*0.01_rk)
    end do
    press = filter%eps*0.999_rk
    pointval(3,4) = press/(euler%isen_coef-1._rk) &
      &           + 0.5_rk/pointval(3,1) &
      &             * (pointval(3,2)**2 + pointval(3,3)**2)


    ! Choose remaining points with random fluctuations of 1% for density and
    ! pressure and 10% for velocity
    ! Some may also violate the tolerance...
    do iPnt = 4, poly_proj%body_2d%nquadPoints
      ! Density larger than eps:
      call random_number(rfac)
      pointval(iPnt, 1) = filter%eps*(1.0_rk+0.01_rk*(rfac-0.1_rk))
      do iVar = 2, nVars-1
        call random_number(rfac)
        pointval(iPnt, iVar) = pointVal(1,iVar)*(1.0_rk + (rfac-0.5_rk)*0.1_rk)
      end do
      call random_number(rfac)
      press = filter%eps*(1.0_rk+0.01_rk*(rfac-0.1_rk))
      pointval(iPnt,4) = press/(euler%isen_coef-1._rk) &
        &              + 0.5_rk/pointval(iPnt,1)       &
        &                * ( pointval(iPnt,2)**2       &
        &                    + pointval(iPnt,3)**2     )
    end do

    ! transform nodal values to modal state
    call ply_poly_project_n2m(me         = poly_proj, &
      &                       dim        = nDim,      &
      &                       nVars      = nVars,     &
      &                       nodal_data = pointval,  &
      &                       modal_data = tmpState   )
    meanval = tmpState(1,:)
    call ply_convertFromOversample(modalCoeffs = tmpState,          &
      &                            poly_proj   = poly_proj,         &
      &                            nDim        = nDim,              &
      &                            state       = state%state(1,:,:) )

    ! run the limiter
    call atl_cons_positivity_preserv_2d( state, mesh, filter, &
      &                                  euler, poly_proj )

    ! convert to nodal
    call ply_convert2oversample(state       = state%state(1,:,:), &
      &                         poly_proj   = poly_proj,          &
      &                         nDim        = 2,                  &
      &                         modalCoeffs = tmpState            )
    call ply_poly_project_m2n(me         = poly_proj, &
      &                       dim        = 2,         &
      &                       nVars      = 4,         &
      &                       nodal_data =pointval,   &
      &                       modal_data =tmpState    )

    ! calculate the pressure
    pressure(:) = (euler%isen_coef-1.0_rk)*( pointVal(:,4) &
                &   - (0.5_rk/pointVal(:,1)) &
                &     *( pointVal(:,2)**2 +  pointVal(:,3)**2 ) )

    minDensity = minval( pointVal(:,1) )
    minPressure = minval( pressure(:) )

    ! Check the change in the mean values ...
    nonCons = maxval( abs( state%state(1,1,:) - meanVal ) )

    res = min( minDensity, minPressure )

  end subroutine

end program atl_cons_positivityFilter_2d_test
