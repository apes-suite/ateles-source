! Copyright (c) 2013-2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2015,2019-2020 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
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

!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program atl_pnt_positivityFilter_2d_test
  use env_module,                     only: rk, fin_env
  use tem_logging_module,             only: logUnit
  use tem_aux_module,                 only: tem_abort
  use tem_general_module,             only: tem_start
  use tem_element_module,             only: eT_fluid
  use atl_solver_param_module,        only: atl_solver_param_type
  use atl_cube_elem_module,           only: atl_cube_elem_type
  use atl_kerneldata_module,          only: atl_statedata_type
  use atl_positivity_preserv_module,  only: atl_positivity_preserv_type
  use atl_stabilize_module,           only: atl_positivity_preserv_2d
  use atl_eqn_euler_module,           only: atl_euler_type
  use atl_modg_2d_scheme_module,      only: atl_modg_2d_scheme_type
  use ply_dof_module,                 only: Q_space
  use ply_dynArray_project_module,    only: ply_prj_init_define, &
    &                                       ply_prj_init_type
  use ply_poly_project_module,        only: ply_poly_project_fillbody, &
    &                                       ply_poly_project_m2n,      &
    &                                       ply_poly_project_n2m,      &
    &                                       ply_poly_project_type
  use ply_oversample_module,          only: ply_convert2oversample
  use ply_prj_header_module,          only: ply_prj_header_type
  use ply_fpt_header_module,          only: ply_fpt_header_define

  implicit none

  integer :: iPower
  real(kind=rk) :: res, newRes
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = params%general      )

  res = 0.0_rk
  do iPower = 1,4
    call atl_check_pntPosLimiter(iPower, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if
  call fin_env()

contains

  subroutine atl_check_pntPosLimiter(power, res)
    !--------------------------------------------------------------------------
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    !--------------------------------------------------------------------------
    type(ply_prj_init_type) :: prj_init
    integer :: maxPolyDegree, iVar, iPoly, basisType
    type(atl_statedata_type) :: state
    type(atl_cube_elem_type) :: mesh
    type(atl_positivity_preserv_type) :: filter
    type(atl_modg_2d_scheme_type) :: modg
!    type(atl_legFpt_type) :: fpt
    type(atl_Euler_type) :: euler
    integer :: nSeeds
    integer, allocatable :: rand_seed(:)
    real(kind=rk) :: rfac
    real(kind=rk), allocatable :: pointVal(:,:), pressure(:)
    real(kind=rk), allocatable :: modalVal(:,:)
    real(kind=rk) :: minDensity, minEnergy, minPressure
    type(ply_poly_project_type) :: poly_proj
    type(ply_prj_header_type) :: poly_proj_header
    !--------------------------------------------------------------------------
    maxPolyDegree =  2**power-1

    ! Setup the mesh info
    mesh%descriptor%nElems = 1
    mesh%descriptor%elem%nElems(eT_fluid) = 1

    ! setup the state
    allocate( state%state(1, (maxPolyDegree+1)**2, 4) )
    call random_seed(size=nSeeds)
    allocate(rand_seed(nSeeds))
    rand_seed = 0
    rand_seed(1) = 8345
    call random_seed(put=rand_seed)
    do iVar = 1, 4
      do iPoly = 1, (maxPolyDegree+1)**2
        call random_number(rfac)
        state%state(1,iPoly, iVar) = rfac
      end do
    end do

    ! setup the filter
    filter%eps = 1.e-09_rk

    ! setup the equation info
    euler%isen_coef = 1.4_rk

    ! set up the projection scheme
    basisType = Q_space
    poly_proj_header%kind = 'fpt'
    call ply_fpt_header_define( me = poly_proj_header%fpt_header, &
      &                         lobattoPoints = .true.,           &
      &                         approx_terms = 18,                &
      &                         striplen = 512                    )

    ! define my poly projection type
    call ply_prj_init_define( me = prj_init,                    &
      &                       header        = poly_proj_header, &
      &                       maxPolyDegree = maxpolydegree,    &
      &                       basisType     = basistype         )
    ! fill the projection body according to the header
    call ply_poly_project_fillbody(me = poly_proj, proj_init=prj_init, scheme_dim=2)
    ! setup the info for modg
    modg%maxPolyDegree = maxPolyDegree
!!VK    call atl_init_legFpt_2D( maxPolyDegree = maxPolyDegree, fpt = fpt, &
!!VK                           & nquadpointsperDir= maxPolydegree+1, &
!!VK                           & nVars = 4, lobattoPoints = .true. )
!!VK
    ! run the limiter
    call atl_positivity_preserv_2d(state, mesh, filter, euler, poly_proj)
    allocate( modalVal( poly_proj%body_2d%oversamp_dofs, 4) )
    call ply_convert2oversample( state       = state%state(1,:,:), &
      &                          poly_proj   = poly_proj,          &
      &                          nDim        = 2,                  &
      &                          modalCoeffs = modalval,           &
      &                          nScalars    = 4                   )

    ! transform to point values
   !!VK allocate( pointVal( (maxPolyDegree+1)**2, 4) )
   !!VK allocate( pressure( (maxPolyDegree+1)**2 ) )
    allocate( pointVal( poly_proj%body_2d%nquadpoints, 4) )
    allocate( pressure( poly_proj%body_2d%nquadpoints ) )
    !!VKcall atl_legToPnt_2D( fpt = fpt, legCoeffs = state%state(1,:,:), &
    !!VK                    & pntVal = pointVal, nVars = 4 )
    call ply_poly_project_m2n( me         = poly_proj, &
      &                        dim        = 2,         &
      &                        nVars      = 4,         &
      &                        nodal_data = pointVal,  &
      &                        modal_data = modalVal   )

    ! calculate the pressure
    pressure(:) = (euler%isen_coef - 1.0_rk) * ( pointVal(:,4) &
      & - (0.5_rk / pointVal(:,1)) * ( pointVal(:,2)**2        &
      & + pointVal(:,3)**2 ) )

    minDensity = minval( pointVal(:,1) )
    minEnergy = minval( pointVal(:,4) )
    minPressure = minval( pressure(:) )

    res = min( minDensity, minEnergy, minPressure)

  end subroutine

end program atl_pnt_positivityFilter_2d_test
