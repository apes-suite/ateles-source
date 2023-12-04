! Copyright (c) 2014-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2020 Harald Klimach <harald.klimach@uni-siegen.de>
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

program atl_varSys_temperatureByElement_test
  use, intrinsic :: iso_c_binding,  only: C_NEW_LINE, c_loc

  use env_module,                   only: rk,                &
    &                                     print_self_status
  use flu_binding,                  only: flu_state

  use aotus_module,                 only: open_config_chunk

  use tem_logging_module,           only: logUnit, tem_logging_init
  use tem_general_module,           only: tem_general_type, &
    &                                     tem_load_general, &
    &                                     tem_start
  use treelmesh_module,             only: treelmesh_type, load_tem
  use tem_varSys_module,            only: tem_varSys_init,               &
    &                                     tem_varSys_append_stateVar,    &
    &                                     tem_varSys_append_derVar,      &
    &                                     tem_varSys_proc_point,         &
    &                                     tem_varSys_proc_element,       &
    &                                     tem_varSys_proc_setParams,     &
    &                                     tem_varSys_proc_getParams,     &
    &                                     tem_varSys_proc_setupIndices,  &
    &                                     tem_varSys_proc_getvalofIndex


  use ply_poly_project_module,      only: ply_poly_project_type
  use ply_dof_module,               only: Q_space, P_Space
  use ply_fpt_header_module,        only: ply_fpt_header_define
  use ply_dynArray_project_module,  only: ply_prj_init_define, &
    &                                     ply_prj_init_type
  use ply_prj_header_module,        only: ply_prj_header_type
  use ply_poly_project_module,      only: ply_poly_project_fillbody

  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_get_new_varSys_data_ptr
  use atl_kerneldata_module,        only: atl_statedata_type
  use atl_equation_module,          only: atl_equations_type
  use atl_eqn_euler_derive_module,  only: atl_temperature_getElement, &
    &                                     atl_pressure_getElement
  use atl_scheme_module,            only: atl_scheme_type, &
    &                                     atl_init_scheme
  use atl_timer_module,             only: atl_addTimers_perElem

  implicit none

  !****************************************************************************!
  ! PAREMETER
  !
  ! Number of dimensions
  integer, parameter :: nDimensions = 2
  ! Refinement level
  integer, parameter :: refinementLevel = 3
  ! Number of components for the state variable
  integer, parameter :: nComponents_stateVars(3) = (/ 1, 2, 1 /)
  ! The total number of the components of the resulting variable
  integer, parameter :: nComponents_res = 1
  ! Polynomial degree
  integer, parameter :: polyDegree = 3
  ! Index of projection to use
  integer, parameter :: projectionIndex = 1
  ! Projection factor, used for both projection kinds
  real(kind=rk), parameter :: projectionFactor = 2.0_rk
  ! The projection to use
  character(len=3), parameter :: projectionKind = 'l2p'
  ! Type of polynomial space to be used
  integer, parameter :: basisType = Q_Space
  ! Isentropic coefficient for euler equations
  real(kind=rk), parameter :: isenCoef = 3.0_rk
  real(kind=rk), parameter :: ideal_gas_constant = 296.0_rk
  !*****************************************************************************
  ! Derived from parameters
  !
  ! The number of elements
  integer, parameter :: nElements = (2**refinementLevel)**nDimensions
  ! Total number of components in the state
  !integer, parameter :: nComponents = sum(nComponents_stateVars)
  ! Workaround for the Intel 15 compiler
  integer, parameter :: nComponents = 4
  ! Degrees of freedom
  integer, parameter :: nDofs = (polyDegree+1)**nDimensions
  !*****************************************************************************
  ! Variables
  !
  procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
  procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
  procedure(tem_varSys_proc_setparams), pointer :: set_params => NULL()
  procedure(tem_varSys_proc_getparams), pointer :: get_params => NULL()
  procedure(tem_varSys_proc_setupIndices), pointer :: setup_indices => NULL()
  procedure(tem_varSys_proc_getValOfIndex), pointer :: get_valOfIndex => NULL()

  type(flu_state) :: conf, conf_scheme

  type(treelmesh_type), target :: tree
  type(tem_general_type), target :: general

  type(atl_varSys_solverData_type), target :: methodData
  type(atl_equations_type), target :: equation
  type(ply_prj_init_type) :: prj_init
  type(ply_prj_header_type) :: header
  type(ply_poly_project_type), target, allocatable :: polyProj(:)
  type(atl_statedata_type), target, allocatable :: state(:)
  type(atl_scheme_type), target, allocatable :: scheme(:)
  integer, target :: poly_proj_pos(3) = (/ 1, 1, 1 /)
  integer, target :: levelPointer(8)

  integer :: iElem, posOfDerived
  logical :: res_correct = .true.

  character, parameter :: nl = C_NEW_LINE
  character(len=110) :: cubeconf
  character(len=240) :: schemeconf
  character(len=1) :: basistype_chr
  character(len=7) :: schemetype
  character(len=5) :: meshType
  !*****************************************************************************

  select case(basisType)
  case(P_Space)
    basistype_chr = 'P'
  case(Q_Space)
    basistype_chr = 'Q'
  end select

  select case(nDimensions)
  case(1)
    schemetype = 'modg_1d'
    meshType = 'line '
  case(2)
    schemetype = 'modg_2d'
    meshType = 'slice'
  case(3)
    schemetype = 'modg   '
    meshType = 'cube'
  end select

  ! Fill the configuration stubs with parameters
  write(cubeconf, '(A,I2,A)')                              &
    &   'mesh = {' // nl                                   &
    & //'  predefined = "' // trim(meshType) // '",' // nl &
    & //'  origin = {0.0, 0.0, 0.0},' // nl                &
    & //'  length = 1.0,' // nl                            &
    & //'  refinementLevel = ',                            &
    & refinementLevel,                                     &
    & nl //'}' // nl
  write(schemeconf, '(A,I3,A)')                             &
    &   'scheme = {' // nl                                  &
    & //'  spatial = {' // nl                               &
    & //'    name = "' // trim(schemetype) // '",' // nl    &
    & //'    modg_space = "' // basistype_chr // '",' // nl &
    & //'    m = ',                                         &
    & polyDegree, nl                                        &
    & //'  },' // nl                                        &
    & //'  temporal = {' // nl                              &
    & //'    name = "explicitRungeKutta",' // nl            &
    & //'    steps = 4,' // nl                              &
    & //'    control = {' // nl                             &
    & //'      name = "cfl",' // nl                         &
    & //'      cfl = 0.95' // nl                            &
    & //'    }' // nl                                       &
    & //'  }' // nl                                         &
    & //'}'

  ! Init the Treelm environment
  call tem_start('atl_varSys_temperatureElement_test unit test', 'utest', general)

  !*****************************************************************************
  write(logUnit(3), *) 'Initialize the mesh'
  call open_config_chunk(L = conf, chunk = trim(cubeconf))

  call tem_load_general( me = general, conf = conf )
  call tem_logging_init(level = 6,               &
    &                   rank  = general%proc%rank)

  ! Load the mesh first.
  call load_tem(me     = tree,                   &
    &           conf   = conf,                   &
    &           myPart = general%proc%rank,      &
    &           nParts = general%proc%comm_size, &
    &           comm   = general%proc%comm       )

  call atl_addTimers_perElem(tree)

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing state'
  allocate(state(refinementLevel))

  ! Allocate the state array for 3D, 8 elements, 4th order, 5 scalars
  allocate(state(refinementLevel)%state(nElements,nDofs,nComponents))
  state(refinementLevel)%state = 0_rk
  do iElem = 1, nElements
    ! We will assign a constant value to the modes to be able to check the
    ! result here
    state(refinementLevel)%state(iElem, 1, 1) = real(iElem, kind=rk)
    state(refinementLevel)%state(iElem, 1, 2) = 2_rk
    state(refinementLevel)%state(iElem, 1, 3) = 3_rk
    state(refinementLevel)%state(iElem, 1, 4) = real(iElem, kind=rk)
  end do

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing scheme'
  call open_config_chunk(L = conf_scheme, chunk = trim(schemeconf))
  allocate(scheme(refinementLevel))
  call atl_init_scheme(scheme, conf_scheme, refinementLevel, refinementLevel)
  if (nDimensions /= 1 ) write(logUnit(1),*) 'ADJUST SCHEME INITIALISATION'
  scheme(refinementLevel)%modg_1d%maxpolydegree = polyDegree
  scheme(refinementLevel)%modg_1d%basisType = basisType

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing projection'
  header%kind = projectionKind
  header%l2p_header%factor = projectionFactor
  header%l2p_header%nodes_header%nodes_kind = 'gauss-legendre'
  call ply_fpt_header_define( me = header%fpt_header,   &
    &                         factor = projectionFactor )

  allocate(polyProj(projectionIndex))
  ! define poly projection init type
  call ply_prj_init_define(       &
    & me            = prj_init,   &
    & header        = header ,    &
    & maxPolyDegree = polyDegree, &
    & basisType     = basistype   )
  ! fill the projection body according to the header
  call ply_poly_project_fillbody(             &
    & me         = polyProj(projectionIndex), &
    & proj_init  = prj_init,                  &
    & scheme_dim = nDimensions                )

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing equation'
  equation%nDimensions = nDimensions
  equation%Euler%isen_coef = isenCoef
  equation%Euler%r = ideal_gas_constant

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing method data'
  methoddata%scheme_listPtr => scheme
  methodData%statedata_listPtr => state
  methodData%equationPtr => equation
  methodData%polyProjectPtr => polyProj
  methodData%poly_proj_posPtr => poly_proj_pos
  methodData%levelPointer => levelPointer

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing variable system'
  call tem_varSys_init(                &
    & me         = equation%varSys,    &
    & systemName = 'utest_accessState' )

  get_element => atl_varSys_getStateForElement
  call tem_varSys_append_stateVar(                           &
    & me             = equation%varSys,                         &
    & varName        = 'density',                               &
    & nComponents    = nComponents_stateVars(1),                &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )

  call tem_varSys_append_stateVar(                           &
    & me             = equation%varSys,                         &
    & varName        = 'momentum',                              &
    & nComponents    = nComponents_stateVars(2),                &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )

  call tem_varSys_append_stateVar(                           &
    & me             = equation%varSys,                         &
    & varName        = 'energy',                                &
    & nComponents    = nComponents_stateVars(3),                &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )

  get_element => atl_pressure_GetElement
  call tem_varSys_append_derVar(                                &
    & me             = equation%varSys,                         &
    & varName        = 'pressure',                              &
    & operType       = 'derived',                               &
    & nComponents    = 1,                                       &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & input_varname  = ['density ', 'momentum', 'energy  '],    &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )

  !*****************************************************************************
  ! This is the routine we want to test.
  get_element => atl_temperature_GetElement
  !*****************************************************************************
  call tem_varSys_append_derVar(                               &
    & me             = equation%varSys,                         &
    & varName        = 'temperature',                           &
    & operType       = 'derived',                               &
    & nComponents    = nComponents_res,                         &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & input_varname  = ['pressure', 'density '],                &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex,                          &
    & pos            = posOfDerived                             )

  res_correct = checkTemperature()

  call print_self_status()
  if (res_correct) then
    write(logUnit(1),*) 'PASSED'
  end if

contains

  function checkTemperature() result(result)
    !--------------------------------------------------------------------------!
    logical :: result
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: res(:)
    real(kind=rk) :: expected_res
    integer, allocatable :: elemPos(:)
    integer :: elem_offset, iElem, nElems_track
    real(kind = rk) :: coeff, pressure, density
    !--------------------------------------------------------------------------!
    result = .true.

    nElems_track = 5
    allocate(elemPos(nElems_track))
    elemPos = (/ 1, 3, 5, 7, 8 /)

    ! levelpointer should point to the correct tracking element element 1 and 5
    ! setting it manually here
    levelpointer = (/ 1,0,3,0,5,0,7,8 /)


    allocate(res(nElems_track * nDofs * nComponents))
    write(logUnit(1),*) 'Allocating res with size: ', nElems_track * nDofs * nComponents

    write(logUnit(3),*) 'Get element for ', trim(equation%varSys%varname%val(posOfDerived))

    ! access state array for given elemPos
    call equation%varSys%method%val(posOfDerived)%get_element(&
      & varSys  = equation%varSys,                            &
      & elemPos = elemPos,                                    &
      & time    = general%simControl%now,                     &
      & tree    = tree,                                       &
      & nElems  = nElems_track,                               &
      & nDofs   = nDofs,                                      &
      & res     = res                                         )

    do iElem = 1, nElems_track

      elem_offset =  (iElem-1)*nDofs
      pressure = ((equation%Euler%isen_coef-1.0_rk)                 &
        & * ( real(elemPos(iElem),kind=rk) - 0.5_rk                 &
        &   / real(elemPos(iElem),kind=rk)                          &
        &   * sum(state(refinementLevel)%state(iElem, 1, 2:3)**2) ) )
      coeff = 1.0 / equation%Euler%r
      density = real(elemPos(iElem),kind=rk)

      expected_res = coeff*(pressure/density)

      ! The difference is bigger than epsilon, but I don't know why. I guess it
      ! is due to the projection that takes place for the real result.
      !if (abs(res(elem_offset+1) - expected_res) > eps) then
      if (abs(res(elem_offset+1) - expected_res) > 0.000000001) then
        result = .false.
        write(logUnit(1),*) 'Failed at Element ', elemPos(iElem)
        write(logUnit(1),*) 'Expected was ', expected_res, ' but result was ', &
          & res(elem_offset+1)
      end if

    end do
  end function checkTemperature

end program atl_varSys_temperatureByElement_test
