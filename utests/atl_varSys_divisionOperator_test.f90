! Copyright (c) 2014-2018, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2018,2020 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> This unit tests the division operator, but also uses it to test the generic
!! routine for modal quantities atl_generic_fromModal_getElement.
!! We test among other things whether the generic routine can handle
!! multi-component input variables as well as multi-component results.
program atl_varSys_divisionOperatorByElement_test
  use, intrinsic :: iso_c_binding,  only: C_NEW_LINE, c_loc, c_ptr

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
  use tem_variable_module,          only: tem_variable_type
  use tem_tools_module,             only: tem_horizontalSpacer
  use tem_operation_var_module,     only: tem_divideVecByScal_forPoint,  &
    &                                     tem_divideVecByScal_fromIndex, &
    &                                     tem_opVar_setupIndices,        &
    &                                     tem_get_new_varSys_data_ptr

  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_varSys_getStateForPoint,   &
    &                                     atl_get_new_varSys_data_ptr,   &
    &                                     atl_varSys_setupStateIndices,  &
    &                                     atl_varSys_getStateValOfIndex, &
    &                                     atl_init_varsys_solverdata
  use atl_kerneldata_module,        only: atl_statedata_type, &
    &                                     atl_kerneldata_type
  use atl_equation_module,          only: atl_equations_type
  use atl_operator_module,          only: atl_op_divideVecByScal_forElement, &
    &                                     atl_set_opVar_getElement
  use ply_poly_project_module,      only: ply_poly_project_type,    &
    &                                     ply_poly_project_fillbody
  use ply_dof_module,               only: P_Space, Q_Space
  use ply_fpt_header_module,        only: ply_fpt_header_define
  use ply_dynArray_project_module,  only: ply_prj_init_define, &
    &                                     ply_prj_init_type
  use ply_prj_header_module,        only: ply_prj_header_type
  use atl_scheme_module,            only: atl_scheme_type, &
    &                                     atl_init_scheme
  use atl_timer_module,             only: atl_addTimers_perElem
  use atl_cube_elem_module,         only: atl_cube_elem_type

  implicit none

  !****************************************************************************!
  ! PAREMETER
  !
  ! Number of dimensions
  integer, parameter :: nDimensions = 2
  ! Refinement level
  integer, parameter :: refinementLevel = 2
  ! Number of components for the state variables
  integer, parameter :: nComponents_stateVars(2) = (/ 2, 1 /)
  ! The total number of the components of the resulting variable
  integer, parameter :: nComponents_res = 2
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
  !*****************************************************************************
  ! Derived from parameters
  !
  ! The number of elements
  integer, parameter :: nElements = (2**refinementLevel)**nDimensions
  ! Total number of components in the state
  !integer, parameter :: nComponents = sum(nComponents_stateVars)
  ! Workaround for Intel 15 compiler
  integer, parameter :: nComponents = 3
  ! Degrees of freedom
  integer, parameter :: nDofs = (polyDegree+1)**nDimensions
  !*****************************************************************************
  ! Variables
  !
  ! Reference for each mesh level to the projection in polyProj that should be
  ! used for the particular level
  integer, target :: poly_proj_pos(nDimensions) = &
    & (/ projectionIndex, projectionIndex /)
  integer, target :: levelPointer(nElements)
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
  type(c_ptr) :: solver_bundle
  type(atl_equations_type), target :: equation
  type(ply_prj_init_type) :: prj_init
  type(ply_prj_header_type) :: header
  type(ply_poly_project_type), target, allocatable :: polyProj(:)
  type(atl_statedata_type), target, allocatable :: state(:)
  type(atl_kerneldata_type), target, allocatable :: kernel(:)
  type(atl_scheme_type), target, allocatable :: scheme(:)
  type(atl_cube_elem_type), target, allocatable :: mesh(:)
  type(atl_varSys_solverData_type), target :: varSys_data

  integer :: iElem, posOfDerived, maxdegree
  logical :: res_correct = .true.
  character, parameter :: nl = C_NEW_LINE

  character(len=110) :: cubeconf
  character(len=240) :: schemeconf
  character(len=1) :: basistype_chr
  character(7) :: schemetype
  character(len=5) :: meshType
  integer :: i
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
    & //'  length = 2.0,' // nl                            &
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
  call tem_start(                                        &
    & 'atl_varSys_divisionOperator_test unit test', &
    & 'utest',                                           &
    & general                                            )

  ! get the level pointer for single level run
  levelpointer = [ (i, i=1,nElements)]
  !*****************************************************************************
  write(logUnit(3), *) 'Initialize the mesh'
  call open_config_chunk(L = conf, chunk = trim(cubeconf))

  call tem_load_general(me = general, conf = conf)
  call tem_logging_init(level = 10, rank  = general%proc%rank)

  call load_tem(me     = tree,                   &
    &           conf   = conf,                   &
    &           myPart = general%proc%rank,      &
    &           nParts = general%proc%comm_size, &
    &           comm   = general%proc%comm       )

  call atl_addTimers_perElem(tree)

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing state'
  !! State for every level
  allocate(state(refinementLevel))
  allocate(kernel(refinementLevel))
  allocate(state(refinementLevel)%state(nElements,nDofs,nComponents))
  state(refinementLevel)%state = 0_rk
  do iElem = 1, nElements
    ! We will assign a constant value to the modes to be able to check the
    ! result here
    state(refinementLevel)%state(iElem, 1, 1) = real(iElem, kind=rk)
    state(refinementLevel)%state(iElem, 1, 2) = real(iElem, kind=rk) * 2_rk
    state(refinementLevel)%state(iElem, 1, 3) = 2_rk
    write(*,*)  state(refinementLevel)%state(iElem,1,:)
  end do

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing scheme'
  call open_config_chunk(L = conf_scheme, chunk = trim(schemeconf))
  allocate(scheme(refinementLevel))
  call atl_init_scheme(scheme, conf_scheme, refinementLevel, refinementLevel)
  if (nDimensions /= 2 ) write(logUnit(1),*) 'ADJUST SCHEME INITIALISATION'
  scheme(refinementLevel)%modg_2d%maxpolydegree = polyDegree
  scheme(refinementLevel)%modg_2d%basisType = basisType

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing projection'
  maxdegree = polyDegree
  header%kind = projectionKind
  header%l2p_header%factor = projectionFactor
  header%l2p_header%nodes_header%nodes_kind = 'gauss-legendre'
  call ply_fpt_header_define( me = header%fpt_header,   &
    &                         factor = projectionFactor )

  allocate(polyProj(projectionIndex))
  ! define poly projection init type
  call ply_prj_init_define(      &
    & me            = prj_init,  &
    & header        = header ,   &
    & maxPolyDegree = maxdegree, &
    & basisType     = basistype  )
  ! fill the projection body according to the header
  call ply_poly_project_fillbody(             &
    & me         = polyProj(projectionIndex), &
    & proj_init  = prj_init,                  &
    & scheme_dim = nDimensions                )

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing equation'
  equation%nDimensions = nDimensions

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing method data'
  methoddata%scheme_listPtr => scheme
  methodData%statedata_listPtr => state
  methodData%equationPtr => equation
  methodData%polyProjectPtr => polyProj
  methodData%poly_proj_posPtr => poly_proj_pos
  methodData%levelPointer => levelPointer

  !*****************************************************************************
  ! levelpointer should point to the correct tracking element element 1 and 5
  ! setting it manually here
  write(logUnit(3),*) 'Initializing Solver data'
  ! Initialize the global instance of atl_varSys_data
  call atl_init_varSys_solverdata(   &
    & this           = varSys_data,  &
    & equation       = equation,     &
    & tree           = tree,         &
    & schemeList     = scheme,       &
    & statedataList  = state,        &
    & kerneldataList = kernel,       &
    & meshList       = mesh,         &
    & polyProjList   = polyProj,     &
    & levelPointer   = levelPointer, &
    & polyProjPos    = poly_proj_pos )

   !*****************************************************************************

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing variable system'
  call tem_varSys_init(                &
    & me         = equation%varSys,    &
    & systemName = 'utest_accessState' )

  get_element => atl_varSys_getStateForElement
  get_point => atl_varSys_getStateForPoint
  setup_indices => atl_varSys_setupStateIndices
  get_valOfIndex => atl_varSys_getStateValOfIndex
  call tem_varSys_append_stateVar(                              &
    & me             = equation%varSys,                         &
    & varName        = 'momentum',                              &
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
    & varName        = 'density',                               &
    & nComponents    = nComponents_stateVars(2),                &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )

  !*****************************************************************************
  ! This is the routine we want to test.
  get_element => atl_op_divideVecByScal_forElement
  get_point => tem_divideVecByScal_forPoint
  setup_indices => tem_opVar_setupIndices
  get_valOfIndex => tem_divideVecByScal_fromIndex

  solver_bundle = atl_get_new_varSys_data_ptr(methodData)
  call tem_varSys_append_derVar(                                   &
    & me             = equation%varSys,                            &
    & varName        = 'velocity',                                 &
    & operType       = 'divide_vector_by_scalar',                  &
    & nComponents    = nComponents_res,                            &
    & method_data    = tem_get_new_varSys_data_ptr(solver_bundle), &
    & input_varname  = ['momentum', 'density '],                   &
    & get_point      = get_point,                                  &
    & get_element    = get_element,                                &
    & set_params     = set_params,                                 &
    & get_params     = get_params,                                 &
    & setup_indices  = setup_indices,                              &
    & get_valOfIndex = get_valOfIndex,                             &
    & pos            = posOfDerived                                )


  !*****************************************************************************
  write(logUnit(3),*) 'Doing the actual test'
  res_correct = checkElementwiseAccess() .and. checkPointwiseAccess() .and. &
    &           checkIndexwiseAccess()
  !*****************************************************************************

  call tem_horizontalSpacer(funit = logUnit(1))
  call print_self_status()
  if (res_correct) then
    write(logUnit(1),*) 'PASSED'
  end if

contains

  function checkPointwiseAccess() result(res)
    !--------------------------------------------------------------------------!
    logical :: res
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: derived(:)
    real(kind=rk) :: pointPos(1,3)
    real(kind=rk) :: expected_res
    !--------------------------------------------------------------------------!
    res = .true.
    pointPos(1,:) = (/ 0.1, 0.1, 0.0 /)

    !! The result is two components long
    allocate(derived(nComponents_res))

    call tem_horizontalSpacer(funit = logUnit(1))
    write(logUnit(3),*) 'Get point for ',             &
      & trim(equation%varSys%varname%val(posOfDerived))

    ! access state array for given elemPos
    call equation%varSys%method%val(posOfDerived)%get_point(&
      & varSys  = equation%varSys,                          &
      & point   = pointPos,                                 &
      & time    = general%simControl%now,                   &
      & tree    = tree,                                     &
      & nPnts   = 1,                                        &
      & res     = derived                                   )

    ! THe result should be 4,5, because the point we are tracking is in element
    ! 9 (the first 8 treeIDs are used for the first refinement level, as the
    ! treeIDs are always given in a 3D manner), and the second variable is
    ! always 2.
    expected_res = 0.5_rk

    if (abs(derived(1) - (expected_res)) > 1.0e-12) then
      res = .false.
      write(logUnit(1),*) 'Failed at point 1'
      write(logUnit(1),*) 'Expected was ', expected_res, ' &
        & but result was ', derived(1)
    end if

  end function checkPointwiseAccess

  function checkElementwiseAccess() result(res)
    !--------------------------------------------------------------------------!
    logical :: res
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: derived(:)
    integer :: iComp, iElem, nElems_track, elem_offset
    integer, allocatable :: elemPos(:)
    real(kind=rk) :: expected_res
    !--------------------------------------------------------------------------!
    res = .true.
    nElems_track = 5
    allocate(elemPos(nElems_track))
    elemPos = (/ 1, 3, 5, 7, 8 /)


    !! The result is two components long
    allocate(derived(nElems_track * nDofs * nComponents_res))

    call tem_horizontalSpacer(funit = logUnit(1))
    write(logUnit(3),*) 'Get element for ',             &
      & trim(equation%varSys%varname%val(posOfDerived))

    ! access state array for given elemPos
    call equation%varSys%method%val(posOfDerived)%get_element(&
      & varSys  = equation%varSys,                        &
      & elemPos = elemPos,                                &
      & time    = general%simControl%now,                 &
      & tree    = tree,                                   &
      & nElems  = nElems_track,                           &
      & nDofs   = nDofs,                                  &
      & res     = derived                                 )

    do iElem = 1, nElems_track

      elem_offset =  (iElem-1)*nDofs*nComponents_res

      expected_res = real(elemPos(iElem),kind=rk)/2_rk

      write(logUnit(4),*) 'iElem ', iElem,     &
        & 'elemPos ', elemPos(iElem),          &
        & ' res(1) ', derived( elem_offset+1), &
        & ' res(2) ', derived( elem_offset+2)

      ! The difference is bigger than epsilon, but I don't know why. I guess it
      ! is due to the projection that takes place for the real result.
      !if (abs(res(elem_offset+1) - expected_res) > eps) then
      do iComp = 1, nComponents_res
        if (abs(derived(elem_offset+iComp) - expected_res*iComp) > 1.0e-12 ) then
          res = .false.
          write(logUnit(1),*) 'Failed at Element ', elemPos(iElem), &
            & ' Component ', iComp
          write(logUnit(1),*) 'Expected was ', expected_res, ' &
            & but result was ', derived(elem_offset+1)
        end if
      end do

    end do

  end function checkElementwiseAccess


  function checkIndexwiseAccess() result(res)
    !--------------------------------------------------------------------------!
    logical :: res
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: derived(:,:)
    real(kind=rk) :: pointPos(5,3)
    real(kind=rk) :: expected_resY(5)
    real(kind=rk) :: expected_resX(5)
    integer, allocatable :: idx(:,:)
    integer :: iLevel, iPnt
    !--------------------------------------------------------------------------!
    call tem_horizontalSpacer(funit = logUnit(1))
    write(logUnit(3),*) "Check Access indexwise"

    res = .true.
    pointPos(1,:) = (/ 0.1, 0.1, 0.0 /)
    pointPos(2,:) = (/ 0.6, 0.6, 0.0 /)
    pointPos(3,:) = (/ 1.1, 1.1, 0.0 /)
    pointPos(4,:) = (/ 1.8, 1.8, 0.0 /)
    pointPos(5,:) = (/ 0.1, 0.1, 0.4 /)

    ! First we need to setup the index
    do iLevel = tree%global%minLevel, tree%global%maxLevel
      allocate(idx(iLevel,5) )

      ! setup indices for the derived variable
      call equation%varSys%method%val(posOfDerived)%setup_indices( &
        &  VarSys = equation%VarSys,   &
        &  point  = pointPos,          &
        &  iLevel = iLevel,            &
        &  tree   = tree,              &
        &  nPnts  = 5,                 &
        &  idx    = idx(iLevel,:)      )
      write(logUnit(3),*) 'Index for derived on level', iLevel, &
        &                 'are=', idx(iLevel,:)
    end do

    !! The result is levelwise and two components long
    allocate(derived(tree%global%minLevel:tree%global%maxLevel,nComponents_res*5))

    write(logUnit(3),*) 'Get value of index for derived variable ', &
      & trim(equation%varSys%varname%val(posOfDerived))

    do iLevel = tree%global%minLevel, tree%global%maxLevel
      ! get dervied variable for given index
      call equation%varSys%method%val(posOfDerived)%get_valOfIndex( &
        & varSys  = equation%varSys,                                &
        & time    = general%simControl%now,                         &
        & iLevel  = iLevel,                                         &
        & idx     = idx(iLevel,:),                                  &
        & nVals   = 5,                                              &
        & res     = derived(iLevel,:)                               )
      write(logUnit(3),*) 'Derived values for level', iLevel, &
        &                 'are=', derived(iLevel,:)

      ! The results are the statevalues at the first and second component,
      ! divided by 2 (density = 2),
      ! the first component is the elemet id
      ! the second coponent is the element id * 2
      ! the points are in elements 1,4,13,16,1
      expected_resX(1) = 0.5_rk
      expected_resX(2) = 2_rk
      expected_resX(3) = 6.5_rk
      expected_resX(4) = 8_rk
      expected_resX(5) = 0.5_rk

      expected_resY(1) = 1_rk
      expected_resY(2) = 4_rk
      expected_resY(3) = 13_rk
      expected_resY(4) = 16_rk
      expected_resY(5) = 1_rk

      do iPnt = 1,5
        ! check X component
        if ( abs(derived(iLevel,((iPnt-1)*ncomponents_res)+1) &
          &  - expected_resX(iPnt)) > 1.0e-12)  then
          res = .false.
          write(logUnit(1),*) 'Failed at point', iPnt
          write(logUnit(1),*) 'Expected was ', expected_resX(iPnt), ' &
            & but result was ', derived(iLevel,((iPnt-1)*ncomponents_res)+1)
        end if
        ! check Y component
        if ( abs(derived(iLevel,((iPnt-1)*ncomponents_res)+ncomponents_res)  &
          &  - expected_resY(iPnt)) > 1.0e-12)  then
          res = .false.
          write(logUnit(1),*) 'Failed at point', iPnt
          write(logUnit(1),*) 'Expected was ', expected_resY(iPnt), '  &
            & but result was ',                                        &
            & derived(iLevel,((iPnt-1)*ncomponents_res)+ncomponents_res)
        end if
      end do
    end do
  end function checkIndexwiseAccess

end program atl_varSys_divisionOperatorByElement_test

