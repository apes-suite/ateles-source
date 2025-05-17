! Copyright (c) 2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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

!> This unit tests the gradient operator
program atl_varSys_gradientOperatorByElement_test
  use, intrinsic :: iso_c_binding,  only: C_NEW_LINE, c_loc

  use env_module,                   only: rk,                &
    &                                     print_self_status, &
    &                                     eps
  use flu_binding,                  only: flu_state

  use aotus_module,                 only: open_config_chunk

  use tem_tools_module,             only: tem_horizontalSpacer
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
    &                                     tem_varSys_proc_getvalofIndex, &
    &                                     tem_varSys_solverData_evalElem_type
  use tem_variable_module,          only: tem_variable_type
  use tem_operation_var_module,     only: tem_varSys_append_operVar

  use ply_poly_project_module,      only: ply_poly_project_type,    &
   &                                      ply_poly_project_fillbody
  use ply_dof_module,               only: Q_space, P_Space
  use ply_fpt_header_module,        only: ply_fpt_header_define
  use ply_dynArray_project_module,  only: ply_prj_init_define, &
    &                                     ply_prj_init_type
  use ply_prj_header_module,        only: ply_prj_header_type
  use ply_modg_basis_module,        only: ply_evalLegendreTensPoly

  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_varSys_getStateForPoint,   &
    &                                     atl_get_new_varSys_data_ptr,   &
    &                                     atl_varSys_setupStateIndices,  &
    &                                     atl_varSys_getStateValOfIndex, &
    &                                     atl_set_stfun_getElement,      &
    &                                     atl_init_varsys_solverdata
  use atl_kerneldata_module,        only: atl_statedata_type, &
    &                                     atl_kerneldata_type
  use atl_equation_module,          only: atl_equations_type
  use atl_operator_module,          only: atl_op_Gradient_forElement, &
    &                                     atl_op_Gradient_forPoint,   &
    &                                     atl_op_Gradient_fromIndex,   &
    &                                     atl_opVar_setupIndices,      & 
    &                                     atl_set_opVar_getElement
  use atl_scheme_module,            only: atl_scheme_type, &
    &                                     atl_init_scheme
  use atl_modg_1d_basis_module,     only: atl_evalLegendreTensPoly1d
  use atl_modg_2d_basis_module,     only: atl_evalLegendreTensPoly2d
  use atl_timer_module,             only: atl_addTimers
  use atl_weights_module,           only: atl_initWeights

  implicit none

  !****************************************************************************!
  ! PAREMETER
  !
  ! Number of dimensions
  integer, parameter :: nDimensions = 3
  ! Refinement level
  integer, parameter :: refinementLevel = 1
  ! length for mesh
  integer, parameter :: meshLength = 4
  ! Number of components for the state variable
  integer, parameter :: nComponents_stateVars(1) = (/ 3 /)
  ! Polynomial degree
  integer, parameter :: polyDegree = 2
  ! Index of projection to use
  integer, parameter :: projectionIndex = 1
  ! Projection factor, used for both projection kinds
  real(kind=rk), parameter :: projectionFactor = 2.0_rk
  ! The projection to use
  character(len=3), parameter :: projectionKind = 'fpt'
  ! Type of polynomial space to be used
  integer, parameter :: basisType = Q_Space
  ! Isentropic coefficient for euler equations
  real(kind=rk), parameter :: isenCoef = 3.0_rk
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
  ! The total number of the components of the resulting variable
  integer, parameter :: nComponents_res = nDimensions*nComponents_stateVars(1)
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
  type(atl_kerneldata_type), target, allocatable :: kernel(:)
  type(atl_scheme_type), target, allocatable :: scheme(:)
  type(atl_cube_elem_type), target, allocatable :: mesh(:)
  integer, target :: poly_proj_pos(1) = (/ 1 /)
  integer, target :: levelPointer(5)

  integer :: iComp, iDof, posOfDerived
  logical :: res_correct = .true.

  character, parameter :: nl = C_NEW_LINE
  character(len=110) :: cubeconf
  character(len=240) :: schemeconf
  character(len=1) :: basistype_chr
  character(len=7) :: schemetype
  character(len=5) :: meshType

  type(tem_varSys_solverData_evalElem_type) :: solverData_evalElem
  type(atl_varSys_solverData_type), target :: varSys_data
  type(tem_variable_type) :: gradMom

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
    & //'  length = 4.0,' // nl                            &
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
  call tem_start(                         &
    & 'atl_varSys_gradientOperator_test', &
    & general                             )

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

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing state'
  allocate(state(refinementLevel))
  allocate(kernel(refinementLevel))
  ! Allocate the state array for 3D, 8 elements, 4th order, 5 scalars
  allocate(state(refinementLevel)%state(nElements,nDofs,nComponents))
  state(refinementLevel)%state = 0_rk
  do iDof = 1, nDofs
    ! We will assign a constant value to the modes to be able to check the
    ! result here
    state(refinementLevel)%state(:, iDof, 1) = real(iDof, kind=rk)
    state(refinementLevel)%state(:, iDof, 2) = real(iDof, kind=rk) * 2_rk
    state(refinementLevel)%state(:, iDof, 3) = real(iDof, kind=rk) * 3_rk
  end do

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing scheme'
  call open_config_chunk(L = conf_scheme, chunk = trim(schemeconf))
  allocate(scheme(refinementLevel))
  call atl_init_scheme(scheme, conf_scheme, refinementLevel, refinementLevel)
  if (nDimensions /= 3 ) write(logUnit(1),*) 'ADJUST SCHEME INITIALISATION'
  scheme(refinementLevel)%modg%maxpolydegree = polyDegree
  scheme(refinementLevel)%modg%basisType = basisType

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

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing method data'
  methoddata%scheme_listPtr => scheme
  methodData%statedata_listPtr => state
  methodData%equationPtr => equation
  methodData%polyProjectPtr => polyProj
  methodData%poly_proj_posPtr => poly_proj_pos
  methodData%levelPointer => levelPointer
  !*****************************************************************************
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
    & nComponents    = 3,                                       &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )

  solverData_evalElem%solver_bundle = c_loc(varSys_data)
  solverData_evalElem%stFun_setter => atl_set_stfun_getElement
  solverData_evalElem%opVar_setter => atl_set_opVar_getElement
  
  !create operation variable

  gradMom%label = 'grad_mom'
  gradMom%nComponents = 9
  gradMom%varType = 'operation'
  gradMom%operType = 'gradient'
  gradMom%input_varName = ['momentum'] 
 
  call tem_varSys_append_operVar( operVar             = gradMom ,           & 
    &                             varSys              = equation%varSys,    &
    &                             pos                 = posOfDerived ,      &
    &                             solverData_evalElem = solverData_evalElem )

  !*****************************************************************************
  write(logUnit(3),*) 'Add timer required in the operation routines'
  call atl_addTimers()
  call atl_initWeights( tree )

  !*****************************************************************************
  call tem_horizontalSpacer(fUnit = logUnit(3))
  write(logUnit(3),*) 'Doing the actual test'
!!  res_correct = checkGradientOperatorByPoint() 
  res_correct = checkGradientOperatorByElement() .and. &
    &           checkGradientOperatorByPoint() .and. &
    &           checkGradientOperatorByIndex()
  !*****************************************************************************

  call print_self_status()
  if (res_correct) then
    write(logUnit(1),*) 'PASSED'
  end if

contains

  function checkGradientOperatorByElement() result(result)
    !--------------------------------------------------------------------------!
    logical :: result
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: res(:)
    real(kind=rk) :: analyticalRes(27,3)
    real(kind=rk) :: expected_res, actual_res
    integer, allocatable :: elemPos(:)
    integer :: elem_offset, iElem, iDim, iDof, nElems_track

    !--------------------------------------------------------------------------!
    write(logUnit(3),*) ''
    write(logUnit(3),*) 'Check Gradient Operator By Element '
    write(logUnit(3),*) ''
    result = .true.

    ! prepare the analytical solutuin
    call get_anl_3d(analyticalRes)

    nElems_track = 2
    allocate(elemPos(nElems_track))
    elemPos = (/ 1, 5 /)

    ! levelpointer should point to the correct tracking element element 1 and 5
    ! setting it manually here
    levelpointer = (/ 1,0,0,0,5 /)

    !! The result is two components long
    allocate(res(nElems_track * nDofs * nComponents_res))

    write(logUnit(3),*) 'Get element for ',             &
      & trim(equation%varSys%varname%val(posOfDerived))

    ! access state array for given elemPos
    call equation%varSys%method%val(posOfDerived)%get_element( &
      & varSys  = equation%varSys,                             &
      & elemPos = elemPos,                                     &
      & time    = general%simControl%now,                      &
      & tree    = tree,                                        &
      & nElems  = nElems_track,                                &
      & nDofs   = nDofs,                                       &
      & res     = res                                          )
    
    do iElem = 1, nElems_track

      elem_offset =  (iElem-1)*nDofs*nComponents_res

      do iDim = 1, nDimensions
        do iDof = 1, nDofs

          expected_res = analyticalRes(iDof, iDim)
          actual_res = res( (iElem-1)*nComponents_res*nDofs &
            & + (iDof-1)*nComponents_res                    &
            & + iDim                                        )

          if (abs(actual_res - expected_res) > eps) then
            result = .false.
            write(logUnit(1),*) 'Failed at Dof ', iDof, &
              & ' Component ', iComp
            write(logUnit(1),*) 'Expected was ', expected_res, ' &
              & but result was ', res(elem_offset+1)
          end if

        end do
      end do
    end do

  end function checkGradientOperatorByElement


  function checkGradientOperatorByPoint() result(result)
    !--------------------------------------------------------------------------!
    logical :: result
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: res(:)
    real(kind=rk), allocatable :: expected_res(:)
    real(kind=rk) :: analyticalRes(27,3)
    real(kind=rk), allocatable :: polyVal(:,:)
    real(kind=rk) :: point(1,3)
    real(kind=rk) :: mxx, mxy, mxz, myx, myy, myz, mzx, mzy, mzz
    real(kind=rk) :: bary(3), refCoord(1,3), halfwidth, extent
    integer :: iDof
    !--------------------------------------------------------------------------!
    write(logUnit(3),*) ''
    write(logUnit(3),*) 'Check Gradient Operator By Point '
    write(logUnit(3),*) ''
    result = .true.

    ! prepare the analytical solutuin
    call get_anl_3d(analyticalRes)

    point(1,:) = (/ 1.0, 1.0, 1.0 /)

    ! levelpointer should point to the correct tracking element element 1 
    levelpointer = (/ 1,0,0,0,5 /)

    ! allcoate array for all gradient at one point
    allocate(res( 1 *  nComponents_res))

    write(logUnit(3),*) 'Get point for ',             &
      & trim(equation%varSys%varname%val(posOfDerived))

    ! access state array for given elemPos
    call equation%varSys%method%val(posOfDerived)%get_point( &
      & varSys  = equation%varSys,                           &
      & point   = point,                                     &
      & time    = general%simControl%now,                    &
      & tree    = tree,                                      &
      & nPnts   = 1,                                         &
      & res     = res                                        )

    ! get the physical coordinates in reference element
    bary(:) = (/ 1.0, 1.0, 1.0 /)
    extent = meshLength/(2**refinementlevel)
    halfwidth = extent/2
    refCoord(1,1) = (point(1,1)- bary(1))/halfwidth
    refCoord(1,2) = (point(1,2)- bary(1))/halfwidth
    refCoord(1,3) = (point(1,3)- bary(1))/halfwidth
    
    ! get the polynomial values ( Li(x,y,z) ) at the point 
    allocate(polyVal(nDofs,1))
    polyVal = 0.0
    select case(nDimensions)
    case(1)
      call atl_evalLegendreTensPoly1d( coords        = refCoord,   &
        &                              ncoords       = 1,          &
        &                              maxPolyDegree = polyDegree, &
        &                              basisType     = basisType,  &
        &                              polyVal       = polyVal     )
    case(2)
      call atl_evalLegendreTensPoly2d( coords        = refCoord,   &
        &                              ncoords       = 1,          &
        &                              maxPolyDegree = polyDegree, &
        &                              basisType     = basisType,  &
        &                              polyVal       = polyVal     )
    case(3)
      call ply_evalLegendreTensPoly( coords        = refCoord,   &
        &                            ncoords       = 1,          &
        &                            maxPolyDegree = polyDegree, &
        &                            basisType     = basisType,  &
        &                            polyVal       = polyVal     )
    end select
    
    !init everything with 0
    mxx=0
    mxy=0
    mxz=0

    do iDof = 1, nDofs
      ! get gradient of first component in x direction --> mxx
      mxx = mxx + polyVal(iDof,1) *analyticalRes(iDof,1)
      ! get gradient of first component in y direction --> mxy
      mxy = mxy + polyVal(iDof,1) *analyticalRes(iDof,2)
      ! get gradient of first component in z direction --> mxz
      mxz = mxz + polyVal(iDof,1) *analyticalRes(iDof,3)
    end do

    ! since state(:, iDof, 1) = iDof.,  state(:, iDof, 2) = iDof * 2_rk
    ! and state(:, iDof, 3) = iDof * 3_rk we can get the y and z component 
    ! via multiplication
    myx = mxx *2    
    myy = mxy *2
    myz = mxz *2
    mzx = mxx *3    
    mzy = mxy *3
    mzz = mxz *3
   
   !! write(logUnit(1),*) 'analytical gradient momentum in x ', mxx, mxy, mxz
   !! write(logUnit(1),*) 'analytical gradient momentum in y ', myx, myy, myz
   !! write(logUnit(1),*) 'analytical gradient momentum in z ', mzx, mzy, mzz

    ! allocate expected res, which should have all gradient at 1 point
    allocate(expected_res( 1 *  nComponents_res))
    ! sort the expected_res accodingly
    expected_res(:) = (/mxx, mxy, mxz, myx, myy, myz, mzx, mzy, mzz /)

    ! check if point value is equal to expected

    if ( any( abs(res - expected_res) > eps ) ) then
      result = .false.
      write(logUnit(1),*) 'Failed: calculated gradient is ', res
      write(logUnit(1),*) '        expected was           ', expected_res
    end if

  end function checkGradientOperatorByPoint


  function checkGradientOperatorByIndex() result(result)
    !--------------------------------------------------------------------------!
    logical :: result
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: res(:,:)
    real(kind=rk), allocatable :: expected_res(:), pnt_res(:)
    real(kind=rk) :: analyticalRes(27,3)
    real(kind=rk), allocatable :: polyVal(:,:)
    real(kind=rk), allocatable :: points(:,:)
    real(kind=rk) :: mxx, mxy, mxz, myx, myy, myz, mzx, mzy, mzz
    integer :: iPnt, nPnts, iLevel, iDof
    real(kind=rk), allocatable :: bary(:,:)
    real(kind=rk) :: refCoord(1,3)
    real(kind=rk) :: halfwidth, extent
    integer, allocatable :: idx(:,:)
    !--------------------------------------------------------------------------!
    write(logUnit(3),*) ''
    write(logUnit(3),*) 'Check Gradient Operator By Index '
    write(logUnit(3),*) ''
    result = .true.

    ! prepare the analytical solutuin
    call get_anl_3d(analyticalRes)
   
    npnts = 3
    allocate ( points(npnts,3) )
    points(1,:) = (/ 1.0, 1.0, 1.0 /)
    points(2,:) = (/ 3.0, 1.0, 1.0 /)
    points(3,:) = (/ 1.0, 3.0, 1.0 /)
    levelpointer = (/ 1,2,3,4,5 /)


    ! First we need to setup the index
    do iLevel = tree%global%minLevel, tree%global%maxLevel
      allocate(idx(iLevel,npnts) )
    
      write(logUnit(3),*) 'Set up indices for ',        &
        & trim(equation%varSys%varname%val(posOfDerived))
      ! setup indices for the derived variable
      call equation%varSys%method%val(posOfDerived)%setup_indices( &
          &  VarSys = equation%VarSys,   &
          &  point  = points,            &
          &  iLevel = iLevel,            &
          &  tree   = tree,              &
          &  nPnts  = npnts,             &
          &  idx    = idx(iLevel,:)      )
      write(logUnit(3),*) 'Index for derived on level', iLevel, &
          &                 'are=', idx(iLevel,:)
    end do

    ! allcoate array for all gradient at one point
    allocate(res (tree%global%minLevel:tree%global%maxLevel, &
       &          nComponents_res*nPnts)                     )

    do iLevel = tree%global%minLevel, tree%global%maxLevel
      write(logUnit(3),*) 'Get value of Index for ',             &
        & trim(equation%varSys%varname%val(posOfDerived))
      write(logUnit(3),*) ''

      call equation%varSys%method%val(posOfDerived)%get_valOfIndex( &
        & varSys  = equation%varSys,                                &
        & time    = general%simControl%now,                         &
        & iLevel  = iLevel,                                         &
        & idx     = idx(iLevel,:),                                  &
        & nVals   = nPnts,                                          &
        & res     = res(iLevel,:)                                   )
  
      write(logUnit(1),*) ''
      write(logUnit(1),*) 'in utest computed gradient ', res(iLevel,:)
      write(logUnit(1),*) ''
    end do!iLevel


    do iLevel = tree%global%minLevel, tree%global%maxLevel
      !! get the expected gradient
      write(logUnit(1),*) ' get the expected gradient'
      allocate( bary(npnts,3) )

      ! allocate expected res, which should have all gradient at 1 point
      allocate(expected_res(nComponents_res))
      allocate(pnt_res(nComponents_res))
      ! get the polynomial values ( Li(x,y,z) ) at the point 
      allocate(polyVal(nDofs,1))


      ! get the physical coordinates in reference element
      ! for points(1,:) = (/ 1.0, 1.0, 1.0 /) bary is
      bary(1,:) = (/ 1.0, 1.0, 1.0 /)
      ! for points(2,:) = (/ 3.0, 1.0, 1.0 /) bary is
      bary(2,:) = (/ 3.0, 1.0, 1.0 /)
      ! for points(3,:) = (/ 1.0, 3.0, 1.0 /) bary is
      bary(3,:) = (/ 1.0, 3.0, 1.0 /)
      ! mesh config is
      extent = meshLength/(2**refinementlevel)
      halfwidth = extent/2
      do iPnt=1,nPnts
        refCoord(1,1) = (points(iPnt,1)- bary(iPnt,1))/halfwidth
        refCoord(1,2) = (points(iPnt,2)- bary(iPnt,2))/halfwidth
        refCoord(1,3) = (points(iPnt,3)- bary(iPnt,3))/halfwidth

        ! init polyVal
        polyVal = 0.0_rk
        select case(nDimensions)
        case(1)
          call atl_evalLegendreTensPoly1d( coords        = refCoord,   &
            &                              ncoords       = 1,          &
            &                              maxPolyDegree = polyDegree, &
            &                              basisType     = basisType,  &
            &                              polyVal       = polyVal     )
        case(2)
          call atl_evalLegendreTensPoly2d( coords        = refCoord,   &
            &                              ncoords       = 1,          &
            &                              maxPolyDegree = polyDegree, &
            &                              basisType     = basisType,  &
            &                              polyVal       = polyVal     )
        case(3)
          call ply_evalLegendreTensPoly( coords        = refCoord,   &
            &                            ncoords       = 1,          &
            &                            maxPolyDegree = polyDegree, &
            &                            basisType     = basisType,  &
            &                            polyVal       = polyVal     )
        end select

        !init everything with 0
        pnt_res = 0.0_rk
        expected_res = 0.0_rk
        mxx=0.0_rk
        mxy=0.0_rk
        mxz=0.0_rk
        myx=0.0_rk
        myy=0.0_rk
        myz=0.0_rk
        mzx=0.0_rk
        mzy=0.0_rk
        mzz=0.0_rk

        do iDof = 1, nDofs
          ! get gradient of first component in x direction --> mxx
          mxx = mxx + polyVal(iDof,1) *analyticalRes(iDof,1)
          ! get gradient of first component in y direction --> mxy
          mxy = mxy + polyVal(iDof,1) *analyticalRes(iDof,2)
          ! get gradient of first component in z direction --> mxz
          mxz = mxz + polyVal(iDof,1) *analyticalRes(iDof,3)
        end do!iDof

        ! since state(:, iDof, 1) = iDof.,  state(:, iDof, 2) = iDof * 2_rk
        ! and state(:, iDof, 3) = iDof * 3_rk we can get the y and z component 
        ! via multiplication
        myx = mxx *2    
        myy = mxy *2
        myz = mxz *2
        mzx = mxx *3    
        mzy = mxy *3
        mzz = mxz *3
   
        write(logUnit(1),*) ''
        write(logUnit(1),*) 'expected gradient momentum in x ', mxx, mxy, mxz
        write(logUnit(1),*) 'expected gradient momentum in y ', myx, myy, myz
        write(logUnit(1),*) 'expected gradient momentum in z ', mzx, mzy, mzz

        ! sort the expected_res accodingly
        expected_res(:) = (/mxx, mxy, mxz, myx, myy, myz, mzx, mzy, mzz /)
        pnt_res(:) = res(iLevel,(iPnt-1)*nComponents_res+1:               &
          &                     (iPnt-1)*nComponents_res +nComponents_res )
        ! check if point value is equal to expected
        if ( any( abs(pnt_res - expected_res) > eps ) ) then
          result = .false.
          write(logUnit(1),*) ''
          write(logUnit(1),*) 'Failed for point' , iPnt, ':'
          write(logUnit(1),*) '  calculated gradient is ', pnt_res
          write(logUnit(1),*) '  expected was           ', expected_res
        end if

      end do!iPnts

      deallocate(expected_res)
      deallocate(res)
    end do! iLevel

  end function checkGradientOperatorByIndex


  subroutine get_anl_3d(anl_diff)
    real(kind=rk), intent(inout)   :: anl_diff(:,:)
    anl_diff(1,1)  =  2.0
    anl_diff(2,1)  =  9.0
    anl_diff(3,1)  =  0.0
    anl_diff(4,1)  =  5.0
    anl_diff(5,1)  = 18.0
    anl_diff(6,1)  =  0.0
    anl_diff(7,1)  =  8.0
    anl_diff(8,1)  = 27.0
    anl_diff(9,1)  =  0.0
    anl_diff(10,1) = 11.0
    anl_diff(11,1) = 36.0
    anl_diff(12,1) =  0.0
    anl_diff(13,1) = 14.0
    anl_diff(14,1) = 45.0
    anl_diff(15,1) =  0.0
    anl_diff(16,1) = 17.0
    anl_diff(17,1) = 54.0
    anl_diff(18,1) =  0.0
    anl_diff(19,1) = 20.0
    anl_diff(20,1) = 63.0
    anl_diff(21,1) =  0.0
    anl_diff(22,1) = 23.0
    anl_diff(23,1) = 72.0
    anl_diff(24,1) =  0.0
    anl_diff(25,1) = 26.0
    anl_diff(26,1) = 81.0
    anl_diff(27,1) =  0.0

    anl_diff(1,2)  =  4.0
    anl_diff(2,2)  =  5.0
    anl_diff(3,2)  =  6.0
    anl_diff(4,2)  = 21.0
    anl_diff(5,2)  = 24.0
    anl_diff(6,2)  = 27.0
    anl_diff(7,2)  =  0.0
    anl_diff(8,2)  =  0.0
    anl_diff(9,2)  =  0.0
    anl_diff(10,2) = 13.0
    anl_diff(11,2) = 14.0
    anl_diff(12,2) = 15.0
    anl_diff(13,2) = 48.0
    anl_diff(14,2) = 51.0
    anl_diff(15,2) = 54.0
    anl_diff(16,2) =  0.0
    anl_diff(17,2) =  0.0
    anl_diff(18,2) =  0.0
    anl_diff(19,2) = 22.0
    anl_diff(20,2) = 23.0
    anl_diff(21,2) = 24.0
    anl_diff(22,2) = 75.0
    anl_diff(23,2) = 78.0
    anl_diff(24,2) = 81.0
    anl_diff(25,2) =  0.0
    anl_diff(26,2) =  0.0
    anl_diff(27,2) =  0.0

    anl_diff(1,3)  = 10.0
    anl_diff(2,3)  = 11.0
    anl_diff(3,3)  = 12.0
    anl_diff(4,3)  = 13.0
    anl_diff(5,3)  = 14.0
    anl_diff(6,3)  = 15.0
    anl_diff(7,3)  = 16.0
    anl_diff(8,3)  = 17.0
    anl_diff(9,3)  = 18.0
    anl_diff(10,3) = 57.0
    anl_diff(11,3) = 60.0
    anl_diff(12,3) = 63.0
    anl_diff(13,3) = 66.0
    anl_diff(14,3) = 69.0
    anl_diff(15,3) = 72.0
    anl_diff(16,3) = 75.0
    anl_diff(17,3) = 78.0
    anl_diff(18,3) = 81.0
    anl_diff(19,3) =  0.0
    anl_diff(20,3) =  0.0
    anl_diff(21,3) =  0.0
    anl_diff(22,3) =  0.0
    anl_diff(23,3) =  0.0
    anl_diff(24,3) =  0.0
    anl_diff(25,3) =  0.0
    anl_diff(26,3) =  0.0
    anl_diff(27,3) =  0.0
  end subroutine get_anl_3d

end program atl_varSys_gradientOperatorByElement_test

