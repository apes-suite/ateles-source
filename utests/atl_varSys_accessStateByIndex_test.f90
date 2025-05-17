! Copyright (c) 2016-2017 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016-2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2018 Harald Klimach <harald.klimach@uni-siegen.de>
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

program atl_varSys_stateVarByIndex_test
  use, intrinsic :: iso_c_binding,  only: C_NEW_LINE, c_loc

  use env_module,                   only: rk, long_k,            &
    &                                     print_self_status, &
    &                                     eps
  use flu_binding,                  only: flu_state

  use aotus_module,                 only: open_config_chunk

  use tem_logging_module,           only: logUnit, tem_logging_init
  use tem_general_module,           only: tem_general_type, &
    &                                     tem_load_general, &
    &                                     tem_start
  use treelmesh_module,             only: treelmesh_type, load_tem
  use tem_varSys_module,            only: tem_varSys_init,               &
    &                                     tem_varSys_append_stateVar,    &
    &                                     tem_varSys_proc_point,         &
    &                                     tem_varSys_proc_element,       &
    &                                     tem_varSys_proc_setParams,     &
    &                                     tem_varSys_proc_getParams,     &
    &                                     tem_varSys_proc_setupIndices,  &
    &                                     tem_varSys_proc_getvalofIndex
  use tem_geometry_module,          only: tem_CoordOfReal, &
    &                                     tem_PosofId, tem_BaryOfId, tem_ElemSize
  use tem_topology_module,          only: tem_IdOfCoord, &
    &                                     tem_levelOf
  use tem_tools_module,             only: tem_horizontalSpacer

  use ply_dof_module,               only: P_Space, Q_Space
  use ply_poly_project_module,      only: ply_poly_project_type

  use atl_varSys_module,            only: atl_varSys_solverData_type,   &
    &                                     atl_varSys_getStateForPoint,  &
    &                                     atl_get_new_varSys_data_ptr,  &
    &                                     atl_varSys_setupStateIndices, &
    &                                     atl_varSys_getStateValOfIndex,&
    &                                     atl_init_varsys_solverdata
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_kerneldata_module,        only: atl_statedata_type, &
    &                                     atl_kerneldata_type
  use atl_equation_module,          only: atl_equations_type
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
  integer, parameter :: refinementLevel = 2
  ! Number of components for the state variable
  integer, parameter :: nComponents_stateVars(3) = (/ 1, 3, 1 /)
  ! Polynomial degree
  integer, parameter :: polyDegree = 3
  ! Type of polynomial space to be used
  integer, parameter :: basisType = Q_Space
  !*****************************************************************************
  ! Derived from parameters
  !
  ! The number of elements
  integer, parameter :: nElements = (2**refinementLevel)**nDimensions
  ! Total number of components in the state
  !integer, parameter :: nComponents = sum(nComponents_stateVars)
  ! workaround for Intel 15
  integer, parameter :: nComponents = 5
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
  type(atl_statedata_type), target, allocatable :: state(:)
  type(atl_kerneldata_type), target, allocatable :: kernel(:)
  type(atl_scheme_type), target, allocatable :: scheme(:)
  integer, allocatable, target   :: levelPointer(:)
  type(atl_varSys_solverData_type), target :: varSys_data
  type(atl_cube_elem_type), target, allocatable :: mesh(:)
  type(ply_poly_project_type), target, allocatable :: polyProj(:)
  integer, target :: poly_proj_pos(1) = (/ 1 /)

  integer :: iElem, iDof, iComponent
  logical :: res_correct = .true.

  character, parameter :: nl = C_NEW_LINE
  character(len=110) :: cubeconf
  character(len=240) :: schemeconf
  character(len=1) :: basistype_chr
  character(len=7) :: schemetype
  character(len=5) :: meshType

  integer :: i, nElems
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
  call tem_start('atl_varSys_accessState_byIndex_test unit test', general)

  !*****************************************************************************
  write(logUnit(3), *) 'Initialize the mesh'
  call open_config_chunk(L = conf, chunk = trim(cubeconf))

  call tem_load_general( me = general, conf = conf )
  call tem_logging_init(level = 6, rank  = general%proc%rank)

  ! Load the mesh first.
  call load_tem(me     = tree,                   &
    &           conf   = conf,                   &
    &           myPart = general%proc%rank,      &
    &           nParts = general%proc%comm_size, &
    &           comm   = general%proc%comm       )
  call atl_addTimers_perElem(tree)  

  ! get the level pointer for single level run 
  nElems = 2**(refinementLevel*nDimensions) ! 2^refinemenut^3 
  levelpointer = [ (i, i=1,nElems)]
  write(*,*) 'levelpointer', levelpointer
  !*****************************************************************************
  write(logUnit(3),*) 'Initializing state'
  !! State for every level
  allocate(state(refinementLevel))
  allocate(kernel(refinementLevel))
  allocate(state(refinementLevel)%state(nElements,nDofs,nComponents))
  state(refinementLevel)%state = 0_rk
  do iComponent = 1, nComponents
    do iDof = 1, nDofs
      do iElem = 1, nElements
        if (iDof .lt. 5) then
          state(refinementLevel)%state(iElem, iDof, iComponent) = iDof
        end if
      end do
    end do
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
  write(logUnit(3),*) 'Initializing Solver data'
  ! Initialize the global instance of atl_varSys_data
  call atl_init_varSys_solverdata(     &
    & this           = varSys_data,    &
    & equation       = equation,       &
    & tree           = tree,           &
    & schemeList     = scheme,         &
    & statedataList  = state,          &
    & kerneldataList = kernel,         &
    & meshList       = mesh,           &
    & polyProjList   = polyProj,       &
    & levelPointer   = levelPointer,   &
    & polyProjPos    = poly_proj_pos   )

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing variable system'
  call tem_varSys_init(                &
    & me         = equation%varSys,    &
    & systemName = 'utest_accessState' )

  !*****************************************************************************
  ! This is the routine we want to test.
  setup_indices => atl_varSys_setupStateIndices
  get_valOfIndex => atl_varSys_getStateValOfIndex
  !*****************************************************************************

  call tem_varSys_append_stateVar(                              &
    & me             = equation%varSys,                         &
    & varName        = 'state1',                                &
    & nComponents    = nComponents_stateVars(1),                &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )
  call tem_varSys_append_stateVar(                              &
    & me             = equation%varSys,                         &
    & varName        = 'state2',                                &
    & nComponents    = nComponents_stateVars(2),                &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )
  call tem_varSys_append_stateVar(                              &
    & me             = equation%varSys,                         &
    & varName        = 'state3',                                &
    & nComponents    = nComponents_stateVars(3),                &
    & method_data    = atl_get_new_varSys_data_ptr(methodData), &
    & get_point      = get_point,                               &
    & get_element    = get_element,                             &
    & set_params     = set_params,                              &
    & get_params     = get_params,                              &
    & setup_indices  = setup_indices,                           &
    & get_valOfIndex = get_valOfIndex                           )

  res_correct = checkIndexwiseAccess()

  call print_self_status()
  if (res_correct) then
    write(*,*) 'PASSED'
  end if

contains

  function checkIndexwiseAccess() result(result)
    !--------------------------------------------------------------------------!
    logical :: result
    !--------------------------------------------------------------------------!
    real(kind=rk), allocatable :: points(:,:)
    real(kind=rk), allocatable :: res(:)
    integer, allocatable :: indices(:,:)
    real(kind=rk) :: analytical_res
    integer :: point_offset
    integer :: iComponent, iPoint, iVar, nComp, nPoints_track
    integer :: iLevel
    !--------------------------------------------------------------------------!
    result = .true.

    nPoints_track = 5
    allocate(points(nPoints_track,3))
    do iPoint = 1, nPoints_track
      points(iPoint,1) = 0.3*iPoint
      points(iPoint,2) = 0.03*iPoint
      points(iPoint,3) = 0.003*iPoint
      write(logUnit(3),*) 'points: ', points(iPoint,:)
    end do

    write(logUnit(3),*) "nstateVars = ", equation%varSys%nStateVars

    call tem_horizontalSpacer(funit = logUnit(1))
    ! Setup the indices for these points
    write(logUnit(3),*) "Setup the indices"
    do iLevel = tree%global%minLevel, tree%global%maxLevel
      allocate(indices(iLevel,nPoints_track))
      write(logUnit(3),*) "iLevel = ", iLevel
      do iVar = 1, equation%varSys%nStateVars
        call equation%varSys%method%val(iVar)%setup_indices( &
          & varSys = equation%VarSys,                        &
          & point  = points,                                 &
          & iLevel = iLevel,                                 &
          & tree   = tree,                                   &
          & nPnts  = nPoints_track,                          &
          & idx    = indices(iLevel,:)                       )
      end do
    end do

    ! get the values for the indices
    call tem_horizontalSpacer(funit = logUnit(1))
    write(logUnit(3),*) "Get the State values of the indices"
    do iLevel = tree%global%minLevel, tree%global%maxLevel
      do iVar = 1, equation%varSys%nStateVars

        call tem_horizontalSpacer(funit = logUnit(1))
        write(logUnit(3),*) 'Get val of index for ', &
          & trim(equation%varSys%varname%val(iVar)), &
          & ' on level = ', iLevel

        nComp = equation%varSys%method%val(iVar)%nComponents
        allocate(res(nPoints_track * nComp))

        write(logUnit(3),*) "nComp of this variable ", nComp
        write(logUnit(3),*) 'Indices are ', indices(iLevel,:)

        ! access state array for given elemPos
        call equation%varSys%method%val(iVar)%get_valOfIndex( &
          & varSys  = equation%varSys,                        &
          & time    = general%simControl%now,                 &
          & iLevel  = iLevel,                                 &
          & idx     = indices(iLevel,:),                      &
          & nVals   = nPoints_track,                          &
          & res     = res                                     )

        write(logUnit(3),*) 'res from get_StatevalOfIndex: ', res

        do iPoint = 1, nPoints_track
          call tem_horizontalSpacer(funit = logUnit(1))

          point_offset = (iPoint-1)*equation%varSys%method%val(iVar)%nComponents

          write(logUnit(3),*) 'iPoint=', iPoint
          write(logUnit(3),*) 'point_offset=', point_offset

          analytical_res = valX(points(iPoint,:))
          write(logUnit(3),*) ' analytical_res=',  analytical_res

          do iComponent = 1, nComp
            ! the epsilon needs to be adjusted to fit the double precision signifivant
            ! digits. So, it needs to me multiplied with significant integer digits
            ! of the result
            if (abs(res(point_offset + iComponent) - analytical_res) &
              & .gt. eps*10**int(log10(abs(analytical_res))+1)) then
              result = .false.
              write(logUnit(1),*) 'Failed at point ', iPoint, &
                & ' Component ', iComponent
              write(logUnit(1),*) 'Expected was ', analytical_res, &
                & ' but result was ', res(point_offset + iComponent)
            end if

          end do !iComp
        end do !iPoint
        deallocate(res)
      end do !iVar
    end do !iLevel
  end function checkIndexwiseAccess

  function valX(points) result (val)
    real(kind=rk) :: points(3)
    real(kind=rk) :: pointX
    real(kind=rk) :: val
    integer :: coord(4)
    integer(kind=long_k) :: treeID
    real(kind=rk) :: bary(3)
    real(kind=rk) :: extent
    real(kind=rk) :: halfwidth

    ! Scale the physical points to element local
    coord =  tem_CoordOfReal(tree, points , tree%global%maxLevel )
    treeId = tem_IdOfCoord(coord)
    bary = tem_BaryOfId(tree,TreeID)
    extent = tem_ElemSize( tree, TreeID)
    halfwidth = extent/2
    pointX = (points(1)-bary(1))/halfwidth

    ! legendre polynomials
    val = 10*pointX**3 + 4.5*pointX**2 - 4*pointX - 0.5

  end function valX

end program atl_varSys_stateVarByIndex_test
