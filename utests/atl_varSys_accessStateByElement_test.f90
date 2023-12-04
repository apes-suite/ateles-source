! Copyright (c) 2014-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

program atl_varSys_stateVarByElement_test
  use, intrinsic :: iso_c_binding,  only: C_NEW_LINE, c_loc

  use env_module,                   only: rk,                &
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

  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_get_new_varSys_data_ptr
  use atl_kerneldata_module,        only: atl_statedata_type
  use atl_equation_module,          only: atl_equations_type
  use atl_timer_module,             only: atl_addTimers_perElem  

  implicit none

  !****************************************************************************!
  ! PAREMETER
  !
  ! Number of dimensions
  integer, parameter :: nDimensions = 3
  ! Refinement level
  integer, parameter :: refinementLevel = 4
  ! Number of components for the state variable
  integer, parameter :: nComponents_stateVars(3) = (/ 1, 3, 1 /)
  ! Polynomial degree
  integer, parameter :: polyDegree = 3
  !*****************************************************************************
  ! Derived from parameters
  !
  ! The number of elements
  integer, parameter :: nElements = (2**refinementLevel)**nDimensions
  ! Total number of components in the state
  !integer, parameter :: nComponents = sum(nComponents_stateVars)
  ! Workaround for intel 15
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

  type(flu_state) :: conf

  type(treelmesh_type), target :: tree
  type(tem_general_type), target :: general

  type(atl_varSys_solverData_type), target :: methodData
  type(atl_equations_type), target :: equation
  type(atl_statedata_type), target, allocatable :: state(:)
  integer, target :: levelPointer(8)

  logical :: res_correct = .true.
  integer :: iElem, iDof, iComponent

  character, parameter :: nl = C_NEW_LINE
  character(len=110) :: cubeconf
  character(len=5) :: meshType
  !*****************************************************************************

  select case(nDimensions)
  case(1)
    meshType = 'line '
  case(2)
    meshType = 'slice'
  case(3)
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

  ! Init the Treelm environment
  call tem_start('atl_varSys_accessState_test unit test', 'utest', general)

  !*****************************************************************************
  write(logUnit(3), *) 'Initialize the mesh'
  call open_config_chunk(L = conf, chunk = trim(cubeconf))

  call tem_load_general( me = general, conf = conf )

  ! Load the mesh first.
  call load_tem( me     = tree,                   &
    &            conf   = conf,                   &
    &            myPart = general%proc%rank,      &
    &            nParts = general%proc%comm_size, &
    &            comm   = general%proc%comm       )

  call atl_addTimers_perElem(tree)  

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing state'
  allocate(state(tree%global%minLevel:tree%global%maxLevel))
  ! Allocate the state array for 3D, 8 elements, 4th order, 5 scalars
  allocate(state(refinementLevel)%state(nElements,nDofs,nComponents))
  do iComponent = 1, nComponents
    do iDof = 1, nDofs
      do iElem = 1, nElements
        state(refinementLevel)%state(iElem, iDof, iComponent) &
          & = iElem * 10000 + iDof * 100 + iComponent
      end do
    end do
  end do

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing method data'
  methodData%statedata_listPtr => state
  methodData%equationPtr => equation
  methodData%levelPointer => levelPointer

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing variable system'
  call tem_varSys_init(                &
    & me         = equation%varSys,    &
    & systemName = 'utest_accessState' )

  !*****************************************************************************
  ! This is the routine we want to test.
  get_element => atl_varSys_getStateForElement
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

  res_correct = checkElementwiseAccess()

  call print_self_status()
  if (res_correct) then
    write(*,*) 'PASSED'
  end if

contains

  function checkElementwiseAccess() result(result)
    !--------------------------------------------------------------------------!
    logical :: result
    !--------------------------------------------------------------------------!
    integer, allocatable :: elemPos(:)
    real(kind=rk), allocatable :: res(:)
    integer :: iVar, iElem, iDof, iComponent, nElems_track
    integer :: nComp, comp_offset, elem_offset
    real :: expected_res
    !--------------------------------------------------------------------------!
    result = .true.

    nElems_track = 5
    allocate(elemPos(nElems_track))
    elemPos = (/ 1, 3, 5, 7, 8 /)

    ! levelpointer should point to the correct tracking element element 1 and 5
    ! setting it manually here
    levelpointer = (/ 1,0,3,0,5,0,7,8 /)

    allocate(res(nElems_track * nComponents * nDofs))

    do iVar = 1, equation%varSys%nStateVars
      write(*,*) 'Get element for ', trim(equation%varSys%varname%val(iVar))
      nComp = equation%varSys%method%val(iVar)%nComponents

      ! access state array for given elemPos
      call equation%varSys%method%val(iVar)%get_element( &
        & varSys  = equation%varSys,                     &
        & elemPos = elemPos,                             &
        & time    = general%simControl%now,              &
        & tree    = tree,                                &
        & nElems  = nElems_track,                        &
        & nDofs   = nDofs,                               &
        & res     = res                                  )

      do iElem = 1, nElems_track
        write(*,*) 'iElem ', iElem, 'elemPos ', elemPos(iElem)

        elem_offset =  (iElem-1)*nComp*nDofs
        do iDof = 1, nDofs
          comp_offset = elem_offset + (iDof-1)*nComp

          write(*,*) 'iDof ', iDof , ' res ',&
            & res( comp_offset+1: &
            &      comp_offset+equation%varSys%method%val(iVar)%nComponents )

          do iComponent = 1, equation%varSys%method%val(iVar)%nComponents

            expected_res = elemPos(iElem) * 10000 + iDof * 100 &
                & + equation%varSys%method%val(iVar)%state_varPos(iComponent)

            if (eps .lt. abs(res(comp_offset+iComponent) - expected_res)) then
              result = .false.
              write(*,*) 'Failed at Element ', elemPos(iElem), ' iDof ', iDof, &
                & ' iComponent ',                                              &
                & equation%varSys%method%val(iVar)%mypos + iComponent
              write(*,*) 'Expected was ', expected_res, ' but result was ', &
                & res(comp_offset+iComponent)
            end if

          end do
        end do
      end do
    end do

  end function checkElementwiseAccess

end program atl_varSys_stateVarByElement_test
