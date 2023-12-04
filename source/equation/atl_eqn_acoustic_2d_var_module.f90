! Copyright (c) 2014, 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014-2015 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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

!> \brief module to configure information about the variables of the acoustic
!! equations
module atl_eqn_acoustic_2d_var_module
  use, intrinsic :: iso_c_binding,  only: c_loc, c_ptr

  use env_module,                   only: rk
  use tem_time_module,              only: tem_time_type
  use tem_varSys_module,            only: tem_varSys_type,              &
    &                                     tem_varSys_init,              &
    &                                     tem_varSys_append_stateVar,   &
    &                                     tem_varSys_proc_element,      &
    &                                     tem_varSys_proc_point,        &
    &                                     tem_varSys_proc_setParams,    &
    &                                     tem_varSys_proc_getParams,    &
    &                                     tem_varSys_proc_setupIndices, &
    &                                     tem_varSys_proc_getValOfIndex
  use tem_varMap_module,            only: tem_possible_variable_type, &
    &                                     init, append
  use tem_grow_array_module,        only: init, append
  use tem_dyn_array_module,         only: init, append

  use ply_poly_project_module,      only: ply_poly_project_type, &
    &                                     assignment(=)

  use atl_varSys_module,            only: atl_varSys_solverData_type,    &
    &                                     atl_varSys_getStateForElement, &
    &                                     atl_varSys_getStateForPoint,   &
    &                                     atl_get_new_varSys_data_ptr,   &
    &                                     atl_varSys_setupStateIndices,  &
    &                                     atl_varSys_getStateValofIndex
  use atl_equation_module,          only: atl_equations_type
  use atl_cube_elem_module,         only: atl_cube_elem_type
  use atl_source_types_module,      only: atl_eqn_sourceMap_type, &
    &                                     atl_source_op_type
  use atl_eqn_acoustic_var_module,  only: atl_eval_massInjection, &
    &                                     atl_eval_force
  use atl_equation_source_module,   only: atl_equation_evaluate_source_modal, &
    &                                     atl_compute_source_interface

  implicit none

  private

  public :: atl_init_acoustic_2d_vars
  public :: atl_append_acoustic_2d_vars
  public :: atl_init_acoustic_2d_sourceTerms


contains


!*******************************************************************************
  !> \brief init the variables for acoustic equation
  subroutine atl_init_acoustic_2d_vars( equation, methodData )
    !---------------------------------------------------------------------------
    !> The equation system
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type) :: methodData
    !---------------------------------------------------------------------------
    ! initialize variable system
    call tem_varSys_init( me = equation%varSys, systemName = 'acoustic_2d' )

    ! Append conservative Variables to variable system
    call atl_append_acoustic_2d_vars(equation, methodData)

  end subroutine atl_init_acoustic_2d_vars
!*******************************************************************************


!*******************************************************************************
  !> \brief append the variables for acoustic simulations
  !! the acoustic equation only has 'primitive' variables or different speaking,
  !! the equation describe the pertubation in primitive variables
  subroutine atl_append_acoustic_2d_vars(equation, methodData)
    !---------------------------------------------------------------------------
    !> The equation type .
    type(atl_equations_type), intent(inout) :: equation
    !> the pointer to the data required for the varsys
    type(atl_varSys_solverData_type), target :: methodData
    procedure(tem_varSys_proc_point), pointer :: get_point => NULL()
    procedure(tem_varSys_proc_element), pointer :: get_element => NULL()
    procedure(tem_varSys_proc_setParams), pointer :: set_params => Null()
    procedure(tem_varSys_proc_getParams), pointer :: get_params => Null()
    procedure(tem_varSys_proc_setupIndices), pointer ::setup_indices => NULL()
    procedure(tem_varSys_proc_getValOfIndex), pointer :: get_ValOfIndex => NULL()
    !---------------------------------------------------------------------------
    allocate(equation%stateVar(2))

    get_element => atl_varSys_getStateForElement
    get_point => atl_varSys_getStateForPoint
    setup_indices => atl_varSys_setupStateIndices
    get_valOfIndex => atl_varSys_getStateValofIndex

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'density',                               &
      & nComponents    = 1,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & pos            = equation%stateVar(1),                    &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

    call tem_varSys_append_stateVar(                              &
      & me             = equation%varSys,                         &
      & varName        = 'velocity',                              &
      & nComponents    = 2,                                       &
      & method_data    = atl_get_new_varSys_data_ptr(methodData), &
      & pos            = equation%stateVar(2),                    &
      & get_point      = get_point,                               &
      & get_element    = get_element,                             &
      & set_params     = set_params,                              &
      & get_params     = get_params,                              &
      & setup_indices  = setup_indices,                           &
      & get_valOfIndex = get_valOfIndex                           )

  end subroutine atl_append_acoustic_2d_vars
!*******************************************************************************


!*******************************************************************************
  !> Init source terms for flow simulations.
  subroutine atl_init_acoustic_2d_sourceTerms(possVars, eval_source)
    !---------------------------------------------------------------------------
    type(tem_possible_variable_type), intent(out) :: possVars
    type(atl_eqn_sourceMap_type), allocatable, intent(out) :: eval_source(:)
    !---------------------------------------------------------------------------
    integer :: pos
    !---------------------------------------------------------------------------

    allocate(eval_source(2))
    call init( me = possVars, length = 2 )

    call append( me          = possVars,        &
      &          varName     = 'massInjection', &
      &          nComponents = 1,               &
      &          pos         = pos              )
    eval_source(pos)%do => atl_eval_source_massInjection_2D

    call append( me          = possVars, &
      &          varName     = 'force',  &
      &          nComponents = 2,        &
      &          pos         = pos       )
    eval_source(pos)%do => atl_eval_source_force_2D

  end subroutine atl_init_acoustic_2d_sourceTerms
!*******************************************************************************


!*******************************************************************************
  subroutine atl_eval_source_massInjection_2D( fun, varSys, time, mesh,    &
    &                                          poly_proj, currentLevel,    &
    &                                          state, material, sourcedata )
    !---------------------------------------------------------------------------

    !> Description of method to update source
    class(atl_source_op_type), intent(in) :: fun

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> Current level mesh information
    type(atl_cube_elem_type), intent(in) :: mesh

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> current level
    integer, intent(in) :: currentLevel

    !> The state in modal space.
    !! This is needed for several source terms that have to be applied to the
    !! current state
    real(kind=rk), intent(in) :: state(:,:,:)

    !> Material description for the complete domain. Used for evaluation of some
    !! source terms.
    real(kind=rk), intent(in) :: material(:)

    !> The source data to update. When all source terms are added to this
    !! buffer, it is applied to the state.
    real(kind=rk), intent(inout) :: sourcedata(:,:,:)
    !---------------------------------------------------------------------------
    procedure(atl_compute_source_interface) , pointer:: evaluate_source
    !---------------------------------------------------------------------------

    ! Set the function pointer for the evaluation of spongeLayer_2d
    evaluate_source => atl_eval_massInjection

    ! Call the common function for updating the sourceData
    call atl_equation_evaluate_source_modal(    &
      & fun          = fun,                     &
      & varSys       = varSys,                  &
      & currentLevel = currentLevel,            &
      & nDim         = 2,                       &
      & time         = time,                    &
      & eval_rhs     = evaluate_source,         &
      & state        = state,                   &
      & poly_proj    = poly_proj,               &
      & polyProjBody = poly_proj%body_2d,       &
      & sourceData   = sourceData               )

  end subroutine atl_eval_source_massInjection_2D
!*******************************************************************************


!*******************************************************************************
  subroutine atl_eval_source_force_2D( fun, varSys, time, mesh, poly_proj, &
    &                                  currentLevel, state, material,      &
    &                                  sourcedata                          )
    !---------------------------------------------------------------------------

    !> Description of method to update source
    class(atl_source_op_type), intent(in) :: fun

    !> variable system
    type(tem_varSys_type), intent(in) :: varSys

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> Current level mesh information
    type(atl_cube_elem_type), intent(in) :: mesh

    !> Parameters for projection
    type(ply_poly_project_type), intent(inout) :: poly_proj

    !> current level
    integer, intent(in) :: currentLevel

    !> The state in modal space.
    !! This is needed for several source terms that have to be applied to the
    !! current state
    real(kind=rk), intent(in) :: state(:,:,:)

    !> Material description for the complete domain. Used for evaluation of some
    !! source terms.
    real(kind=rk), intent(in) :: material(:)

    !> The source data to update. When all source terms are added to this
    !! buffer, it is applied to the state.
    real(kind=rk), intent(inout) :: sourcedata(:,:,:)
    !---------------------------------------------------------------------------
    procedure(atl_compute_source_interface) , pointer:: evaluate_source
    !---------------------------------------------------------------------------
    ! Set the function pointer for the evaluation of spongeLayer_2d
    evaluate_source => atl_eval_force

    ! Call the common function for updating the sourceData
    call atl_equation_evaluate_source_modal(    &
      & fun          = fun,                     &
      & varSys       = varSys,                  &
      & currentLevel = currentLevel,            &
      & nDim         = 2,                       &
      & time         = time,                    &
      & eval_rhs     = evaluate_source,         &
      & state        = state,                   &
      & poly_proj    = poly_proj,               &
      & polyProjBody = poly_proj%body_2d,       &
      & sourceData   = sourceData               )

  end subroutine atl_eval_source_force_2D
!*******************************************************************************

end module atl_eqn_acoustic_2d_var_module
