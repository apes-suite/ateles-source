! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Harald Klimach <harald.klimach@uni-siegen.de>
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

program atl_dumpEquationTable_test

  use, intrinsic :: iso_c_binding,  only: C_NEW_LINE

  use env_module,                   only: newunit,  &
    &                                     labelLen, &
    &                                     init_env

  use flu_binding,                  only: flu_state
  use aotus_module,                 only: open_config_chunk
  use aot_out_module,               only: aot_out_type, &
    &                                     aot_out_open
  use aot_table_module,             only: aot_table_open, &
    &                                     aot_table_close

  use tem_float_module,             only: operator(.feq.)
  use tem_logging_module,           only: logUnit
  use tem_general_module,           only: tem_general_type, &
    &                                     tem_load_general
  use tem_varSys_module,            only: tem_varSys_init
  use tem_variable_module,          only: tem_variable_type, &
    &                                     tem_variable_load
  use tem_varMap_module,            only: tem_variable_loadMapping
  use tem_derived_module,           only: tem_varSys_append_luaVar
  use tem_stringKeyValuePair_module,          &
    & only: grw_stringKeyValuePairArray_type, &
    &       tem_stringKeyValuePair_type,      &
    &       init,                             &
    &       operator(==)

  use atl_varSys_module,            only: atl_varSys_solverData_type
  use atl_equation_module,          only: atl_equations_type
  use atl_equation_init_module,     only: atl_init_equation, atl_eqn_write
  use atl_source_types_module,      only: atl_init_source_type
  use atl_materialIni_module,       only: atl_append_newMaterialVars
  use atl_materialPrp_module,       only: atl_init_material_type

  implicit none

  ! ***************************************************************************!
  ! PAREMETER
  !
  character, parameter :: nl = C_NEW_LINE
  ! ***************************************************************************!


  ! ****************************************************************************
  ! Variables
  !
  type(aot_out_type) :: out_conf
  type(flu_state) :: eqnconf, varconf, emptyconf, readconf
  type(tem_general_type) :: general
  type(tem_variable_type), allocatable :: userVars(:)

  type(atl_equations_type), target :: equation, equation_read
  type(atl_varSys_solverData_type) :: methodData, methodData_read
  type(atl_init_source_type) :: initSource, initSource_read
  type(atl_init_material_type) :: initMaterial, initMaterial_read
  type(tem_stringKeyValuePair_type) :: strKvpA, strKvpB

  character(len=1000) :: eqnconf_str, varconf_str, file_str, line_str
  character(len=labelLen), allocatable :: materials(:)
  integer :: funit, iMat, tblHandle, iError
  integer, allocatable :: vError(:)
  logical :: isEqual = .true.
  ! ****************************************************************************

  call init_env()
  write(varconf_str, '(A)' )                                  &
    &    'variable = {'                                 // nl &
    & // '  {'                                          // nl &
    & // '    name = "global_euler_characteristic",'    // nl &
    & // '    ncomponents = 1,'                         // nl &
    & // '    vartype = "st_fun",'                      // nl &
    & // '    st_fun = { const = 0.0 }'                 // nl &
    & // '  },'                                         // nl &
    & // '  {'                                          // nl &
    & // '    name = "global_euler_relax_velocity",'    // nl &
    & // '    ncomponents = 3,'                         // nl &
    & // '    vartype = "st_fun",'                      // nl &
    & // '    st_fun = {'                               // nl &
    & // '      const = { 0.0, 0.0, 0.0 }'              // nl &
    & // '    }'                                        // nl &
    & // '  },'                                         // nl &
    & // '  {'                                          // nl &
    & // '    name = "global_euler_relax_temperature",' // nl &
    & // '    ncomponents = 1,'                         // nl &
    & // '    vartype = "st_fun",'                      // nl &
    & // '    st_fun = { const = 0.0 }'                 // nl &
    & // '  }'                                          // nl &
    & // '}'

  write(eqnconf_str, *)                                                   &
    &    'equation = {'                                             // nl &
    & // '  name = "euler",'                                        // nl &
    & // '  therm_cond = 2.555e-02,'                                // nl &
    & // '  isen_coef = 1.4,'                                       // nl &
    & // '  r = 296.0,'                                             // nl &
    & // '  cv = 740.0,'                                            // nl &
    & // '  material = {'                                           // nl &
    & // '    characteristic = "global_euler_characteristic",'      // nl &
    & // '    relax_velocity = "global_euler_relax_velocity",'      // nl &
    & // '    relax_temperature = "global_euler_relax_temperature"' // nl &
    & // '  }'                                                      // nl &
    & // '}'

  !*****************************************************************************
  write(logUnit(3), *) 'Initialize general data structures'
  call open_config_chunk( L = emptyconf, chunk = '' )

  call tem_load_general( me = general, conf = emptyconf )

  !*****************************************************************************
  write(logUnit(3), *) 'Initializing equation'
  call open_config_chunk( L = eqnconf, chunk = trim(eqnconf_str) )
  call atl_init_equation(          &
    & equation     = equation,     &
    & conf         = eqnconf,      &
    & varSys_data  = methodData,   &
    & initSource   = initSource,   &
    & initMaterial = initMaterial  )
  methodData%equationPtr => equation

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing variable system'
  call tem_varSys_init( me         = equation%varSys,    &
    &                   systemName = 'utest_accessState' )

  call open_config_chunk( L = varconf, chunk = trim(varconf_str) )
  call tem_variable_load( me     = userVars, &
    &                     conf   = varconf,  &
    &                     vError = vError    )

  call tem_varSys_append_luaVar( luaVar     = userVars,          &
    &                            varSys     = equation%varSys,   &
    &                            st_funList = equation%stFunList )

  !*****************************************************************************
  write(logUnit(3),*) 'Initializing material'
  call open_config_chunk( L = eqnconf, chunk = trim(eqnconf_str) )

  ! The following lines got extarcted from atl_materialIni_module. This is
  ! mostly to avoid initializaing the whole solver just to have an initialized
  ! variable systen and equation type, which we want to dump later on.
  call aot_table_open( L       = eqnconf,   &
    &                  thandle = tblHandle, &
    &                  key     = "equation" )
  if( tblHandle > 0 ) then
    ! Initilized keyValuePair dict for material
    call init( me = initMaterial%materialDict )

    call tem_variable_loadMapping(                 &
      & possVars = initMaterial%poss_materialVars, &
      & conf     = eqnconf,                        &
      & parent   = tblHandle,                      &
      & key      = "material",                     &
      & varSys   = equation%varSys,                &
      & varDict  = initMaterial%materialDict       )
    call aot_table_close( eqnconf, tblHandle )
  end if
  call atl_append_newMaterialVars(                   &
    & varSys       = equation%varSys,                &
    & varSys_data  = methodData,                     &
    & poss_matVars = initMaterial%poss_materialVars, &
    & materialFun  = equation%material,              &
    & variables    = initMaterial%materialDict       )
  write(logUnit(3),*) 'Number of materials: ', equation%material%nMat

  ! Store the positions of variables acting as data providers for materials
  ! in a sepparate list.
  allocate(materials(initMaterial%poss_materialVars%varName%nVals))
  do iMat=1, initMaterial%poss_materialVars%varName%nVals
    materials(iMat) = equation%varSys                                  &
      &                       %varname                                 &
      &                       %val( equation%material%matvar_pos(iMat) )
  end do

  !*****************************************************************************
  write(logUnit(3),*) 'Writing output'

  funit = newunit()
  open(unit=funit, status='scratch')
  call aot_out_open( put_conf = out_conf, &
    &                outUnit  = funit     )

  call atl_eqn_write( equation = equation, &
    &                 pconf    = out_conf  )

  rewind(funit)

  !*****************************************************************************
  write(logUnit(3),*) 'Reading output'
  file_str = ''
  do
    read(funit,'(A)', iostat=iError) line_str
    if (iError /= 0) exit
    file_str = trim(file_str) // adjustl(trim(line_str))
  end do
  close(funit)


  call open_config_chunk( L     = readconf, &
    &                     chunk = file_str  )
  call atl_init_equation(               &
    & equation     = equation_read,     &
    & conf         = readconf,          &
    & varSys_data  = methodData_read,   &
    & initSource   = initSource_read,   &
    & initMaterial = initMaterial_read  )
  !*****************************************************************************
  write(logUnit(3),*) 'Initializing variable system for reading'
  call tem_varSys_init( me         = equation_read%varSys, &
    &                   systemName = 'utest_accessState'   )
  call tem_variable_load( me     = userVars, &
    &                     conf   = varconf,  &
    &                     vError = vError    )

  call tem_varSys_append_luaVar( luaVar     = userVars,               &
    &                            varSys     = equation_read%varSys,   &
    &                            st_funList = equation_read%stFunList )

  call aot_table_open( L       = readconf,   &
    &                  thandle = tblHandle, &
    &                  key     = "equation" )
  if( tblHandle > 0 ) then
    ! Initilized keyValuePair dict for material
    call init( me = initMaterial_read%materialDict )

    call tem_variable_loadMapping(                      &
      & possVars = initMaterial_read%poss_materialVars, &
      & conf     = readconf,                            &
      & parent   = tblHandle,                           &
      & key      = "material",                          &
      & varSys   = equation_read%varSys,                &
      & varDict  = initMaterial_read%materialDict       )
    call aot_table_close( eqnconf, tblHandle )
  end if
  call atl_append_newMaterialVars(                        &
    & varSys       = equation_read%varSys,                &
    & varSys_data  = methodData,                          &
    & poss_matVars = initMaterial_read%poss_materialVars, &
    & materialFun  = equation_read%material,              &
    & variables    = initMaterial_read%materialDict       )
  write(logUnit(3),*) 'Number of materials: ', equation_read%material%nMat

  isEqual = isEqual .and. equation%eq_kind == equation_read%eq_kind
  isEqual = isEqual .and. &
    & (equation%euler%isen_coef .feq. equation_read%euler%isen_coef)
  isEqual = isEqual .and. (equation%euler%r .feq. equation_read%euler%r)
  isEqual = isEqual .and. (equation%euler%cv .feq. equation_read%euler%cv)
  isEqual = isEqual .and. &
    & (equation%euler%porosity .feq. equation_read%euler%porosity)
  isEqual = isEqual .and. &
    & (equation%euler%viscous_permeability .feq. &
    & equation_read%euler%viscous_permeability)
  isEqual = isEqual .and. &
    & (equation%euler%thermal_permeability .feq. &
    & equation_read%euler%thermal_permeability)

  isEqual = equation_read%material%nMat == equation%material%nMat
  do iMat = 1, initMaterial_read%materialDict%nVals
    strKvpA = initMaterial_read%materialDict%val(iMat)
    strKvpB = initMaterial%materialDict%val(iMat)
    isEqual = isEqual .and. strKvpA == strKvpB
  end do

  if (isEqual) then
    write(*,*) 'PASSED'
  end if

end program atl_dumpEquationTable_test
