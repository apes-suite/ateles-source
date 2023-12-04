! Copyright (c) 2013 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
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

?? include 'treelm/source/deriveMacros.inc'

!> Routines to derive quantities from the state in the Maxwell equation system.
module atl_eqn_maxwell_derive_module
  use, intrinsic :: iso_c_binding,  only: c_f_pointer
  use env_module,                   only: rk

  use tem_varSys_module,            only: tem_varSys_type,            &
    &                                     tem_varSys_op_type
  use tem_time_module,              only: tem_time_type
  use treelmesh_module,             only: treelmesh_type

  use atl_equation_module,          only: atl_equations_type

  implicit none

  private

  public :: atl_xFrom3D_getElement, atl_xFrom3D_getPoint
  public :: atl_yFrom3D_getElement, atl_yFrom3D_getPoint
  public :: atl_zFrom3D_getElement, atl_zFrom3D_getPoint
  public :: atl_eqn_maxwell_cons2prim
  public :: atl_eqn_maxwell_prim2cons

contains

  !> Convert primitive varibales to conservative variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_maxwell_prim2cons(equation, instate, outstate, material )
    ! ------------------------------------------------------------------------ !
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout)         :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional, intent(in) :: material(:,:)
    ! ------------------------------------------------------------------------ !

    if( size(material,1).eq.1 ) then
      if(present(outstate)) then
        ! displacement field, i.e. D from
        outstate(:,1) = instate(:,1) * material(1,1)
        outstate(:,2) = instate(:,2) * material(1,1)
        outstate(:,3) = instate(:,3) * material(1,1)

        ! magnetic field, i.e. B from H
        outstate(:,4) = instate(:,4) * material(1,2)
        outstate(:,5) = instate(:,5) * material(1,2)
        outstate(:,6) = instate(:,6) * material(1,2)
      else
        ! displacement field, i.e. D from
        instate(:,1) = instate(:,1) * material(1,1)
        instate(:,2) = instate(:,2) * material(1,1)
        instate(:,3) = instate(:,3) * material(1,1)

        ! magnetic field, i.e. B from H
        instate(:,4) = instate(:,4) * material(1,2)
        instate(:,5) = instate(:,5) * material(1,2)
        instate(:,6) = instate(:,6) * material(1,2)
      end if
    else
      if(present(outstate)) then
        ! displacement field, i.e. D from
        outstate(:,1) = instate(:,1) * material(:,1)
        outstate(:,2) = instate(:,2) * material(:,1)
        outstate(:,3) = instate(:,3) * material(:,1)

        ! magnetic field, i.e. B from H
        outstate(:,4) = instate(:,4) * material(:,2)
        outstate(:,5) = instate(:,5) * material(:,2)
        outstate(:,6) = instate(:,6) * material(:,2)
      else
        ! displacement field, i.e. D from
        instate(:,1) = instate(:,1) * material(:,1)
        instate(:,2) = instate(:,2) * material(:,1)
        instate(:,3) = instate(:,3) * material(:,1)

        ! magnetic field, i.e. B from H
        instate(:,4) = instate(:,4) * material(:,2)
        instate(:,5) = instate(:,5) * material(:,2)
        instate(:,6) = instate(:,6) * material(:,2)
      end if
    end if

  end subroutine atl_eqn_maxwell_prim2cons


  !> Convert conservative to primitive variables.
  !!
  !! The interface has to comply to the abstract interface
  !! atl_equation_module#eqn_var_trafo "eqn_var_trafo".
  subroutine atl_eqn_maxwell_cons2prim(equation, instate, outstate, material )
    ! --------------------------------------------------------------------------!
    !> Description of the equation system.
    class(atl_equations_type), intent(in) :: equation

    !> Primitive variables to convert.
    real(kind=rk), intent(inout) :: instate(:,:)

    !> Converted variables.
    real(kind=rk), optional, intent(out) :: outstate(:,:)

    !> The material information.
    real(kind=rk), optional,  intent(in) :: material(:,:)
    ! --------------------------------------------------------------------------!

    if(size(material,1).eq.1) then
      if(present(outstate)) then
        ! electric field, i.e. E from D
        outstate(:,1) = instate(:,1) / material(1,1)
        outstate(:,2) = instate(:,2) / material(1,1)
        outstate(:,3) = instate(:,3) / material(1,1)

        ! magnetizing field, i.e. H from B
        outstate(:,4) = instate(:,4) / material(1,2)
        outstate(:,5) = instate(:,5) / material(1,2)
        outstate(:,6) = instate(:,6) / material(1,2)
      else
        ! electric field, i.e. E from D
        instate(:,1) = instate(:,1) / material(1,1)
        instate(:,2) = instate(:,2) / material(1,1)
        instate(:,3) = instate(:,3) / material(1,1)

        ! magnetizing field, i.e. H from B
        instate(:,4) = instate(:,4) / material(1,2)
        instate(:,5) = instate(:,5) / material(1,2)
        instate(:,6) = instate(:,6) / material(1,2)
      end if
    else
      if(present(outstate)) then
        ! electric field, i.e. E from D
        outstate(:,1) = instate(:,1) / material(:,1)
        outstate(:,2) = instate(:,2) / material(:,1)
        outstate(:,3) = instate(:,3) / material(:,1)

        ! magnetizing field, i.e. H from B
        outstate(:,4) = instate(:,4) / material(:,2)
        outstate(:,5) = instate(:,5) / material(:,2)
        outstate(:,6) = instate(:,6) / material(:,2)
      else
        ! electric field, i.e. E from D
        instate(:,1) = instate(:,1) / material(:,1)
        instate(:,2) = instate(:,2) / material(:,1)
        instate(:,3) = instate(:,3) / material(:,1)

        ! magnetizing field, i.e. H from B
        instate(:,4) = instate(:,4) / material(:,2)
        instate(:,5) = instate(:,5) / material(:,2)
        instate(:,6) = instate(:,6) / material(:,2)
      end if
    end if
  end subroutine atl_eqn_maxwell_cons2prim


  subroutine atl_xFrom3D_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nPnts*3)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = input                                    )

    do iPoint = 1, nPnts
      res?IDXPNT?(1,iPoint,1) = input?IDXPNT?(1,iPoint,3)
    end do

  end subroutine atl_xFrom3D_getPoint

  subroutine atl_xFrom3D_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                     nDofs, res                          )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nDofs*nElems*3)
    integer :: iElem, iDof
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = input                                      )

    !! We only want to have one component, thus we don't need to loop over them
    do iElem = 1, nElems
      do iDof = 1, nDofs
        res?IDXELEM?(1, iDof, iElem, 1, nDofs) =    &
          & input?IDXELEM?(1, iDof, iElem, 3, nDofs)
      end do
    end do

  end subroutine atl_xFrom3D_getElement

  subroutine atl_yFrom3D_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nPnts*3)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = input                                    )

    do iPoint = 1, nPnts
      res?IDXPNT?(1,iPoint,1) = input?IDXPNT?(2,iPoint,3)
    end do

  end subroutine atl_yFrom3D_getPoint

  subroutine atl_yFrom3D_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                     nDofs, res                          )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nDofs*nElems*3)
    integer :: iElem, iDof
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = input                                      )

    !! We only want to have one component, thus we don't need to loop over them
    do iElem = 1, nElems
      do iDof = 1, nDofs
        res?IDXELEM?(1, iDof, iElem, 1, nDofs) =    &
          & input?IDXELEM?(2, iDof, iElem, 3, nDofs)
      end do
    end do

  end subroutine atl_yFrom3D_getElement

  subroutine atl_zFrom3D_getPoint(fun, varsys, point, time,tree, nPnts, res )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the point time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> Three-dimensional coordinates at which the variable should be
    !! evaluated. Only useful for variables provided as space-time functions.
    real(kind=rk), intent(in) :: point(:,:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nPnts

    !> Resulting values for the requested variable.
    !!
    !! Dimension: n requested entries x nComponents of this variable
    !! Access: (iElem-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nPnts*3)
    integer :: iPoint
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_point( &
      & varSys  = varSys,                                  &
      & point   = point,                                   &
      & time    = time,                                    &
      & tree    = tree,                                    &
      & nPnts   = nPnts,                                   &
      & res     = input                                    )

    do iPoint = 1, nPnts
      res?IDXPNT?(1,iPoint,1) = input?IDXPNT?(3,iPoint,3)
    end do

  end subroutine atl_zFrom3D_getPoint

  subroutine atl_zFrom3D_getElement(fun, varsys, elempos, time, tree, nElems, &
    &                                     nDofs, res                          )
    ! ---------------------------------------------------------------------------
    !> Description of the method to obtain the variables, here some preset
    !! values might be stored, like the space time function to use or the
    !! required variables.
    class(tem_varSys_op_type), intent(in) :: fun

    !> The variable system to obtain the variable from.
    type(tem_varSys_type), intent(in) :: varSys

    !> TreeID of the element to get the variable for.
    integer, intent(in) :: elempos(:)

    !> Point in time at which to evaluate the variable.
    type(tem_time_type), intent(in)  :: time

    !> global treelm mesh info
    type(treelmesh_type), intent(in) :: tree

    !> Number of values to obtain for this variable (vectorized access).
    integer, intent(in) :: nElems

    !> Number of degrees of freedom within an element.
    integer, intent(in) :: nDofs

    !> Resulting values for the requested variable.
    !!
    !! Linearized array dimension:
    !! (n requested entries) x (nComponents of this variable)
    !! x (nDegrees of freedom)
    !! Access: (iElem-1)*fun%nComponents*nDofs +
    !!         (iDof-1)*fun%nComponents + iComp
    real(kind=rk), intent(out) :: res(:)
    ! ---------------------------------------------------------------------------
    real(kind=rk) :: input(nDofs*nElems*3)
    integer :: iElem, iDof
    ! ---------------------------------------------------------------------------

    call varSys%method%val(fun%input_varPos(1))%get_element( &
      & varSys  = varSys,                                    &
      & elemPos = elemPos,                                   &
      & time    = time,                                      &
      & tree    = tree,                                      &
      & nElems  = nElems,                                    &
      & nDofs   = nDofs,                                     &
      & res     = input                                      )

    !! We only want to have one component, thus we don't need to loop over them
    do iElem = 1, nElems
      do iDof = 1, nDofs
        res?IDXELEM?(1, iDof, iElem, 1, nDofs) =    &
          & input?IDXELEM?(3, iDof, iElem, 3, nDofs)
      end do
    end do

  end subroutine atl_zFrom3D_getElement

end module atl_eqn_maxwell_derive_module
