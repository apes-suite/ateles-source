! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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

?? include "ply_dof_module.inc"
!> Module for routines and datatypes of Modal Discontinuous Galerkin (MODG)
!! scheme for the LinearEuler equation. This scheme is a spectral scheme for linear, purley hyperbolic
!! partial differential equation systems.
module atl_modg_LinearEuler_kernel_module
  use env_module,               only: rk

  use ply_poly_project_module,  only: ply_poly_project_type, assignment(=)
  use ply_dof_module,           only: Q_space, P_space

  use atl_equation_module,      only: atl_equations_type
  use atl_facedata_module,      only: atl_facedata_type
  use atl_cube_elem_module,     only: atl_cube_elem_type
  use atl_scheme_module,        only: atl_scheme_type, atl_modg_scheme_type
!  use atl_LinearEuler_numflux_module,  only: atl_LinearEuler_numflux
  use atl_penalization_module,  only: atl_penalizationData_type
  use atl_materialPrp_module,   only: atl_material_type

  implicit none
  private

  public :: atl_modg_LinearEuler_numflux, atl_modg_LinearEuler_physFlux


contains


  ! ****************************************************************************
  !> Calculate the physical flux for the MODG scheme and
  !! Linearized euler equation.
  subroutine atl_modg_LinearEuler_physFlux( equation, res, state, iElem, iDir, &
    &                                       penalizationData, poly_proj,       &
    &                                       material, nodal_data,              &
    &                                       nodal_gradData, nodal_res,         &
    &                                       elemLength,  scheme_min,           &
    &                                       scheme_current                     )
    ! --------------------------------------------------------------------------
    !> The equation system we are working with
    type(atl_equations_type), intent(in) :: equation
    !> The result in the modal form
    real(kind=rk), intent(inout)     :: res(:,:)
    !> The state in the modal form
    real(kind=rk), intent(in), optional :: state(:,:)
    !> The current element index
    integer, intent(in) :: iElem
    !> The current direction
    integer, intent(in) :: iDir
    !> The Penalization data
    type(atl_penalizationData_type), intent(inout) :: penalizationData
    !> The projection datatype for the projection information
    type(ply_poly_project_type), intent(inout) :: poly_proj
    !> The material information
    type(atl_material_type), intent(inout) :: material
    !> The state data in the nodal form
    real(kind=rk), intent(in), optional :: nodal_data(:,:)
    real(kind=rk), intent(in), optional :: nodal_GradData(:,:,:)
    !> The result in the nodal form
    real(kind=rk), intent(inout)     :: nodal_res(:,:)
    !> The length of the current element
    real(kind=rk), intent(in) :: ElemLength
    !> The scheme information of the min level (This is needed for the temp
    ! buffer array for evaluating the physical fluxes )
    type(atl_scheme_type), intent(inout) :: scheme_min
    !> Information about the current level
    type(atl_scheme_type), intent(inout) :: scheme_current
    ! -------------------------------------------------------------------------!
    real(kind=rk) :: v0, d0, p0, inv_d0
    integer :: iDof
    integer :: nDofs
    ! Rotation indices for physical flux calculation
    integer :: rot1, rot2, rot3, rot4, rot5
    real(kind=rk) :: tmp1, tmp2
    ! -------------------------------------------------------------------------!

    ! get the rotation for the physical flux calculation
    rot1 = equation%varRotation(iDir)%varTransformIndices(1)
    rot2 = equation%varRotation(iDir)%varTransformIndices(2)
    rot3 = equation%varRotation(iDir)%varTransformIndices(3)
    rot4 = equation%varRotation(iDir)%varTransformIndices(4)
    rot5 = equation%varRotation(iDir)%varTransformIndices(5)

    v0 = equation%LinearEuler%velocity_0(iDir)
    d0 = equation%LinearEuler%density_0
    p0 = equation%LinearEuler%pressure_0 * equation%LinearEuler%isen_coef
    inv_d0 = 1.0_rk / equation%LinearEuler%density_0
    nDofs = poly_proj%body_3d%ndofs

    do iDof = 1, nDofs
      tmp1 = state(iDof,rot2)
      tmp2 = state(iDof,rot5)
      res( iDof,rot1 ) = v0 * state(iDof,rot1) + d0 * tmp1
      res( iDof,rot2 ) = v0 * tmp1 + tmp2 * inv_d0
      res( iDof,rot3 ) = v0 * state(iDof,rot3)
      res( iDof,rot4 ) = v0 * state(iDof,rot4)
      res( iDof,rot5 ) = v0 * tmp2 + p0 * tmp1
    end do

  end subroutine atl_modg_LinearEuler_physFlux
  ! ****************************************************************************


  ! ****************************************************************************
  !> Calculate the numerical flux for LinearEuler equation and MODG scheme
  subroutine atl_modg_LinearEuler_numFlux( mesh, equation, facedata, scheme )
    ! --------------------------------------------------------------------------
    !> The mesh you are working with.
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation you solve.
    type(atl_equations_type), intent(in) :: equation
    !> The face representation of the state.
    type(atl_facedata_type), intent(inout) :: facedata
    !> Parameters of the modal dg scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    ! --------------------------------------------------------------------------
    integer :: iDir, nFaceDofs
    ! --------------------------------------------------------------------------

    ! Numerical flux for faces in all 3 spatial face directions
    select case(scheme%basisType)
      case(Q_space)
        nFaceDofs = (scheme%maxPolyDegree+1)**2
      case(P_space)
?? copy :: getDofsPTens2D(scheme%maxPolyDegree, nFaceDofs)
    end select


    ! Calculate the numerical fluxes for the faces in all 3 spatial face
    ! directions
    do iDir = 1,3
      call equation%LinearEuler%dir_proc(iDir)%numFlux(               &
        &  nSides = size(mesh%faces%faces(iDir)%computeFace%leftPos), &
        &  nFaceDofs = nFaceDofs,                                     &
        &  faceRep = facedata%faceRep(iDir)%dat,                      &
        &  faceFlux = facedata%faceFlux(iDir)%dat,                    &
        &  leftPos = mesh%faces%faces(iDir)%computeFace%leftPos,      &
        &  rightPos = mesh%faces%faces(iDir)%computeFace%rightPos,    &
        &  var = equation%varRotation(iDir)%varTransformIndices(1:5), &
        &  LinearEuler = equation%LinearEuler,                        &
        &  iDir = iDir                                                )
    end do


  end subroutine atl_modg_LinearEuler_numFlux
  ! ****************************************************************************

end module atl_modg_LinearEuler_kernel_module
