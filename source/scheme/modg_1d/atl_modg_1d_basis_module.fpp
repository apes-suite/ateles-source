! Copyright (c) 2013 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
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

?? include "ply_dof_module.inc"
!> author: Jens Zudrop
!! Routines and datatypes related to the modal basis functions of the
!! modal discontinuous Galerkin scheme.
module atl_modg_1d_basis_module

  use env_module,            only: rk
  use ply_dof_module,        only: Q_space, &
    &                              P_space

  implicit none

  private

  public :: atl_evalLegendreTensPoly1d

contains


  subroutine atl_evalLegendreTensPoly1d( coords, nCoords, maxPolyDegree, &
    &                                    basisType, polyVal              )
    !> Array of coordinates (on the reference element) to evaluate the tensor
    !! product polynomials at. First dimension is nCoord, second is 3 for x,y,z
    !! component.
    real(kind=rk), intent(in) :: coords(:,:)

    !> The number of coordinates to evaluate the polynomials at.
    integer, intent(in) :: nCoords

    !> The maximum polynomial degree of the MODG scheme.
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: basisType

    !> The polynomial values. First dimension is the number of tensor product
    !! polynomials and the second dimension is the number of points, i.e.
    !! nCoords.
    real(kind=rk), allocatable, intent(out) :: polyVal(:,:)
    ! --------------------------------------------------------------------------
    integer :: iAns
    real(kind=rk) :: n_q
    ! --------------------------------------------------------------------------

    ! allocate the output array

    allocate( polyVal( (maxPolyDegree+1), nCoords) )

    ! Evaluate the Legendre polynomials per direction:
    ! ... first Legendere polynmoial is constant
    polyVal(1,:) = 1.0_rk
    if (basisType == Q_space .or. basisType == P_space) then
      if (maxPolyDegree > 0) then
        ! ... second Legendere polynmoial is identity
        polyVal(2,:) = coords(:,1)
        ! ... higher order polynomials are build recursively
        do iAns = 3, maxPolyDegree+1
          n_q = 1.0_rk / real(iAns-1,kind=rk)
          ! x recursion
          polyVal(iAns,:) = ( (2*(iAns-1)-1)*coords(:,1)*polyVal(iAns-1,:) &
            &                 - ((iAns-1)-1)*polyVal(iAns-2,:) )*n_q
        end do
      end if
    end if

  end subroutine atl_evalLegendreTensPoly1D


end module atl_modg_1d_basis_module
