! Copyright (c) 2013 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016-2017, 2019-2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2014 Rishabh Chandola <rishabh.chandola@student.uni-siegen.de>
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
!! \author{Jens Zudrop}
!> Routines and datatypes related to the modal basis functions of the
!! modal discontinuous Galerkin scheme.
module atl_modg_2d_basis_module

  use env_module,                  only: rk
  use ply_dof_module,              only: Q_space, P_space

  implicit none
  private

  public :: atl_evalLegendreTensPoly2d


contains


  subroutine atl_evalLegendreTensPoly2d( coords, nCoords, maxPolyDegree, &
    &                                    basisType, polyVal              )
    !> Array of coordinates (on the reference element) to evaluate the tensor
    !! product polynomials at. First dimension is nCoord, second is 2 for x,y
    !! component.
    real(kind=rk), intent(in) :: coords(:,:)

    !> The number of coordinates to evaluate the polynomials at.
    integer, intent(in) :: nCoords

    !> The maximum polynomail degree of the MODG scheme.
    integer, intent(in) :: maxPolyDegree
    integer, intent(in) :: basisType

    !> The polynomial values. First dimension is the number of tensor product
    !! polynomials and the second dimension is the number of points, i.e.
    !! nCoords.
    real(kind=rk), allocatable, intent(out) :: polyVal(:,:)
    ! ---------------------------------------------------------------------------
    real(kind=rk), allocatable :: polyValX(:,:), polyValY(:,:)
    integer :: iAnsX, iAnsY, iAns, ansPos, ansPosMax
    real(kind=rk) :: n_q
    ! ---------------------------------------------------------------------------

    ! allocate the output array
    select case(basisType)
      case(Q_space)
        allocate( polyVal( (maxPolyDegree+1)**2 ,nCoords) )
      case(P_space)
        allocate( polyVal((maxPolydegree+1)*(maxPolydegree+2)/2, nCoords ) )
    end select

    allocate( polyValX( (maxPolyDegree+1) ,nCoords) )
    allocate( polyValY( (maxPolyDegree+1) ,nCoords) )
!    allocate( polyValZ( (maxPolyDegree+1) ,nCoords) )

    ! Evaluate the Legendre polynomials per direction:
    ! ... first Legendere polynmoial is constant
    polyValX(1,:) = 1.0_rk
    polyValY(1,:) = 1.0_rk
    if(maxPolyDegree > 0) then
      ! ... second Legendere polynmoial is identity
      polyValX(2,:) = coords(:,1)
      polyValY(2,:) = coords(:,2)
      ! ... higher order polynomials are build recursively
      do iAns = 3, maxPolyDegree+1
        n_q = 1.0_rk / real(iAns-1,kind=rk)
        ! x recursion
        polyValX(iAns,:) = ( (2*(iAns-1)-1)*coords(:,1)*polyValX(iAns-1,:) &
          &                 - ((iAns-1)-1)*polyValX(iAns-2,:) )*n_q
        ! y recursion
        polyValY(iAns,:) = ( (2*(iAns-1)-1)*coords(:,2)*polyValY(iAns-1,:) &
          &                 - ((iAns-1)-1)*polyValY(iAns-2,:) )*n_q
      end do
    end if

    ! Now, build the complete point value.

    select case(basisType)
      case(Q_space)
        do iAnsX = 1, maxPolyDegree+1
          do iAnsY = 1, maxPolyDegree+1
              ! get the position of this ansatz function combination.
?? copy :: posOfModgCoeffQTens2D(iAnsX, iAnsY, maxPolyDegree, ansPos)
              polyVal(ansPos, :) = polyValX(iAnsX,:) * polyValY(iAnsY,:)
          end do
        end do
      case(P_space)
        iAnsX = 1
        iAnsY = 1
?? copy :: getDofsPTens2D(maxPolyDegree, ansPosMax)
        do ansPos = 1, ansPosMax
          polyVal(ansPos, :) = polyValX(iAnsX,:) * polyValY(iAnsY,:)
?? copy :: nextModgCoeffPTens2D(iAnsX, iAnsY)
        end do
    end select


  end subroutine atl_evalLegendreTensPoly2D


end module atl_modg_2d_basis_module
