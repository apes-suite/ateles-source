! Copyright (c) 2015, 2018 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Girresser <tobias.girresser@student.uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
module atl_project_physflux_module
  use env_module,                only: rk
  ! ateles
  use atl_cube_elem_module,      only: atl_cube_elem_type
  use atl_equation_module,       only: atl_equations_type
  use atl_kerneldata_module,     only: atl_kerneldata_type
  use atl_scheme_module,         only: atl_modg_scheme_type
  use atl_modg_2d_scheme_module, only: atl_modg_2d_scheme_type
  ! polynomials
  use ply_dof_module,            only: Q_space, P_space
  use ply_modg_basis_module,     only: ply_scalProdDualLegDiff

  implicit none

  private

  public :: atl_modg_project_physflux_testfunc, &
    &       atl_modg_2d_project_physFlux_testfunc


contains


  !> Subroutine to project modal representations of physical flux, numerical flux
  !! and source terms onto test functions.
  subroutine atl_modg_project_PhysFlux_testFunc( mesh, equation, kerneldata, &
    &                                           scheme, iDir, dl_prod,       &
    &                                           dirVec, iElem, state_der     )
    ! --------------------------------------------------------------------------
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The direction
    integer, intent(in) :: iDir
    !> The parameters of the MODG scheme
    type(atl_modg_scheme_type), intent(in) :: scheme
    !> stored scalar products of the testfunction and ansatz function
    real(kind=rk), intent(in) :: dl_prod(2, scheme%maxPolyDegree+1)
    !> vector for direction indicators
    integer, intent(in) :: dirVec(3,3)
    integer, intent(in) :: iElem
    real(kind=rk), intent(in)  :: state_der(kerneldata%nDofs,equation%varSys%nScalars)
    ! --------------------------------------------------------------------------

    ! Projection of the physical flux
      select case(scheme%basisType)
      case(Q_space)

        select case(equation%varSys%nScalars)
        case(6)
          select case(iDir)
          case(1)
            call modg_prj_pFlux1_Q_6( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(2)
            call modg_prj_pFlux2_Q_6( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(3)
            call modg_prj_pFlux3_Q_6( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          end select
        case(5)
          select case(iDir)
          case(1)
            call modg_prj_pFlux1_Q_5( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(2)
            call modg_prj_pFlux2_Q_5( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          case(3)
            call modg_prj_pFlux3_Q_5( maxPolyDegree = scheme%maxPolyDegree, &
              &                       length = mesh%length, &
              &                       dl_prod = dl_prod, &
              &                       state = kerneldata%state_der, &
              &                       iElem  = iElem,              &
              &                       state_der = state_der        )
          end select
        case default
          call modg_project_physFlux_Q( nScalars = equation%varSys%nScalars, &
            &                           maxPolyDegree = scheme%maxPolyDegree, &
            &                           length = mesh%length, &
            &                           dl_prod = dl_prod, &
            &                           state = kerneldata%state_der, &
            &                           dirVec = dirVec(:,iDir),     &
            &                           iElem  = iElem,              &
            &                           state_der = state_der        )
        end select

      case(P_space)
        call modg_project_physFlux_P( nScalars = equation%varSys%nScalars, &
          &                           maxPolyDegree = scheme%maxPolyDegree, &
          &                           length = mesh%length, &
          &                           dl_prod = dl_prod, &
          &                           iElem = iElem, &
          &                           state = kerneldata%state_der, &
          &                           dirVec = dirVec(:,iDir), &
          &                           state_der = state_der                  )
      end select


  end subroutine atl_modg_project_PhysFlux_testFunc


  !> Projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_project_physFlux_Q( nScalars, maxPolyDegree, length, state, &
    &                                 dl_prod, dirVec, iElem, state_der )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> ordering of xyz for current direction
    integer, intent(in) :: dirVec(3)
    integer, intent(in) :: iElem
    !> The state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj, scalProd1(maxPolyDegree)
    integer :: testPos
    integer :: iTest2, iTest1, iTest3, iTestVec(3), min2mpd
    integer :: iAnsVec(3)
    integer :: iVar
    integer :: ansPos(4)
    real(kind=rk) :: scalProd(4)
    integer :: jk
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**2

    do iTest1=2,maxPolyDegree+1
      scalProd1(iTest1-1) = ply_scalProdDualLegDiff(iTest1-1, iTest1) &
        &                   * jacobiDetStiffProj
    end do

    min2mpd = min(2, maxPolyDegree+1)

    ! unrolled loop

    do iTest3 = 1, min2mpd

      do jk=1,min2mpd*maxPolyDegree
        iTest2 = mod(jk-1,min2mpd) + 1
        iTest1 = (jk-1)/min2mpd + 2

        ! one entry

        iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens( iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), maxPolyDegree, testPos )

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(1))
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(ansPos(1),iVar) * scalProd(1)
        end do

      end do

      do jk=1,(maxPolyDegree-1)*maxPolyDegree
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 3

        ! two entries

        iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens( iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), maxPolyDegree, testPos )


        iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(1))
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(1,iTest2) &
          &           * dl_prod(2,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(2))
        scalProd(2) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(ansPos(1),iVar) * scalProd(1) &
            & + state_der(ansPos(2),iVar) * scalProd(2)
        end do

      end do

    end do

    do iTest3 = 3, maxPolyDegree+1

      do jk=1,min2mpd*maxPolyDegree
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 1

        ! two entries

        iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens( iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), maxPolyDegree, testPos )


        iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(1))
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(1,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(2))
        scalProd(2) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar) &
            & + state_der(ansPos(1),iVar) * scalProd(1) &
            & + state_der(ansPos(2),iVar) * scalProd(2)
        end do
      end do

      do jk=1,maxPolyDegree*(maxPolyDegree-1)
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 3

        ! four entries

        iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens(iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), maxPolyDegree, testPos)


        iAnsVec = (/iTest1-1, iTest2-2, iTest3-2/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(1))
        scalProd(1) =   scalProd1(iTest1-1) &
          &           * dl_prod(1,iTest2) &
          &           * dl_prod(1,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(2))
        scalProd(2) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(1,iTest3)

        iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(3))
        scalProd(3) =   scalProd1(iTest1-1) &
          &           * dl_prod(1,iTest2) &
          &           * dl_prod(2,iTest3)

        iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffQTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), maxPolyDegree, ansPos(4))
        scalProd(4) =   scalProd1(iTest1-1) &
          &           * dl_prod(2,iTest2) &
          &           * dl_prod(2,iTest3)

        do iVar=1,nScalars
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(ansPos(1),iVar) * scalProd(1) &
            & + state_der(ansPos(2),iVar) * scalProd(2) &
            & + state_der(ansPos(3),iVar) * scalProd(3) &
            & + state_der(ansPos(4),iVar) * scalProd(4)
        end do

      end do

    end do


  end subroutine modg_project_physFlux_Q


  !> Projection of the physical flux onto the testfunctions, with unrolled loops
  !! => fewer loop-overhead/instructions, but more "random" memory accesses
  !! MZ: perhaps this version is faster for low order (or always, depending on the machine?)
  subroutine modg_project_physFlux_P( nScalars, iElem, maxPolyDegree, length, &
    &                                 state, dl_prod, dirVec, state_der       )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> ordering of xyz for current direction
    integer, intent(in) :: dirVec(3)
    !> The state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    !> The element index
    integer, intent(in) :: iElem
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj, scalProd1(maxPolyDegree)
    integer :: testPos
    integer :: iTest2, iTest1, iTest3, iTestVec(3), min2mpd
    integer :: iAnsVec(3)
    integer :: iVar
    integer :: ansPos(4)
    real(kind=rk) :: scalProd(4)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! This is the stiffness term!
    !
    ! We have cubiic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**2

    do iTest1=2,maxPolyDegree+1
      scalProd1(iTest1-1) = ply_scalProdDualLegDiff(iTest1-1, iTest1) &
        &                   * jacobiDetStiffProj
    end do

    min2mpd = min(2, maxPolyDegree+1)

    ! unrolled loop

    do iTest3 = 1, min2mpd
      do iTest2 = 1, min(2, maxPolyDegree+1 - (iTest3-1))
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! one entry

          iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens( iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), testPos)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(1))
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1)
          end do

        end do
      end do

      do iTest2 = 3, maxPolyDegree+1 - (iTest3-1)
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! two entries

          iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens( iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), testPos )


          iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(1) )
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(1,iTest2) &
            &           * dl_prod(2,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(2) )
          scalProd(2) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1) &
              & + state_der(ansPos(2),iVar) * scalProd(2)
          end do

        end do
      end do
    end do

    do iTest3 = 3, maxPolyDegree+1
      do iTest2 = 1, min(2, maxPolyDegree+1 - (iTest3-1))
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! two entries

          iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens( iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), testPos )


          iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(1) )
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(1,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(2) )
          scalProd(2) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1) &
              & + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      end do

      do iTest2 = 3, maxPolyDegree+1 - (iTest3-1)
        do iTest1 = 2, maxPolyDegree+1 - (iTest3-1) - (iTest2-1)

          ! four entries

          iTestVec = (/iTest1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens( iTestVec(dirVec(1)), iTestVec(dirVec(2)), iTestVec(dirVec(3)), testPos )


          iAnsVec = (/iTest1-1, iTest2-2, iTest3-2/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(1) )
          scalProd(1) =   scalProd1(iTest1-1) &
            &           * dl_prod(1,iTest2) &
            &           * dl_prod(1,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3-2/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(2) )
          scalProd(2) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(1,iTest3)

          iAnsVec = (/iTest1-1, iTest2-2, iTest3/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(3) )
          scalProd(3) =   scalProd1(iTest1-1) &
            &           * dl_prod(1,iTest2) &
            &           * dl_prod(2,iTest3)

          iAnsVec = (/iTest1-1, iTest2, iTest3/)
?? copy :: posOfModgCoeffPTens(iAnsVec(dirVec(1)), iAnsVec(dirVec(2)), iAnsVec(dirVec(3)), ansPos(4) )
          scalProd(4) =   scalProd1(iTest1-1) &
            &           * dl_prod(2,iTest2) &
            &           * dl_prod(2,iTest3)

          do iVar=1,nScalars
            state(iElem,testPos,iVar) &
              & = state(iElem,testPos,iVar) &
              & + state_der(ansPos(1),iVar) * scalProd(1) &
              & + state_der(ansPos(2),iVar) * scalProd(2) &
              & + state_der(ansPos(3),iVar) * scalProd(3) &
              & + state_der(ansPos(4),iVar) * scalProd(4)
          end do

        end do
      end do
    end do

  end subroutine modg_project_physFlux_P


  !> Subroutine to project modal representations of physical flux, numerical 
  !! flux and source terms onto test functions.
  subroutine atl_modg_2d_project_physFlux_testFunc( mesh, equation,        &
    &                                               kerneldata, iElem,     &
    &                                               dl_prod, iDir, scheme, &
    &                                               state_data, ndofs      )
    ! --------------------------------------------------------------------------
    !> Descritption of the cubical elements in the mesh
    type(atl_cube_elem_type), intent(in) :: mesh
    !> The equation description.
    type(atl_equations_type), intent(in) :: equation
    !> The data of the kernel. Holds the physical fluxes.
    type(atl_kerneldata_type), intent(inout) :: kerneldata
    !> The parameters of the MODG scheme
    type(atl_modg_2d_scheme_type), intent(in) :: scheme
    !> The total degrees of freedom
    integer, intent(in) :: nDofs
    !> The direction
    integer, intent(in) :: iDir
    !> The element index
    integer, intent(in) :: iElem
    real(kind=rk) , intent(in) :: dl_prod(2, scheme%maxPolyDegree+1)
    !> The physical fluxes that needs to be projected
    real(kind=rk), intent(in)  :: state_data(nDofs,equation%varSys%nScalars)
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Iterate over all the elements and do the following:
    ! 1. Project physical fluxes (3 directions) to test functions (stiffness
    !    terms).
    !    Attention: physical x flux: state_der(:,:,2,:)
    !    Attention: physical y flux: state_der(:,:,3,:)
    ! 2. Project numerical fluxes (4 faces) to test functions.
    ! 3. Project source terms to test functions.
    ! Attention: write the projections to state_der(:,:,1,:) because inverse
    ! of mass matrix will be applied to these entries.

    select case(scheme%basisType)
    case(Q_space)
      if (iDir == 1) then
        ! Projection of the physical flux
        ! ... x direction
        call modg_2d_project_physFluxX_Q(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & state_der     = state_data                )
      else
        ! ... y direction
        call modg_2d_project_physFluxY_Q(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & state_der     = state_data                )
      endif

    case(P_space)
      if (iDir == 1) then
        ! Projection of the physical flux
        ! ... x direction
        call modg_2d_project_physFluxX_P(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & nDofs         = nDofs,                    &
          & state_der     = state_data                )
      else
        ! ... y direction
        call modg_2d_project_physFluxY_P(             &
          & nScalars      = equation%varSys%nScalars, &
          & maxPolyDegree = scheme%maxPolyDegree,     &
          & length        = mesh%length,              &
          & dl_prod       = dl_prod,                  &
          & iElem         = iElem,                    &
          & state         = kerneldata%state_der,     &
          & nDofs         = nDofs,                    &
          & state_der     = state_data                )
      endif
    end select

  end subroutine atl_modg_2d_project_physFlux_testfunc


  !> Projection of the physical flux in x direction onto the testfunctions.
  subroutine modg_2d_project_physFluxX_Q( nScalars, maxPolyDegree, length, &
    &                                     dl_prod, state, iElem, state_der )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der((maxPolyDegree+1)**2,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdX(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestY
    integer :: iAnsX
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test
    ! functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    ! for x direction (x test function differentiated)
    ! get the relevant indices for the ansatz function
    do iAnsX=1,maxPolyDegree
      scalProdX(iAnsX) = ply_scalProdDualLegDiff(iAnsX, iAnsX+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    YLoop: do iTestY = 1, maxPolyDegree+1

?? copy :: posOfModgCoeffQTens(1, iTestY, 1, maxPolyDegree, testPos)
      do iVar=var_lb,var_ub
        state(iElem,testPos,iVar) = 0.0_rk
      end do

      if (iTestY > 2) then
        ! Need to add two terms
        do iAnsX = 1, maxPolyDegree
          ! the position of the current test functions
?? copy :: posOfModgCoeffQTens(iAnsX+1, iTestY,1, maxPolyDegree, testPos)
          scalProd(1) = dl_prod(1, iTestY) * scalProdX(iAnsX)
          scalProd(2) = dl_prod(2, iTestY) * scalProdX(iAnsX)
?? copy :: posOfModgCoeffQTens(iAnsX, iTestY-2 , 1, maxPolyDegree, ansPos(1))
?? copy :: posOfModgCoeffQTens(iAnsX, iTestY, 1, maxPolyDegree, ansPos(2))
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsX = 1, maxPolyDegree
          ! the position of the current test functions
?? copy :: posOfModgCoeffQTens(iAnsX+1, iTestY, 1, maxPolyDegree, testPos)
          scalProd(1) = dl_prod(2, iTestY) * scalProdX(iAnsX)
?? copy :: posOfModgCoeffQTens(iAnsX, iTestY, 1, maxPolyDegree, ansPos(1))
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                     &
              & = state(iElem,testPos,iVar)               &
              &   + state_der(ansPos(1),iVar) * scalProd(1)
          end do

        end do
      end if
    end do YLoop

  end subroutine modg_2d_project_physFluxX_Q


  !> Projection of the physical flux in y direction onto the testfunctions.
  subroutine modg_2d_project_physFluxY_Q( nScalars, maxPolyDegree, length,  &
    &                                     dl_prod, state , ielem, state_der )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der((maxPolyDegree+1)**2,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdY(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestX
    integer :: iAnsY
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test
    ! functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    !  for y direction (y test function differentiated)
    do iAnsY=1,maxPolyDegree
      scalProdY(iAnsY) = ply_scalProdDualLegDiff(iAnsY, iAnsY+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    XLoop: do iTestX = 1, maxPolyDegree+1

      if (iTestX > 2) then
        ! Need to add two terms
        do iAnsY = 1, maxPolyDegree
?? copy :: posOfModgCoeffQTens(iTestX, iAnsY+1, 1, maxPolyDegree, testPos)
          scalProd(1) = dl_prod(1, iTestX) * scalProdY(iAnsY)
          scalProd(2) = dl_prod(2, iTestX) * scalProdY(iAnsY)
?? copy :: posOfModgCoeffQTens(iTestX-2, iAnsY, 1, maxPolyDegree, ansPos(1))
?? copy :: posOfModgCoeffQTens(iTestX, iAnsY, 1, maxPolyDegree, ansPos(2))
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsY = 1, maxPolyDegree
?? copy :: posOfModgCoeffQTens(iTestX, iAnsY+1, 1, maxPolyDegree, testPos)
          scalProd(1) = dl_prod(2, iTestX) * scalProdY(iAnsY)
?? copy :: posOfModgCoeffQTens(iTestX, iAnsY, 1, maxPolyDegree, ansPos(1))
          do iVar = var_lb,var_ub
            state(iElem,testPos,iVar)                    &
              &  = state(iElem,testPos,iVar)             &
              &  + state_der(ansPos(1),iVar) * scalProd(1)
          end do
        end do
      end if
    end do XLoop

  end subroutine modg_2d_project_physFluxY_Q


  !> Projection of the physical flux in x direction onto the testfunctions.
  subroutine modg_2d_project_physFluxX_P( nScalars, maxPolyDegree, length, &
    &                                     dl_prod, state, iElem, nDofs,    &
    &                                     state_der                        )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> Number of degrees of freedom
    integer, intent(in) :: nDofs
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der(nDofs,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdX(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestY
    integer :: iAnsX
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! This is the stiffness term!
    !
    ! We have cubiic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.

    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    ! for x direction (x test function differentiated)
    ! get the relevant indices for the ansatz function
    do iAnsX=1,maxPolyDegree
      scalProdX(iAnsX) = ply_scalProdDualLegDiff(iAnsX, iAnsX+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    YLoop: do iTestY = 1, maxPolyDegree+1

?? copy :: posOfModgCoeffPTens2D(1, iTestY, testPos)
      do iVar=var_lb,var_ub
        state(iElem,testPos,iVar) = 0.0_rk
      end do

      if (iTestY > 2) then
        ! Need to add two terms
        do iAnsX = 1, maxPolyDegree-(iTestY-1)
          ! the position of the current test functions
?? copy :: posOfModgCoeffPTens2D(iAnsX+1, iTestY, testPos)
          scalProd(1) = dl_prod(1, iTestY) * scalProdX(iAnsX)
          scalProd(2) = dl_prod(2, iTestY) * scalProdX(iAnsX)
?? copy :: posOfModgCoeffPTens2D(iAnsX, iTestY-2, ansPos(1))
?? copy :: posOfModgCoeffPTens2D(iAnsX, iTestY, ansPos(2))
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsX = 1, maxPolyDegree - (iTestY-1)
          ! the position of the current test functions
?? copy :: posOfModgCoeffPTens2D(iAnsX+1, iTestY, testPos)
          scalProd(1) = dl_prod(2, iTestY) * scalProdX(iAnsX)
?? copy :: posOfModgCoeffPTens2D(iAnsX, iTestY, ansPos(1))
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                     &
              & = state(iElem,testPos,iVar)               &
              &   + state_der(ansPos(1),iVar) * scalProd(1)
          end do

        end do
      end if
    end do YLoop

  end subroutine modg_2d_project_physFluxX_P


  !> Projection of the physical flux in y direction onto the testfunctions.
  subroutine modg_2d_project_physFluxY_P( nScalars, maxPolyDegree, length, &
    &                                     dl_prod, state , iElem, nDofs,   &
    &                                     state_der                        )
    ! --------------------------------------------------------------------------
    !> The number of scalar variables in your equation system.
    integer, intent(in) :: nScalars
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    !> The element index
    integer, intent(in) :: iElem
    !> Number of degrees of freedom
    integer, intent(in) :: nDofs
    !> The state data for the element
    real(kind=rk), intent(in)  :: state_der(nDofs,nScalars)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj
    real(kind=rk) :: scalProdY(maxPolyDegree)
    integer :: ansPos(2), testPos
    integer :: iTestX
    integer :: iAnsY
    integer :: iVar
    integer :: var_lb, var_ub
    real(kind=rk) :: scalProd(2)
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test
    ! functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**1

    var_lb = lbound(state,3)
    var_ub = ubound(state,3)

    !  for y direction (y test function differentiated)
    do iAnsY=1,maxPolyDegree
      scalProdY(iAnsY) = ply_scalProdDualLegDiff(iAnsY, iAnsY+1) &
        &                * jacobiDetStiffProj
    end do

    ! now, project onto all test functions
    XLoop: do iTestX = 1, maxPolyDegree+1

      if (iTestX > 2) then
        ! Need to add two terms
        do iAnsY = 1, maxPolyDegree - (iTestX-1)
?? copy :: posOfModgCoeffPTens2D(iTestX, iAnsY+1, testPos)
          scalProd(1) = dl_prod(1, iTestX) * scalProdY(iAnsY)
          scalProd(2) = dl_prod(2, iTestX) * scalProdY(iAnsY)
?? copy :: posOfModgCoeffPTens2D(iTestX-2, iAnsY, ansPos(1))
?? copy :: posOfModgCoeffPTens2D(iTestX, iAnsY, ansPos(2))
          do iVar=var_lb,var_ub
            state(iElem,testPos,iVar)                      &
              &  = state(iElem,testPos,iVar)               &
              &  + state_der(ansPos(1),iVar) * scalProd(1) &
              &  + state_der(ansPos(2),iVar) * scalProd(2)
          end do
        end do
      else
        ! Need to add one term
        do iAnsY = 1, maxPolyDegree - (iTestX-1)
?? copy :: posOfModgCoeffPTens2D(iTestX, iAnsY+1, testPos)
          scalProd(1) = dl_prod(2, iTestX) * scalProdY(iAnsY)
?? copy :: posOfModgCoeffPTens2D(iTestX, iAnsY, ansPos(1))
          do iVar = var_lb,var_ub
            state(iElem,testPos,iVar)                    &
              &  = state(iElem,testPos,iVar)             &
              &  + state_der(ansPos(1),iVar) * scalProd(1)
          end do
        end do
      end if
    end do XLoop

  end subroutine modg_2d_project_physFluxY_P


?? text :: prj_pflux_q(x, y, z, dir, s)
  !> Projection of the physical flux onto the testfunctions, with unrolled loops
  subroutine modg_prj_pFlux?dir?_Q_?s?( maxPolyDegree, length, state, dl_prod, &
    &                                 iElem, state_der )
    ! --------------------------------------------------------------------------
    !> The maximal polynomial degree in each spatial direction.
    integer, intent(in) :: maxPolyDegree
    !> The length of the cubes.
    real(kind=rk), intent(in) :: length
    !> The state to alter.
    real(kind=rk), intent(inout) :: state(:,:,:)
    !> Precomputed dual Legendre products:
    real(kind=rk), intent(in) :: dl_prod(2, maxPolyDegree+1)
    integer, intent(in) :: iElem
    !> The state to be used to project the physical fluxes
    real(kind=rk), intent(in)  :: state_der(:,:)
    ! --------------------------------------------------------------------------
    real(kind=rk) :: jacobiDetStiffProj, prescal !!(maxPolyDegree)
    integer :: testPos
    integer :: iTest2, iTest1, iTest3, min2mpd
    integer :: iAns1, iAns2, iAns3
    integer :: iVar
    integer :: ansPos1, anspos2, anspos3, anspos4
    real(kind=rk) :: scalProd1, scalprod2, scalprod3, scalprod4
    integer :: jk
    ! --------------------------------------------------------------------------

    ! Jacobi determinant for pojections of the physical fluxes onto the test functions.
    ! This is the stiffness term!
    !
    ! We have cubic elements, so the determinant of the jacobian of the mapping
    ! from reference element to physical element is the same everywhere.
    ! Please notice, that the mapping of the element itself is usually
    ! (mesh%length/2.0)**3, but the derivative in the volume integral
    ! gives an additional prefactor of 2.0/mesh%length and therefore
    ! the following is the correct scaling factor for the volume integrals.
    jacobiDetStiffProj = (0.5_rk*length)**2

    !!do iTest1=2,maxPolyDegree+1
    !!  prescal(iTest1-1) = ply_scalProdDualLegDiff(iTest1-1, iTest1) &
    !!    &                 * jacobiDetStiffProj
    !!end do
    ! We only consider the non-zeroes here, and these are all 2.
    prescal = 2*jacobiDetStiffProj

    min2mpd = min(2, maxPolyDegree+1)

    ! unrolled loop

    do iTest3 = 1, min2mpd

      !$NEC ivdep
      do jk=1,min2mpd*maxPolyDegree
        iTest2 = mod(jk-1,min2mpd) + 1
        iTest1 = (jk-1)/min2mpd + 2

        ! one entry

        testpos = itest?x?                                      &
          &      + ( ( itest?y?-1)                             &
          &      + (itest?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)

        iAns1 = iTest1-1
        iAns2 = iTest2
        iAns3 = iTest3
        anspos1 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,iTest2) * dl_prod(2,iTest3)

        !$NEC unroll(?s?)
        do iVar=1,?s?
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(anspos1,iVar) * scalprod1
        end do

      end do

      !$NEC ivdep
      do jk=1,(maxPolyDegree-1)*maxPolyDegree
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 3

        ! two entries

        testpos = itest?x?                                      &
          &      + ( ( itest?y?-1)                             &
          &      + (itest?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAns1 = iTest1-1
        iAns2 = iTest2-2
        iAns3 = iTest3
        anspos1 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,iTest2) * dl_prod(2,iTest3)

        iAns1 = iTest1-1
        iAns2 = iTest2
        iAns3 = iTest3
        anspos2 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,iTest2) * dl_prod(2,iTest3)

        !$NEC unroll(?s?)
        do iVar=1,?s?
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(anspos1,iVar) * scalprod1 &
            & + state_der(anspos2,iVar) * scalprod2
        end do

      end do

    end do

    do iTest3 = 3, maxPolyDegree+1

      !$NEC ivdep
      do jk=1,min2mpd*maxPolyDegree
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 1

        ! two entries

        testpos = itest?x?                                      &
          &      + ( ( itest?y?-1)                             &
          &      + (itest?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAns1 = iTest1-1
        iAns2 = iTest2
        iAns3 = iTest3-2
        anspos1 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(2,iTest2) * dl_prod(1,iTest3)

        iAns1 = iTest1-1
        iAns2 = iTest2
        iAns3 = iTest3
        anspos2 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,iTest2) * dl_prod(2,iTest3)

        !$NEC unroll(?s?)
        do iVar=1,?s?
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar) &
            & + state_der(anspos1,iVar) * scalprod1 &
            & + state_der(anspos2,iVar) * scalprod2
        end do
      end do

      !$NEC ivdep
      do jk=1,maxPolyDegree*(maxPolyDegree-1)
        iTest1 = mod(jk-1,maxPolyDegree) + 2
        iTest2 = (jk-1)/maxPolyDegree + 3

        ! four entries

        testpos = itest?x?                                      &
          &      + ( ( itest?y?-1)                             &
          &      + (itest?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)


        iAns1 = iTest1-1
        iAns2 = iTest2-2
        iAns3 = iTest3-2
        anspos1 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod1 = prescal * dl_prod(1,iTest2) * dl_prod(1,iTest3)

        iAns1 = iTest1-1
        iAns2 = iTest2
        iAns3 = iTest3-2
        anspos2 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod2 = prescal * dl_prod(2,iTest2) * dl_prod(1,iTest3)

        iAns1 = iTest1-1
        iAns2 = iTest2-2
        iAns3 = iTest3
        anspos3 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod3 = prescal * dl_prod(1,iTest2) * dl_prod(2,iTest3)

        iAns1 = iTest1-1
        iAns2 = iTest2
        iAns3 = iTest3
        anspos4 = ians?x?                                      &
          &      + ( ( ians?y?-1)                             &
          &      + (ians?z?-1)*(maxpolydegree+1))*(maxpolydegree+1)
        scalprod4 = prescal * dl_prod(2,iTest2) * dl_prod(2,iTest3)

        !$NEC unroll(?s?)
        do iVar=1,?s?
          state(iElem,testPos,iVar) &
            & = state(iElem,testPos,iVar)  &
            & + state_der(anspos1,iVar) * scalprod1 &
            & + state_der(anspos2,iVar) * scalprod2 &
            & + state_der(anspos3,iVar) * scalprod3 &
            & + state_der(anspos4,iVar) * scalprod4
        end do

      end do

    end do


  end subroutine modg_prj_pFlux?dir?_Q_?s?
?? end text

  !!dirVec(:,1) = [ 1,2,3 ]
  !!dirVec(:,2) = [ 2,1,3 ]
  !!dirVec(:,3) = [ 2,3,1 ]

  !> X direction for 6 scalars
?? copy :: prj_pflux_q(1,2,3,1,6)

  !> Y direction for 6 scalars
?? copy :: prj_pflux_q(2,1,3,2,6)

  !> Z direction for 6 scalars
?? copy :: prj_pflux_q(2,3,1,3,6)

  !> X direction for 5 scalars
?? copy :: prj_pflux_q(1,2,3,1,5)

  !> Y direction for 5 scalars
?? copy :: prj_pflux_q(2,1,3,2,5)

  !> Z direction for 5 scalars
?? copy :: prj_pflux_q(2,3,1,3,5)

end module atl_project_physflux_module
