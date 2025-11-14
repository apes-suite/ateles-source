! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2013-2014, 2016-2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
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

?? include "ply_dof_module.inc"
!> Unit test to check functionallity of atl_modg_scaledTransposedInvMassMatrix_Q
program atl_modg_localPrediction
  use env_module,              only: rk, fin_env
  use tem_logging_module,      only: logUnit
  use tem_aux_module,          only: tem_abort
  use atl_modg_kernel_module,  only:               &
    & atl_modg_scaledTransposedInvMassMatrix_Q,    &
    & atl_modg_scaledTransposedInvMassMatrix_P,    &
    & atl_modg_scaledTransposedProject_physFlux_Q, &
    & atl_modg_scaledTransposedProject_physFlux_P
  use tem_general_module,      only: tem_start
  use atl_solver_param_module, only: atl_solver_param_type

  implicit none

  integer :: polyDegree
  real(kind=rk) :: res
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            general  = params%general      )

  do polyDegree = 1, 7
    write(logUnit(10),*) '---- CHECKING scaledTransposedInvMassMatrix_Q FOR degree', polyDegree
    call checkMass_Q(polyDegree, res)
    if( res .gt. 1.e-8 ) then
      write(logUnit(10),*) 'ERROR, res', res
      call tem_abort()
    end if
    write(logUnit(10),*) '---- DONE'
  end do

  do polyDegree = 1, 7
    write(logUnit(10),*) '---- CHECKING scaledTransposedProject_physFlux_Q FOR degree', polyDegree
    call checkProjection_Q(polyDegree, res)
    if( res .gt. 1.e-8 ) then
      write(logUnit(10),*) 'ERROR, res', res
      call tem_abort()
    end if
    write(logUnit(10),*) '---- DONE'
  end do

  do polyDegree = 1, 7
    write(logUnit(10),*) '---- CHECKING scaledTransposedInvMassMatrix_P FOR degree', polyDegree
    call checkMass_P(polyDegree, res)
    if( res .gt. 1.e-8 ) then
      write(logUnit(10),*) 'ERROR, res', res
      call tem_abort()
    end if
    write(logUnit(10),*) '---- DONE'
  end do

  do polyDegree = 1, 7
    write(logUnit(10),*) '---- CHECKING scaledTransposedProject_physFlux_P FOR degree', polyDegree
    call checkProjection_P(polyDegree, res)
    if( res .gt. 1.e-8 ) then
      write(logUnit(10),*) 'ERROR, res', res
      call tem_abort()
    end if
    write(logUnit(10),*) '---- DONE'
  end do


  write(logUnit(1),*) 'PASSED'
  call fin_env()

contains

  subroutine checkMass_Q(degree, res)
    integer, intent(in) :: degree
    real(kind=rk), intent(out) :: res

    real(kind=rk), allocatable :: state(:,:,:,:)
    integer :: i, j, k, pos
    integer :: i_, j_, k_, pos_
    real(kind=rk) :: factor

    ! alloc
    allocate(state(1,(degree+1)**3,4,1))

    ! initialize
    do k = 1, degree+1
      do j = 1, degree+1
        do i = 1, degree+1

?? copy :: posOfModgCoeffQTens(i, j, k, degree, pos)
          state(:,pos,2,:) = pos

        end do
      end do
    end do


    ! multiply with scaled mass matrix
    state(:,:,:,:) = 0.0_rk
    do i = 1, degree+1
      do i_ = i, min(degree+1,i+2),2
        do j = 1, degree+1
          do j_ = j, min(degree+1,j+2),2
            do k = 1, degree+1
              do k_ = k, min(degree+1,k+2),2

?? copy :: posOfModgCoeffQTens(i, j, k, degree, pos)
?? copy :: posOfModgCoeffQTens(i_, j_, k_, degree, pos_)

                factor = delta(i_,i) * delta(j_,j) * delta(k_,k)

                state(:,pos,1,:) = state(:,pos,1,:) + factor*state(:,pos_,2,:)
              end do
            end do
          end do
        end do
      end do
    end do

    call atl_modg_scaledTransposedInvMassMatrix_Q(nTotalElems = 1, &
      &                                           nDofs = (degree+1)**3, &
      &                                           nScalars = 1, &
      &                                           maxPolyDegree = degree, &
      &                                           nElems = 1, &
      &                                           state = state )

    res = sqrt(sum((state(1,:,1,1)-state(1,:,2,1))**2))


  end subroutine checkMass_Q


  subroutine checkProjection_Q(degree, res)
    integer, intent(in) :: degree
    real(kind=rk), intent(out) :: res

    real(kind=rk), allocatable :: state(:,:,:,:)
    integer :: i, j, k, pos
    integer :: i_, j_, k_, pos_
    real(kind=rk) :: factor

    ! alloc
    allocate(state(1,(degree+1)**3,4,1))

    ! initialize
    state(:,:,:,:) = 0.0_rk
    do k = 1, degree+1
      do j = 1, degree+1
        do i = 1, degree+1

?? copy :: posOfModgCoeffQTens(i, j, k, degree, pos)
          state(:,pos,2,:) = pos

        end do
      end do
    end do


    ! multiply with scaled transposed stiffness matrix
    do i = 1, degree
      i_ = i+1
        do j = 1, degree+1
          do j_ = j, min(degree+1,j+2),2
            do k = 1, degree+1
              do k_ = k, min(degree+1,k+2),2

?? copy :: posOfModgCoeffQTens(i, j, k, degree, pos)
?? copy :: posOfModgCoeffQTens(i_, j_, k_, degree, pos_)

                factor = - (2*i-1) * delta(j_,j) * delta(k_,k)

                state(:,pos,3,:) = state(:,pos,3,:) + factor*state(:,pos_,2,:)
              end do
            end do
          end do
        end do

    end do

    call atl_modg_scaledTransposedProject_physFlux_Q( nScalars = 1,           &
                                                    & maxPolyDegree = degree, &
                                                    & length = 2.0_rk,        &
                                                    & nElems = 1,             &
                                                    & state = state,          &
                                                    & iDir = 1,               &
                                                    & dirVec = (/1, 2, 3/)    )


    res = sqrt(sum((state(1,:,1,1)-state(1,:,3,1))**2))


  end subroutine checkProjection_Q


  subroutine checkMass_P(degree, res)
    integer, intent(in) :: degree
    real(kind=rk), intent(out) :: res

    real(kind=rk), allocatable :: state(:,:,:,:)
    integer :: i, j, k, pos, dofs
    integer :: i_, j_, k_, pos_
    real(kind=rk) :: factor

    ! alloc
?? copy :: getDofsPTens(degree, dofs)
    allocate(state(1,dofs,4,1))

    ! initialize
    do k = 1, degree+1, 1
      do j = 1, degree+1 - (k-1), 1
        do i = 1, degree+1 - (k-1) - (j-1), 1

?? copy :: posOfModgCoeffPTens(i,j,k,pos)
          state(:,pos,2,:) = pos

        end do
      end do
    end do


    ! multiply with scaled mass matrix
    state(:,:,:,:) = 0.0_rk
    do i = 1, degree+1, 1
      do i_ = i, min(degree+1,i+2),2
        do j = 1, degree+1 - (i-1), 1
          do j_ = j, min(degree+1 - (i_-1),j+2),2
            do k = 1, degree+1 - (i-1) - (j-1), 1
              do k_ = k, min(degree+1 - (i_-1) - (j_-1),k+2),2

?? copy :: posOfModgCoeffPTens(i,j,k,pos)
?? copy :: posOfModgCoeffPTens(i_,j_,k_,pos_)

                factor = delta(i_,i) * delta(j_,j) * delta(k_,k)
                state(:,pos,1,:) = state(:,pos,1,:) + factor*state(:,pos_,2,:)
              end do
            end do
          end do
        end do
      end do
    end do

    call atl_modg_scaledTransposedInvMassMatrix_P(nTotalElems = 1, &
      &                                           nDofs = size(state,2), &
      &                                           nScalars = 1, &
      &                                           maxPolyDegree = degree, &
      &                                           nElems = 1, &
      &                                           state = state )

    res = sqrt(sum((state(1,:,1,1)-state(1,:,2,1))**2))


  end subroutine checkMass_P


  subroutine checkProjection_P(degree, res)
    integer, intent(in) :: degree
    real(kind=rk), intent(out) :: res

    real(kind=rk), allocatable :: state(:,:,:,:)
    integer :: i, j, k, pos, dofs
    integer :: i_, j_, k_, pos_
    real(kind=rk) :: factor

    ! alloc
?? copy :: getDofsPTens(degree, dofs)
    allocate(state(1,dofs,4,1))
    ! initialize
    state(:,:,:,:) = 0.0_rk
    do k = 1, degree+1, 1
      do j = 1, degree+1 - (k-1), 1
        do i = 1, degree+1 - (k-1) - (j-1), 1

?? copy :: posOfModgCoeffPTens(i,j,k,pos)
          state(:,pos,2,:) = pos

        end do
      end do
    end do


    ! multiply with scaled transposed stiffness matrix
    do i = 1, degree, 1
      i_ = i+1
        do j = 1, degree+1 - (i-1), 1
          do j_ = j, min(degree+1 - (i_-1),j+2),2
            do k = 1, degree+1 - (i-1) - (j-1), 1
              do k_ = k, min(degree+1 - (i_-1) - (j_-1),k+2),2

?? copy :: posOfModgCoeffPTens(i,j,k,pos)
?? copy :: posOfModgCoeffPTens(i_,j_,k_,pos_)

                factor = - (2*i-1) * delta(j_,j) * delta(k_,k)
                state(:,pos,3,:) = state(:,pos,3,:) + factor*state(:,pos_,2,:)
              end do
            end do
          end do
        end do

    end do

    call atl_modg_scaledTransposedProject_physFlux_P( nScalars = 1,           &
                                                    & maxPolyDegree = degree, &
                                                    & length = 2.0_rk,        &
                                                    & nElems = 1,             &
                                                    & state = state,          &
                                                    & iDir = 1,               &
                                                    & dirVec = (/1, 2, 3/)    )


    res = sqrt(sum((state(1,:,1,1)-state(1,:,3,1))**2))


  end subroutine checkProjection_P



  function delta(i,j)
    integer, intent(in) :: i,j
    real :: delta

    if( i .eq. j ) then
      delta = 1.0_rk
    elseif( i-2 .eq. j ) then
      delta = -1.0_rk
    else
      delta = 0.0_rk
    end if
  end function delta

end program atl_modg_localPrediction
