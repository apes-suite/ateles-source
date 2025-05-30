! Copyright (c) 2013-2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2013-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
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

! This is the unit test for the projection module.
program atl_project_2d_test
  use env_module,                   only: rk, fin_env
  use tem_logging_module,           only: logUnit
  use ply_dof_module,               only: Q_space
  use ply_dynArray_project_module,  only: ply_prj_init_define, &
    &                                     ply_prj_init_type
  use ply_poly_project_module,      only: ply_poly_project_fillbody, &
    &                                     ply_poly_project_m2n,&
    &                                     ply_poly_project_n2m, &
    &                                     ply_poly_project_type
  use ply_prj_header_module,        only: ply_prj_header_type
  use tem_general_module,           only: tem_start
  use atl_solver_param_module,      only: atl_solver_param_type
implicit none

  real(kind=rk) :: res, newRes
  integer :: power
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            general  = params%general      )
  res = 0.0_rk

  !>\todo Check reading of the projection configuration.
  !!      The input has to be configured on the fly!

  !>\todo Put various projections into a projection descriptor.

  !>\todo Check those projections by doing m2n and n2m, where
  !!      input and output should be the same within certain
  !!      bounds.

  ! check l2p Q-Space
  do power = 1,7
    write(logUnit(10),*) '---------------------------   CHECKING CHEB->LEG L2P Q-SPACE TRAFO FOR ', 2**power
    call check_l2p_qspace_2d(power, newRes)
    if (newRes.gt.res) then
      res = newRes
    end if
    write(logUnit(10),*) '--------------------------- DONE'
  end do

  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if
  call fin_env()

contains

  subroutine check_l2p_qspace_2d(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    !---------------------------------------------
    type(ply_prj_header_type) :: header
    type(ply_poly_project_type)   :: me
    type(ply_prj_init_type)   :: prj_init
    real(kind=rk),allocatable     :: nodal_data(:,:)
    real(kind=rk), allocatable    :: modal_data(:,:)
    real(kind=rk), allocatable    :: oversamp_modal(:,:)
    real(kind=rk), allocatable    :: ref_modes(:)
    !-----------for init
    integer ::basisType, maxdegree, i
    !-----------for oversamp
    integer :: iDegX, iDegY, dof, dofOverSamp

    basisType = Q_space
    maxdegree = power
    header%kind = 'l2p'
    header%l2p_header%nodes_header%nodes_kind = 'gauss-legendre'
    header%l2p_header%factor = 2.0_rk

    ! define my poly projection type
    call ply_prj_init_define(me= prj_init,             &
      &                          header = header ,         &
      &                          maxPolyDegree= maxdegree, &
      &                          basisType = basistype )
    ! fill the projection body according to the header
    call ply_poly_project_fillbody(me, proj_init=prj_init,scheme_dim=2)

    allocate(ref_modes(1:me%body_2d%ndofs) )
    allocate(modal_data(1:me%body_2d%ndofs,1) )
    allocate(oversamp_modal(1:me%body_2d%oversamp_dofs,1) )
    allocate(nodal_data(1:me%body_2d%nQuadPoints,1) )

    do i=1, me%body_2d%ndofs
       ref_modes(i)=1.0/real(i, kind=rk)
    end do

    modal_data(:,:) = 0.0_rk
    oversamp_modal(:,:) = 0.0_rk
    ! oversampling of modes
    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1)
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)
        oversamp_modal(dofOverSamp,1) = ref_modes(dof)
      end do
    end do

    call ply_poly_project_m2n(me = me ,                &
      &                       dim = 2 ,                &
      &                       nVars = 1,               &
      &                       nodal_data=nodal_data,   &
      &                       modal_data=oversamp_modal)
    call ply_poly_project_n2m(me = me,                  &
      &                       dim = 2 ,                 &
      &                       nVars = 1,                &
      &                       nodal_data=nodal_data,    &
      &                       modal_data= oversamp_modal)

    do iDegX = 1, maxdegree+1
      do iDegY = 1, maxdegree+1
        dof = 1 + (iDegX-1) + (iDegY-1)*(maxdegree+1)
        dofOverSamp = 1 + (iDegX-1) + (iDegY-1)*(me%oversamp_degree+1)
        modal_data(dof,:) = oversamp_modal(dofOverSamp,:)
      end do
    end do

    res= maxval( abs(ref_modes-modal_data(:,1)) )

    write(*,*) 'power=', power
    write(*,*) 'res=', res
  end subroutine  check_l2p_qspace_2d

end program atl_project_2d_test
