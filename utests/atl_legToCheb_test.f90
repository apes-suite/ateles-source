! Copyright (c) 2012, 2014 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2013-2014, 2017 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2014, 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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

!> Unit test to check functionallity of fast polynomial transformations.
!! \author{Jens Zudrop}
program atl_legToCheb_test
  use env_module,               only: rk, fin_env
  use tem_logging_module,       only: logUnit
  use tem_aux_module,           only: tem_abort
  use ply_polyBaseExc_module,   only: ply_fpt_init, ply_fpt_exec_striped, &
                                    & ply_legToCheb_param, &
                                    & ply_trafo_params_type
  use tem_general_module,       only: tem_start
  use atl_solver_param_module,  only: atl_solver_param_type 

  implicit none

  integer :: power
  real(kind=rk) :: res, newRes
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            general  = params%general      )

  ! Check the transformation for different number of polynomial coefficients.
  res = 0.0_rk
  do power = 1,7
    write(logUnit(1),*) '---------------------------   CHECKING LEG->CHEB TRAFO FOR ', 2**power
    call check(power, newRes)
    if(newRes.gt.res) then
      res = newRes
    end if
  end do

  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  end if 

  call fin_env()


contains


  subroutine check(power, res)
    integer, intent(in) :: power
    real(kind=rk), intent(out) :: res
    integer :: maxPolyDegree
    type(ply_trafo_params_type) :: fptParams, fptParamsRef
    real(kind=rk), allocatable :: legCoeffs(:), chebCoeffs(:), refResult(:)

    ! Define the maximal polynomial degree we want to calculate the
    ! bases exchange for.
    maxPolyDegree = 2**power - 1  ! maxPolyDegree+1 has to be a power of 2

    ! Create the Legendre expansion coefficients
    allocate(legCoeffs(1:maxPolyDegree+1)) 
    legCoeffs(:) = 1.0_rk
    !write(logUnit(1),*) 'Input as Legendre coefficients:'
    !write(*,*) legCoeffs(0:maxPolyDegree)

    ! Define the reference results
    call ply_fpt_init(n                = maxPolyDegree + 1,   &
      &               params           = fptParamsRef,        &
      &               trafo            = ply_legToCheb_param, &
      &               striplen         = 1,                   &
      &               blocksize        = maxPolyDegree + 1,   &
      &               subblockingWidth = 2 * power            )
    write(logUnit(1),*) 'Calculating reference result ... '
    allocate(refResult(1 : maxPolyDegree + 1)) 
    call ply_fpt_exec_striped( alph    = legCoeffs,   &
      &                        gam     = refResult,   &
      &                        nIndeps = 1,           &
      &                        params  = fptParamsRef )

    ! The reference results are ...
    !write(logUnit(1),*) 'Results for Chebyshev coefficients: '
    !write(*,*) refResult(0:maxPolyDegree)

    ! Initialize the matrix for Legendre to Chebyshev bases exchange.
    write(logUnit(1),*) 'Initializing Leg->Cheb base exchange ... '
    call ply_fpt_init(n        = maxPolyDegree + 1,  &
      &               striplen = 1,                  &
      &               params   = fptParams,          &
      &               trafo    = ply_legToCheb_param )
    write(logUnit(1),*) '... finished.'
    
    ! Now, execute a bases transformation with a given vector. Noitce that 
    ! gam will be allocated inside the routine.
    write(logUnit(1),*) 'Calculating Leg->Cheb base exchange ... '
    allocate(chebCoeffs(1:maxPolyDegree+1)) 
    legCoeffs(:) = 1.0_rk
    call ply_fpt_exec_striped( alph    = legCoeffs,  &
      &                        gam     = chebCoeffs, &
      &                        nIndeps = 1,          &
      &                        params  = fptParams   )
    write(logUnit(1),*) '... finished.'

    ! The Chebyshev weights are ...
    !write(logUnit(1),*) 'Results for Chebyshev coefficients: '
    !write(*,*) chebCoeffs(0:maxPolyDegree)

    ! Check the error for each Chebyshev coefficient
    write(*,*) 'Maximum error for Cheb. coeff ', maxloc(abs(refResult(:)-chebCoeffs(:))), &
              & ' is: ', maxval(abs(refResult(:)-chebCoeffs(:)))

    res = maxval(abs(refResult(:)-chebCoeffs(:)))
 
  end subroutine 


end program atl_legToCheb_test
