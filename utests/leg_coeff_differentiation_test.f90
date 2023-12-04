! Copyright (c) 2014 Nikhil Anand <nikhil.anand@uni-siegen.de>
! Copyright (c) 2014 Timo Stentenbach
! Copyright (c) 2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017, 2020 Peter Vitt <peter.vitt2@uni-siegen.de>
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

! This is the unit test for checking the differentiation of legendre coeffs.
!
program atl_leg_coeff_diff_test
  use env_module,               only: rk, fin_env
  use tem_logging_module,       only: logUnit
  use tem_general_module,       only: tem_start
  use ply_leg_diff_module,      only: ply_calcDiff_leg,    &
   &                                  ply_calcDiff_leg_2d, &
   &                                  ply_calcDiff_leg_1d
  use atl_solver_param_module,  only: atl_solver_param_type
implicit none

  real(kind=rk) :: res, newRes
  integer :: dimen
  type(atl_solver_param_type) :: params

  ! Init the Treelm environment, needed to init the log Unit
  call tem_start(codeName = 'Ateles unit test', &
    &            version  = 'utest',            &
    &            general  = params%general      )
  res = 0.0_rk


  do dimen = 1,3
    write(logUnit(1),*) '-------- CHECKING LEGCOEFF  DIFF FOR DIM = ',dimen
    call check_diff(dimen, newRes)
    if (newRes.gt.res) then
      res = newRes
    end if
    write(logUnit(1),*) '--------------------------- DONE'
  end do

  !>\todo If everything worked fine, write PASSED on the very last line of output, to
  !!      indicate a successful run of the unit test:
  if(res.lt.1e-08) then
    write(logUnit(1),*) 'PASSED'
  else
    write(logUnit(1),*) "FAILED error=", res
  end if
  call fin_env()

contains

  subroutine check_diff(dimen, res)
    integer, intent(in) :: dimen
    real(kind=rk), intent(out) :: res
    !---------------------------------------------
    real(kind=rk), allocatable    :: legCoeffs(:,:)
    real(kind=rk), allocatable    :: legCoeffsDiff(:,:,:), anl_diff(:,:)
    real(kind=rk)  ::  elemLength
    integer ::  maxPolyDegree, nVars, i
    res = 0.0

    if (dimen ==3) then
      allocate(legcoeffs(27,1))
      allocate(legcoeffsDiff(27,1,3))
      allocate(anl_diff(27,3))
      call get_coeff(legCoeffs)
      call get_anl_3d(anl_diff)

      maxPolyDegree = 2
      nVars = 1
      elemLength = 2.0

      call ply_calcDiff_leg( legCoeffs = legCoeffs,         &
        &                    legCoeffsDiff = legCoeffsDiff, &
        &                    maxPolyDegree = maxPolyDegree, &
        &                    nVars         = nVars,         &
        &                    elemLength    = elemLength     )

       do i = 1, size(legcoeffs,1)
         res = res + legCoeffsDiff(i,1,1) - anl_diff(i,1)
         res = res + legCoeffsDiff(i,1,2) - anl_diff(i,2)
         res = res + legCoeffsDiff(i,1,3) - anl_diff(i,3)
       enddo

    elseif (dimen ==2) then
      allocate(legcoeffs(16,1))
      allocate(legcoeffsDiff(16,1,2))
      allocate(anl_diff(16,2))
      call get_coeff(legCoeffs)
      call get_anl_2d(anl_diff)

      maxPolyDegree = 3
      nVars = 1
      elemLength = 2.0

      call ply_calcDiff_leg_2d( legCoeffs     = legCoeffs,     &
        &                       legCoeffsDiff = legCoeffsDiff, &
        &                       maxPolyDegree = maxPolyDegree, &
        &                       nVars         = nVars,         &
        &                       elemLength    = elemLength     )

       do i = 1, size(legcoeffs,1)
         res = res + legCoeffsDiff(i,1,1) - anl_diff(i,1)
         res = res + legCoeffsDiff(i,1,2) - anl_diff(i,2)
       enddo

    else

      allocate(legcoeffs(4,1))
      allocate(legCoeffsDiff(4,1,1))
      allocate(anl_diff(4,1))
      call get_coeff(legCoeffs)
      call get_anl_1d(anl_diff)

      maxPolyDegree = 3
      nVars = 1
      elemLength = 2.0

      call ply_calcDiff_leg_1d( legCoeffs     = legCoeffs,            &
        &                       legCoeffsDiff = legCoeffsDiff(:,:,1), &
        &                       maxPolyDegree = maxPolyDegree,        &
        &                       elemLength    = elemLength            )

       do i = 1, size(legcoeffs,1)
         res = res + legCoeffsDiff(i,1,1) - anl_diff(i,1)
       enddo

    endif

    res = abs(res)

  end subroutine check_diff

  subroutine get_coeff(legCoeffs)
    real(kind=rk), intent(inout)    :: legCoeffs(:,:)
    integer coeffs, iCoeff
    coeffs = size(legCoeffs,1)
    do iCoeff = 1, coeffs
      legCoeffs(iCoeff,1) = iCoeff
    enddo
  end subroutine get_coeff

  subroutine get_anl_3d(anl_diff)
    real(kind=rk), intent(inout)   :: anl_diff(:,:)
    anl_diff(1,1)  =  2.0
    anl_diff(2,1)  =  9.0
    anl_diff(3,1)  =  0.0
    anl_diff(4,1)  =  5.0
    anl_diff(5,1)  = 18.0
    anl_diff(6,1)  =  0.0
    anl_diff(7,1)  =  8.0
    anl_diff(8,1)  = 27.0
    anl_diff(9,1)  =  0.0
    anl_diff(10,1) = 11.0
    anl_diff(11,1) = 36.0
    anl_diff(12,1) =  0.0
    anl_diff(13,1) = 14.0
    anl_diff(14,1) = 45.0
    anl_diff(15,1) =  0.0
    anl_diff(16,1) = 17.0
    anl_diff(17,1) = 54.0
    anl_diff(18,1) =  0.0
    anl_diff(19,1) = 20.0
    anl_diff(20,1) = 63.0
    anl_diff(21,1) =  0.0
    anl_diff(22,1) = 23.0
    anl_diff(23,1) = 72.0
    anl_diff(24,1) =  0.0
    anl_diff(25,1) = 26.0
    anl_diff(26,1) = 81.0
    anl_diff(27,1) =  0.0

    anl_diff(1,2)  =  4.0
    anl_diff(2,2)  =  5.0
    anl_diff(3,2)  =  6.0
    anl_diff(4,2)  = 21.0
    anl_diff(5,2)  = 24.0
    anl_diff(6,2)  = 27.0
    anl_diff(7,2)  =  0.0
    anl_diff(8,2)  =  0.0
    anl_diff(9,2)  =  0.0
    anl_diff(10,2) = 13.0
    anl_diff(11,2) = 14.0
    anl_diff(12,2) = 15.0
    anl_diff(13,2) = 48.0
    anl_diff(14,2) = 51.0
    anl_diff(15,2) = 54.0
    anl_diff(16,2) =  0.0
    anl_diff(17,2) =  0.0
    anl_diff(18,2) =  0.0
    anl_diff(19,2) = 22.0
    anl_diff(20,2) = 23.0
    anl_diff(21,2) = 24.0
    anl_diff(22,2) = 75.0
    anl_diff(23,2) = 78.0
    anl_diff(24,2) = 81.0
    anl_diff(25,2) =  0.0
    anl_diff(26,2) =  0.0
    anl_diff(27,2) =  0.0

    anl_diff(1,3)  = 10.0
    anl_diff(2,3)  = 11.0
    anl_diff(3,3)  = 12.0
    anl_diff(4,3)  = 13.0
    anl_diff(5,3)  = 14.0
    anl_diff(6,3)  = 15.0
    anl_diff(7,3)  = 16.0
    anl_diff(8,3)  = 17.0
    anl_diff(9,3)  = 18.0
    anl_diff(10,3) = 57.0
    anl_diff(11,3) = 60.0
    anl_diff(12,3) = 63.0
    anl_diff(13,3) = 66.0
    anl_diff(14,3) = 69.0
    anl_diff(15,3) = 72.0
    anl_diff(16,3) = 75.0
    anl_diff(17,3) = 78.0
    anl_diff(18,3) = 81.0
    anl_diff(19,3) =  0.0
    anl_diff(20,3) =  0.0
    anl_diff(21,3) =  0.0
    anl_diff(22,3) =  0.0
    anl_diff(23,3) =  0.0
    anl_diff(24,3) =  0.0
    anl_diff(25,3) =  0.0
    anl_diff(26,3) =  0.0
    anl_diff(27,3) =  0.0
  end subroutine get_anl_3d

  subroutine get_anl_2d(anl_diff)
    real(kind=rk), intent(inout)   :: anl_diff(:,:)
    anl_diff(1,1)  =  6.0
    anl_diff(2,1)  =  9.0
    anl_diff(3,1)  = 20.0
    anl_diff(4,1)  =  0.0
    anl_diff(5,1)  = 14.0
    anl_diff(6,1)  = 21.0
    anl_diff(7,1)  = 40.0
    anl_diff(8,1)  =  0.0
    anl_diff(9,1)  = 22.0
    anl_diff(10,1) = 33.0
    anl_diff(11,1) = 60.0
    anl_diff(12,1) =  0.0
    anl_diff(13,1) = 30.0
    anl_diff(14,1) = 45.0
    anl_diff(15,1) = 80.0
    anl_diff(16,1) =  0.0

    anl_diff(1,2)  = 18.0
    anl_diff(2,2)  = 20.0
    anl_diff(3,2)  = 22.0
    anl_diff(4,2)  = 24.0
    anl_diff(5,2)  = 27.0
    anl_diff(6,2)  = 30.0
    anl_diff(7,2)  = 33.0
    anl_diff(8,2)  = 36.0
    anl_diff(9,2)  = 65.0
    anl_diff(10,2) = 70.0
    anl_diff(11,2) = 75.0
    anl_diff(12,2) = 80.0
    anl_diff(13,2) =  0.0
    anl_diff(14,2) =  0.0
    anl_diff(15,2) =  0.0
    anl_diff(16,2) =  0.0
  end subroutine get_anl_2d

  subroutine get_anl_1d(anl_diff)
    real(kind=rk), intent(inout)   :: anl_diff(:,:)
    anl_diff(1,1)  =  6.0
    anl_diff(2,1)  =  9.0
    anl_diff(3,1)  = 20.0
    anl_diff(4,1)  =  0.0
  end subroutine get_anl_1d
end program atl_leg_coeff_diff_test
