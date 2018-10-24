module simplex
    use sort_n_swap
    implicit none

    ! FOR COMPILATION WITH F2PY COPY DECLARATIONS INTO SUBROUTINES

    ! parameters
    ! reflection, expansion, contraction, shrinking
    real(8), parameter :: a = 1.d0, g = 2.d0, q = .5d0, s = .5d0
    ! dimension
    integer, parameter :: n = 2

    ! indices
    ! integer, dimension(3) :: idx
    ! values fn(x)
    real(8), dimension(n+1) :: z
    ! middle, reflected, expanded, contracted
    real(8), dimension(n) :: m, r, e, c
    ! reflected value
    real(8) :: fr
    ! loop index
    integer :: i

contains
    subroutine step(x, fn)
        ! arguments
        real(8), dimension(n, n+1), intent(inout) :: x
        ! declaring fn two times is required for compilation with f2py
        real(8) :: fn
        external :: fn

        ! calculate and sort values fn(x)
        ! ONLY FOR size(x) == 3
        z = [( fn(x(:, i)), i=1,n+1 )]
        call sort3r(x, z)
        ! argsort doesn't quite work yet...
        ! idx = argsort3(z)
        ! z = z(idx)
        ! x = x(:,idx)

        ! calculate middle, omitting x(:,n+1)
        m = sum(x(:,1:n), dim=1) / n

        ! calculate reflection
        r = m + a*(m - x(:, n+1))
        fr = fn(r)
        if (z(1) <= fr .and. fr <= z(n)) then
            x(:, n+1) = r
            return
        end if

        ! calculate expansion
        if (fr < z(1)) then
            e = m + g*(r - m)
            if (fn(e) < fr) then
                x(:, n+1) = e
            else
                x(:, n+1) = r
            end if
            return
        end if

        ! calculate contraction
        ! here it is known, that fn(r) >= fn(x(n))
        c = m + q*(x(:, n+1) - m)
        if (fn(c) < z(n+1)) then
            x(:, n+1) = c
            return
        end if

        ! if all else fails: shrinking
        do i=2, n+1
            x(:, i) = x(:, 1) + s*(x(:, i) - x(:, 1))
        end do
        return
    end subroutine step

    subroutine find(x, fn, tol, maxiter, conv)
        ! arguments
        real(8), dimension(n, n+1), intent(inout) :: x
        ! declaring fn two times is required for compilation with f2py
        real(8) :: fn
        external :: fn
        real(8), intent(in) :: tol
        integer, intent(in) :: maxiter
        logical, intent(out) :: conv
        integer :: step
        conv = .false.
        step = 0

        100 step = step + 1
        if (all(std(x, n) < tol)) then
            conv = .true.
            return
        else if (step > maxiter) then
            return
        end if

        ! calculate and sort values fn(x)
        ! ONLY FOR size(x) == 3
        z = [( fn(x(:, i)), i=1,n+1 )]
        call sort3r(x, z)

        ! calculate middle, omitting x(:,n+1)
        m = sum(x(:,1:n), dim=1) / n

        ! calculate reflection
        r = m + a*(m - x(:, n+1))
        fr = fn(r)
        if (z(1) <= fr .and. fr <= z(n)) then
            x(:, n+1) = r
            goto 100
        end if

        ! calculate expansion
        if (fr < z(1)) then
            e = m + g*(r - m)
            if (fn(e) < fr) then
                x(:, n+1) = e
            else
                x(:, n+1) = r
            end if
            goto 100
        end if

        ! calculate contraction
        ! here it is known, that fn(r) >= fn(x(n))
        c = m + q*(x(:, n+1) - m)
        if (fn(c) < z(n+1)) then
            x(:, n+1) = c
            goto 100
        end if

        ! if all else fails: shrinking
        do i=2, n+1
            x(:, i) = x(:, 1) + s*(x(:, i) - x(:, 1))
        end do
        goto 100
    end subroutine find

end module simplex


program test
    use simplex, only: step, find
    implicit none
    ! integer :: i
    real(8), dimension(2, 3) :: x = reshape([-1., -3.5, 0., -4.5, 1., -3.5], &
                                            [2, 3])
    logical :: conv

    100 format(2(3(3X, F9.6)/))

    call find(x, himmelblau, 1d-7, 100, conv)
    print 100, transpose(x)
    print *, 'convergence: ', conv

    ! print 100, transpose(x)
    ! do i=1, 20
    !     call step(x, himmelblau)
    !     print 100, transpose(x)
    ! end do

contains
    real(8) pure &
    function himmelblau(x)
        real(8), dimension(2), intent(in) :: x
        himmelblau = (x(1)**2+x(2)-11d0)**2 + (x(1)+x(2)**2-7d0)**2
    end function himmelblau

end program
