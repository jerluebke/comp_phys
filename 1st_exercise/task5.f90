module simplex
    implicit none

contains
    subroutine step(x, fn)
        ! parameters
        ! reflection, expansion, contraction, shrinking
        real(8), parameter :: a = 1.d0, g = 2.d0, q = .5d0, s = .5d0
        ! dimension
        integer, parameter :: n = 2

        ! arguments
        real(8), dimension(n, n+1), intent(inout) :: x
        real(8) :: fn
        external :: fn

        ! values fn(x)
        real(8), dimension(n+1) :: z
        ! middle, reflected, expanded, contracted
        real(8), dimension(n) :: m, r, e, c
        ! reflected value
        real(8) :: fr
        ! loop index
        integer :: i

        ! calculate and sort values fn(x)
        ! ONLY FOR size(x) == 3
        z = [( fn(x(:, i)), i=1,n+1 )]
        call sort3r(x, z)

        ! calculate middle, omitting x(:,n+1)
        m = sum(x(:,1:n), dim=2) / n

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


    subroutine sort3r(x, z)
        real(8), dimension(3) :: x, z
        if (z(1) > z(2)) then
            call swapr(x(1), x(2))
            call swapr(z(1), z(2))
        end if
        if (z(2) > z(3)) then
            call swapr(x(2), x(3))
            call swapr(z(2), z(3))
        end if
        if (z(1) > z(2)) then
            call swapr(x(1), x(2))
            call swapr(z(1), z(2))
        end if
    end subroutine sort3r

    subroutine swapr(a, b)
        real(8) :: a, b, t
        t = a
        a = b
        b = t
    end subroutine swapr
end module simplex

program test
    use simplex, only: step
    implicit none
    integer :: i
    real(8), dimension(2, 3) :: x = reshape([-1.9, -1.2, -0.7, -1.9, -1.5, -1.8], &
                                            [2, 3])

    100 format(2(3(3X, F4.1)/))

    print 100, x

    do i=1, 20
        call step(x, himmelblau)
        print 100, x
    end do

contains
    real(8) pure &
    function himmelblau(x)
        real(8), dimension(2), intent(in) :: x
        himmelblau = (x(1)**2+x(2)-11d0)**2 + (x(1)+x(2)**2-7d0)**2
    end function himmelblau

end program
