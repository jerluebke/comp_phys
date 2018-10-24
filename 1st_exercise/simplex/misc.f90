module sort_n_swap
    implicit none
contains
    function std(x, n)
        real(8), dimension(:,:), intent(in) :: x
        integer, intent(in) :: n
        real(8), dimension(n) :: std
        real(8), dimension(n) :: mean
        integer :: s, i
        s = size(x, 2)
        do i=1,n
            mean(i) = sum(x(i,:)) / s
            std(i) = sqrt(sum((x(i,:) - mean(i))**2)) / s
        end do
    end function std

    function argsort3(x) result(idx)
        real(8), dimension(3), intent(in) :: x
        integer, dimension(3) :: idx
        idx = [1, 2, 3]
        if (x(1) > x(2)) call swapi(idx(1), idx(2))
        if (x(2) > x(3)) call swapi(idx(2), idx(3))
        if (x(1) > x(2)) call swapi(idx(1), idx(2))
    end function

    subroutine sort3r(x, z)
        real(8), dimension(3) :: z
        real(8), dimension(2, 3) :: x
        if (z(1) > z(2)) then
            call swapr(x(1,1), x(1,2))
            call swapr(x(2,1), x(2,2))
            call swapr(z(1), z(2))
        end if
        if (z(2) > z(3)) then
            call swapr(x(1,2), x(1,3))
            call swapr(x(2,2), x(2,3))
            call swapr(z(2), z(3))
        end if
        if (z(1) > z(2)) then
            call swapr(x(1,1), x(1,2))
            call swapr(x(2,1), x(2,2))
            call swapr(z(1), z(2))
        end if
    end subroutine sort3r

    subroutine swapr(a, b)
        real(8) :: a, b, t
        t = a
        a = b
        b = t
    end subroutine swapr

    subroutine swapi(a, b)
        integer :: a, b, t
        t = a
        a = b
        b = t
    end subroutine

end module sort_n_swap
