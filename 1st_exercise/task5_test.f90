module test
    use ogpf
    use simplex, only: step
    implicit none

    ! procedure(interf_func), pointer :: fnptr => null()

    abstract interface
        real(8) pure function test_func(x)
            real(8), dimension(2), intent(in) :: x
        end function

        ! real(8) pure function interf_func(x)
        !     real(8), intent(in) :: x
        ! end function
    end interface

contains
    subroutine do_test(func, xdom, ydom, sim, tol, maxiter)
        ! function to test
        procedure(test_func) :: func
        ! domains
        real(8), dimension(2), intent(in) :: xdom, ydom
        ! simplex vertices to start with
        real(8), dimension(3), intent(inout) :: sim
        ! convergence tolerance
        real(8), intent(in) :: tol
        ! max iterations
        integer, intent(in) :: maxiter

        integer, parameter :: N = 100

        type(gpf) :: gp
        real(8), dimension(N) :: xl, yl
        real(8), allocatable :: x(:,:), y(:,:), z(:,:)

        ! fnptr => func

        xl = linspace(xdom(1), xdom(2), N)
        yl = linspace(ydom(1), ydom(2), N)
        call meshgrid(x, y, xl, yl)
        allocate( z(size(x, dim=1), size(x, dim=2)) )
        ! z = elem_func(x, y)
        z = func([x, y])

        call gp%contour(x, y, z)

    end subroutine do_test


    ! real(8) elemental &
    ! function elem_func(x, y) result(res)
    !     real(8), intent(in) :: x, y
    !     if (associated(fnptr)) res = fnptr([x, y])
    ! end function elem_func

end module test


program main
    use test, only: do_test
    implicit none

    real(8), parameter :: PI = 4*atan(1d0)
    real(8), parameter :: eps = epsilon(1d0)
    integer :: res


    call do_test(rastrigin, [-5d0, 5d0], [-5d0, 5d0], [1, 2, 3], eps, 100)


contains
    ! global minimum: f(0, 0) = 0
    ! domain: -5.12 <= x(i) <= 5.12
    real(8) pure &
    function rastrigin(x)
        real(8), dimension(2), intent(in) :: x
        real(8), parameter :: a = 10d0
        rastrigin = 2d0*A + sum(x**2 - A*cos(2d0*PI*x))
    end function rastrigin

    ! global minimum: f(0, 0) = 0
    ! domain: -5 <= x(i) <= 5
    real(8) pure &
    function ackley(x)
        real(8), dimension(2), intent(in) :: x
        ackley = -20d0 * exp( -0.2d0 * sqrt(sum(x**2)/2d0) ) &
                -exp( sum(cos(2d0*PI*x))/2d0 ) + exp(1d0) + 20d0
    end function ackley

    ! global minimum: f(1, 1) = 0
    ! domain: -INF <= x(i) <= INF
    real(8) pure &
    function rosenbrock(x)
        real(8), dimension(2), intent(in) :: x
        real(8), parameter :: a = 1d0, b = 100d0
        rosenbrock = (a - x(1))**2 + b*(x(2) - x(1)**2)**2
    end function rosenbrock

    ! global minimum: f(0, -1) = 3
    ! domain: -2 <= x(i) <= 2
    real(8) pure &
    function goldstein_prince(x)
        real(8), dimension(2), intent(in) :: x
        goldstein_prince = ( 1d0 + (sum(x)+1)**2                            &
            * (19d0 - 14d0*x(1) + 3d0*x(1)**2 - 14d0*x(2) + 6d0*x(1)*x(2)   &
                + 3d0*x(2)**2) )                                            & 
            * ( 30d0 + (2d0*x(1) - 3d0*x(2))**2                             &
            * (18d0 - 32d0*x(1) + 12d0*x(1)**2 + 48d0*x(2) - 36d0*x(1)*x(2) &
                + 27d0*x(2)**2) )
    end function goldstein_prince

    ! global minimum: f(-10, 1) = 0
    ! domain: -15 <= x <= -5, -3 <= y <= 3
    real(8) pure &
    function bukin_n6(x)
        real(8), dimension(2), intent(in) :: x
        bukin_n6 = 100d0*sqrt(abs(x(2)-0.01d0*x(1)**2)) + 0.01d0*abs(x(1)+10)
    end function bukin_n6

    ! global minimum: f(1, 1) = 0
    ! -10 <= x(i) <= 10
    real(8) pure &
    function levi_n13(x)
        real(8), dimension(2), intent(in) :: x
        levi_n13 = sin(3d0*PI*x(1))**2 + (x(1)-1d0)**2*(1d0+sin(3d0*PI*x(2))) &
            + (x(2)-1)**2*(1+sin(2*PI*x(2)))
    end function levi_n13

    ! global minimum:
    !   f(3, 2)                 \
    !   f(-2.805118, 3.131312)  |   = 0
    !   f(-3.779310, -3.283186) |
    !   f(3.584428, -1.848126)  /
    !
    ! domain: -5 <= x(i) <= 5
    real(8) pure &
    function himmelblau(x)
        real(8), dimension(2), intent(in) :: x
        himmelblau = (x(1)**2+x(2)-11d0)**2 + (x(1)+x(2)**2-7d0)**2
    end function himmelblau

    ! global minimum: f(PI, PI) = -1
    ! domain: -100 <= x(i) <= 100
    real(8) pure &
    function easom(x)
        real(8), dimension(2), intent(in) :: x
        easom = -product(cos(x)) * exp(-sum((x-PI)**2))
    end function easom

    ! global minimum:
    !     f(1.34941, -1.34941)  \
    !     f(1.34941, 1.34941)   |   = -2.06261
    !     f(-1.34941, 1.34941)  |
    !     f(-1.34941, -1.34941) /
    !
    ! domain: -10 <= x(i) <= 10
    real(8) pure &
    function cross_in_tray(x)
        real(8), dimension(2), intent(in) :: x
        cross_in_tray = -0.0001d0*( abs(product(sin(x)) &
            * exp( abs(100d0-sqrt(sum(x**2))/PI) ))+1d0)**(0.1d0)
    end function cross_in_tray

    ! global minimum: f(512, 404.2319) = -959.6407
    ! domain: -512 <= x(i) <= 512
    real(8) pure &
    function eggholder(x)
        real(8), dimension(2), intent(in) :: x
        eggholder = -(x(2)+47d0) * sin(sqrt(abs(x(1)/2d0+x(2)+47d0)))   &
            - x(1) * sin(sqrt(abs(x(1)-x(2)+47d0)))
    end function eggholder

    ! global minimum:
    !     f(8.05502, 9.66459)   \
    !     f(-8.05502, 9.66459)  |   = -19.2085
    !     f(8.05502, -9.66459)  |
    !     f(-8.05502, -9.66459) /
    !
    ! domain: -10 <= x(i) <= 10
    real(8) pure &
    function hoelder_table(x)
        real(8), dimension(2), intent(in) :: x
        hoelder_table = -abs(sin(x(1))*cos(x(2))*exp(abs(1-sqrt(sum(x**2))/PI)))
    end function hoelder_table

    ! global minimum: -78.33234 < f(-2.903534, -2.903534) < -78.33232
    ! domain: -5 <= x(i) <= 5
    real(8) pure &
    function styblinski_tang(x)
        real(8), dimension(2), intent(in) :: x
        styblinski_tang = sum(x**4 - 16d0*x**2 + 5d0*x) / 2d0
    end function styblinski_tang

end program
