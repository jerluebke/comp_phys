module integrator
    implicit none

contains
    function verlet_step(f, x0, x1, dt) result(xn)
        interface
            function f(x) result(ddx)
                real(8), dimension(:), intent(in) :: x
                real(8), dimension(size(x)) :: ddx
            end function
        end interface
        real(8), dimension(:), intent(in) :: x0, x1
        real(8), intent(in) :: dt
        real(8), dimension(size(x0)) :: xn
        xn = 2d0 * x1 - x0 + dt**2 * f(x1)
        return
    end function
end module integrator


program main
    use ogpf
    use integrator, only: verlet_step
    implicit none

    type(gpf) :: gp
    integer, parameter :: steps = 9733
    real(8), parameter :: mu = 2d0
    real(8), dimension(2) :: x0, x1
    real(8), dimension(2, steps) :: xn
    integer :: i
    real(8) :: dt, t
    x0 = [1d0, 1d0]
    x1 = [1.05d0, 1d0]
    dt = .05d0
    t = 0d0

    ! call gp%animation_start(1)
    call gp%axis([-2d0, 2d0, -2d0, 2d0])

    do i=1, steps
        xn(:,i) = verlet_step(f, x0, x1, dt)
        x0 = x1
        x1 = xn(:,i)
        t = t + dt
        ! if (mod(100, i) == 0) then
        !     call gp%plot(xn(1,:i), xn(2,:i), 'w lines lc "blue" lw 1')
        !     ! call gp%plot(xn(1,i:i), xn(2,i:i), 'w points ps 2 lc "red"')
        ! end if
    end do

    call gp%plot(xn(1,:), xn(2,:), 'w lines lc "blue" lw 1')
    ! call gp%animation_show()

contains
    function f(x) result(ddx)
        real(8), dimension(:), intent(in) :: x
        real(8), dimension(size(x)) :: ddx
        ddx = -mu * x / sum(x**2)**(3d0/2d0)
        return
    end function
end program
