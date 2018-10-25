module n_bodies
    implicit none
contains
    function solve_all(pos0, pos1, m, dt, dist, steps) result(traj)
        real(8), dimension(:,:), intent(in) :: pos0, pos1
        real(8), dimension(size(pos0,2)), intent(inout) :: m
        real(8), intent(in) :: dt, dist
        integer, intent(in) :: steps
        real(8), dimension(size(pos0,1),size(pos0,2),steps) :: traj
        real(8), parameter :: eps = epsilon(1e0)
        integer :: i, n
        n = size(pos0, 2)
        traj(:,:,1) = pos0
        traj(:,:,2) = pos1
        do i=3,steps
            traj(:,:,i) = 2d0*traj(:,:,i-1)-traj(:,:,i-2)+dt*f(traj(:,:,i-1))
        end do

    contains
        function f(x) result(ddx)
            real(8), dimension(:,:), intent(in) :: x
            real(8), dimension(size(x,1),size(x,2)) :: ddx
            real(8), dimension(size(x,1),size(x,2)-1) :: r
            real(8), dimension(size(x,2)-1) :: r_sq
            integer :: i, j, k
            do i=1,n
                if (isnan(x(1,i))) cycle
                do j=1,i-1
                    r(:,j) = x(:,j) - x(:,i)
                    r_sq(j) = sum(r(:,j)**2)
                    call collision(r_sq(j), m(j)+m(i), m(i), x(1,j))
                end do
                do j=i,n-1
                    r(:,j) = x(:,j+1) - x(:,i)
                    r_sq(j) = sum(r(:,j)**2)
                    call collision(r_sq(j), m(j)+m(i), m(i), x(1,j+1))
                end do
                ddx(:,i) = sum([(m(k)*r(:,k)/(r_sq(k)+eps)**(3d0/2d0), k=1,n-1)])
            end do
        end function

        subroutine collision(r, m, mi, xj)
            real(8) :: r, m, mi, xj
            real(8) :: i = -1.
            if (r < m*dist) then
                mi = m
                xj = sqrt(i)
            end if
        end subroutine
    end function
end module n_bodies


program test
    use n_bodies, only: solve_all
    implicit none
    integer, parameter :: steps = 100
    real(8), dimension(2,3) :: q0, q1
    real(8), dimension(3) :: m
    real(8), dimension(2,3,steps) :: traj
    q0 = reshape([-1d0, 4.5d0, 0d0, -4d0, -4d0, 8d0], [2, 3])
    q1 = q0 + reshape([.05d0, 0d0, -.05d0, 0d0, .02d0, -.02d0], [2, 3])
    m = [4.5d0, 2.8d0, 3.9d0]
    100 format(2(3(3X, F6.3)/))
    traj = solve_all(q0, q1, m, .05d0, .1d0, steps)
end program
