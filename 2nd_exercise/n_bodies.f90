module n_bodies
    implicit none
contains
    function solve_all(pos0, pos1, m, dt, dist, steps) result(traj)
        real(8), dimension(:,:), intent(in) :: pos0, pos1
        real(8), dimension(size(pos0,2)), intent(in) :: m
        real(8), intent(in) :: dt, dist
        integer, intent(in) :: steps
        real(8), dimension(size(pos0,1),size(pos0,2),steps) :: traj
        integer :: i, n
        n = size(pos0, 2)
        traj(:,:,1) = pos0
        traj(:,:,2) = pos1
        do i=3,steps
            traj(:,:,i) = 2d0*traj(:,:,i-1)-traj(:,:,i-2)+dt*f(traj(:,:,1))
        end do

    contains
        function f(x) result(ddx)
            real(8), dimension(:,:), intent(in) :: x
            real(8), dimension(size(x,1),size(x,2)) :: ddx
            real(8), dimension(size(x,1),size(x,2)-1) :: r
            real(8), dimension(size(x,2)-1) :: r_sq
            integer :: i, j
            do i=1,n
                if (isnan(x(1,i))) cycle
                do j=1,i-1
                    r(:,j) = x(:,j) - x(:,i)
                    r_sq(j) = sum(r(:,j)**2)
                    call collision(r_sq(j), m(j)+m(i), m(i), x(1,j))
                end do
                do j=i+1,n
                    r(:,j) = x(:,j) - x(:,i)
                    r_sq(j) = sum(r(:,j)**2)
                    call collision(r_sq(j), m(j)+m(i), m(i), x(1,j))
                end do
                ddx(:,i) = sum([(m(j)*r(:,j)/r_sq(j)**(3d0/2d0), j=1,n-1)])
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
