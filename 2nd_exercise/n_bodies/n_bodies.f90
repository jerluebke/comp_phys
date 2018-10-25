module n_bodies
    implicit none
contains
    function solve_all(q0, q1, m, dt, dist, steps) result(traj)
        real(8), dimension(:,:), intent(in)             :: q0, q1
        real(8), dimension(size(q0,2)), intent(inout)   :: m
        real(8), intent(in)                             :: dt, dist
        integer, intent(in)                             :: steps
        real(8), dimension(size(q0,1),size(q0,2),steps) :: traj
        real(8), dimension(size(q0,1),size(q0,2))       :: r
        real(8), dimension(size(q0,2))                  :: r_sq
        integer, dimension(size(q0,2))                  :: idx
        integer :: i, n
        n = size(q0, 2)
        idx = [(i,i=1,n)]
        traj(:,:,1) = q0
        traj(:,:,2) = q1
        do i=3,steps
            traj(:,:,i) = 2d0*traj(:,:,i-1)-traj(:,:,i-2)+dt**2*f(traj(:,:,i-1))
        end do

    contains
        function f(x) result(ddx)
            real(8), dimension(:,:), intent(in) :: x
            real(8), dimension(size(x,1),size(x,2)) :: ddx
            integer :: i, j, k
            ! do i=1,n
            !     do j=1,n
            !         if (j == i) then
            !             r(:,idx(j)) = 0
            !             cycle
            !         end if
            !         r(:,idx(j)) = x(:,idx(j)) - x(:,idx(i))
            !         r_sq(idx(j)) = sum(r(:,idx(j))**2)
            !         if (r_sq(idx(j)) < (m(idx(i))+m(idx(j))*dist)) then
            !             m(idx(i)) = m(idx(i))+m(idx(j))
            !             idx(j:) = idx(j:)+1
            !             n = n-1
            !         end if
            !     end do
            !     ddx(:,idx(i)) = sum([(m(idx(k))*r(:,idx(k))/ &
            !                         (r_sq(idx(k))+dist)**(3d0/2d0), k=1,n)])
            ! end do

            do i=1,n
                do j=1,n
                    r(:,j) = x(:,j) - x(:,i)
                    r_sq(j) = sum(r(:,j)**2)
                end do
                ddx(:,i) = sum([(m(k)*r(:,k)/(r_sq(k)+dist)**(3d0/2d0), k=1,n)])
            end do
        end function
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
    traj = solve_all(q0, q1, m, .05d0, .01d0, steps)
    print *, 'done'
end program
