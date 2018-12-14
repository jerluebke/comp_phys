program heat
    implicit none
    integer, parameter :: mx=256, my=256
    real(kind=8) :: u(0:mx+1,0:my+1)
    real(kind=8) :: unp1(0:mx+1,0:my+1)
    integer :: i,j,ic
    real(kind=8) :: x,y,xb,xe,yb,ye,dx,dy
    real(kind=8) :: ax,ay,pi,time,dt
    ax = 1
    ay = 1
    pi = 4.*atan(1.d0)
    dt = 0.00025 ! critical value
    dt = 0.0002
    xb = 0.d0
    xe = 1
    yb = 0.d0
    ye = 1
    dx = (xe-xb)/mx
    dy = (ye-yb)/my
    dt = 0.2*min(dx,dy)**2

    u=0.
    do j = -1,my+1
        y = yb + j*dy
        do i = -1,mx+1
            x = xb + i*dx
            if ((x-0.5)**2+(y-0.5)**2 .lt. 0.2) then
                u(i,j)= 1
            end if
        end do
    end do
    time = 0.
    ic = 0
    ! call VTKout(u,mx,my)

    ! c --> simple Euler step
    do while (time .le. 0.1)
        ! c--> heat equation
        unp1(1:mx,1:my) = u(1:mx,1:my) &
                +dt*((u(2:mx+1,1:my)-2.*u(1:mx,1:my) &
                +u(0:mx-1,1:my))/dx**2 &
                +(u(1:mx,2:my+1)-2.*u(1:mx,1:my) &
                +u(1:mx,0:my-1))/dy**2)
        call boundary(unp1,mx,my)
        u = unp1
        time = time+dt
        ic = ic+1
        write(*,*) 'time = ',time
        if (ic > 10) then
            ic = 0
            ! call VTKout(u,mx,my)
        end if
    end do
    return
end program heat


subroutine boundary(u,mx,my)
    implicit none
    integer :: mx,my
    real(kind=8) :: u(0:mx+1,0:my+1)
    integer :: i,j
    do i=1,mx
        u(i,0   ) = u(i,my)
        u(i,my+1) = u(i,1 )
    end do
    do j=0,my+1
        u(0,j   ) = u(mx,j)
        u(mx+1,j) = u(1,j)
    end do
    return
end subroutine boundary
