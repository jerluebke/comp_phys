program heat
    implicit none

    interface
	subroutine outputVtk(u,dimX,nproc,me,outputDirectory)
	    implicit none
	    integer :: dimX(2),nproc(2),me(2)
	    real(kind=8), dimension(:,:) :: u
	    integer :: i,mx,my,mz,nx,ny,nz
	    integer :: ib(2),ie(2)
	    integer, save :: ic=1
	    character(len=256) :: outputDirectory
	end subroutine

    end interface

    include "mpif.h"

    integer, parameter :: mx=128, my=128
    real(kind=8) :: u(0:mx+1,0:my+1)
    real(kind=8) :: u_tmp(0:mx+1,0:my+1)
    integer :: i,j
    real(kind=8) :: x,y,xb,xe,yb,ye,xb_loc,yb_loc,dx,dy, &
		    pi,t,te,dt,output_interval

    integer :: nproc(2) = [2,2], dimX(2) = [mx,my], &
	ierror, iprocs, myRank, &
	mpi_comm2d, mpi_coords(2), mpi_neighbours(4)
    logical :: periods(2), reorder

    character(len=256) :: outputDirectory = "heat_output"

    call MPI_Init(ierror)
    call MPI_Comm_Size(MPI_COMM_WORLD, iprocs, ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, myRank, ierror)

    periods = .true.
    reorder = .false.

    call MPI_Cart_create(MPI_COMM_WORLD, 2, nproc, periods, reorder, mpi_comm2d, ierror)
    call MPI_Cart_coords(mpi_comm2d, myRank, 2, mpi_coords, ierror)
    call MPI_Cart_shift(mpi_comm2d, 0, 1, mpi_neighbours(1), mpi_neighbours(2), ierror)
    call MPI_Cart_shift(mpi_comm2d, 1, 1, mpi_neighbours(3), mpi_neighbours(4), ierror)

    pi = 4.*atan(1.d0)

    ! setup
    xb = 0.d0
    xe = 1d0
    yb = 0.d0
    ye = 1d0
    dx = (xe-xb)/(mx*nproc(1))
    dy = (ye-yb)/(mx*nproc(2))

    xb_loc = xb + dx*mx*mpi_coords(1)
    yb_loc = yb + dy*my*mpi_coords(2)

    dt = 0.2*min(dx,dy)**2

    t = 0.
    te = .1

    ! initial condition
    u=0.
    do j = -1,my+1
        y = yb_loc + j*dy
        do i = -1,mx+1
            x = xb_loc + i*dx
            if ((x-0.5)**2+(y-0.5)**2 .lt. 0.2) then
                u(i,j)= 1
            end if
        end do
    end do

    ! create directories for output
    if (mpi_coords(1) + mpi_coords(2) == 0) then
	call execute_command_line('mkdir -p ' // trim(outputDirectory) // '/vti')
    end if

    ! output of initial condition
    call outputVtk(u,dimX,nproc,mpi_coords,outputDirectory)

    output_interval = te/100.
    if (mpi_coords(1) + mpi_coords(2) == 0) write(*,*) 'dt = ', dt

    ! Euler steps
    do while (t < te)
        ! heat equation
        u_tmp(1:mx,1:my) = u(1:mx,1:my) &
                +dt*((u(2:mx+1,1:my)-2.*u(1:mx,1:my) &
                +u(0:mx-1,1:my))/dx**2 &
                +(u(1:mx,2:my+1)-2.*u(1:mx,1:my) &
                +u(1:mx,0:my-1))/dy**2)
	call boundary(u_tmp,dimX,mpi_neighbours,nproc,mpi_coords)
        u = u_tmp
        t = t+dt

	if (mpi_coords(1) + mpi_coords(2) == 0) write(*,*) 'time = ',t

        if (mod(t,output_interval) < dt) then
	    call outputVtk(u,dimX,nproc,mpi_coords,outputDirectory)
        end if
    end do

    call MPI_Finalize()

end program heat


subroutine boundary(u,dimX,mpi_neighbours,nproc,me)
    implicit none

    include "mpif.h"

    integer :: dimX(2), mpi_neighbours(4), nproc(2), me(2)
    real(kind=8) :: u(0:dimX(1)+1,0:dimX(2)+1)

    integer :: buffer_size,ierror

    real(kind=8), dimension(:), allocatable :: mpi_buffer_in_xl, mpi_buffer_in_xu, &
	mpi_buffer_in_yl, mpi_buffer_in_yu, &
	mpi_buffer_out_xl, mpi_buffer_out_xu, &
	mpi_buffer_out_yl, mpi_buffer_out_yu

    integer :: mpi_requests(8), mpi_status(MPI_STATUS_SIZE, 8), mpi_tags(4), mpi_err

    call MPI_Barrier(MPI_COMM_WORLD,ierror)

    allocate(mpi_buffer_in_xl(0:dimX(2)+1))
    allocate(mpi_buffer_in_xu(0:dimX(2)+1))
    allocate(mpi_buffer_in_yl(0:dimX(1)+1))
    allocate(mpi_buffer_in_yu(0:dimX(1)+1))

    allocate(mpi_buffer_out_xl(0:dimX(2)+1))
    allocate(mpi_buffer_out_xu(0:dimX(2)+1))
    allocate(mpi_buffer_out_yl(0:dimX(1)+1))
    allocate(mpi_buffer_out_yu(0:dimX(1)+1))

    mpi_buffer_out_xl(:) = u(1,:)
    mpi_buffer_out_xu(:) = u(dimX(1),:)
    mpi_buffer_out_yl(:) = !
    mpi_buffer_out_yu(:) = !

    mpi_requests = MPI_REQUEST_NULL
    mpi_tags = (/1, 2, 3, 4/)

    ! x lower boundary receive
    buffer_size = size(mpi_buffer_in_xl)
    call MPI_Irecv(mpi_buffer_in_xl, buffer_size, MPI_REAL8,&
	mpi_neighbours(1), mpi_tags(1), MPI_COMM_WORLD, &
	mpi_requests(1), mpi_err)

    ! x upper boundary receive
    buffer_size = size(mpi_buffer_in_xu)
    call MPI_Irecv(mpi_buffer_in_xu, buffer_size, MPI_REAL8, &
	mpi_neighbours(2), mpi_tags(2), MPI_COMM_WORLD, &
	mpi_requests(2), mpi_err)

    ! y lower boundary receive
    buffer_size = size(mpi_buffer_in_yl)
    call MPI_Irecv(mpi_buffer_in_yl, buffer_size, MPI_REAL8, &
	mpi_neighbours(3), mpi_tags(3), MPI_COMM_WORLD, &
	mpi_requests(3), mpi_err)

    ! y upper boundary receive
    buffer_size = size(mpi_buffer_in_yu)
    call MPI_Irecv(mpi_buffer_in_yu, buffer_size, MPI_REAL8, &
	mpi_neighbours(4), mpi_tags(4), MPI_COMM_WORLD, &
	mpi_requests(4), mpi_err)

    ! x lower boundary send
    buffer_size = size(mpi_buffer_out_xl)
    call MPI_Isend(mpi_buffer_out_xl, buffer_size, MPI_REAL8,  ) ! TODO

    ! x upper boundary send
    buffer_size = size(mpi_buffer_out_xu)
    call MPI_Isend( ) ! TODO

    ! y lower boundary send
    buffer_size = size(mpi_buffer_out_yl)
    call MPI_Isend( ) ! TODO

    ! y upper boundary send
    buffer_size = size(mpi_buffer_out_yu)
    call MPI_Isend( ) ! TODO

    call MPI_Waitall(8, mpi_requests, mpi_status, mpi_err)

    ! TODO
    u(0,:) = mpi_buffer_in_xl(:)
     = mpi_buffer_in_xu(:)
     = mpi_buffer_in_yl(:)
     = mpi_buffer_in_yu(:)

    call MPI_Barrier(MPI_COMM_WORLD,ierror)

    deallocate(mpi_buffer_in_xl)
    deallocate(mpi_buffer_in_xu)
    deallocate(mpi_buffer_in_yl)
    deallocate(mpi_buffer_in_yu)
				
    deallocate(mpi_buffer_out_xl)
    deallocate(mpi_buffer_out_xu)
    deallocate(mpi_buffer_out_yl)
    deallocate(mpi_buffer_out_yu)

end subroutine


subroutine outputVtk(u,dimX,nproc,me,outputDirectory)
    implicit none

    integer :: dimX(2),nproc(2),me(2)
    real(kind=8), dimension(:,:) :: u
    integer :: i,mx,my,mz,nx,ny,nz
    integer :: ib(2),ie(2)
    integer, save :: ic=1

    character(len=256) :: outputDirectory

    ib(:) = me(:)*(dimX(:))
    ie(:) = ib(:) + dimX(:)

    call write_vti_cells(sngl(u(1:dimX(1),1:dimX(2))), &
	me(1), me(2), ic, &
	trim(outputDirectory)//"/vti"//CHAR(0), &
	ib(1), ie(1), ib(2), ie(2))

    mx = dimX(1)
    my = dimX(2)
    nx = nproc(1)*mx
    ny = nproc(2)*my

    if (me(1) + me(2) == 0) then
	call write_pvti_cells(trim(outputDirectory)//CHAR(0), ic, nx, ny, mx, my)
    end if

    ic = ic+1
end subroutine

