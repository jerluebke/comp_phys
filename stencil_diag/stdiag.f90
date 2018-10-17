module diag_from_stencil_module
    contains
        function stdiag(n, s, err)
            implicit none
            ! dimension of returned stdiag
            integer(4), intent(in) :: n
            ! stencil vector
            real(8), dimension(:), intent(in) :: s
            ! errornumber
            integer(4), intent(out) :: err

            ! result
            real(8), dimension(n, n) :: stdiag
            real(8), dimension(n*n) :: tmp
            ! size of s, middle index, size of tmp
            integer(4) :: l, m, nn
            ! loop indices
            integer(4) :: i, k

            stdiag = 0
            l = size(s)
            ! checking conditions
            if (mod(l, 2) == 0) then
                err = -1
                return
            else if (n <= l/2) then
                err = -2
                return
            end if
            err = 0
            ! if s contains one element, return diag(s(1))
            if (l == 1) then
                do i=1, n
                    stdiag(i, i) = s(1)
                end do
                return
            end if

            tmp = 0
            nn = size(tmp)
            m = (l-1) / 2
            k = 0
            ! fill upper and lower elements of stdiag with slices of s
            do i = 1, n*m, n
                tmp(i : i+m+k) = s(m+1-k : l)           ! upper 
                tmp(nn-i-(m+k) : nn-i) = s(0 : m+1+k)   ! lower
                k = k + 1       ! extend slice of s by 1 for next iter
            end do
            ! fill rest of stdiag with copies of s along the diagonal
            do i = n*m+1, n*(n-m), n+1
                tmp(i:i+l) = s
            end do

            stdiag = reshape(tmp, [n, n])
        end function stdiag

end module diag_from_stencil_module


program test
    use diag_from_stencil_module, only: stdiag

    implicit none

    character(len=16) :: fmtstr
    integer(4) :: n
    integer(4) :: err
    real(8), dimension(5) :: s
    real(8), dimension(10, 10) :: a
    real(8), dimension(1) :: e
    real(8), dimension(10, 10) :: b
    n = 10

    write(fmtstr, '(A, I4, A, F3.1, A)') '(', n, '(2X,F', 4.0, '))'
    fmtstr = '(10(2X, F4.0))'
    s = [-1, -8, 0, 8, 1]
    a = stdiag(n, s, err)
    e = [6]
    b = stdiag(n, e, err)
    write(*, *) 'a = '
    write(*, fmtstr) a
    write(*, *) ''
    write(*, *) 'b = '
    write(*, fmtstr) b
end program
