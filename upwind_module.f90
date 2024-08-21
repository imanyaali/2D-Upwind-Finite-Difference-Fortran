module upwind_module
    implicit none
contains

    subroutine Upwind_Full(U, M, N, dt, dx, dy, n, Cp, Cn)
        implicit none
        integer, intent(in) :: M, N, n
        real(8), intent(in) :: dt, dx, dy
        real(8), dimension(:,:), intent(inout) :: U
        real(8), dimension(:), intent(in) :: Cp, Cn

        integer :: i, j, ip, ie, iw, in, is
        real(8) :: FX, FY

        ! Central Nodes
        do i = 2, M
            do j = 2, N
                ip = (j-1)*(M+1) + i
                ie = ip + 1
                iw = ip - 1
                in = ip + M + 1
                is = ip - M - 1
                call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
                call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
                U(ip, n+1) = U(ip, n) - (FX + FY)
            end do
        end do

        ! Up B.C.
        do i = 2, M
            ip = N*(M+1) + i
            ie = ip + 1
            iw = ip - 1
            in = (M+1) + i
            is = ip - M - 1
            call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
            call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
            U(ip, n+1) = U(ip, n) - (FX + FY)
        end do

        ! Right B.C.
        do j = 2, N
            ip = j*(M+1)
            ie = (j-1)*(M+1) + 2
            iw = ip - 1
            in = ip + M + 1
            is = ip - M - 1
            call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
            call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
            U(ip, n+1) = U(ip, n) - (FX + FY)
        end do

        ! Down B.C.
        do i = 2, M
            ip = i
            ie = ip + 1
            iw = ip - 1
            in = ip + M + 1
            is = (N-1)*(M+1) + i
            call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
            call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
            U(ip, n+1) = U(ip, n) - (FX + FY)
        end do

        ! Left B.C.
        do j = 2, N
            ip = (j-1)*(M+1) + 1
            ie = ip + 1
            iw = j*(M+1) - 1
            in = ip + M + 1
            is = ip - M - 1
            call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
            call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
            U(ip, n+1) = U(ip, n) - (FX + FY)
        end do

        ! Left & Down Corner
        i = 1
        j = 1
        ip = (j-1)*(M+1) + i
        ie = ip + 1
        iw = j*(M+1) - 1
        in = ip + M + 1
        is = (N-1)*(M+1) + i
        call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
        call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
        U(ip, n+1) = U(ip, n) - (FX + FY)

        ! Right & Down Corner
        i = M+1
        j = 1
        ip = (j-1)*(M+1) + i
        ie = (j-1)*(M+1) + 2
        iw = ip - 1
        in = ip + M + 1
        is = (N-1)*(M+1) + i
        call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
        call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
        U(ip, n+1) = U(ip, n) - (FX + FY)

        ! Left & Up Corner
        i = 1
        j = N+1
        ip = (j-1)*(M+1) + i
        ie = ip + 1
        iw = j*(M+1) - 1
        in = (M+1) + i
        is = ip - M - 1
        call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
        call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
        U(ip, n+1) = U(ip, n) - (FX + FY)

        ! Right & Up Corner
        i = M+1
        j = N+1
        ip = (j-1)*(M+1) + i
        ie = (j-1)*(M+1) + 2
        iw = ip - 1
        in = (M+1) + i
        is = ip - M - 1
        call Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
        call Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
        U(ip, n+1) = U(ip, n) - (FX + FY)

    end subroutine Upwind_Full

    subroutine Upwind_X(ip, ie, iw, U, dt, dx, Cp, Cn, n, FX)
        implicit none
        integer, intent(in) :: ip, ie, iw, n
        real(8), intent(in) :: dt, dx
        real(8), dimension(:,:), intent(in) :: U
        real(8), dimension(:), intent(in) :: Cp, Cn
        real(8), intent(out) :: FX

        ! Implement the upwind scheme for the X direction here
        FX = (Cp(ip)*U(ie, n) - Cn(ip)*U(iw, n))*dt/dx
    end subroutine Upwind_X

    subroutine Upwind_Y(ip, in, is, U, dt, dy, Cp, Cn, n, FY)
        implicit none
        integer, intent(in) :: ip, in, is, n
        real(8), intent(in) :: dt, dy
        real(8), dimension(:,:), intent(in) :: U
        real(8), dimension(:), intent(in) :: Cp, Cn
        real(8), intent(out) :: FY

        ! Implement the upwind scheme for the Y direction here
        FY = (Cp(ip)*U(in, n) - Cn(ip)*U(is, n))*dt/dy
    end subroutine Upwind_Y

end module upwind_module
