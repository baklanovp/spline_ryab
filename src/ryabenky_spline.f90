module ryabenky_spline
    implicit none
    private
    public :: bilinear_interp, trilinear_interp, quadlinear_interp

    contains

    ! ====== Универсальный бинарный поиск индекса в неравномерной сетке ======
    pure integer function search_index(x0, grid, n)
        implicit none
        real(8), intent(in) :: x0, grid(:)
        integer, intent(in) :: n

        integer :: low, high, mid

        low = 1
        high = n - 1

        ! Обработка выхода за границы
        if (x0 <= grid(1)) then
            search_index = 1
            return
        end if

        if (x0 >= grid(n)) then
            search_index = n - 1
            return
        end if

        ! Бинарный поиск
        do while (high > low + 1)
            mid = (low + high) / 2
            if (x0 >= grid(mid)) then
            low = mid
            else
            high = mid
            end if
        end do

        search_index = low
    end function search_index


    ! ====== Билинейная интерполяция (2D) ======
    pure real(8) function bilinear_interp(x0, y0, x, y, f, nx, ny)
        implicit none
        real(8), intent(in) :: x0, y0, x(:), y(:), f(:,:)
        integer, intent(in) :: nx, ny
        integer :: i, j
        real(8) :: x1, x2, y1, y2, Q11, Q12, Q21, Q22, denom

        i = max(1, min(nx-1, search_index(x0, x, nx)))
        j = max(1, min(ny-1, search_index(y0, y, ny)))

        x1 = x(i)
        x2 = x(i+1)
        y1 = y(j)
        y2 = y(j+1)
        Q11 = f(j,i)
        Q21 = f(j,i+1)
        Q12 = f(j+1,i)
        Q22 = f(j+1,i+1)

        denom = (x2 - x1) * (y2 - y1)
        bilinear_interp = (Q11*(x2 - x0)*(y2 - y0) + &
            Q21*(x0 - x1)*(y2 - y0) + &
            Q12*(x2 - x0)*(y0 - y1) + &
            Q22*(x0 - x1)*(y0 - y1)) / denom
    end function bilinear_interp


    ! ====== Трилинейная интерполяция (3D) ======
    real(8) function trilinear_interp(x0, y0, z0, x, y, z, f, nx, ny, nz)
        implicit none
        real(8), intent(in) :: x0, y0, z0, x(:), y(:), z(:), f(:,:,:)
        integer, intent(in) :: nx, ny, nz
        integer :: i, j, k
        real(8) :: xd, yd, zd
        real(8) :: c00, c01, c10, c11, c0, c1

        i = max(1, min(nx-1, search_index(x0, x, nx)))
        j = max(1, min(ny-1, search_index(y0, y, ny)))
        k = max(1, min(nz-1, search_index(z0, z, nz)))

        xd = (x0 - x(i)) / (x(i+1) - x(i))
        yd = (y0 - y(j)) / (y(j+1) - y(j))
        zd = (z0 - z(k)) / (z(k+1) - z(k))

        c00 = f(k,j,i)*(1.0d0 - xd) + f(k,j,i+1)*xd
        c01 = f(k,j+1,i)*(1.0d0 - xd) + f(k,j+1,i+1)*xd
        c10 = f(k+1,j,i)*(1.0d0 - xd) + f(k+1,j,i+1)*xd
        c11 = f(k+1,j+1,i)*(1.0d0 - xd) + f(k+1,j+1,i+1)*xd

        c0 = c00*(1.0d0 - yd) + c01*yd
        c1 = c10*(1.0d0 - yd) + c11*yd
        trilinear_interp = c0*(1.0d0 - zd) + c1*zd
    end function trilinear_interp


    ! ====== Вспомогательная: локальная трилинейная интерполяция по блоку 2x2x2 ======
    real(8) function trilinear_sub(xd, yd, zd, block)
        implicit none
        real(8), intent(in) :: xd, yd, zd
        real(8), intent(in) :: block(2,2,2)
        real(8) :: c00, c01, c10, c11, c0, c1

        c00 = block(1,1,1)*(1.0d0 - xd) + block(1,1,2)*xd
        c01 = block(1,2,1)*(1.0d0 - xd) + block(1,2,2)*xd
        c10 = block(2,1,1)*(1.0d0 - xd) + block(2,1,2)*xd
        c11 = block(2,2,1)*(1.0d0 - xd) + block(2,2,2)*xd

        c0 = c00*(1.0d0 - yd) + c01*yd
        c1 = c10*(1.0d0 - yd) + c11*yd
        trilinear_sub = c0*(1.0d0 - zd) + c1*zd
    end function trilinear_sub


    ! ====== Четырёхлинейная интерполяция (4D) ======
    real(8) function quadlinear_interp(x0, y0, z0, w0, x, y, z, w, f, nx, ny, nz, nw)
        implicit none
        real(8), intent(in) :: x0, y0, z0, w0, x(:), y(:), z(:), w(:)
        real(8), intent(in) :: f(:,:,:,:)
        integer, intent(in) :: nx, ny, nz, nw
        integer :: i, j, k, l
        real(8) :: xd, yd, zd, wd
        real(8) :: c0, c1
        real(8), dimension(2,2,2) :: block

        i = max(1, min(nx-1, search_index(x0, x, nx)))
        j = max(1, min(ny-1, search_index(y0, y, ny)))
        k = max(1, min(nz-1, search_index(z0, z, nz)))
        l = max(1, min(nw-1, search_index(w0, w, nw)))

        xd = (x0 - x(i)) / (x(i+1) - x(i))
        yd = (y0 - y(j)) / (y(j+1) - y(j))
        zd = (z0 - z(k)) / (z(k+1) - z(k))
        wd = (w0 - w(l)) / (w(l+1) - w(l))

        block = f(l,k:k+1,j:j+1,i:i+1)
        c0 = trilinear_sub(xd, yd, zd, block)

        block = f(l+1,k:k+1,j:j+1,i:i+1)
        c1 = trilinear_sub(xd, yd, zd, block)

        quadlinear_interp = (1.0d0 - wd) * c0 + wd * c1
    end function quadlinear_interp

end module ryabenky_spline
