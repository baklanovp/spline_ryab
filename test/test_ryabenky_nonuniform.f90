program test_ryabenky_nonuniform
   use ryabenky_spline
   implicit none

   integer, parameter :: nx=10, ny=10, nz=10, nw=10
   real(8), dimension(nx) :: x
   real(8), dimension(ny) :: y
   real(8), dimension(nz) :: z
   real(8), dimension(nw) :: w

   real(8), dimension(ny,nx) :: f2d
   real(8), dimension(nz,ny,nx) :: f3d
   real(8), dimension(nw,nz,ny,nx) :: f4d

   real(8) :: x0, y0, z0, w0
   real(8) :: val2d, val3d, val4d
   real(8) :: exact2d, exact3d, exact4d
   integer :: i, j, k, l

   ! ====== Построение НЕРАВНОМЕРНЫХ сеток ======
   do i = 1, nx
      x(i) = log(1.0d0 + 9.0d0 * (i - 1) / (nx - 1)) ! log-scale: [0, ~2.3]
   end do
   do j = 1, ny
      y(j) = (j - 1)**2 * 1.0d0 / (ny - 1)**2 ! квадратичная сетка: [0, 1]
   end do
   do k = 1, nz
      z(k) = sin((k - 1) * 3.1415926d0 / (2*(nz - 1))) ! синусная сетка: [0, 1]
   end do
   do l = 1, nw
      w(l) = sqrt((l - 1) * 1.0d0 / (nw - 1)) ! корневая сетка: [0, 1]
   end do

   ! ====== Значения функции ======
   do j = 1, ny
      do i = 1, nx
         f2d(j,i) = sin(x(i)) * cos(y(j))
      end do
   end do

   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            f3d(k,j,i) = sin(x(i)) * cos(y(j)) * exp(z(k))
         end do
      end do
   end do

   do l = 1, nw
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               f4d(l,k,j,i) = sin(x(i)) * cos(y(j)) * exp(z(k)) * sqrt(w(l) + 1.0d0)
            end do
         end do
      end do
   end do

   ! ====== Точка внутри области ======
   x0 = 0.75d0 ! лог(x+1) примерно 0.56
   y0 = 0.3d0 ! допустим, внутри
   z0 = 0.7d0
   w0 = 0.5d0

   ! ====== Интерполяция ======
   val2d = bilinear_interp(x0, y0, x, y, f2d, nx, ny)
   val3d = trilinear_interp(x0, y0, z0, x, y, z, f3d, nx, ny, nz)
   val4d = quadlinear_interp(x0, y0, z0, w0, x, y, z, w, f4d, nx, ny, nz, nw)

   ! ====== Точное значение функции ======
   exact2d = sin(x0) * cos(y0)
   exact3d = exact2d * exp(z0)
   exact4d = exact3d * sqrt(w0 + 1.0d0)

   ! ====== Вывод ======
   print *, "== НЕРАВНОМЕРНЫЕ СЕТКИ =="
   print *, "Точка: x=", x0, " y=", y0, " z=", z0, " w=", w0
   print *
   print *, "2D: Интерполяция = ", val2d, " | Точно = ", exact2d, " | Ошибка = ", abs(val2d - exact2d)
   print *, "3D: Интерполяция = ", val3d, " | Точно = ", exact3d, " | Ошибка = ", abs(val3d - exact3d)
   print *, "4D: Интерполяция = ", val4d, " | Точно = ", exact4d, " | Ошибка = ", abs(val4d - exact4d)
   
end program test_ryabenky_nonuniform
