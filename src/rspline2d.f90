module rspline2d
    use kinds,   only: dp, alloc1d

    implicit none
    private

    public :: spline2d_type

    character(len=*), parameter  :: mdl_name = 'rspline2d'
    integer, parameter, public :: p_dim_2d = 2 ! dimension of the tables

    type spline2d_type
        integer :: ndim = p_dim_2d   ! dimension of the tables
        real(dp) :: V_2d(-1:2,p_dim_2d)  ! local volume
        real(dp) :: f_2d(-1:2,-1:2)  ! function values in the local volume
        integer, dimension(p_dim_2d) :: n_cur 
        real(dp), dimension(:), allocatable :: x_tab
        real(dp), dimension(:), allocatable :: y_tab
        real(dp), dimension(:,:), allocatable :: funcTab
        integer :: n_x, n_y
        contains
            procedure :: init => spline2d_init
            procedure :: value => spline2d_value
            procedure :: check_value => spline2d_check_value
            procedure :: destroy => spline2d_destroy
    end type spline2d_type


    contains


    subroutine spline2d_init(this, x_tab, y_tab, funcTab) 
        class(spline2d_type), intent(inout)  :: this
        real(8), dimension(:), intent(in) :: x_tab, y_tab
        real(8), dimension(:,:) :: funcTab
        character(len=*), parameter ::  subrtn_name = 'spline2d_init', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        integer :: ierr               
        ! TODO add checking

        this%n_x = size(x_tab)
        this%n_y = size(y_tab)

        call alloc1d('x_tab',  this%n_x, this%x_tab, path=fullPathSubrtn)
        this%x_tab = x_tab
        call alloc1d('y_tab',  this%n_y, this%y_tab, path=fullPathSubrtn)
        this%y_tab = y_tab

        allocate(this%funcTab(this%n_x,this%n_y), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a, 3i4)') fullPathSubrtn, &
               ' Not enough memory for funcTab where n_x,n_y =', this%n_x,this%n_y;
            error stop 666;
        endif
        this%funcTab = funcTab

    end subroutine spline2d_init


    pure subroutine spline2d_destroy(this)
        implicit none
        class(spline2d_type), intent(inout)  :: this
        this%n_x    = 0
        this%n_y    = 0
        if (allocated(this%x_tab))      deallocate(this%x_tab)
        if (allocated(this%y_tab))      deallocate(this%y_tab)
        if (allocated(this%funcTab))    deallocate(this%funcTab)
    end subroutine spline2d_destroy


    subroutine spline2d_check_value(this, point, ierr)
        class(spline2d_type), intent(inout)  :: this
        real(8), intent(in) :: point(p_dim_2d)
        integer, intent(in) :: ierr

        associate(n_x=>this%n_x, n_y=>this%n_y)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab)
        if(ierr == 10)then
            print*,'variables are out of range'
            print*,'X',x_tab(2),point(1),x_tab(n_x-1)
            read*
            stop
        end if

        if(ierr == 20)then
            print*,'variables are out of range'
            print*,'Y',y_tab(2),point(2),y_tab(n_y-1)
            read*
            stop
        end if
        endassociate
        endassociate
    end subroutine spline2d_check_value


    real(dp) function spline2d_value(this, point, ierr) result(res)
        ! use ryabmod
        implicit none
        class(spline2d_type), intent(inout)  :: this
        real(8), intent(in) :: point(p_dim_2d)
        integer, intent(out) :: ierr

        logical :: reload
        logical :: first_run=.true.

        ierr = 0

        associate(n_x=>this%n_x, n_y=>this%n_y)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab)
        associate(n_cur=>this%n_cur, V_2d=>this%V_2d, f_2d=>this%f_2d, funcTab=>this%funcTab)

        !------/checking if current position is not out of the table's ranges/------
        if(point(1) < x_tab(2) .or. point(1) > x_tab(n_x-1))then
            ierr = 10
            return            
        end if

        if(point(2) < y_tab(2) .or. point(2) > y_tab(n_y-1))then
            ierr = 20
            return            
        end if

        reload=.false.

        if(first_run)then
            first_run=.false.
            reload=.true.
        end if

        !---------/check if we're at the old xyz box (from the previous call)/----------------
        if(reload .or. point(1) < x_tab(n_cur(1)) .or. point(1) > x_tab(n_cur(1)+1))then
            reload=.true.
            n_cur(1)=minloc(point(1)-x_tab,mask=point(1)-x_tab >= 0.d0,dim=1)
        end if

        if(reload .or. point(2) < y_tab(n_cur(2)) .or. point(2) > y_tab(n_cur(2)+1))then
            reload=.true.
            n_cur(2)=minloc(point(2)-y_tab,mask=point(2)-y_tab >= 0.d0,dim=1)
        end if

        if(reload)then
            V_2d(:,1) = x_tab(n_cur(1)-1:n_cur(1)+2)
            V_2d(:,2) = y_tab(n_cur(2)-1:n_cur(2)+2)
            f_2d = funcTab( n_cur(1)-1:n_cur(1)+2, n_cur(2)-1:n_cur(2)+2 )
        end if

        call ryab_2d(V_2d, f_2d, point,res)

        endassociate
        endassociate
        endassociate
        return
    endfunction spline2d_value


    pure subroutine ryab_2d(V_2d, f_2d, Y, res)
        ! To calculate two-dimensional Ryabenkii spline with P=2,s=1
        implicit none
        real(dp), dimension(-1:2,p_dim_2d), intent(in) :: V_2d
        real(dp), dimension(-1:2,-1:2), intent(in) :: f_2d

        real(8), dimension(:), intent(in):: Y
        real(8), intent(out) :: res
        real(dp) :: Q(0:3,size(Y))
        real(dp) :: T(size(Y))
        INTEGER I,J,K,VIN(size(Y)),VBASE(size(Y))
        !---------------------------------
        Q(0,:) = 1.d0

        Q(1,:) = Y-V_2d(-1,:)

        Q(2,:) = 0.5d0*(Y-V_2d(-1,:)) * (Y-V_2d(0,:))

        T = (Y-V_2d(0,:)) / (V_2d(1,:)-V_2d(0,:))
        Q(3,:) = 0.5d0*(V_2d(1,:)-V_2d(0,:))**2 * (V_2d(2,:)-V_2d(-1,:)) * T**3 * (T-1) * (1-2.d0/3.d0*T)

        VBASE = -1
        res = 0.D0

        do i=0,3
            VIN(1)=I
            do j=0,3
                VIN(2)=J
                res = res+Q(I,1)*Q(J,2)*delta2(V_2d, f_2d, VIN,VBASE)
            end do
        end do        
    end subroutine ryab_2d
    !=======================================================
    !*******************************************************

    pure recursive function delta2(V_2d, f_2d, VIN,VBASE) RESULT(res)
        ! To calculate 2d Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        IMPLICIT NONE
        real(dp), dimension(-1:2,p_dim_2d), intent(in) :: V_2d
        real(dp), dimension(-1:2,-1:2), intent(in) :: f_2d
        INTEGER, dimension(p_dim_2d), intent(in) :: VIN, VBASE

        real(dp) :: res
        INTEGER K
        INTEGER V1(p_dim_2d),V2(p_dim_2d)
        !--------------------------------------

        K = MAXLOC(VIN,DIM=1)
        if(VIN(k) == 0)then
            res = f_2d(VBASE(1),VBASE(2))
            return
        end if

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        res = VIN(K)*(delta2(V_2d, f_2d, V1,V2)-delta2(V_2d, f_2d, V1,VBASE)) / (V_2d(VBASE(K)+VIN(K),K)-V_2d(VBASE(K),K))
        return
    end function

end module rspline2d

