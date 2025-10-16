module rspline2d
    use kinds,   only: dp, alloc1d

    implicit none
    private

    public :: spline2d_type

    character(len=*), parameter  :: mdl_name = 'rspline2d'
    integer, parameter, public :: p_dim_2d = 2 ! dimension of the tables

    type spline2d_type
        logical :: first_run
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
        real(8), dimension(:,:), intent(in) :: funcTab
        character(len=*), parameter ::  subrtn_name = 'spline2d_init', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        integer :: ierr               
        ! TODO add checking
        this%first_run = .true.

        this%n_x = size(x_tab)
        this%n_y = size(y_tab)

        call alloc1d('x_tab',  this%n_x, this%x_tab, path=fullPathSubrtn)
        this%x_tab = x_tab
        call alloc1d('y_tab',  this%n_y, this%y_tab, path=fullPathSubrtn)
        this%y_tab = y_tab

        this%n_cur = 1 !  todo check?

        if (size(funcTab,1) /= this%n_x) then
            write(*, '(a,2(a,i4))') fullPathSubrtn, ' Size(x_tab)=', this%n_x,' is not equal  size(funcTab,1)= ', size(funcTab,1);
            error stop 666;
        endif

        if (size(funcTab,2) /= this%n_y) then
            write(*, '(a,2(a,i4))') fullPathSubrtn, ' Size(y_tab)=', this%n_y,' is not equal  size(funcTab,2)= ', size(funcTab,2);
            error stop 666;
        endif

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
        character(len=*), parameter ::  subrtn_name = 'spline2d_check_value', &
                            fullPathSubrtn = mdl_name//'.'//subrtn_name

        associate(n_x=>this%n_x, n_y=>this%n_y)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab)
        if(ierr == 10)then
            write(*,'(A,100(1pe12.4))') 'x_tab : ',x_tab
            print*,fullPathSubrtn//': variables are out of range'
            print*,'x_tab(2),point(1),x_tab(n_x-2): ',x_tab(2),point(1),x_tab(n_x-2)
            ! read*
            ! stop 888
        endif

        if(ierr == 20)then
            write(*,'(A,100(1pe12.4))') 'y_tab : ',y_tab
            print*, fullPathSubrtn//': variables are out of range'
            print*,'y_tab(2),point(2),y_tab(n_y-2): ',y_tab(2),point(2),y_tab(n_y-2)
            ! read*
            ! stop 888
        endif
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

        ierr = 0

        associate(n_x=>this%n_x, n_y=>this%n_y)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab)
        associate(n_cur=>this%n_cur, V_2d=>this%V_2d, f_2d=>this%f_2d, funcTab=>this%funcTab)

        !------/checking if current position is not out of the table's ranges/------
        if(point(1) < x_tab(2) .or. point(1) >= x_tab(n_x-1))then
            ierr = 10
            return            
        endif

        if(point(2) < y_tab(2) .or. point(2) >= y_tab(n_y-1))then
            ierr = 20
            return            
        endif

        reload=.false.

        if ( this%first_run ) then
            this%first_run=.false.
            reload=.true.
        endif

        !---------/check if we're at the old xyz box (from the previous call)/----------------
        if(reload .or. point(1) < x_tab(n_cur(1)) .or. point(1) > x_tab(n_cur(1)+1))then
            reload=.true.
            n_cur(1) = minloc(point(1)-x_tab,mask=point(1)-x_tab >= 0.d0,dim=1)
            ! n_cur(1) = findloc(point(1)-x_tab >= 0._dp, .TRUE., dim=1)
        endif

        if(reload .or. point(2) < y_tab(n_cur(2)) .or. point(2) > y_tab(n_cur(2)+1))then
            reload=.true.
            n_cur(2) = minloc(point(2)-y_tab,mask=point(2)-y_tab >= 0.d0,dim=1)
        endif

        if(reload)then
            V_2d(:,1) = x_tab(n_cur(1)-1:n_cur(1)+2)
            V_2d(:,2) = y_tab(n_cur(2)-1:n_cur(2)+2)
            f_2d = funcTab( n_cur(1)-1:n_cur(1)+2, n_cur(2)-1:n_cur(2)+2 )
        endif

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
        integer, dimension(size(Y)) :: VIN, VBASE
        integer :: I, J
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
            enddo
        enddo        
    end subroutine ryab_2d
    !=======================================================
    !*******************************************************

    pure recursive function delta2(V_2d, f_2d, VIN,VBASE) result(res)
        ! To calculate 2d Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        implicit none
        real(dp), dimension(-1:2,p_dim_2d), intent(in) :: V_2d
        real(dp), dimension(-1:2,-1:2), intent(in) :: f_2d
        integer, dimension(p_dim_2d), intent(in) :: VIN, VBASE

        real(dp) :: res
        integer :: K
        integer, dimension(p_dim_2d) ::  V1, V2

        K = maxloc(VIN,DIM=1)
        if(VIN(k) == 0) then
            res = f_2d(VBASE(1),VBASE(2))
            return
        endif

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        res = VIN(K)*(delta2(V_2d, f_2d, V1,V2)-delta2(V_2d, f_2d, V1,VBASE)) / (V_2d(VBASE(K)+VIN(K),K)-V_2d(VBASE(K),K))
        return
    end function

end module rspline2d

