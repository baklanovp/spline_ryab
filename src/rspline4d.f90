module rspline4d
    use kinds,   only: dp, alloc1d

    implicit none
    private

    public :: spline4d_type

    character(len=*), parameter  :: mdl_name = 'rspline4d'
    integer, parameter, public :: p_dim_4d = 4 ! dimension of the tables

    type spline4d_type
        integer :: ndim = p_dim_4d   ! dimension of the tables
        real(dp) :: V_4d(-1:2,p_dim_4d)  ! local volume
        real(dp) :: f_4d(-1:2,-1:2,-1:2,-1:2)  ! function values in the local volume
        integer, dimension(p_dim_4d) :: n_cur 
        real(dp), dimension(:), allocatable :: x_tab
        real(dp), dimension(:), allocatable :: y_tab
        real(dp), dimension(:), allocatable :: z_tab
        real(dp), dimension(:), allocatable :: w_tab
        real(dp), dimension(:,:,:,:), allocatable :: funcTab
        integer :: n_x, n_y, n_z, n_w
        contains
            procedure :: init => spline4d_init
            procedure :: value => spline4d_value
            procedure :: check_value => spline4d_check_value
            procedure :: destroy => spline4d_destroy
    end type spline4d_type


    contains


    subroutine spline4d_init(this, x_tab, y_tab, z_tab, w_tab, funcTab) 
        class(spline4d_type), intent(inout)  :: this
        real(8), dimension(:), intent(in) :: x_tab, y_tab, z_tab, w_tab
        real(8), dimension(:,:,:,:) :: funcTab
        character(len=*), parameter ::  subrtn_name = 'spline4d_init', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        integer :: ierr               
        ! TODO add checking

        this%n_x = size(x_tab)
        this%n_y = size(y_tab)
        this%n_z = size(z_tab)
        this%n_w = size(w_tab)

        call alloc1d('x_tab', this%n_x, this%x_tab, path=fullPathSubrtn)
        this%x_tab = x_tab
        call alloc1d('y_tab', this%n_y, this%y_tab, path=fullPathSubrtn)
        this%y_tab = y_tab
        call alloc1d('z_tab', this%n_z, this%z_tab, path=fullPathSubrtn)
        this%z_tab = z_tab
        call alloc1d('w_tab', this%n_w, this%w_tab, path=fullPathSubrtn)
        this%w_tab = w_tab

        allocate(this%funcTab(this%n_x,this%n_y,this%n_z,this%n_w), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a, 4i4)') fullPathSubrtn, ' Not enough memory for funcTab where n_x,n_y,n_z,n_w =', &
                this%n_x, this%n_y, this%n_z, this%n_w;
            error stop 666;
        endif
        this%funcTab = funcTab

    end subroutine spline4d_init


    pure subroutine spline4d_destroy(this)
        implicit none
        class(spline4d_type), intent(inout)  :: this
        this%n_x    = 0
        this%n_y    = 0
        this%n_z    = 0
        this%n_w    = 0
        if (allocated(this%x_tab))      deallocate(this%x_tab)
        if (allocated(this%y_tab))      deallocate(this%y_tab)
        if (allocated(this%z_tab))      deallocate(this%z_tab)
        if (allocated(this%w_tab))      deallocate(this%w_tab)
        if (allocated(this%funcTab))      deallocate(this%funcTab)
    end subroutine spline4d_destroy


    subroutine spline4d_check_value(this, point, ierr)
        class(spline4d_type), intent(inout)  :: this
        real(8), intent(in) :: point(p_dim_4d)
        integer, intent(in) :: ierr

        associate(n_x=>this%n_x, n_y=>this%n_y, n_z=>this%n_z, n_w=>this%n_w)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab, z_tab=>this%z_tab, w_tab=>this%w_tab)
        if(ierr == 10)then
            print*,'X-variables are out of range'
            print*,'X',x_tab(2),point(1),x_tab(n_x-1)
            read*
            stop
        end if

        if(ierr == 20)then
            print*,'Y-variables  are out of range'
            print*,'Y',y_tab(2),point(2),y_tab(n_y-1)
            read*
            stop
        end if

        if(ierr == 30)then
            print*,'Z-variables 3 are out of range'
            print*,'Z',z_tab(2),point(3),z_tab(n_z-1)
            read*
            stop
        end if

        if(ierr == 40)then
            print*,'W-variables 4 are out of range'
            print*,'w_tab(2)= ',w_tab(2),' point(4)= ',point(4),' w_tab(n_z-1)= ',w_tab(n_w-1)
            read*
            stop
        end if

        endassociate
        endassociate
    end subroutine spline4d_check_value


    real(dp) function spline4d_value(this, point, ierr) result(res)
        ! use ryabmod
        implicit none
        class(spline4d_type), intent(inout)  :: this
        real(8), intent(in) :: point(p_dim_4d)
        integer, intent(out) :: ierr

        logical :: reload
        logical :: first_run=.true.

        ierr = 0

        associate(n_x=>this%n_x, n_y=>this%n_y, n_z=>this%n_z, n_w=>this%n_w)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab, z_tab=>this%z_tab, w_tab=>this%w_tab)
        associate(n_cur=>this%n_cur, V_4d=>this%V_4d, f_4d=>this%f_4d, funcTab=>this%funcTab)

        !------/checking if current position is not out of the table's ranges/------
        if(point(1) < x_tab(2) .or. point(1) > x_tab(n_x-1))then
            ierr = 10
            return            
        end if

        if(point(2) < y_tab(2) .or. point(2) > y_tab(n_y-1))then
            ierr = 20
            return            
        end if

        if(point(3) < z_tab(2) .or. point(3) > z_tab(n_z-1))then
            ierr = 30
            return            
        end if

        if(point(4) < w_tab(2) .or. point(4) > w_tab(n_w-1))then
            ierr = 40
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

        if(reload .or. point(3) < z_tab(n_cur(3)) .or. point(3) > z_tab(n_cur(3)+1))then
            reload=.true.
            n_cur(3)=minloc(point(3)-z_tab,mask=point(3)-z_tab >= 0.d0,dim=1)
        end if

        if(reload .or. point(4) < w_tab(n_cur(4)) .or. point(4) > w_tab(n_cur(4)+1))then
            reload=.true.
            n_cur(4) = minloc(point(4)-w_tab,mask=point(4)-w_tab >= 0.d0,dim=1)
        end if        

        if(reload)then
            V_4d(:,1) = x_tab(n_cur(1)-1:n_cur(1)+2)
            V_4d(:,2) = y_tab(n_cur(2)-1:n_cur(2)+2)
            V_4d(:,3) = z_tab(n_cur(3)-1:n_cur(3)+2)
            V_4d(:,4) = w_tab(n_cur(4)-1:n_cur(4)+2)
            f_4d = funcTab( n_cur(1)-1:n_cur(1)+2, n_cur(2)-1:n_cur(2)+2, n_cur(3)-1:n_cur(3)+2, n_cur(4)-1:n_cur(4)+2 )
        end if

        call ryab_4d(V_4d, f_4d, point,res)

        endassociate
        endassociate
        endassociate
        return
    endfunction spline4d_value


    pure subroutine ryab_4d(V_4d, f_4d, Y, RES)
        ! To calculate two-dimensional Ryabenkii spline with P=2,s=1
        implicit none
        real(dp), dimension(-1:2,p_dim_4d), intent(in) :: V_4d
        real(dp), dimension(-1:2,-1:2,-1:2,-1:2), intent(in) :: f_4d

        real(8), dimension(:), intent(in):: Y
        real(8), intent(out) :: res
        real(dp) :: Q(0:3,size(Y))
        real(dp) :: T(size(Y))
        INTEGER I,J,K, l
        integer, dimension(size(Y)) :: VIN, VBASE
        !---------------------------------
        Q(0,:) = 1.d0

        Q(1,:) = Y-V_4d(-1,:)

        Q(2,:) = 0.5d0*(Y-V_4d(-1,:)) * (Y-V_4d(0,:))

        T = (Y-V_4d(0,:)) / (V_4d(1,:)-V_4d(0,:))

        Q(3,:) = 0.5d0*(V_4d(1,:)-V_4d(0,:))**2 * (V_4d(2,:)-V_4d(-1,:)) * T**3 * (T-1) * (1-2.d0/3.d0*T)

        VBASE = -1
        RES = 0.D0

        do I=0,3
            VIN(1)=I
            do J=0,3
                VIN(2)=J
                do K=0,3
                    VIN(3)=K
                    do L=0,3
                        VIN(4)=L
                        RES = RES+Q(I,1)*Q(J,2)*Q(K,3)*Q(L,4)*DELTA4D(V_4d, f_4d, VIN,VBASE)
                    enddo
                enddo
            enddo
        enddo
    end subroutine ryab_4d
    !=======================================================
    !*******************************************************

    pure recursive function DELTA4D(V_4d, f_4d, VIN,VBASE) RESULT(RES)
        ! To calculate 4d Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        IMPLICIT NONE
        real(dp), dimension(-1:2,p_dim_4d), intent(in) :: V_4d
        real(dp), dimension(-1:2,-1:2,-1:2,-1:2), intent(in) :: f_4d
        INTEGER, dimension(p_dim_4d), intent(in) :: VIN, VBASE

        real(dp) :: RES
        INTEGER K
        INTEGER V1(p_dim_4d),V2(p_dim_4d)
        !--------------------------------------

        K = MAXLOC(VIN,DIM=1)
        IF(VIN(K).EQ.0)THEN
            RES = f_4d(VBASE(1),VBASE(2),VBASE(3),VBASE(4))
            RETURN
        END IF

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        RES = VIN(K)*(DELTA4D(V_4d, f_4d, V1,V2)-DELTA4D(V_4d, f_4d, V1,VBASE)) / (V_4d(VBASE(K)+VIN(K),K)-V_4d(VBASE(K),K))
        return
    end function


end module rspline4d

