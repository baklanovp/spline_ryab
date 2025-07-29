module rspline3d
    use kinds,   only: dp, alloc1d

    implicit none
    private

    public :: spline3d_type

    character(len=*), parameter  :: mdl_name = 'rspline3d'
    integer, parameter, public :: p_dim_3d = 3 ! dimension of the tables

    type spline3d_type
        integer :: ndim = p_dim_3d   ! dimension of the tables
        real(dp) :: V_3d(-1:2,p_dim_3d)  ! local volume
        real(dp) :: f_3d(-1:2,-1:2,-1:2)  ! function values in the local volume
        integer, dimension(p_dim_3d) :: n_cur 
        real(dp), dimension(:), allocatable :: x_tab
        real(dp), dimension(:), allocatable :: y_tab
        real(dp), dimension(:), allocatable :: z_tab
        real(dp), dimension(:,:,:), allocatable :: funcTab
        integer :: n_x, n_y, n_z
        contains
            procedure :: init => spline3d_init
            procedure :: value => spline3d_value
            procedure :: check_value => spline3d_check_value
            procedure :: destroy => spline3d_destroy
    end type spline3d_type


    contains


    subroutine spline3d_init(this, x_tab, y_tab, z_tab, funcTab) 
        class(spline3d_type), intent(inout)  :: this
        real(8), dimension(:), intent(in) :: x_tab, y_tab, z_tab
        real(8), dimension(:,:,:) :: funcTab
        character(len=*), parameter ::  subrtn_name = 'spline3d_init', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        integer :: ierr               
        ! TODO add checking

        this%n_x = size(x_tab)
        this%n_y = size(y_tab)
        this%n_z = size(z_tab)

        call alloc1d('x_tab',  this%n_x, this%x_tab, path=fullPathSubrtn)
        this%x_tab = x_tab
        call alloc1d('y_tab',  this%n_y, this%y_tab, path=fullPathSubrtn)
        this%y_tab = y_tab
        call alloc1d('z_tab',  this%n_z, this%z_tab, path=fullPathSubrtn)
        this%z_tab = z_tab

        allocate(this%funcTab(this%n_x,this%n_y,this%n_z), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a, 3i4)') fullPathSubrtn, &
               ' Not enough memory for funcTab where n_x,n_y,n_z =', this%n_x,this%n_y,this%n_z;
            error stop 666;
        endif
        this%funcTab = funcTab

    end subroutine spline3d_init


    pure subroutine spline3d_destroy(this)
        implicit none
        class(spline3d_type), intent(inout)  :: this
        this%n_x    = 0
        this%n_y    = 0
        this%n_z    = 0
        if (allocated(this%x_tab))      deallocate(this%x_tab)
        if (allocated(this%y_tab))      deallocate(this%y_tab)
        if (allocated(this%z_tab))      deallocate(this%z_tab)
        if (allocated(this%funcTab))    deallocate(this%funcTab)
    end subroutine spline3d_destroy


    subroutine spline3d_check_value(this, point, ierr)
        class(spline3d_type), intent(inout)  :: this
        real(8), intent(in) :: point(p_dim_3d)
        integer, intent(in) :: ierr

        associate(n_x=>this%n_x, n_y=>this%n_y, n_z=>this%n_z)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab, z_tab=>this%z_tab)
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

        if(ierr == 30)then
            print*,'variables are out of range'
            print*,'res',z_tab(2),point(3),z_tab(n_z-1)
            read*
            stop
        end if

        endassociate
        endassociate
    end subroutine spline3d_check_value


    real(dp) function spline3d_value(this, point, ierr) result(res)
        ! use ryabmod
        implicit none
        class(spline3d_type), intent(inout)  :: this
        real(8), intent(in) :: point(p_dim_3d)
        integer, intent(out) :: ierr

        logical :: reload
        logical :: first_run=.true.

        ierr = 0

        associate(n_x=>this%n_x, n_y=>this%n_y, n_z=>this%n_z)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab, z_tab=>this%z_tab)
        associate(n_cur=>this%n_cur, V_3d=>this%V_3d, f_3d=>this%f_3d, funcTab=>this%funcTab)

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

        if(reload)then
            V_3d(:,1) = x_tab(n_cur(1)-1:n_cur(1)+2)
            V_3d(:,2) = y_tab(n_cur(2)-1:n_cur(2)+2)
            V_3d(:,3) = z_tab(n_cur(3)-1:n_cur(3)+2)
            f_3d = funcTab( n_cur(1)-1:n_cur(1)+2, n_cur(2)-1:n_cur(2)+2, n_cur(3)-1:n_cur(3)+2 )
        end if

        call ryab_3d(V_3d, f_3d, point,res)

        endassociate
        endassociate
        endassociate
        return
    endfunction spline3d_value


    pure subroutine ryab_3d(V_3d, f_3d, Y, RES)
        ! To calculate two-dimensional Ryabenkii spline with P=2,s=1
        implicit none
        real(dp), dimension(-1:2,p_dim_3d), intent(in) :: V_3d
        real(dp), dimension(-1:2,-1:2,-1:2), intent(in) :: f_3d

        real(8), dimension(:), intent(in):: Y
        real(8), intent(out) :: res
        real(dp) :: Q(0:3,size(Y))
        real(dp) :: T(size(Y))
        INTEGER I,J,K,VIN(size(Y)),VBASE(size(Y))
        !---------------------------------
        Q(0,:) = 1.d0

        Q(1,:) = Y-V_3d(-1,:)

        Q(2,:) = 0.5d0*(Y-V_3d(-1,:)) * (Y-V_3d(0,:))

        T = (Y-V_3d(0,:)) / (V_3d(1,:)-V_3d(0,:))
        Q(3,:) = 0.5d0*(V_3d(1,:)-V_3d(0,:))**2 * (V_3d(2,:)-V_3d(-1,:)) * T**3 * (T-1) * (1-2.d0/3.d0*T)

        VBASE = -1
        RES = 0.D0

        DO I=0,3
            VIN(1)=I
            DO J=0,3
                VIN(2)=J
                DO K=0,3
                VIN(3)=K
                RES=RES+Q(I,1)*Q(J,2)*Q(K,3)*DELTA3(V_3d, f_3d, VIN,VBASE)
                END DO
            END DO
        END DO
    end subroutine ryab_3d
    !=======================================================
    !*******************************************************

    pure recursive function delta3(V_3d, f_3d, VIN,VBASE) RESULT(RES)
        ! To calculate 3D Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        IMPLICIT NONE
        real(dp), dimension(-1:2,p_dim_3d), intent(in) :: V_3d
        real(dp), dimension(-1:2,-1:2,-1:2), intent(in) :: f_3d
        INTEGER, dimension(p_dim_3d), intent(in) :: VIN, VBASE

        real(dp) :: RES
        INTEGER K
        INTEGER V1(p_dim_3d),V2(p_dim_3d)
        !--------------------------------------

        K = MAXLOC(VIN,DIM=1)
        IF(VIN(K).EQ.0)THEN
            RES = f_3d(VBASE(1),VBASE(2),VBASE(3))
            RETURN
        END IF

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        RES = VIN(K)*(DELTA3(V_3d, f_3d, V1,V2)-DELTA3(V_3d, f_3d, V1,VBASE)) / (V_3d(VBASE(K)+VIN(K),K)-V_3d(VBASE(K),K))
        return
    end function

end module rspline3d

