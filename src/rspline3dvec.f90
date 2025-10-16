module rspline3dvec
    use kinds,   only: dp, alloc1d
    use array_expand,  only: expand_3dvec, expand_1d    

    implicit none
    private

    public :: spline3dvec_type

    character(len=*), parameter  :: mdl_name = 'rspline3dvec'
    integer, parameter, public :: p_dim_3d = 3 ! dimension of the tables

    real(dp), parameter :: p_val_max = HUGE ( 0.d0 ) / 2
    type spline3dvec_type
        integer :: ndim = p_dim_3d   ! dimension of the tables
        real(dp) :: V_3d(-1:2,p_dim_3d)  ! local volume
        real(dp), allocatable :: f_3d(:,:,:,:)  ! function values in the local volume
        integer, dimension(p_dim_3d) :: n_cur 
        real(dp), dimension(:), allocatable :: x_tab
        real(dp), dimension(:), allocatable :: y_tab
        real(dp), dimension(:), allocatable :: z_tab
        real(dp), dimension(:,:,:,:), allocatable :: funcTab
        integer :: n_vec
        integer :: n_x, n_y, n_z
        logical :: is_cache = .false.
        integer, dimension(:,:,:,:,:,:), allocatable :: cache_idx
        integer :: cache_length
        integer :: cache_pos
        integer :: cache_counter_reset
        real(dp), dimension(:,:), allocatable :: cache_delta3
        contains
            procedure :: init => spline3d_init_vec
            procedure :: value => spline3d_vec
            procedure :: check_value => spline3d_check_value
            procedure :: destroy => spline3d_destroy
    end type spline3dvec_type

    contains

    subroutine spline3d_init_vec(this, x_tab, y_tab, z_tab, funcTab, cache_length, is_expand) 
        class(spline3dvec_type), intent(inout)  :: this
        real(dp), dimension(:), intent(in) :: x_tab, y_tab, z_tab
        real(dp), dimension(:,:,:,:), intent(in) :: funcTab
        integer, optional, intent(in) :: cache_length
        logical, optional, intent(in) :: is_expand
        character(len=*), parameter ::  subrtn_name = 'spline3d_init_vec', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name
        integer :: ierr               
        logical :: is_cache, is_expand_

        is_expand_ = .false.
        if ( present(is_expand) ) then
            is_expand_ = is_expand
        endif
        
        is_cache = .false.
        if (present(cache_length)) is_cache = cache_length > 0

        this%n_vec = size(funcTab,1)        
        
        call spline3d_init_grid(this, x_tab, y_tab, z_tab, is_expand_) 
        
        if ( is_expand) then
            call expand_3dvec(funcTab, this%funcTab)
        else
            allocate(this%funcTab(this%n_vec, this%n_x, this%n_y, this%n_z), STAT=ierr);
            if (ierr /= 0) then
                write(*, '(2a,4i10)') fullPathSubrtn, &
                ' Not enough memory for funcTabExp where nvec, lx, ly, lz =', this%n_vec, this%n_x, this%n_y, this%n_z;
                error stop 666;
            endif                    
            ! Same as original if no expansion is needed
            ! write(*, '(2a, 4i4)') fullPathSubrtn, ' No expansion needed for funcTabExp where nvec, lx, ly, lz =', nvec, lx, ly, lz; 
            this%funcTab(:,:,:,:) = funcTab(:,:,:,:)
        endif

        ! Check sizes
        if (size(this%funcTab,2) /= this%n_x) then
            write(*, '(a,2(a,i4))') fullPathSubrtn, ' Size(x_tab)=', this%n_x,' is not equal  size(this%funcTab,2)= ', size(this%funcTab,2);
            error stop 666;
        endif
        if (size(this%funcTab,3) /= this%n_y) then
            write(*, '(a,2(a,i4))') fullPathSubrtn, ' Size(y_tab)=', this%n_y,' is not equal  size(this%funcTab,3)= ', size(this%funcTab,3);
            error stop 666;
        endif        
        if (size(this%funcTab,4) /= this%n_z) then
            write(*, '(a,2(a,i4))') fullPathSubrtn, ' Size(z_tab)=', this%n_z,' is not equal  size(this%funcTab,4)= ', size(this%funcTab,4);
            error stop 666;
        endif

        allocate(this%f_3d(this%n_vec,-1:2,-1:2,-1:2), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a, 4i4)') fullPathSubrtn, &
               ' Not enough memory for f_3d where n_vec =', this%n_vec;
            error stop 666;
        endif
        this%f_3d = 0.

        if ( is_cache ) then;
            call spline3d_init_cache(this, cache_length)
        endif

    end subroutine spline3d_init_vec


    subroutine spline3d_init_grid(this, x_tab, y_tab, z_tab, is_expand) 
        class(spline3dvec_type), intent(inout)  :: this
        real(dp), dimension(:), intent(in) :: x_tab, y_tab, z_tab
        logical, intent(in) :: is_expand
        character(len=*), parameter ::  subrtn_name = 'spline3d_init_grid', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        this%n_cur = 1
        if (is_expand) then
            call expand_1d('this%x_tab', x_tab, this%x_tab)
            call expand_1d('this%y_tab', y_tab, this%y_tab)
            call expand_1d('this%z_tab', z_tab, this%z_tab)
        else
            call alloc1d('x_tab', size(x_tab), this%x_tab, path=fullPathSubrtn)
            this%x_tab = x_tab
            call alloc1d('y_tab', size(y_tab), this%y_tab, path=fullPathSubrtn)
            this%y_tab = y_tab
            call alloc1d('z_tab', size(z_tab), this%z_tab, path=fullPathSubrtn)
            this%z_tab = z_tab
        endif

        this%n_x = size(this%x_tab)
        this%n_y = size(this%y_tab)
        this%n_z = size(this%z_tab)
        
    end subroutine spline3d_init_grid


    subroutine spline3d_init_cache(this, clength)
        class(spline3dvec_type), intent(inout)  :: this
        integer, intent(in) :: clength
        character(len=*), parameter ::  subrtn_name = 'spline3d_init_cache', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name
        integer :: ierr

        this%cache_length = clength
        this%is_cache = .true.;
        allocate(this%cache_delta3(this%n_vec,this%cache_length), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a,i5,i10)') fullPathSubrtn, &
            ' Not enough memory for cache_delta3 where this%n_vec,cache_length =', this%n_vec,this%cache_length;
            error stop 666;
        endif

        allocate(this%cache_idx(0:3,0:3,0:3,this%n_x,this%n_y,this%n_z), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a, 3i4)') fullPathSubrtn, &
            ' Not enough memory for cache_idx where 4*4*4* n_x,n_y,n_z =', this%n_x,this%n_y,this%n_z;
            error stop 666;
        endif

        this%cache_pos = 0
        this%cache_counter_reset = 0
        this%cache_idx = 0
        this%cache_delta3 = 0. !p_val_max
    endsubroutine spline3d_init_cache


    pure subroutine spline3d_destroy(this)
        implicit none
        class(spline3dvec_type), intent(inout)  :: this
        this%n_x    = 0
        this%n_y    = 0
        this%n_z    = 0
        if (allocated(this%x_tab))      deallocate(this%x_tab)
        if (allocated(this%y_tab))      deallocate(this%y_tab)
        if (allocated(this%z_tab))      deallocate(this%z_tab)
        if (allocated(this%funcTab))    deallocate(this%funcTab)
        if (allocated(this%cache_delta3)) deallocate(this%cache_delta3)
        if (allocated(this%cache_idx))    deallocate(this%cache_idx)
    end subroutine spline3d_destroy


    subroutine spline3d_check_value(this, point, ierr)
        class(spline3dvec_type), intent(inout)  :: this
        real(dp), intent(in) :: point(p_dim_3d)
        integer, intent(in) :: ierr
        character(len=*), parameter ::  subrtn_name = 'spline3d_check_value', &
                            fullPathSubrtn = mdl_name//'.'//subrtn_name
        associate(n_x=>this%n_x, n_y=>this%n_y, n_z=>this%n_z)
        associate(x_tab=>this%x_tab, y_tab=>this%y_tab, z_tab=>this%z_tab)
        if(ierr == 10)then
            write(*,'(A,100(1pe12.4))') 'x_tab : ',x_tab
            print*, fullPathSubrtn//': X variables are out of range'
            print*,'x_tab(2),point(1),x_tab(n_x-1): ',x_tab(2),point(1),x_tab(n_x-1)
            ! stop 880
        elseif(ierr == 20)then
            write(*,'(A,100(1pe12.4))') 'y_tab : ',y_tab
            print*,fullPathSubrtn//': Y variables are out of range'
            print*,'y_tab(2),point(2),y_tab(n_y-1): ',y_tab(2),point(2),y_tab(n_y-1)
            ! stop 882
        elseif(ierr == 30)then
            write(*,'(A,100(1pe12.4))') 'z_tab : ',z_tab
            print*,fullPathSubrtn//': Z variables are out of range'
            print*,'z_tab(2),point(3),z_tab(n_z-1): ',z_tab(2),point(3),z_tab(n_z-1)
            ! stop 884
        else
            print*,fullPathSubrtn//': unknown error'
            print*,'ierr= ',ierr
            ! stop 886
        end if

        endassociate
        endassociate
    end subroutine spline3d_check_value


    function spline3d_vec(this, point, ierr) result(res)
        ! use ryabmod
        implicit none
        class(spline3dvec_type), intent(inout)  :: this
        real(dp), intent(in) :: point(p_dim_3d)
        integer, intent(out) :: ierr
        real(dp), dimension(this%n_vec) :: res

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
            f_3d = funcTab(:, n_cur(1)-1:n_cur(1)+2, n_cur(2)-1:n_cur(2)+2, n_cur(3)-1:n_cur(3)+2 )
        end if

        if ( this%is_cache) then
            call ryab_3d_cache(this, V_3d, f_3d, point,res)
        else 
            call ryab_3d(this%n_vec,V_3d, f_3d, point,res)
        endif

        endassociate
        endassociate
        endassociate
        return
    endfunction spline3d_vec


    subroutine ryab_3d_cache(this, V_3d, f_3d, Y, RES)
        class(spline3dvec_type), intent(inout)  :: this
        ! To calculate two-dimensional Ryabenkii spline with P=2,s=1
        real(dp), dimension(-1:2,p_dim_3d), intent(in) :: V_3d
        real(dp), dimension(this%n_vec,-1:2,-1:2,-1:2), intent(in) :: f_3d

        real(dp), dimension(p_dim_3d), intent(in):: Y
        real(dp), dimension(this%n_vec), intent(out) :: res

        real(dp) :: Q(0:3,p_dim_3d)
        real(dp) :: T(p_dim_3d)
        INTEGER I,J,K,VIN(p_dim_3d),VBASE(p_dim_3d)
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
                    RES = RES + Q(I,1)*Q(J,2)*Q(K,3)*delta3cache(this, V_3d, f_3d, VIN,VBASE)
                END DO
            END DO
        END DO
    endsubroutine ryab_3d_cache



    pure subroutine ryab_3d(n_vec, V_3d, f_3d, Y, RES)
        ! To calculate two-dimensional Ryabenkii spline with P=2,s=1
        implicit none
        integer, intent(in) :: n_vec
        real(dp), dimension(-1:2,p_dim_3d), intent(in) :: V_3d
        real(dp), dimension(n_vec,-1:2,-1:2,-1:2), intent(in) :: f_3d

        real(dp), dimension(p_dim_3d), intent(in):: Y
        real(dp), dimension(n_vec), intent(out) :: res

        real(dp) :: Q(0:3,p_dim_3d)
        real(dp), dimension(p_dim_3d) :: T
        integer, dimension(p_dim_3d) :: VIN, VBASE
        integer :: i,j,k
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
                    RES = RES + Q(I,1)*Q(J,2)*Q(K,3)*delta3(n_vec, V_3d, f_3d, VIN,VBASE)
                END DO
            END DO
        END DO
    end subroutine ryab_3d
    !=======================================================
    !*******************************************************

    function delta3cache(this,V_3d, f_3d, VIN,VBASE) result(res)
        ! Cache delta3
        class(spline3dvec_type), intent(inout)  :: this
        real(dp), dimension(-1:2,p_dim_3d), intent(in) :: V_3d
        real(dp), dimension(this%n_vec,-1:2,-1:2,-1:2), intent(in) :: f_3d
        integer, dimension(p_dim_3d), intent(in) :: VIN, VBASE
        real(dp), dimension(this%n_vec) :: res
        integer :: idx

        associate(n_cur=>this%n_cur, n_vec=>this%n_vec)
        idx = this%cache_idx(VIN(1),VIN(2),VIN(3),n_cur(1), n_cur(2), n_cur(3))
        if ( idx > 0 ) then
            res(1:n_vec) = this%cache_delta3(1:n_vec,idx)
            ! write(*,*) 'res(:)= ', res(:)
            ! ! write(*,*) 'this%cache_delta3(:,idx)= ', this%cache_delta3(1,idx)
            ! write(*,*) 'idx= ',idx, 'this%cache_delta3(1,idx)= ', this%cache_delta3(1,idx)
            ! stop 
        else
            res = delta3(n_vec, V_3d, f_3d, VIN,VBASE)
            if (this%cache_pos == this%cache_length) then  ! clean cache
                call cache_reset(this)
            endif
            this%cache_pos = this%cache_pos + 1
            this%cache_idx(VIN(1),VIN(2),VIN(3),n_cur(1), n_cur(2), n_cur(3)) = this%cache_pos
            this%cache_delta3(1:n_vec,this%cache_pos) = res(1:n_vec)
            ! write(*,'(5x,1I3,2e12.4)') this%cache_pos, this%cache_delta3(1,this%cache_pos),this%cache_delta3(n_vec,this%cache_pos)
        endif        
        ! write(*,'(8I3,2e12.4)') idx, this%cache_pos, VIN(1),VIN(2),VIN(3),n_cur(1), n_cur(2), n_cur(3),res(1), res(n_vec)
        endassociate

    end function
        

    subroutine cache_reset(this)
        class(spline3dvec_type), intent(inout)  :: this
        this%cache_pos = 0
        this%cache_idx = 0
        this%cache_delta3 = 0. 
        this%cache_counter_reset = this%cache_counter_reset + 1
    endsubroutine cache_reset


    pure recursive function delta3(n_vec, V_3d, f_3d, VIN,VBASE) result(res)
        ! To calculate 3D Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        IMPLICIT NONE
        integer, intent(in) :: n_vec
        real(dp), dimension(-1:2,p_dim_3d), intent(in) :: V_3d
        real(dp), dimension(n_vec,-1:2,-1:2,-1:2), intent(in) :: f_3d
        INTEGER, dimension(p_dim_3d), intent(in) :: VIN, VBASE

        real(dp), dimension(n_vec) :: res
        INTEGER K
        INTEGER V1(p_dim_3d),V2(p_dim_3d)
        !--------------------------------------

        K = MAXLOC(VIN,DIM=1)
        IF(VIN(K).EQ.0)THEN
            RES = f_3d(:,VBASE(1),VBASE(2),VBASE(3))
            RETURN
        END IF

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        RES = VIN(K)*(delta3(n_vec, V_3d, f_3d, V1,V2)-delta3(n_vec, V_3d, f_3d, V1,VBASE)) / (V_3d(VBASE(K)+VIN(K),K)-V_3d(VBASE(K),K))
        return
    end function

end module rspline3dvec

