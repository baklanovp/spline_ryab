module rspline2dvec
    use kinds,   only: dp, alloc1d
    use array_expand,  only: expand_2dvec, expand_1d

    implicit none
    private

    public :: spline2dvec_type

    character(len=*), parameter  :: mdl_name = 'rspline2dvec'
    integer, parameter, public :: p_dim_2d = 2 ! dimension of the tables

    type spline2dvec_type
        logical :: first_run
        integer :: ndim = p_dim_2d   ! dimension of the tables
        real(dp) :: V_2d(-1:2,p_dim_2d)  ! local volume
        real(dp), allocatable :: f_2d(:,:,:)  ! function values in the local volume
        integer, dimension(p_dim_2d) :: n_cur 
        real(dp), dimension(:), allocatable :: x_tab
        real(dp), dimension(:), allocatable :: y_tab
        real(dp), dimension(:,:,:), allocatable :: funcTab
        integer :: n_vec
        integer :: n_x, n_y
        logical :: is_cache = .false.
        integer, dimension(:,:,:,:), allocatable :: cache_idx
        integer :: cache_length
        integer :: cache_pos
        integer :: cache_counter_reset
        real(dp), dimension(:,:), allocatable :: cache_delta2
        contains
            procedure :: init => spline2d_init
            procedure :: value => spline2d_vec
            procedure :: check_value => spline2d_check_value
            procedure :: destroy => spline2d_destroy
    end type spline2dvec_type


    contains


    subroutine spline2d_init(this, x_tab, y_tab, funcTab, cache_length, is_expand) 
        class(spline2dvec_type), intent(inout)  :: this
        real(8), dimension(:), intent(in) :: x_tab, y_tab
        real(8), dimension(:,:,:), intent(in) :: funcTab
        integer, optional, intent(in) :: cache_length
        logical, intent(in), optional :: is_expand
        character(len=*), parameter ::  subrtn_name = 'spline2d_init', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        integer :: ierr                
        logical :: is_cache, is_expand_
        
        is_expand_ = .false.
        if ( present(is_expand) ) is_expand_ = is_expand
        is_cache = .false.
        if (present(cache_length)) is_cache = cache_length > 0

        if (is_expand_) then
            call expand_1d('this%x_tab', x_tab, this%x_tab)
            call expand_1d('this%y_tab', y_tab, this%y_tab)
        else
            call alloc1d('x_tab', size(x_tab), this%x_tab, path=fullPathSubrtn)
            this%x_tab = x_tab
            call alloc1d('y_tab', size(y_tab), this%y_tab, path=fullPathSubrtn)
            this%y_tab = y_tab
        endif

        this%n_vec = size(funcTab,1)
        call expand_2dvec(funcTab, this%funcTab)

        this%first_run = .true.

        this%n_x = size(x_tab)
        this%n_y = size(y_tab)

        this%n_cur = 1 !  todo check?

        if (size(funcTab,2) /= this%n_x) then
            write(*, '(4a, i4)') fullPathSubrtn, ' Size(x_tab)=', this%n_x,' is not equal  size(funcTab,2)= ', size(funcTab,2);
            error stop 666;
        endif

        if (size(funcTab,3) /= this%n_y) then
            write(*, '(4a, i4)') fullPathSubrtn, ' Size(y_tab)=', this%n_y,' is not equal  size(funcTab,3)= ', size(funcTab,3);
            error stop 666;
        endif

        if (is_expand_) then
            call expand_2dvec(funcTab, this%funcTab)
        else
            allocate(this%funcTab(this%n_vec,this%n_x,this%n_y), STAT=ierr);
            if (ierr /= 0) then
                write(*, '(2a, 3i4)') fullPathSubrtn, &
                ' Not enough memory for funcTab where n_vec, n_x,n_y =', this%n_vec,this%n_x,this%n_y;
                error stop 666;
            endif
            this%funcTab = funcTab
        endif

        ! Check sizes
        if (size(this%funcTab,2) /= this%n_x) then
            write(*, '(4a, i4)') fullPathSubrtn, ' Size(x_tab)=', this%n_x,' is not equal  size(this%funcTab,2)= ', size(this%funcTab,2);
            error stop 666;
        endif
        if (size(this%funcTab,3) /= this%n_y) then
            write(*, '(4a, i4)') fullPathSubrtn, ' Size(y_tab)=', this%n_y,' is not equal  size(this%funcTab,3)= ', size(this%funcTab,3);
            error stop 666;
        endif
        

        allocate(this%f_2d(this%n_vec,-1:2,-1:2), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a, 4i4)') fullPathSubrtn, &
               ' Not enough memory for f_2d where n_vec =', this%n_vec;
            error stop 666;
        endif
        this%f_2d = 0.

        if ( this%is_cache ) then;            
            call spline2d_init_cache(this, cache_length)
        endif   

    end subroutine spline2d_init

    subroutine spline2d_init_cache(this, clength)
        class(spline2dvec_type), intent(inout)  :: this
        integer, intent(in) :: clength
        character(len=*), parameter ::  subrtn_name = 'spline2d_init_cache', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name
        integer :: ierr

        this%cache_length = clength
        this%is_cache = .true.;
        allocate(this%cache_delta2(this%n_vec,this%cache_length), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a,i5,i10)') fullPathSubrtn, &
            ' Not enough memory for cache_delta2 where this%n_vec,cache_length =', this%n_vec,this%cache_length;
            error stop 666;
        endif

        allocate(this%cache_idx(0:3,0:3,this%n_x,this%n_y), STAT=ierr);
        if (ierr /= 0) then
            write(*, '(2a, 3i4)') fullPathSubrtn, &
            ' Not enough memory for cache_idx where 4*4* n_x,n_y =', this%n_x,this%n_y;
            error stop 666;
        endif

        this%cache_pos = 0
        this%cache_counter_reset = 0
        this%cache_idx = 0
        this%cache_delta2 = 0.
    endsubroutine spline2d_init_cache

    pure subroutine spline2d_destroy(this)
        implicit none
        class(spline2dvec_type), intent(inout)  :: this
        this%n_x    = 0
        this%n_y    = 0
        if (allocated(this%x_tab))      deallocate(this%x_tab)
        if (allocated(this%y_tab))      deallocate(this%y_tab)
        if (allocated(this%funcTab))    deallocate(this%funcTab)
        if (allocated(this%cache_delta2)) deallocate(this%cache_delta2)
        if (allocated(this%cache_idx))    deallocate(this%cache_idx)
    end subroutine spline2d_destroy


    subroutine spline2d_check_value(this, point, ierr)
        class(spline2dvec_type), intent(inout)  :: this
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
            ! stop 880
        endif

        if(ierr == 20)then
            write(*,'(A,100(1pe12.4))') 'y_tab : ',y_tab
            print*, fullPathSubrtn//': variables are out of range'
            print*,'y_tab(2),point(2),y_tab(n_y-2): ',y_tab(2),point(2),y_tab(n_y-2)
            ! read*
            ! stop 882
        endif
        endassociate
        endassociate
    end subroutine spline2d_check_value


    function spline2d_vec(this, point, ierr) result(res)
        class(spline2dvec_type), intent(inout)  :: this
        real(8), intent(in) :: point(p_dim_2d)
        integer, intent(out) :: ierr
        real(dp), dimension(this%n_vec) :: res

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
            f_2d = funcTab(:, n_cur(1)-1:n_cur(1)+2, n_cur(2)-1:n_cur(2)+2 )
        endif

        
        if ( this%is_cache) then
            call ryab_2d_cache(this,  V_2d, f_2d, point,res)
        else 
            call ryab_2d(this%n_vec, V_2d, f_2d, point,res)
        endif
        endassociate
        endassociate
        endassociate
        return
    endfunction spline2d_vec



    subroutine ryab_2d_cache(this, V_2d, f_2d, Y, res)
        class(spline2dvec_type), intent(inout)  :: this
        ! To calculate two-dimensional Ryabenkii spline with P=2,s=1
        real(dp), dimension(-1:2,p_dim_2d), intent(in) :: V_2d
        real(dp), dimension(this%n_vec,-1:2,-1:2,-1:2), intent(in) :: f_2d

        real(dp), dimension(p_dim_2d), intent(in):: Y
        real(dp), dimension(this%n_vec), intent(out) :: res

        real(dp) :: Q(0:3,p_dim_2d)
        real(dp) :: T(p_dim_2d)
        INTEGER I,J,K,VIN(p_dim_2d),VBASE(p_dim_2d)
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
                res = res+Q(I,1)*Q(J,2)*delta2cache(this,V_2d, f_2d, VIN,VBASE)
            enddo
        enddo   
    endsubroutine ryab_2d_cache


    pure subroutine ryab_2d(n_vec, V_2d, f_2d, Y, res)
        ! To calculate two-dimensional Ryabenkii spline with P=2,s=1
        integer, intent(in) :: n_vec
        real(dp), dimension(-1:2,p_dim_2d), intent(in) :: V_2d
        real(dp), dimension(-1:2,-1:2), intent(in) :: f_2d
        real(8), dimension(p_dim_2d), intent(in):: Y
        real(dp), dimension(n_vec), intent(out) :: res

        real(dp) :: Q(0:3,p_dim_2d)
        real(dp) :: T(p_dim_2d)
        integer, dimension(p_dim_2d) :: VIN, VBASE
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
                res = res+Q(I,1)*Q(J,2)*delta2(n_vec,V_2d, f_2d, VIN,VBASE)
            enddo
        enddo        
    end subroutine ryab_2d
    !=======================================================
    !*******************************************************

    function delta2cache(this,V_2d, f_2d, VIN,VBASE) result(res)
        ! Cache delta2
        class(spline2dvec_type), intent(inout)  :: this
        real(dp), dimension(-1:2,p_dim_2d), intent(in) :: V_2d
        real(dp), dimension(this%n_vec,-1:2,-1:2), intent(in) :: f_2d
        integer, dimension(p_dim_2d), intent(in) :: VIN, VBASE
        real(dp), dimension(this%n_vec) :: res
        integer :: idx

        associate(n_cur=>this%n_cur, n_vec=>this%n_vec)
        idx = this%cache_idx(VIN(1),VIN(2),n_cur(1),n_cur(2))
        if ( idx > 0 ) then
            res(1:n_vec) = this%cache_delta2(1:n_vec,idx)
        else
            res = delta2(n_vec, V_2d, f_2d, VIN,VBASE)
            if (this%cache_pos == this%cache_length) then  ! clean cache
                call cache_reset(this)
            endif
            this%cache_pos = this%cache_pos + 1
            this%cache_idx(VIN(1),VIN(2),n_cur(1), n_cur(2)) = this%cache_pos
            this%cache_delta2(1:n_vec,this%cache_pos) = res(1:n_vec)
            ! write(*,'(5x,1I3,2e12.4)') this%cache_pos, this%cache_delta2(1,this%cache_pos),this%cache_delta2(n_vec,this%cache_pos)
        endif        
        ! write(*,'(8I3,2e12.4)') idx, this%cache_pos, VIN(1),VIN(2),VIN(3),n_cur(1), n_cur(2), n_cur(3),res(1), res(n_vec)
        endassociate

    end function


    subroutine cache_reset(this)
        class(spline2dvec_type), intent(inout)  :: this
        this%cache_pos = 0
        this%cache_idx = 0
        this%cache_delta2 = 0. 
        this%cache_counter_reset = this%cache_counter_reset + 1
    endsubroutine cache_reset


    pure recursive function delta2(n_vec, V_2d, f_2d, VIN,VBASE) result(res)
        ! To calculate 2d Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        integer, intent(in) :: n_vec
        real(dp), dimension(-1:2,p_dim_2d), intent(in) :: V_2d
        real(dp), dimension(n_vec,-1:2,-1:2), intent(in) :: f_2d
        integer, dimension(p_dim_2d), intent(in) :: VIN, VBASE

        real(dp),dimension(n_vec) :: res
        integer :: K
        integer, dimension(p_dim_2d) ::  V1, V2

        K = maxloc(VIN,DIM=1)
        if(VIN(k) == 0) then
            res = f_2d(:,VBASE(1),VBASE(2))
            return
        endif

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        res = VIN(K)*(delta2(n_vec, V_2d, f_2d, V1,V2)-delta2(n_vec, V_2d, f_2d, V1,VBASE)) / (V_2d(VBASE(K)+VIN(K),K)-V_2d(VBASE(K),K))
        return
    end function

end module rspline2dvec

