module array_expand
    use kinds,   only: dp, alloc1d

    implicit none
    private

    public :: expand_1d, expand_2dvec, expand_3dvec

    character(len=*), parameter  :: mdl_name = 'array_expand'


contains
    
    subroutine expand_2dvec(funcTab, funcTabExp)
        ! Expand funcTab by 1 point on each side if is_expand = .true.
        real(dp), dimension(:,:,:), intent(in) :: funcTab
        character(len=*), parameter ::  subrtn_name = 'expand_2dvec', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        real(dp), dimension(:,:,:), allocatable, intent(out) :: funcTabExp
        integer :: nvec, lx, ly
        integer :: ierr

        nvec = size(funcTab,1)            
        lx = size(funcTab,2)
        ly = size(funcTab,3)

        if (allocated(funcTabExp)) deallocate(funcTabExp);

        allocate(funcTabExp(nvec, lx+2, ly+2), STAT=ierr); ! Add 2 points (one on each side)
        if (ierr /= 0) then
            write(*, '(2a,4i10)') fullPathSubrtn, &
            ' Not enough memory for funcTabExp where nvec, lx, ly =', nvec, lx+2, ly+2;
            error stop 666;
        endif
        funcTabExp = -99 ! Initialize to some wrong value to catch errors

        ! Fill the expanded array
        funcTabExp(:,2:lx+1,2:ly+1) = funcTab(:,1:lx,1:ly)  
        !! Points
        funcTabExp(:,1,1)       = funcTab(:,1,1)
        funcTabExp(:,lx+2,1)    = funcTab(:,lx,1)
        funcTabExp(:,1,ly+2)    = funcTab(:,1,ly)
        funcTabExp(:,lx+2,ly+2) = funcTab(:,lx,ly)
        ! X
        funcTabExp(:,2:lx+1,1)    = funcTab(:,1:lx,1) 
        funcTabExp(:,2:lx+1,ly+2) = funcTab(:,1:lx,ly)
        ! Y
        funcTabExp(:,1,2:ly+1)       = funcTab(:,1,1:ly)
        funcTabExp(:,lx+2,2:ly+1)    = funcTab(:,lx,1:ly)        
        ! write(*, '(2a, 4i4)') fullPathSubrtn, ' Expansion done for funcTabExp where nvec, lx, ly =', nvec, lx+2, ly+2;
    end subroutine expand_2dvec


    subroutine expand_3dvec(funcTab, funcTabExp)
        ! Expand funcTab by 1 point on each side if is_expand = .true.
        real(dp), dimension(:,:,:,:), intent(in) :: funcTab
        character(len=*), parameter ::  subrtn_name = 'expand_3dvec', &
                    fullPathSubrtn = mdl_name//'.'//subrtn_name

        real(dp), dimension(:,:,:,:), allocatable, intent(out) :: funcTabExp
        integer :: nvec, lx, ly, lz
        integer :: ierr
        ! integer :: i, j, k

        nvec = size(funcTab,1)            
        lx = size(funcTab,2)
        ly = size(funcTab,3)
        lz = size(funcTab,4)

        if (allocated(funcTabExp)) deallocate(funcTabExp);

        allocate(funcTabExp(nvec, lx+2, ly+2, lz+2), STAT=ierr); ! Add 2 points (one on each side)
        if (ierr /= 0) then
            write(*, '(2a,4i10)') fullPathSubrtn, &
            ' Not enough memory for funcTabExp where nvec, lx, ly, lz =', nvec, lx+2, ly+2, lz+2;
            error stop 666;
        endif
        funcTabExp = -99 ! Initialize to some wrong value to catch errors

        ! Fill the expanded array
        funcTabExp(:,2:lx+1,2:ly+1,2:lz+1) = funcTab(:,:,:,:)

        !! Points
        ! XY plane: 1
        funcTabExp(:,1,1,1) = funcTab(:,1,1,1)
        funcTabExp(:,lx+2,1,1) = funcTab(:,lx,1,1)
        funcTabExp(:,1,ly+2,1) = funcTab(:,1,ly,1)
        funcTabExp(:,lx+2,ly+2,1) = funcTab(:,lx,ly,1)
        ! XY plane: lz+2
        funcTabExp(:,1,1,lz+2) = funcTab(:,1,1,lz)
        funcTabExp(:,lx+2,1,lz+2) = funcTab(:,lx,1,lz)
        funcTabExp(:,1,ly+2,lz+2) = funcTab(:,1,ly,lz)
        funcTabExp(:,lx+2,ly+2,lz+2) = funcTab(:,lx,ly,lz)

        !! Faces
        funcTabExp(:,1,2:ly+1,2:lz+1) = funcTab(:,1,1:ly,1:lz)  ! XY
        funcTabExp(:,lx+2,2:ly+1,2:lz+1) = funcTab(:,lx,1:ly,1:lz)  ! XY
        funcTabExp(:,2:lx+1,1,2:lz+1) = funcTab(:,1:lx,1,1:lz)  ! Y
        funcTabExp(:,2:lx+1,ly+2,2:lz+1) = funcTab(:,1:lx,ly,1:lz)  ! Y
        funcTabExp(:,2:lx+1,2:ly+1,1) = funcTab(:,1:lx,1:ly,1)  ! Z
        funcTabExp(:,2:lx+1,2:ly+1,lz+2) = funcTab(:,1:lx,1:ly,lz)  ! Z

        !! Edges
        ! Z
        funcTabExp(:,1,1,2:lz+1)       = funcTab(:,1,1,1:lz)  
        funcTabExp(:,lx+2,1,2:lz+1)    = funcTab(:,lx,1,1:lz) 
        funcTabExp(:,1,ly+2,2:lz+1)    = funcTab(:,1,ly,1:lz) 
        funcTabExp(:,lx+2,ly+2,2:lz+1) = funcTab(:,lx,ly,1:lz)
        ! X
        funcTabExp(:,2:lx+1,1,1)       = funcTab(:,1:lx,1,1) 
        funcTabExp(:,2:lx+1,1,lz+2)    = funcTab(:,1:lx,1,lz)
        funcTabExp(:,2:lx+1,ly+2,1)    = funcTab(:,1:lx,ly,1)
        funcTabExp(:,2:lx+1,ly+2,lz+2) = funcTab(:,1:lx,ly,lz)
        ! Y
        funcTabExp(:,1,2:ly+1,1)       = funcTab(:,1,1:ly,1)
        funcTabExp(:,1,2:ly+1,lz+2)    = funcTab(:,1,1:ly,lz)
        funcTabExp(:,lx+2,2:ly+1,1)    = funcTab(:,lx,1:ly,1)
        funcTabExp(:,lx+2,2:ly+1,lz+2) = funcTab(:,lx,1:ly,lz)    
        
        ! write(*, '(2a, 4i4)') fullPathSubrtn, ' Expansion done for funcTabExp where nvec, lx, ly, lz =', nvec, lx+2, ly+2, lz+2;
    end subroutine expand_3dvec


    subroutine expand_1d(nm, x, x_exp)
        character(len=*), intent(in) :: nm
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), allocatable, intent(out) :: x_exp
        integer :: l

        l = size(x)
        call alloc1d(nm, l+2, x_exp)
        x_exp(1) = x(1)
        x_exp(2:l+1) = x(1:l)
        x_exp(l+2) = x(l)
        ! write(*, '(3a,2i10)') 'expand_1d: ', nm, ' lx, lx_exp =', l, size(x_exp);
    end subroutine expand_1d


end module array_expand