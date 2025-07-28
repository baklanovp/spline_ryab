module ryabmod
    use rspline3d, only: spline3d_type, p_dim_3d

    implicit none
    private

    public :: spline3d_type, p_dim_3d

    character(len=*), parameter  :: mdl_name = 'ryabmod'

end module ryabmod