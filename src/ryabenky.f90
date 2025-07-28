module ryabmod
    use rspline3d, only: spline3d_type, p_dim_3d
    use rspline4d, only: spline4d_type, p_dim_4d

    implicit none
    private

    public :: spline3d_type, p_dim_3d
    public :: spline4d_type, p_dim_4d


    character(len=*), parameter, private  :: mdl_name = 'ryabmod'

end module ryabmod