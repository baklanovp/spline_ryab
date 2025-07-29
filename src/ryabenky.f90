module ryabmod
    use rspline2d, only: spline2d_type, p_dim_2d
    use rspline3d, only: spline3d_type, p_dim_3d
    use rspline4d, only: spline4d_type, p_dim_4d

    implicit none
    private

    public :: spline2d_type, p_dim_2d
    public :: spline3d_type, p_dim_3d
    public :: spline4d_type, p_dim_4d


    character(len=*), parameter, private  :: mdl_name = 'ryabmod'

end module ryabmod