module ryabmod
    use rspline2d, only: spline2d_type, p_dim_2d
    use rspline3d, only: spline3d_type, p_dim_3d
    use rspline2dvec, only: spline2dvec_type
    use rspline3dvec, only: spline3dvec_type
    use rspline4d, only: spline4d_type, p_dim_4d    
    use array_expand, only: expand_1d, expand_2dvec, expand_3dvec

    implicit none
    private

    public :: spline2d_type, spline2dvec_type, p_dim_2d
    public :: spline3d_type, spline3dvec_type, p_dim_3d
    public :: spline4d_type, p_dim_4d
    public :: expand_1d, expand_2dvec, expand_3dvec


    character(len=*), parameter, private  :: mdl_name = 'ryabmod'

end module ryabmod