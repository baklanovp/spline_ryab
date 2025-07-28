MODULE RYABMOD

    implicit none

    public :: rspline3d

    INTEGER, PARAMETER :: p_Ndim_3d = 3 ! dimension of the tables
    REAL(8) V_3d(-1:2,p_Ndim_3d)
    REAL(8) f_3d(-1:2,-1:2,-1:2)

    !====================================================
    !====================================================

   contains

    real(8) function rspline3d(n_x, n_y, n_z, x_tab, y_tab, z_tab, funcTab, point) result(res)
        ! use ryabmod
        implicit none
        real(8), dimension(:), intent(in) :: x_tab, y_tab, z_tab
        real(8), dimension(:,:,:) :: funcTab
        integer, intent(in)  :: n_x, n_y, n_z

        real(8), intent(in) :: point(p_Ndim_3d)
        integer, save :: n_cur(p_Ndim_3d)
        logical :: reload,first_run=.true.

        !------/checking if current position is not out of the table's ranges/------

        if(point(1).lt.x_tab(2).or.point(1).gt.x_tab(n_x-1))then
            print*,'variables are out of range'
            print*,'X',x_tab(2),point(1),x_tab(n_x-1)
            read*
            stop
        end if

        if(point(2).lt.y_tab(2).or.point(2).gt.y_tab(n_y-1))then
            print*,'variables are out of range'
            print*,'Y',y_tab(2),point(2),y_tab(n_y-1)
            read*
            stop
        end if

        if(point(3).lt.z_tab(2).or.point(3).gt.z_tab(n_z-1))then
            print*,'variables are out of range'
            print*,'res',z_tab(2),point(3),z_tab(n_z-1)
            read*
            stop
        end if

        !---------/check if we're at the old xyz box (from the previous call)/----------------

        reload=.false.

        if(first_run)then
            first_run=.false.
            reload=.true.
        end if

        if(reload.or.point(1).lt.x_tab(n_cur(1)).or.point(1).gt.x_tab(n_cur(1)+1))then
            reload=.true.
            n_cur(1)=minloc(point(1)-x_tab,mask=point(1)-x_tab.ge.0.d0,dim=1)
        end if

        if(reload.or.point(2).lt.y_tab(n_cur(2)).or.point(2).gt.y_tab(n_cur(2)+1))then
            reload=.true.
            n_cur(2)=minloc(point(2)-y_tab,mask=point(2)-y_tab.ge.0.d0,dim=1)
        end if

        if(reload.or.point(3).lt.z_tab(n_cur(3)).or.point(3).gt.z_tab(n_cur(3)+1))then
            reload=.true.
            n_cur(3)=minloc(point(3)-z_tab,mask=point(3)-z_tab.ge.0.d0,dim=1)
        end if

        if(reload)then
            V_3d(:,1) = x_tab(n_cur(1)-1:n_cur(1)+2)
            V_3d(:,2) = y_tab(n_cur(2)-1:n_cur(2)+2)
            V_3d(:,3) = z_tab(n_cur(3)-1:n_cur(3)+2)
            f_3d = funcTab(n_cur(1)-1:n_cur(1)+2,n_cur(2)-1:n_cur(2)+2,n_cur(3)-1:n_cur(3)+2)
        end if

        call ryab_3d(point,res)

        return

    endfunction rspline3d


    pure subroutine ryab_3d(Y,RES)
        ! To calculate two-dimensional Ryabenkii spline
        ! with P=2,s=1
        ! USE RYABMOD, only: p_Ndim_3d, X
        implicit none
        real(8), dimension(:), intent(in):: Y
        real(8), intent(out) :: res
        REAL(8) Q(0:3,size(Y))
        REAL(8) T(size(Y))
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
                RES=RES+Q(I,1)*Q(J,2)*Q(K,3)*DELTA3(VIN,VBASE)
                END DO
            END DO
        END DO

    END SUBROUTINE ryab_3d
    !=======================================================
    !*******************************************************

    pure RECURSIVE FUNCTION DELTA3(VIN,VBASE) RESULT(RES)
        ! To calculate 3D Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        IMPLICIT NONE
        INTEGER, dimension(p_Ndim_3d), intent(in) :: VIN, VBASE

        REAL(8) RES
        INTEGER K
        INTEGER V1(p_Ndim_3d),V2(p_Ndim_3d)
        !--------------------------------------

        K=MAXLOC(VIN,DIM=1)

        IF(VIN(K).EQ.0)THEN
            RES = f_3d(VBASE(1),VBASE(2),VBASE(3))
            RETURN
        END IF

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        RES = VIN(K)*(DELTA3(V1,V2)-DELTA3(V1,VBASE)) / (V_3d(VBASE(K)+VIN(K),K)-V_3d(VBASE(K),K))

        RETURN
    END FUNCTION

END MODULE

