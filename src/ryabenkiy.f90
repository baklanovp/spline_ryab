MODULE RYABMOD
    use datas

    ! implicit none

    public :: FUNCINTERP

   INTEGER, PARAMETER :: NDIM=3 ! dimension of the tables
   REAL(8) X(-1:2,NDIM)
   ! REAL(8) F(-1:2,-1:2) ! two-dimensional
   REAL(8) F(-1:2,-1:2,-1:2)

    !====================================================
    !====================================================

   contains

    real(8) function FUNCINTERP(point) result(z)
        ! use ryabmod
        implicit none
        real(8), intent(in) :: point(NDIM)
        integer, save :: n_cur(NDIM)
        logical :: reload,first_run=.true.

        !------/checking if current position is not out of the table's ranges/------

        if(point(1).lt.TpTab(2).or.point(1).gt.TpTab(n_tp-1))then
            print*,'variables are out of range'
            print*,'T',TpTab(2),point(1),TpTab(n_tp-1)
            read*
            stop
        end if

        if(point(2).lt.RhoTab(2).or.point(2).gt.RhoTab(n_rho-1))then
            print*,'variables are out of range'
            print*,'Rho',RhoTab(2),point(2),RhoTab(n_rho-1)
            read*
            stop
        end if

        if(point(3).lt.lnTimeTab(2).or.point(3).gt.lnTimeTab(n_times-1))then
            print*,'variables are out of range'
            print*,'Rho',lnTimeTab(2),point(3),lnTimeTab(n_times-1)
            read*
            stop
        end if

        !---------/check if we're at the old xyz box (from the previous call)/----------------

        reload=.false.

        if(first_run)then
            first_run=.false.
            reload=.true.
        end if

        if(reload.or.point(1).lt.TpTab(n_cur(1)).or.point(1).gt.TpTab(n_cur(1)+1))then
            reload=.true.
            n_cur(1)=minloc(point(1)-TpTab,mask=point(1)-TpTab.ge.0.d0,dim=1)
        end if

        if(reload.or.point(2).lt.RhoTab(n_cur(2)).or.point(2).gt.RhoTab(n_cur(2)+1))then
            reload=.true.
            n_cur(2)=minloc(point(2)-RhoTab,mask=point(2)-RhoTab.ge.0.d0,dim=1)
        end if

        if(reload.or.point(3).lt.lnTimeTab(n_cur(3)).or.point(3).gt.lnTimeTab(n_cur(3)+1))then
            reload=.true.
            n_cur(3)=minloc(point(3)-lnTimeTab,mask=point(3)-lnTimeTab.ge.0.d0,dim=1)
        end if

        if(reload)then
            X(:,1)=TpTab(n_cur(1)-1:n_cur(1)+2)
            X(:,2)=RhoTab(n_cur(2)-1:n_cur(2)+2)
            X(:,3)=lnTimeTab(n_cur(3)-1:n_cur(3)+2)
            F=arr_dump(n_cur(1)-1:n_cur(1)+2,n_cur(2)-1:n_cur(2)+2,n_cur(3)-1:n_cur(3)+2)
        end if

        call RYAB3(point,z)

        return

    endfunction FUNCINTERP


    SUBROUTINE RYAB3(Y,RES)
        ! To calculate two-dimensional Ryabenkii spline
        ! with P=2,s=1
        ! USE RYABMOD, only: NDIM, X
        IMPLICIT NONE
        REAL(8), INTENT(IN):: Y(NDIM)
        REAL(8), INTENT(OUT) :: RES
        REAL(8) Q(0:3,NDIM)
        REAL(8) T(NDIM)
        INTEGER I,J,K,VIN(NDIM),VBASE(NDIM)
        !---------------------------------

        Q(0,:)=1.d0

        Q(1,:)=Y-X(-1,:)

        Q(2,:)=0.5d0*(Y-X(-1,:))*(Y-X(0,:))

        T=(Y-X(0,:))/(X(1,:)-X(0,:))
        Q(3,:)=0.5d0*(X(1,:)-X(0,:))**2*(X(2,:)-X(-1,:))*T**3*(T-1)*(1-2.d0/3.d0*T)

        VBASE=-1
        RES=0.D0

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

    END SUBROUTINE RYAB3
    !=======================================================
    !*******************************************************

    RECURSIVE FUNCTION DELTA3(VIN,VBASE) RESULT(RES)
        ! To calculate 3D Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
        ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
        ! USE RYABMOD, only: NDIM, X, F
        IMPLICIT NONE
        REAL(8) RES
        INTEGER K,I,J,M,N
        INTEGER VIN(NDIM),VBASE(NDIM)
        INTEGER V1(NDIM),V2(NDIM)
        !--------------------------------------

        K=MAXLOC(VIN,DIM=1)

        IF(VIN(K).EQ.0)THEN
            RES = F(VBASE(1),VBASE(2),VBASE(3))
            RETURN
        END IF

        V1 = VIN
        V1(K) = V1(K)-1
        V2 = VBASE
        V2(K) = V2(K)+1

        RES = VIN(K)*(DELTA3(V1,V2)-DELTA3(V1,VBASE)) / (X(VBASE(K)+VIN(K),K)-X(VBASE(K),K))

        RETURN
    END FUNCTION

END MODULE

