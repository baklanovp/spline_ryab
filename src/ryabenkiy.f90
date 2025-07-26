MODULE RYABMOD
   INTEGER, PARAMETER :: NDIM=3 ! dimension of the tables
   REAL(8) X(-1:2,NDIM)
   ! REAL(8) F(-1:2,-1:2) ! two-dimensional
   REAL(8) F(-1:2,-1:2,-1:2)
END MODULE


!====================================================
!====================================================

SUBROUTINE RYAB3(Y,RES)
    ! To calculate two-dimensional Ryabenkii spline
    ! with P=2,s=1
    USE RYABMOD, only: NDIM, X
    IMPLICIT NONE
    REAL(8), INTENT(IN):: Y(NDIM)
    REAL(8), INTENT(OUT) :: RES
    REAL(8) Q(0:3,NDIM)
    REAL(8) T(NDIM)
    REAL(8) DELTA3
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

END SUBROUTINE
!=======================================================
!*******************************************************

RECURSIVE FUNCTION DELTA3(VIN,VBASE) RESULT(RES)
    ! To calculate 3D Delta_X^I*Delta_Y^J*Delta_Z^K F_{M,N,P}
    ! VIN=(/I,J,K/), VBASE=(/M,N,P/)
    USE RYABMOD, only: NDIM, X, F
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
