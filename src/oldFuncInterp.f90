real(8) function FUNCINTERP(point) result(z)
    use datas
    use ryabmod
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

end function
