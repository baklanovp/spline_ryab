program MAIN
    use kinds,            only: dp
    implicit none
    

    call test_spline3d;


contains
    
  
  subroutine test_spline3d
    use ryabmod, only: spline3d_type, p_dim_3d

    type(spline3d_type) :: rspline
    real(dp) s, per(p_dim_3d)
    real :: t1, t2

    ! integer K
    integer :: c
    integer :: ierr
    real(dp), dimension(:), allocatable :: TpTab, RhoTab, lnTimeTab
    real(dp), dimension(:,:,:), allocatable :: arr_dump
    integer :: n_tp, n_rho, n_times

    call load_3d(n_tp, n_rho, n_times, TpTab, RhoTab, lnTimeTab, arr_dump)

    open(11,file='result.dat')

    ! do K=1,n_tp
    !   write(11,'(10es23.15)') TpTab(K),arr_dump(K,5,3)
    ! end do

    ! close(11)
    ! stop


    call rspline%init(TpTab, RhoTab, lnTimeTab, arr_dump)

    c = 0
    s=0.d0
    call cpu_time(t1)
    
    do while(s.le.1.d0)
      c = c+1
      per(1)=TpTab(2)+s*(TpTab(n_tp-1)-TpTab(2))
      per(2)=RhoTab(5) !RhoTab(2)+s*(RhoTab(n_rho-1)-RhoTab(2))
      per(3)=lnTimeTab(3) !lnTimeTab(2)+s*(lnTimeTab(n_times-1)-lnTimeTab(2))

      write(11,'(10es23.15)') per(1), rspline%value(per, ierr) 
      if (ierr > 0) then
        call rspline%check_value(per, ierr)
      endif
      ! print*,s
      !  read*
      s=s+1.d-4
    end do

    close(11)

    ! Code segment to be timed
    call cpu_time(t2)

    write(*,*) 'Time taken for ', c, ' calls of rspline3d: ', t2 - t1, ' seconds.'

    stop
  endsubroutine test_spline3d
    

  subroutine  load_3d(n_tp, n_rho, n_times, TpTab, RhoTab, lnTimeTab, arr_dump)
    real(dp), dimension(:), allocatable, intent(out) :: TpTab, RhoTab, lnTimeTab
    real(dp), dimension(:,:,:), allocatable, intent(out) :: arr_dump
    integer, intent(out) :: n_tp, n_rho, n_times

    character(30) fname
    integer :: ierr, ui

    fname = 'data/neM20Ni01Z002.3d.dump'

    write(*,"(A,A)") 'Loading from file: ', trim(fname);

    open(newunit=ui, file=trim(fname), status='unknown', form='formatted', IOSTAT=ierr);
    read(ui,*) n_tp, n_rho, n_times
    write(*,*) ' n_tp, n_rho, n_times: ', n_tp, n_rho, n_times

    allocate(TpTab(n_tp),RhoTab(n_rho),lnTimeTab(n_times))
    allocate(arr_dump(n_tp,n_rho,n_times))

    read(ui,*) TpTab
    read(ui,*) RhoTab
    read(ui,*) lnTimeTab
    read(ui,*) arr_dump

    close(ui)

    write(*,"(A,A/)") ' Show data ';


    write(*,*) 'TpTab: ', TpTab
    !		read*
    write(*,*) 'RhoTab: ', RhoTab
    !		read*
    write(*,*) 'lnTimeTab: ', lnTimeTab
    !		read*
    write(*,*) 'arr_dump: ', arr_dump
    ! read*

  end subroutine  load_3d

end program
