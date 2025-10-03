program main_value
    use kinds,            only: dp
    use cla, only: cla_init, cla_register, cla_get, cla_help, cla_int, cla_flag, cla_key_present;

    implicit none
    
    character(len=*), parameter  :: mdl_name = 'main_value'

    character(len=*), parameter  :: dir_data = '../data'
    logical :: is_3d, is_4d;
    logical :: is_cache;

    call args_init(is_3d, is_4d, is_cache)

    if (is_3d) then
      call test_spline3d(is_cache);
      call test_spline3d_vec(is_cache);
    endif
    
    if (is_4d) then
      call test_spline4d;
    endif

contains
    
    
  subroutine test_spline4d
    use ryabmod, only: spline4d_type, p_dim_4d

    character(len=*), parameter ::  subrtn_name = 'test_spline4d', &
                      fullPathSubrtn = mdl_name//'.'//subrtn_name
    type(spline4d_type) :: rspline
    real(dp) s, per(p_dim_4d)
    real :: t1, t2

    integer :: i, j, k, c, ui
    integer :: ierr
    real(dp), dimension(:), allocatable :: TpTab, RhoTab, lnTimeTab, LcoordTab
    real(dp), dimension(:,:,:,:), allocatable :: arr_dump
    integer :: n_tp, n_rho, n_coord, n_times

    write(*,*) '---  RUN test: ', fullPathSubrtn

    call load_4d(n_tp, n_rho, n_coord, n_times, TpTab, RhoTab, LcoordTab, lnTimeTab, arr_dump)

    
    open(newunit=ui, file='input_4d.dat', status='unknown', form='formatted', IOSTAT=ierr);

    do i=1,n_tp
      do j=1,n_rho
        write(ui,'(10es23.15)') TpTab(i), RhoTab(j), arr_dump(i,j,6,3)
      end do
    end do

    close(ui)
    ! stop

    open(newunit=ui, file='result_4d.dat', status='unknown', form='formatted', IOSTAT=ierr);
    call rspline%init(TpTab, RhoTab, LcoordTab, lnTimeTab, arr_dump)

    c = 0
    s=0.d0
    call cpu_time(t1)
    
    do while(s.le.1.d0)
      c = c+1
      ! per(1) = TpTab(5) 
      per(1) = TpTab(2) + s*(TpTab(n_tp-1)-TpTab(2))
      per(2) = RhoTab(2) + s*(RhoTab(n_rho-1)-RhoTab(2))
      ! per(2) = RhoTab(n_rho-1) - s*(RhoTab(n_rho-1)-RhoTab(2))
      
      per(3) = LcoordTab(6)
      per(4) = lnTimeTab(3)
      write(ui,'(10es23.15)') (per(i),i=1,2), rspline%value(per, ierr)

      ! per(1) = TpTab(2)+s*(TpTab(n_tp-1)-TpTab(2))
      ! per(2) = RhoTab(2)+s*(RhoTab(n_rho-1)-RhoTab(2))
      ! per(3) = LcoordTab(2)+s*(LcoordTab(n_coord-1)-LcoordTab(2))
      ! per(4) = lnTimeTab(2)+s*(lnTimeTab(n_times-1)-lnTimeTab(2))
      ! write(ui,'(10es23.15)') (per(i),i=1,4), rspline%value(per, ierr) 

      if (ierr > 0) then
        call rspline%check_value(per, ierr)
      endif
      ! print*,s
      !  read*
      s=s+1.d-3
    end do

    close(ui)
    call rspline%destroy()
    ! Code segment to be timed
    call cpu_time(t2)

    write(*,*) 'Time taken for ', c, ' calls of rspline4d: ', t2 - t1, ' seconds.'
  endsubroutine test_spline4d


  subroutine test_spline3d(is_cache)
    use ryabmod, only: spline3d_type, p_dim_3d
    logical, intent(in) :: is_cache

    character(len=*), parameter ::  subrtn_name = 'test_spline3d', &
                      fullPathSubrtn = mdl_name//'.'//subrtn_name
    type(spline3d_type) :: rspline
    real(dp) s, per(p_dim_3d)
    real :: t1, t2

    ! integer K
    integer :: c
    integer :: ierr
    real(dp), dimension(:), allocatable :: TpTab, RhoTab, lnTimeTab
    real(dp), dimension(:,:,:), allocatable :: arr_dump
    integer :: n_tp, n_rho, n_times

    write(*,*) '---  RUN test: ', fullPathSubrtn, ' is_cache= ',is_cache

    call load_3d(n_tp, n_rho, n_times, TpTab, RhoTab, lnTimeTab, arr_dump)

    open(11,file='result_3d.dat')

    ! do K=1,n_tp
    !   write(11,'(10es23.15)') TpTab(K),arr_dump(K,5,3)
    ! end do

    ! close(11)
    ! stop


    call rspline%init(TpTab, RhoTab, lnTimeTab, arr_dump, is_cache)

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
    call rspline%destroy()
    ! Code segment to be timed
    call cpu_time(t2)

    write(*,*) 'Time taken for ', c, ' calls of rspline3d: ', t2 - t1, ' seconds.'

  endsubroutine test_spline3d
    

  subroutine test_spline3d_vec(is_cache)
    use ryabmod, only: spline3dvec_type, p_dim_3d
    logical, intent(in) :: is_cache
    character(len=*), parameter ::  subrtn_name = 'test_spline3d_vec', &
                      fullPathSubrtn = mdl_name//'.'//subrtn_name
    type(spline3dvec_type) :: rspline
    real(dp) s, per(p_dim_3d)
    real :: t1, t2
    integer, parameter :: p_nfreq = 11
    real(dp), parameter :: p_wlmax = 50000._dp;! in A
    real(dp), parameter :: p_wlmin = 1._dp;! in A
    real(dp), parameter :: p_c = 2.9979245800e+10_dp

    ! integer K
    integer :: c
    integer :: ierr
    real(dp), dimension(:), allocatable :: TpTab, RhoTab, lnTimeTab
    real(dp), dimension(:,:,:), allocatable :: arr_dump
    real(dp), dimension(:,:,:,:), allocatable :: funcTab

    integer :: i, n_tp, n_rho, n_times
    real(dp), dimension(p_nfreq) :: fr, res;
    real(dp) :: basis

    write(*,*) '---  RUN test: ', fullPathSubrtn, ' is_cache= ',is_cache

    call load_3d(n_tp, n_rho, n_times, TpTab, RhoTab, lnTimeTab, arr_dump)

    allocate(funcTab(p_nfreq,n_tp,n_rho,n_times))

    ! заглушка
    do i = 1, p_nfreq
      funcTab(i,:,:,:) = arr_dump(:,:,:) 
    enddo

    open(11,file='result_3d_vec.dat')
    open(21,file='result_3d_vecall.dat')


    fr(1) = p_c * 1.d+08 / p_wlmax;
    fr(p_nfreq) = p_c * 1.d+08 / p_wlmin;
    basis = (p_wlmax/p_wlmin)**(1.D0/(DBLE(p_nfreq)));
    !  Geometric PROGRESSION
    do i = 2, p_nfreq
        fr(i) = fr(i-1)*basis;
    enddo


    call rspline%init(TpTab, RhoTab, lnTimeTab, funcTab, is_cache)

    c = 0
    s=0.d0
    call cpu_time(t1)
    
    do while(s.le.1.d0)
      c = c+1
      per(1)=TpTab(2)+s*(TpTab(n_tp-1)-TpTab(2))
      per(2)=RhoTab(5) !RhoTab(2)+s*(RhoTab(n_rho-1)-RhoTab(2))
      per(3)=lnTimeTab(3) !lnTimeTab(2)+s*(lnTimeTab(n_times-1)-lnTimeTab(2))

      res = rspline%value(per, ierr)
      write(11,'(10es23.15)') per(1), res(1)
      ! write(11,'(10es23.15)') per(1), per(2), per(3), res(1)
      write(21,'(100es23.15)') per(1), res
      if (ierr > 0) then
        call rspline%check_value(per, ierr)
      endif
      ! print*,s
      !  read*
      s=s+1.d-4
    end do

    close(11)
    close(21)

    call rspline%destroy()
    ! Code segment to be timed
    call cpu_time(t2)

    write(*,*) 'Time taken for ', c, ' calls of rspline3d: ', t2 - t1, ' seconds.'

    stop
  endsubroutine test_spline3d_vec
    

  subroutine  load_4d(n_tp, n_rho, n_coord, n_times, TpTab, RhoTab, LcoordTab, lnTimeTab, arr_dump)
    real(dp), dimension(:), allocatable, intent(out) :: TpTab, RhoTab, LcoordTab, lnTimeTab
    real(dp), dimension(:,:,:,:), allocatable, intent(out) :: arr_dump
    integer, intent(out) :: n_tp, n_rho, n_coord, n_times

    character(30) fname
    integer :: ierr, ui

    fname = dir_data//'/neM20Ni01Z002.4d.dump'

    write(*,"(A,A)") 'Loading from file: ', trim(fname);

    open(newunit=ui, file=trim(fname), status='unknown', form='formatted', IOSTAT=ierr);
    read(ui,*) n_tp, n_rho, n_coord, n_times
    write(*,*) ' n_tp, n_rho, n_coord, n_times: ', n_tp, n_rho, n_coord, n_times

    allocate(TpTab(n_tp))
    allocate(RhoTab(n_rho))
    allocate(LcoordTab(n_coord))
    allocate(lnTimeTab(n_times))
    allocate(arr_dump(n_tp,n_rho,n_coord,n_times))

    read(ui,*) TpTab
    read(ui,*) RhoTab
    read(ui,*) LcoordTab
    read(ui,*) lnTimeTab

    read(ui,*) arr_dump
    
    close(ui)

    write(*,"(A,A/)") ' Show data ';


    write(*,*) 'TpTab: ', TpTab
    !		read*
    write(*,*) 'RhoTab: ', RhoTab
    ! read*
    write(*,*) 'LcoordTab: ', LcoordTab
    !		read*
    write(*,*) 'lnTimeTab: ', lnTimeTab
    ! read*
    write(*,*) 'arr_dump: ', arr_dump
    ! read*

  end subroutine  load_4d


  subroutine  load_3d(n_tp, n_rho, n_times, TpTab, RhoTab, lnTimeTab, arr_dump)
    real(dp), dimension(:), allocatable, intent(out) :: TpTab, RhoTab, lnTimeTab
    real(dp), dimension(:,:,:), allocatable, intent(out) :: arr_dump
    integer, intent(out) :: n_tp, n_rho, n_times

    character(30) fname
    integer :: ierr, ui

    fname = dir_data//'/neM20Ni01Z002.3d.dump'

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


  subroutine args_init(is_3d, is_4d, is_cache)
    use cla, only: cla_init, cla_register, cla_get, cla_help, cla_int, cla_flag, cla_key_present;

    logical, intent(out) :: is_3d, is_4d, is_cache;

    !  Init commang arguments
    call cla_init();
    call cla_register('-3d', 'Run 3d test',  cla_flag, 'f');
    call cla_register('-4d', 'Run 4d test',  cla_flag, 'f');
    call cla_register('--cache', 'Cached results',  cla_flag, 'f');
    
    call cla_register('-h',  'Print this help',  cla_flag, 'f');

    is_3d = cla_key_present('-3d');
    is_4d = cla_key_present('-4d');
    is_cache = cla_key_present('--cache');

    if( cla_key_present('-h') ) then;
        call cla_help();
        stop;
    endif;   

  end subroutine args_init

end program main_value
