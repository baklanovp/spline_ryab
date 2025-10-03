module cla
  
  implicit none


  ! Command Line Arguments

  ! A key-value parser for the commandline

 ! integer(kind=4), public, parameter :: int_kind = 4
 ! integer(kind=4), public, parameter :: real_kind = 4
 ! integer(kind=4), public, parameter :: ptr_kind = 8

  integer, parameter  :: isp = SELECTED_INT_KIND (8)
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
  
  
  integer(kind=4), public, parameter :: STRLEN = 120
  integer(kind=4), public, parameter :: XSTRLEN = 256
  
  integer(kind=isp), parameter :: &
       cla_int   = 1, &
       cla_float = 2, &
       cla_char  = 3, &
       cla_xchar = 4, & ! NOT IMPLEMENTED
       cla_logical=5, &
       cla_flag  = 6

  character(len=STRLEN), dimension(6) :: cla_kindstr
  character(len=STRLEN), private :: cla_empty
  character(len=STRLEN), dimension(6) :: cla_true_str
 
  type, private :: cla_t
     character(len=STRLEN)  :: key
     character(len=XSTRLEN) :: description
     integer(kind=isp) :: kind
     character(len=STRLEN)  :: default
  end type cla_t

  type(cla_t), private, dimension(:), pointer :: cla_registry
  integer(kind=isp), private :: cla_num

  interface cla_get
     module procedure &
          cla_get_float_r4, &
          cla_get_float_r8, &
          cla_get_int_i4, &
          cla_get_int_i8, &
          cla_get_char, &
          cla_get_logical
  end interface

  contains

    subroutine cla_init
      ! Allocate a zero size registry, just so that it gets
      ! associated.
      cla_num = 0
      allocate(cla_registry(0))
      cla_kindstr(cla_int)     = 'integer'
      cla_kindstr(cla_float)   = 'float'
      cla_kindstr(cla_char)    = 'character'
      cla_kindstr(cla_xchar)   = 'xcharacter' !NOT IMPLEMENTED
      cla_kindstr(cla_logical) = 'logical'
      cla_kindstr(cla_flag)    = 'flag'
      cla_empty='THIS_IS_THE_EMPTY_STRING'
      cla_true_str(1)='true'
      cla_true_str(2)='on'
      cla_true_str(3)='1'
      cla_true_str(4)='t'
      cla_true_str(5)='T'
      cla_true_str(6)='.true.'
    end subroutine cla_init

    subroutine cla_register(key,description,kkind,default)
      character(len=*)  :: key
      character(len=*) :: description
      integer(kind=isp) :: kkind
      character(len=*)  :: default
      type(cla_t), dimension(:), pointer :: cla_registry_tmp
      integer(kind=isp) :: i

      ! This is a dumb way to increase the size of the
      ! registry of command line arguments, but there
      ! should not be so many arguments that either speed
      ! or memory is an issue.
      allocate(cla_registry_tmp(cla_num+1))
      do i=1,cla_num
         cla_registry_tmp(i)%key         = cla_registry(i)%key
         cla_registry_tmp(i)%description = cla_registry(i)%description
         cla_registry_tmp(i)%kind        = cla_registry(i)%kind
         cla_registry_tmp(i)%default     = cla_registry(i)%default
         if (index(trim(key),' ') /= 0) then
            print *,'Error: attempt to register cla key'
            print *,'containing a space character.'
            stop 1
         end if
         if (cla_str_eq(trim(cla_registry(i)%key),trim(key))) then
            print *,'Error: attempt to register cla key which'
            print *,'has already been registered'
            print *,'           i: ',i
            print *,'         key: ',trim(key)
            print *,' conflicting key: ',trim(cla_registry(i)%key)
            print *,' description: ',trim(cla_registry(i)%description)
            print *,'        kind: ',trim(cla_kindstr(cla_registry(i)%kind))
            print *,' default val: ',trim(cla_registry(i)%default)
            stop 2
         end if
      end do
      cla_num = cla_num + 1
      deallocate(cla_registry)
      allocate(cla_registry(cla_num))
      do i=1,cla_num-1
         cla_registry(i)%key         = cla_registry_tmp(i)%key
         cla_registry(i)%description = cla_registry_tmp(i)%description
         cla_registry(i)%kind        = cla_registry_tmp(i)%kind
         cla_registry(i)%default     = cla_registry_tmp(i)%default
      end do
      i = cla_num
      cla_registry(i)%key         = key
      cla_registry(i)%description = description
      cla_registry(i)%kind        = kkind
      cla_registry(i)%default     = default
      deallocate(cla_registry_tmp)
    end subroutine cla_register

    subroutine cla_show
      integer(kind=isp) :: i
      if (cla_num == 0) then
         print *,'There are no command line arguments registered.'
         return
      end if
      print *,'General usage:'
      print *,'  command -[key] [value] | -[flag]'
      print *,'  The key/value pairs must be matched if they appear.'
      print *,'  Key/value pairs and flags may be in any order.'
      print *,' '
      print *,'The following command line arguments are registered for use:'
      do i=1,cla_num
         print *,'---------- i: ',i
         print *,'         key: ',trim(cla_registry(i)%key)
         print *,' description: ',trim(cla_registry(i)%description)
         print *,'        kind: ',trim(cla_kindstr(cla_registry(i)%kind))
         print *,'     default: ',trim(cla_registry(i)%default)
      end do

      print *,' '
      print *,'Also, -?, -h, -H, -help, --help, and --usage are recognized.'
      print *,' '
    end subroutine cla_show

    subroutine cla_help
      integer(kind=isp) :: i
      if (cla_num == 0) then
         print *,'There are no command line arguments registered.'
         return
      end if
      print *,'General usage:'
      print *,'  command -[key] [value] | -[flag]'
      print *,' '
      print *,'The following command line arguments may be used:'
      do i=1,cla_num
         if (cla_registry(i)%kind == cla_flag) then
            write(6,'(i2,a,a)')i,'. ',trim(cla_registry(i)%description)
            print *,'                ',trim(cla_registry(i)%key)
         else
            write(6,'(i2,a,a)')i,'. ',trim(cla_registry(i)%description)
            print *,'                ',trim(cla_registry(i)%key), &
                 ' [',trim(cla_registry(i)%default),']'
         endif
      end do

      print *,' '
      print *,'Also, -?, -h, -H, -help, --help, and --usage are recognized.'
      print *,' '
    end subroutine cla_help

    pure integer function cla_eq(str1,str2)
      implicit none
      character(*), intent(in) :: str1, str2
      cla_eq = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
    end function cla_eq

    pure logical function cla_str_eq(str1,str2)
      implicit none
      character(*), intent(in) :: str1, str2
      integer :: str_test
      str_test = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
      cla_str_eq = .false.
      if (str_test /= 0) cla_str_eq = .true.
    end function cla_str_eq

    subroutine cla_validate
      implicit none
      character(len=STRLEN) :: arg
      character(len=STRLEN)  :: value, key
      integer(kind=isp) :: ncla, k, kk
!      integer :: command_argument_count
!      external command_argument_count

      ncla = command_argument_count()
      
      if (ncla == 0) return
      
      ! First check for -?, -h, -H, -help, or --help flags.
      call get_command_argument(1,arg)
      key = trim(arg)
      if (cla_str_eq(trim(key),'-h')      .or. &
          cla_str_eq(trim(key),'-?')      .or. &
          cla_str_eq(trim(key),'-H')      .or. &
          cla_str_eq(trim(key),'-help')   .or. &
          cla_str_eq(trim(key),'--help')  .or. &
          cla_str_eq(trim(key),'--usage')      &
          ) then
         call cla_help
         stop
      endif

!      if (mod(ncla,2) /= 0) then
!         print *,'From cla_validate, ncla = ',ncla
!         print *,'Usage:'
!         print *,'  command -[key] [value]'
!         print *,'  The key/value pairs must be matched.'
!         print *,'Did you forget to put a space between a key/value pair?'
!         stop
!      endif
      
      k=1
!      do k=1,ncla,2
      do while (k <= ncla)
         call get_command_argument(k,arg)
         
         key = trim(arg)
         value = cla_empty
         do kk=1,cla_num
            if (cla_str_eq(trim(cla_registry(kk)%key),trim(key))) then
               value = trim(cla_registry(kk)%default)
               if (cla_registry(kk)%kind == cla_flag) then
                  k = k+1
               else
                  k = k+2
                  if ( k == ncla+2 ) then
                     ! I wonder if there is any easy way to verify the type of the
                     ! input argument? Or, I could check that the values do not
                     ! begin with a minus sign.
                     print *,'Error: A key/value pair is incomplete.'
                     print *,'Check that key/value pairs are correctly matched.'
                     print *,' '
                     call cla_show
                     print *,' '
                     print *,'Error: A key/value pair is incomplete.'
                     stop     
                  endif
               endif
               exit
            end if
         end do ! kk=1,cla_num loop over possible cla
         if (cla_eq(trim(value),trim(cla_empty)) /= 0) then
            print *,'Error: You used an unknown command line argument: ',trim(key)
            print *,'Are your key/value pairs matched?'
            print *,'Do not use -h, --help, etc. These are reserved flags.'
            print *,' '
            call cla_show
            print *,' '
            print *,'Error: You used an unknown command line argument: ',trim(key)
            stop
         endif
      enddo ! k=1,ncla,2 loop over used cla
    end subroutine cla_validate
    
    logical function cla_key_present(key)
      implicit none
      character(len=STRLEN) :: arg
      character(len=*)  :: key
      character(len=STRLEN)  :: value
      integer(kind=isp) :: ncla, k, kk
!      integer :: command_argument_count
!      external command_argument_count

      !     Loop over the command line arguments to assign to
      !     value.
      !     Note that no error is reported if the key was NOT
      !     registered, but it is present on the command line.
 
      cla_key_present = .false.

!      print *,'Calling cla_key_present with key = ',trim(key)
      value = trim(cla_empty)
      do kk=1,cla_num
         ! must test for exact match, not just substring
         if (cla_str_eq(trim(cla_registry(kk)%key),trim(key)) ) then
            value = trim(cla_registry(kk)%default)
            exit
         end if
      end do
      
      if (index(trim(value),trim(cla_empty)) /= 0) then
         print *,'Error: You used an unknown command line argument: ',trim(key)
         call cla_show
         stop 5
      endif
      
      ncla = command_argument_count()
!      print *,'Found ncla = ',ncla

!      if (mod(ncla,2) /= 0) then
!         print *,'Usage:'
!         print *,'  command -[key] [value]'
!         print *,'  The key/value pairs must be matched.'
!         stop
!      endif
      
      if (ncla == 0) return
      
      do k=1,ncla
         call get_command_argument(k,arg)
         ! test for exact match
         if (cla_str_eq(trim(arg),trim(key))) then
            !          print *,'    Found: ',trim(key),' = ',trim(value),' on command line.'
            cla_key_present = .true.
            return
         endif
      enddo
      
    end function cla_key_present

    subroutine cla_get_char(key,value)
      implicit none
      character(len=*), intent(in)  :: key
      character(len=STRLEN), intent(out)  :: value
      character(len=STRLEN) :: arg
      integer(kind=isp) :: ncla, k, kkind
!      integer :: command_argument_count
!      external command_argument_count

      !     Loop over the command line arguments to assign to
      !     value.
      !     Note that no error is reported if the key was NOT
      !     registered, but it is present on the command line.
      
      kkind = 0
      value = trim(cla_empty)
      do k=1,cla_num
         ! must test for exact match, not just substring
         if (cla_str_eq(trim(cla_registry(k)%key),trim(key))) then
            value = trim(cla_registry(k)%default)
            kkind = cla_registry(k)%kind
         end if
      end do
      
      if (index(trim(value),trim(cla_empty)) /= 0) then
         print *,'Error: You used and unknown command line argument: ',trim(key)
         call cla_show
         stop 5
      endif
      
      ncla = command_argument_count()
!      if (mod(ncla,2) /= 0) then
!         print *,'From cla_get_char, ncla = ',ncla
!         print *,'Usage:'
!         print *,'  command -[key] [value]'
!         print *,'  The key/value pairs must be matched.'
!         stop
!      endif
      
      if (ncla == 0) return
      
      do k=1,ncla
         call get_command_argument(k,arg)
         ! test for exact match
         if (cla_str_eq(trim(arg),trim(key))) then
            if (kkind == cla_flag) then
               value = 't'  ! '
               return
            else
               call get_command_argument(k+1,arg)
               value = trim(arg)
                  ! print *,'    Found: ',trim(key),' = ',trim(value),' on command line.'
               return
            endif
         endif
      enddo
      
    end subroutine cla_get_char
    
!  subroutine cla_get_float(key,float_value)
!    implicit none
!    character(len=*)  :: key
!    character(len=STRLEN)  :: value
!    real(kind=sp)   :: float_value
!    integer(kind=isp) :: ncla, k
!    
!    call cla_get_char(key,value)
!    if (index(trim(value),trim(cla_empty)) == 0) read(value,*)float_value
!  end subroutine cla_get_float

  subroutine cla_get_float_r4(key,float_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    real(kind=sp)           :: float_value
   !  integer(kind=isp) :: ncla, k
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*)float_value
  end subroutine cla_get_float_r4

  subroutine cla_get_float_r8(key,float_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    real(kind=dp)           :: float_value
   !  integer(kind=isp) :: ncla, k
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*)float_value
  end subroutine cla_get_float_r8

!  subroutine cla_get_int(key,int_value)
!    implicit none
!    character(len=*)  :: key
!    character(len=STRLEN)  :: value
!    integer(kind=isp)   :: int_value
!    integer(kind=isp) :: ncla, k
!    
!    call cla_get_char(key,value)
!    if (index(trim(value),trim(cla_empty)) == 0) read(value,*)int_value
!  end subroutine cla_get_int

  subroutine cla_get_int_i4(key,int_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    integer(kind=4)        :: int_value
   !  integer(kind=isp) :: ncla, k
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*)int_value
  end subroutine cla_get_int_i4

  subroutine cla_get_int_i8(key,int_value)
    implicit none
    character(len=*)       :: key
    character(len=STRLEN)  :: value
    integer(kind=dp)        :: int_value
   !  integer(kind=isp) :: ncla, k
    
    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) read(value,*)int_value
  end subroutine cla_get_int_i8

  subroutine cla_get_logical(key,logical_value)
    implicit none
    character(len=*)  :: key
    character(len=STRLEN)  :: value
    logical :: logical_value
    integer(kind=isp) :: k
    
    logical_value = .false.

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) then
       do k=1,6
          if (index(trim(value),trim(cla_true_str(k))) /= 0) then
             logical_value = .true.
          endif
       end do
    end if
  end subroutine cla_get_logical

  subroutine cla_get_flag(key,logical_value)
    implicit none
    character(len=*)  :: key
    character(len=STRLEN)  :: value
    logical :: logical_value
    integer(kind=isp) ::  k
    
    logical_value = .false.

    call cla_get_char(key,value)
    if (index(trim(value),trim(cla_empty)) == 0) then
       do k=1,6
          if (index(trim(value),trim(cla_true_str(k))) /= 0) then
             logical_value = .true.
          endif
       end do
    end if
  end subroutine cla_get_flag


end module cla
