module kinds
     use,intrinsic :: iso_fortran_env, only : int16, int32, real32, real64
     ! use iso_fortran_env, only: int32, int64, sp=>real32, dp=>real64

     implicit none

     private
     public :: sp, dp, qp, isp, idp, ip, i4, i8
     public :: kinds_write_info

     ! integer, parameter :: sp = selected_real_kind(p=4)  ! standart
     ! integer, parameter :: sp = selected_real_kind(p=5)  ! same as MESA
     ! integer, parameter :: dp = selected_real_kind(p=15)
     integer, parameter :: sp = real32 ! kind(1.)  
     integer, parameter :: dp = real64 ! kind(1.d0)
     ! integer, parameter :: sp = SELECTED_REAL_KIND ( 6, 30 )
     ! integer, parameter :: dp = SELECTED_REAL_KIND ( 14, 200 )
     integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp)) 
     integer, parameter :: i4 = int16 !  selected_int_kind(9)
     integer, parameter :: i8 = int32 ! selected_int_kind(14)
     integer, parameter :: isp = i4 !SELECTED_INT_KIND (8)
     integer, parameter :: idp = i8 !SELECTED_INT_KIND (18)
     integer, parameter :: ip = idp ! selected_int_kind(14)

     CONTAINS

     subroutine kinds_write_info (iw)
          integer, intent(in) :: iw

          write( iw, '( /, T2, A )' ) 'DATA TYPE INFORMATION:'

          write( iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E15.8) )' ) &
               'REAL: Data type name:', 'default', '      Kind value:', KIND ( 0.d0 ), &
               '      Precision:', PRECISION ( 0.d0 ), &
               '      Smallest non-negligible quantity relative to 1:', &
               EPSILON ( 0.d0 ), &
               '      Smallest positive number:', TINY ( 0.d0 ), &
               '      Largest representable number:', HUGE ( 0.d0 )
          write( iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E15.8) )' ) &
               'REAL: Data type name:', 'dp', '      Kind value:', KIND ( 0.0_dp ), &
               '      Precision:', PRECISION ( 0.0_dp ), &
               '      Smallest non-negligible quantity relative to 1:', &
               EPSILON ( 0.0_dp ), &
               '      Smallest positive number:', TINY ( 0.0_dp ), &
               '      Largest representable number:', HUGE ( 0.0_dp )
          write( iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E15.8) )' ) &
               '      Data type name:', 'sp', '      Kind value:', KIND ( 0.0_sp ), &
               '      Precision:', PRECISION ( 0.0_sp ), &
               '      Smallest non-negligible quantity relative to 1:', &
               EPSILON ( 0.0_sp ), &
               '      Smallest positive number:', TINY ( 0.0_sp ), &
               '      Largest representable number:', HUGE ( 0.0_sp )
          write( iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E32.16E4) )' ) &
               '      Data type name:', 'qp', '      Kind value:', KIND ( 0.0_qp ), &
               '      Precision:', PRECISION ( 0.0_qp ), &
               '      Smallest non-negligible quantity relative to 1:', &
               EPSILON ( 0.0_qp ), &
               '      Smallest positive number:', TINY ( 0.0_qp ), &
               '      Largest representable number:', HUGE ( 0.0_qp )
          write( iw, '( /,T2,A,T72,A,4(/,T2,A,T61,I20) )' ) &
               'integer: Data type name:', '(default)', '         Kind value:', &
               KIND ( 0 ), &
               '         Bit size:', BIT_SIZE ( 0 ), &
               '         Largest representable number:', HUGE ( 0 )
          write( iw, '( /,T2,A,T72,A,4(/,T2,A,T61,I20) )' ) &
               '      Data type name:', 'isp', '      Kind value:', KIND ( 0_isp ), &
               '         Bit size:', BIT_SIZE ( 0_isp ), &
               '      Largest representable number:', HUGE ( 0_isp )
          write( iw, '( /,T2,A,T72,A,/,T2,A,T75,I6,/ )' ) &
               'LOGICAL: Data type name:', '(default)', &
               '         Kind value:', KIND ( .TRUE. )
          write( iw, '( /,T2,A,T72,A,/,T2,A,T75,I6,/ )' ) &
               'CHARACTER: Data type name:', '(default)', &
               '           Kind value:', KIND ( 'C' )

     end subroutine kinds_write_info

end module kinds