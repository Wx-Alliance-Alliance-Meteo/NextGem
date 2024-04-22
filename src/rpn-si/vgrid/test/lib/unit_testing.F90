module Unit_Testing

  private
  public :: ut_report                                   !Report test status

  character(len=64), parameter :: REPORT_FILE_NAME='test_report.txt' !File to create containing report
  integer, parameter :: REPORT_FILE_UNIT=99             !File unit of report file


  interface ut_report
     module procedure report_stat
     module procedure report_logical
  end interface

contains
  
  subroutine report_stat(status,message)
    ! Report on the exit status of a test
    integer, intent(in) :: status                       !Test status
    character(len=*), intent(in), optional :: message   !Message to report on failure

    ! Local variables
    character(len=1024) :: my_message

    ! Set default values
    my_message = ''
    if (present(message)) my_message = ': '//message

    ! Open output file
    open(unit=REPORT_FILE_UNIT,file=trim(REPORT_FILE_NAME))

    ! Check status and report
    if (status == 0) then
       write(REPORT_FILE_UNIT,*) 'ok'
    else
       write(REPORT_FILE_UNIT,*) 'failed'//trim(my_message)
    endif

    ! Close output file
    close(REPORT_FILE_UNIT)
    
    if (status /= 0) then
      error stop 1
    endif

    ! End of subprogram
    return
  end subroutine report_stat

  subroutine report_logical(logic,message)
    ! Report on the exit status of a test
    logical, intent(in) :: logic                        !Assertion result
    character(len=*), intent(in), optional :: message   !Message to report on failure

    ! Local variables
    character(len=1024) :: my_message

    ! Set default values
    my_message = ''
    if (present(message)) my_message = ': '//message

    ! Open output file
    open(unit=REPORT_FILE_UNIT,file=trim(REPORT_FILE_NAME))

    ! Check status and report
    if (logic) then
       write(REPORT_FILE_UNIT,*) 'ok'
    else
       write(REPORT_FILE_UNIT,*) 'failed'//trim(my_message)
    endif

    ! Close output file
    close(REPORT_FILE_UNIT)

    if (.not. logic) then
      error stop 1
    endif

   ! End of subprogram
    return
  end subroutine report_logical

end module Unit_Testing
