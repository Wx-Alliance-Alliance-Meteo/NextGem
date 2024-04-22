module app
    use, intrinsic :: iso_c_binding
    implicit none

    enum, bind(C) 
       enumerator :: APP_VERBATIM=-1, APP_ALWAYS=0, APP_FATAL=1, APP_SYSTEM=2, APP_ERROR=3, APP_WARNING=4, APP_INFO=5,          &
          APP_TRIVIAL=6, APP_DEBUG=7, APP_EXTRA=8, APP_QUIET=9
       enumerator :: APP_MAIN=0, APP_LIBRMN=1, APP_LIBFST=2, APP_LIBWB=3, APP_LIBGMM=4, APP_LIBVGRID=5, APP_LIBINTERPV=6,       &
          APP_LIBGEOREF=7, APP_LIBRPNMPI=8, APP_LIBIRIS=9, APP_LIBIO=10, APP_LIBMDLUTIL=11, APP_LIBGEMDYN=12, APP_LIBRPNPHY=13, &
          APP_LIBMIDAS=14, APP_LIBEER=15, APP_LIBTDPACK=16
       enumerator :: APP_MASTER=0, APP_THREAD=1
    end enum
    
    type(C_PTR) :: app_ptr          !Global (opaque) app structure pointer
    integer :: app_status           !To recuperate application status
    character(len=4096) :: app_msg  !String to write output messages     
    character(len=*) , parameter :: EOL = char(13)//char(11)
 
    interface

    ! 
    ! Bindings using C adapters
    !

!   TApp *App_Init(int Type,char* Name,char* Version,char* Desc,char* Stamp);
    type(C_PTR) FUNCTION app_init4fortran(type,name,version,desc,stamp) BIND(C,name="App_Init")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: type
        character(C_CHAR), dimension(*) :: name
        character(C_CHAR), dimension(*) :: version
        character(C_CHAR), dimension(*) :: desc
        character(C_CHAR), dimension(*) :: stamp
    end FUNCTION

!   void  App_Free(void);
    SUBROUTINE app_free() BIND(C, name="App_Free")
        use, intrinsic :: iso_c_binding
    end SUBROUTINE

!   void App_LibRegister(char *Lib,char *Version) {
    SUBROUTINE app_libregister(lib,version) BIND(C, name="App_LibRegister")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: lib
        character(C_CHAR), dimension(*) :: version
    end SUBROUTINE

!   void  App_Start(void);
    SUBROUTINE app_start() BIND(C, name="App_Start")
        use, intrinsic :: iso_c_binding
    end SUBROUTINE

!   int   App_End(int Status);
    integer(C_INT) FUNCTION app_end(status) BIND(C, name="App_End")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: status
    end FUNCTION

!    void  App_Log(TApp_LogLevel Level,const char *Format,...);
    SUBROUTINE app_log4fortran(level,msg) BIND(C, name="App_Log4Fortran")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: level
        character(kind=C_CHAR), dimension(*), intent(in) :: msg
    end SUBROUTINE

!    void  Lib_Log(TApp_Lib Lib,TApp_LogLevel Level,const char *Format,...);
    SUBROUTINE lib_log4Fortran(lib,level,msg) BIND(C, name="Lib_Log4Fortran")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: level
        integer(C_INT), value :: lib
        character(kind=C_CHAR), dimension(*), intent(in) :: msg
    end SUBROUTINE
    
!   void  App_Progress(float Percent,const char *Format,...);

!   int   App_LogLevel(char *Level);
    integer(C_INT) FUNCTION app_loglevel(level) BIND(C, name="App_LogLevel")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=C_CHAR), dimension(*), intent(in):: level
    end FUNCTION

    !   int App_LogLevelNo(TApp_LogLevel Val) {
    integer(C_INT) FUNCTION app_loglevelno(levelno) BIND(C, name="App_LogLevelNo")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: levelno
    end FUNCTION

    !   int Lib_LogLevel(TApp_Lib Lib,char *Val) {
    integer(C_INT) FUNCTION lib_loglevel(lib,level) BIND(C, name="Lib_LogLevel")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: lib
        character(kind=C_CHAR), dimension(*), intent(in) :: level
    end FUNCTION

    !   int Lib_LogLevelNo(TApp_Lib Lib,TApp_LogLevel Val) {
    integer(C_INT) FUNCTION lib_loglevelno(lib,levelno) BIND(C, name="Lib_LogLevelNo")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: lib
        integer(C_INT), value :: levelno
    end FUNCTION

    !    int   App_ParseArgs(TApp_Arg *AArgs,int argc,char *argv[],int Flags);
!    int   App_ParseInput(void *Def,char *File,TApp_InputParseProc *ParseProc);

!   int   App_ParseBool(char *Param,char *Value,char *Var);
    logical(C_BOOL) FUNCTION app_parsebool(param,value,var) BIND(C, name="App_ParseBool")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=C_CHAR), dimension(*), intent(in) :: param
        character(kind=C_CHAR), dimension(*), intent(in) :: value
        type(C_PTR), intent(out) :: var
    end FUNCTION

!   int   App_ParseDate(char *Param,char *Value,time_t *Var);
    integer(C_INT) FUNCTION app_parsedate(param,value,var) BIND(C, name="App_ParseDate")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=C_CHAR), dimension(*), intent(in) :: param
        character(kind=C_CHAR), dimension(*), intent(in) :: value
        integer(C_LONG), intent(out) :: var
    end FUNCTION

!   int   App_ParseDateSplit(char *Param,char *Value,int *Year,int *Month,int *Day,int *Hour,int *Min);
    integer(C_INT) FUNCTION app_parsedatesplit(param,value,year,month,day,hour,min) BIND(C, name="App_ParseDateSplit")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=C_CHAR), dimension(*), intent(in) :: param
        character(kind=C_CHAR), dimension(*), intent(in) :: value
        integer(C_INT), intent(out) :: year
        integer(C_INT), intent(out) :: month
        integer(C_INT), intent(out) :: day
        integer(C_INT), intent(out) :: hour
        integer(C_INT), intent(out) :: min
    end FUNCTION

!   int   App_ParseCoords(char *Param,char *Value,double *Lat,double *Lon,int Index);
    integer(C_INT) FUNCTION app_parsecoords(param,value,lat,lon,index) BIND(C, name="App_ParseCoords")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=C_CHAR), dimension(*), intent(in) :: param
        character(kind=C_CHAR), dimension(*), intent(in) :: value
        integer(C_INT), intent(out) :: lat
        integer(C_INT), intent(out) :: lon
        integer(C_INT), value :: index
    end FUNCTION

!   void  App_SeedInit(void);
    SUBROUTINE app_seedinit() BIND(C, name="App_SeedInit")
        use, intrinsic :: iso_c_binding
    end SUBROUTINE

!   char* App_ErrorGet(void);
    type(C_PTR) FUNCTION app_errorget() BIND(C, name="App_ErrorGet")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   int   App_ThreadPlace(void);
    integer(C_INT) FUNCTION app_threadplace() BIND(C, name="App_ThreadPlace")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   void  App_Trap(int Signal);
    SUBROUTINE app_trap(signal) BIND(C, name="App_Trap")
        use, intrinsic :: iso_c_binding
        integer(C_INT), value :: signal
    end SUBROUTINE

!   int   App_IsDone(void); 
    logical(C_BOOL) FUNCTION app_isdone() BIND(C, name="App_IsDone")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   int   App_IsMPI(void);
    logical(C_BOOL) FUNCTION app_ismpi() BIND(C, name="App_IsMPI")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   int   App_IsOMP(void);
    logical(C_BOOL) FUNCTION app_isomp() BIND(C, name="App_IsOMP")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   int   App_IsSingleNode(void);
    logical(C_BOOL) FUNCTION app_issinglenode() BIND(C, name="App_IsSingleNode")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   int   App_IsAloneNode(void);
    logical(C_BOOL) FUNCTION app_isalonenode() BIND(C, name="App_IsAloneNode")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   int   App_NodeGroup(void);
    integer(C_INT) FUNCTION app_nodegroup() BIND(C, name="App_NodeGroup")
        use, intrinsic :: iso_c_binding
    end FUNCTION

!   void App_SetMPIComm(MPI_Comm Comm);
    SUBROUTINE app_setmpicomm(comm) BIND(C, name="App_SetMPIComm")
        use, intrinsic :: iso_c_binding
        integer(C_INT), value :: comm
    end SUBROUTINE
   
end interface

contains
    type(C_PTR) FUNCTION app_init(type,name,version,desc,stamp)
        use, intrinsic :: iso_c_binding
        implicit none
        integer(C_INT), value :: type
        character(len=*) :: name
        character(len=*) :: version
        character(len=*) :: desc
        character(len=*) :: stamp
        app_init = app_init4fortran(type,name//achar(0),version//achar(0),desc//achar(0),stamp//achar(0))
    end FUNCTION
    
    SUBROUTINE app_log(level,msg)
        use, intrinsic :: iso_c_binding
        implicit none
        integer :: level
        integer :: i
        character(len=*) :: msg
        character(len=4097) :: c_msg
               
        i=len_trim(msg)
        c_msg=msg
        c_msg(i+1:i+1)=C_NULL_CHAR
        call app_log4fortran(level,c_msg)
    end SUBROUTINE

    SUBROUTINE lib_log(lib,level,msg)
        use, intrinsic :: iso_c_binding
        implicit none
        integer :: level
        integer :: lib
        integer :: i
        character(len=*) :: msg
        character(len=4097) :: c_msg

        i=len_trim(msg)
        c_msg=msg
        c_msg(i+1:i+1)=C_NULL_CHAR
        call lib_log4fortran(lib,level,c_msg)
    end SUBROUTINE
end module
