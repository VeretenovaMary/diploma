MODULE module_Settings
! <!>  - mark places for check

    !! Real precision definitions:
    integer, parameter:: r16p = selected_real_kind(33,4931) ! 33  digits, range $[\pm 10^{-4931}  ,\pm 10^{+4931}   -1]$
    integer, parameter:: r8p  = selected_real_kind(15,307)  ! 15  digits, range $[\pm 10^{-307}~~ ,\pm 10^{+307}~~  -1]$
    integer, parameter:: r4p  = selected_real_kind(6,37)    ! 6~~~digits, range $[\pm 10^{-37}~~~~,\pm 10^{+37}~~~~ -1]$
    integer, parameter:: r_p  = r8p                         ! default real precision
    !! Integer precision definitions:
    integer, parameter:: i8p  = selected_int_kind(18)       ! range $[-2^{63} ,+2^{63}  -1]$
    integer, parameter:: i4p  = selected_int_kind(9)        ! range $[-2^{31} ,+2^{31}  -1]$
    integer, parameter:: i2p  = selected_int_kind(4)        ! range $[-2^{15} ,+2^{15}  -1]$
    integer, parameter:: i1p  = selected_int_kind(2)        ! range $[-2^{7}~~,+2^{7}~~ -1]$
    integer, parameter:: i_p  = i4p                         ! default integer precision

    real(R8P), parameter :: rdpprm = 1.0e0_r8p ! example of using precision
    character(128) :: OutputDir
    character(8)   :: vtkFmtOut
    character(8)   :: OSType
    integer        :: NumThreads

!     files units:
!    integer, parameter :: cfgUnit = 201

    ! physical measure units and dimensions
    real(R8P) :: dimLen
    real(R8P) :: dimVel
    real(R8P) :: dimDen
    real(R8P) :: dimTem

END MODULE module_Settings


subroutine set_base_settings(cfgName)
    use module_Settings
    use module_TxtRead
    implicit none
    character(*) :: cfgName
    integer      :: ios
    integer(i4p) :: GetUnit,cfgUnit

    ! set OutputDir, vtkFmtOut
    cfgUnit = GetUnit()
    open(cfgUnit,file=cfgName)
    ! read atom names
    !call txt_read_value(cfgUnit,"OutputDir",OutputDir,ios)
    call txt_read_value(cfgUnit,"vtkFmtOut",vtkFmtOut,ios)
    call txt_read_value(cfgUnit,"NumThreads",NumThreads,ios)

    call txt_read_value(cfgUnit,"dimLen",dimLen,ios)
    call txt_read_value(cfgUnit,"dimDen",dimDen,ios)
    call txt_read_value(cfgUnit,"dimTem",dimTem,ios)
    call txt_read_value(cfgUnit,"dimVel",dimVel,ios)

    close(cfgUnit)

!    call SetOSType()
!    write(*,*) "OSType: ",trim(OSType)
!    pause

    return
end subroutine set_base_settings

! - - - - - - - - - - - - - - - - - - - - -
subroutine SetOSType()
    use module_Settings
    implicit none
    ! local:
    integer        :: getcwd,ios
    character(128) :: fpath

    ios = getcwd(fpath)
    if(index(fpath,"\").gt.0) then
        OSType = "_WIN32"
    elseif(index(fpath,"/").gt.0) then
        OSType = "__linux"
    else
        OSType = "undef"
    end if

end subroutine

! - - - - - - - - - - - - - - - - - - - - -
function GetUnit() result(Free_Unit)
    use module_Settings
    !--------------------------------------------------------------------------------------------------------------------------------
    !!The GetUnit function is used for getting a free logic unit. The users of \LIBVTKIO does not know which is
    !!the logical unit: \LIBVTKIO handels this information without boring the users. The logical unit used is safe-free: if the
    !!program calling \LIBVTKIO has others logical units used \LIBVTKIO will never use these units, but will choice one that is free.
    !--------------------------------------------------------------------------------------------------------------------------------

    implicit none

    !--------------------------------------------------------------------------------------------------------------------------------
    integer(i4p):: Free_Unit ! free logic unit
    integer(i4p):: n1        ! counter
    integer(i4p):: ios       ! inquiring flag
    logical(4)::   lopen     ! inquiring flag
    !--------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------
    !!The following is the code snippet of GetUnit function: the units 0, 5, 6, 9 and all non-free units are discarded.
    !!
    !(\doc)codesnippet
    Free_Unit = -1_i4p                                      ! initializing free logic unit
    n1=1_i4p                                                ! initializing counter
    do
        if ((n1/=5_i4p).AND.(n1/=6_i4p).AND.(n1/=9_i4p)) then
            inquire (unit=n1,opened=lopen,iostat=ios)           ! verify logic units
            if (ios==0_i4p) then
                if (.NOT.lopen) then
                    Free_Unit = n1                                  ! assignment of free logic
                    return
                endif
            endif
        endif
        n1=n1+1_i4p                                           ! updating counter
    enddo
    return
    !(doc/)codesnippet
    !!GetUnit function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope to public.
    !--------------------------------------------------------------------------------------------------------------------------------
endfunction GetUnit
