! The program based on strings module from:
! http://gbenthien.net/strings/index.php
module module_TxtRead
! <!>  - mark places for check

    integer, parameter :: N_ch_max  = 2048               ! maximal number of characters in dimension
    integer, parameter :: N_str_ch  = 128                ! maximal number of characters in structure argument
    integer, parameter :: N_str_arg = 128                ! maximal number of arguments in structure
    private :: N_ch_max
    private :: N_str_ch
    private :: N_str_arg

    ! Real kinds
    integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
    integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

    ! Integer kinds
    integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
    integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

    !Complex kinds
    integer, parameter :: kc4 = kr4                            ! single precision complex
    integer, parameter :: kc8 = kr8                            ! double precision complex


!    private :: value_str,value_dr,value_sr,value_di,value_si
!    private :: write_dr,write_sr,write_di,write_si
!    private :: writeq_dr,writeq_sr,writeq_di,writeq_si

    !private :: txt_findvar
    private :: txt_read_arr_di,txt_read_arr_si
    private :: txt_read_arr_dr,txt_read_arr_sr,txt_read_arr_str
    private :: txt_read_val_di,txt_read_val_si
    private :: txt_read_val_dr,txt_read_val_sr,txt_read_val_str

    private :: txt_read_struct_arr_di,txt_read_struct_arr_si
    private :: txt_read_struct_arr_dr,txt_read_struct_arr_sr
    private :: txt_read_struct_arr_str
    private :: txt_read_struct_val_di,txt_read_struct_val_si
    private :: txt_read_struct_val_dr,txt_read_struct_val_sr
    private :: txt_read_struct_val_str
    private :: txt_read_struct_data_string

    private :: txt_findvar

    public :: txt_read_array
    public :: txt_read_value
    public :: txt_read_struct_value
    public :: txt_read_struct_array
    public :: slash_bck2frw
    public :: slash_frw2bck
    public :: slash_win2lin
    public :: slash_lin2win
    public :: slash_end_del
    public :: slash_end_set
!    public :: readline

    interface txt_read_value  ! Generic operator for converting a number string to a
        ! number. Calling syntax is 'call value(numstring,number,ios)'
        ! where 'numstring' is a number string and 'number' is a
        ! real number or an integer (single or double precision).
        module procedure txt_read_val_str
        module procedure txt_read_val_dr
        module procedure txt_read_val_sr
        module procedure txt_read_val_di
        module procedure txt_read_val_si
    end interface

    interface txt_read_array  ! Generic operator for converting a number string to a
        ! number. Calling syntax is 'call value(numstring,number,ios)'
        ! where 'numstring' is a number string and 'number' is a
        ! real number or an integer (single or double precision).
        module procedure txt_read_arr_str
        module procedure txt_read_arr_dr
        module procedure txt_read_arr_sr
        module procedure txt_read_arr_di
        module procedure txt_read_arr_si
    end interface

    interface txt_read_struct_value  ! Generic operator for converting a number string to a
        ! number. Calling syntax is 'call value(numstring,number,ios)'
        ! where 'numstring' is a number string and 'number' is a
        ! real number or an integer (single or double precision).
        module procedure txt_read_struct_val_str
        module procedure txt_read_struct_val_dr
        module procedure txt_read_struct_val_sr
        module procedure txt_read_struct_val_di
        module procedure txt_read_struct_val_si
    end interface

    interface txt_read_struct_array  ! Generic operator for converting a number string to a
        ! number. Calling syntax is 'call value(numstring,number,ios)'
        ! where 'numstring' is a number string and 'number' is a
        ! real number or an integer (single or double precision).
        module procedure txt_read_struct_arr_str
        module procedure txt_read_struct_arr_dr
        module procedure txt_read_struct_arr_sr
        module procedure txt_read_struct_arr_di
        module procedure txt_read_struct_arr_si
    end interface

    interface txt_value  ! Generic operator for converting a number string to a
        ! number. Calling syntax is 'call value(numstring,number,ios)'
        ! where 'numstring' is a number string and 'number' is a
        ! real number or an integer (single or double precision).
        module procedure value_str
        module procedure value_dr
        module procedure value_sr
        module procedure value_di
        module procedure value_si
    end interface

    interface writenum  ! Generic  interface for writing a number to a string. The
        ! number is left justified in the string. The calling syntax
        ! is 'call writenum(number,string,format)' where 'number' is
        ! a real number or an integer, 'string' is a character string
        ! containing the result, and 'format' is the format desired,
        ! e.g., 'e15.6' or 'i5'.
        module procedure write_dr
        module procedure write_sr
        module procedure write_di
        module procedure write_si
    end interface

    interface writeq  ! Generic interface equating a name to a numerical value. The
        ! calling syntax is 'call writeq(unit,name,value,format)' where
        ! unit is the integer output unit number, 'name' is the variable
        ! name, 'value' is the real or integer value of the variable,
        ! and 'format' is the format of the value. The result written to
        ! the output unit has the form <name> = <value>.
        module procedure writeq_dr
        module procedure writeq_sr
        module procedure writeq_di
        module procedure writeq_si
    end interface


    !**********************************************************************

contains

    !**********************************************************************
    subroutine slash_end_set(fpath)
        implicit none
        character(*) :: fpath
        integer      :: ipos

        ipos = len(trim(fpath))
        if(.not.ANY(['/','\'].eq.fpath(ipos:ipos))) then
            fpath(ipos+1:ipos+1) = '\'
        endif
        return
    end subroutine slash_end_set

    !**********************************************************************
    subroutine slash_end_del(fpath)
        implicit none
        character(*) :: fpath
        integer      :: ipos

        ipos = len(trim(fpath))
        if(ANY(['/','\'].eq.fpath(ipos:ipos))) then
            fpath(ipos:ipos) = ' '
        endif
        return
    end subroutine slash_end_del

    !**********************************************************************
    subroutine slash_bck2frw(fpath)
        implicit none
        character(*) :: fpath
        integer      :: ipos

        do ipos = 1,len(fpath)
            if(fpath(ipos:ipos).eq.'\') then
                fpath(ipos:ipos) = '/'
            endif
        enddo

        ipos = index(fpath,"//")
        do while(ipos.gt.0)
            fpath = fpath(1:ipos)//trim(fpath(ipos+2:))
            ipos = index(fpath,"//")
        enddo

        return
    end subroutine slash_bck2frw

    !**********************************************************************
    subroutine slash_frw2bck(fpath)
        implicit none
        character(*) :: fpath
        integer      :: ipos

        do ipos = 1,len(fpath)
            if(fpath(ipos:ipos).eq.'/') then
                fpath(ipos:ipos) = '\'
            endif
        enddo

        ipos = index(fpath,"\\")
        do while(ipos.gt.0)
            fpath = fpath(1:ipos)//trim(fpath(ipos+2:))
            ipos = index(fpath,"\\")
        enddo

        return
    end subroutine slash_frw2bck

    !**********************************************************************
    subroutine slash_win2lin(fpath)
        implicit none
        character(*) :: fpath
        integer      :: ipos

        do ipos = 1,len(fpath)
            if(fpath(ipos:ipos).eq.'\') then
                fpath(ipos:ipos) = '/'
            endif
        enddo

        ipos = index(fpath,"//")
        do while(ipos.gt.0)
            fpath = fpath(1:ipos)//trim(fpath(ipos+2:))
            ipos = index(fpath,"//")
        enddo

        return
    end subroutine slash_win2lin

    !**********************************************************************
    subroutine slash_lin2win(fpath)
        implicit none
        character(*) :: fpath
        integer      :: ipos

        do ipos = 1,len(fpath)
            if(fpath(ipos:ipos).eq.'/') then
                fpath(ipos:ipos) = '\'
            endif
        enddo

        ipos = index(fpath,"\\")
        do while(ipos.gt.0)
            fpath = fpath(1:ipos)//trim(fpath(ipos+2:))
            ipos = index(fpath,"\\")
        enddo

        return
    end subroutine slash_lin2win

    !**********************************************************************
    subroutine slash_correct(fpath,ostype)
        implicit none
        character(*)  :: ostype
        character(*)  :: fpath
        character(16) :: sbuf

        sbuf = trim(ostype)
        sbuf = lowercase(sbuf)
        call slash_lin2win(fpath)
        if(index(sbuf,'win').gt.0) then
          call slash_lin2win(fpath)
        elseif(index(sbuf,'lin').gt.0) then
          call slash_win2lin(fpath)
        else
          call slash_lin2win(fpath)
        endif

        return
    end subroutine


    !**********************************************************************
    subroutine txt_read_arr_str(cfgUnit,arrTit,array,num,ios)
        implicit none
        integer                    :: cfgUnit
        character(*)               :: arrTit
        character(*), dimension(:) :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        character(N_ch_max), dimension(size(array)) :: strdim
        character(N_ch_max)    :: str
        integer                :: i,ipos,imatch

        call txt_findvar(cfgUnit,arrTit,str)
        ipos = 1
        imatch = 0
        call match(str,ipos,imatch)
        str = str(ipos+1:imatch-1)
        call parse(str,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_arr_str

    !**********************************************************************
    subroutine txt_read_arr_dr(cfgUnit,arrTit,array,num,ios)
        implicit none
        integer                    :: cfgUnit
        character(*)               :: arrTit
        real(kr8), dimension(:)    :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        character(N_ch_max), dimension(size(array)) :: strdim
        character(N_ch_max)   :: str
        integer          :: i,ipos,imatch

        call txt_findvar(cfgUnit,arrTit,str)
        ipos = 1
        imatch = 0
        call match(str,ipos,imatch)
        str = str(ipos+1:imatch-1)
        call parse(str,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_arr_dr

    !**********************************************************************

    subroutine txt_read_arr_sr(cfgUnit,arrTit,array,num,ios)
        implicit none
        integer                    :: cfgUnit
        character(*)               :: arrTit
        real(kr4), dimension(:)    :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        character(N_ch_max), dimension(size(array)) :: strdim
        character(N_ch_max)   :: str
        integer          :: i,ipos,imatch

        call txt_findvar(cfgUnit,arrTit,str)
        ipos = 1
        imatch = 0
        call match(str,ipos,imatch)
        str = str(ipos+1:imatch-1)
        call parse(str,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_arr_sr

    !**********************************************************************

    subroutine txt_read_arr_di(cfgUnit,arrTit,array,num,ios)
        implicit none
        integer                    :: cfgUnit
        character(*)               :: arrTit
        integer(ki8), dimension(:) :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        character(N_ch_max), dimension(size(array)) :: strdim
        character(N_ch_max)   :: str
        integer          :: i,ipos,imatch

        call txt_findvar(cfgUnit,arrTit,str)
        ipos = 1
        imatch = 0
        call match(str,ipos,imatch)
        str = str(ipos+1:imatch-1)
        call parse(str,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_arr_di

    !**********************************************************************

    subroutine txt_read_arr_si(cfgUnit,arrTit,array,num,ios)
        implicit none
        integer                    :: cfgUnit
        character(*)               :: arrTit
        integer(ki4), dimension(:) :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        character(N_ch_max), dimension(size(array)) :: strdim
        character(N_ch_max)   :: str
        integer          :: i,ipos,imatch

        call txt_findvar(cfgUnit,arrTit,str)
        ipos = 1
        imatch = 0
        call match(str,ipos,imatch)
        str = str(ipos+1:imatch-1)
        call parse(str,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_arr_si

    !**********************************************************************

    subroutine txt_read_val_str(cfgUnit,varTit,varval,ios)
        implicit none
        integer :: cfgUnit
        character(*)   :: varTit
        character(*)   :: varval
        integer        :: ios
        character(N_ch_max)   :: str

        call txt_findvar(cfgUnit,varTit,str)
        call txt_value(str,varval,ios)
        return
    end subroutine txt_read_val_str

    !**********************************************************************

    subroutine txt_read_val_dr(cfgUnit,varTit,varval,ios)
        implicit none
        integer :: cfgUnit
        character(*)   :: varTit
        real(kr8)      :: varval
        integer        :: ios
        character(N_ch_max)   :: str

        call txt_findvar(cfgUnit,varTit,str)
        call txt_value(str,varval,ios)

        return
    end subroutine txt_read_val_dr

    !**********************************************************************

    subroutine txt_read_val_sr(cfgUnit,varTit,varval,ios)
        implicit none
        integer :: cfgUnit
        character(*)   :: varTit
        real(kr4)      :: varval
        integer        :: ios
        character(N_ch_max)   :: str

        call txt_findvar(cfgUnit,varTit,str)
        call txt_value(str,varval,ios)

        return
    end subroutine txt_read_val_sr

    !**********************************************************************

    subroutine txt_read_val_di(cfgUnit,varTit,varval,ios)
        implicit none
        integer :: cfgUnit
        character(*)   :: varTit
        integer(ki8)   :: varval
        integer        :: ios
        character(N_ch_max)   :: str

        call txt_findvar(cfgUnit,varTit,str)
        call txt_value(str,varval,ios)

        return
    end subroutine txt_read_val_di

    !**********************************************************************

    subroutine txt_read_val_si(cfgUnit,varTit,varval,ios)
        implicit none
        integer :: cfgUnit
        character(*)   :: varTit
        integer(ki4)   :: varval
        integer        :: ios
        character(N_ch_max)   :: str

        call txt_findvar(cfgUnit,varTit,str)
        call txt_value(str,varval,ios)

        return
    end subroutine txt_read_val_si

    !**********************************************************************

    subroutine txt_findvar(nunitr,varTit,valstr)
        implicit none

        integer        :: nunitr
        character(*)   :: varTit
        character(*)   :: valstr
        integer        :: ios
        character(N_ch_max) :: bufstr,namestr

        rewind(nunitr)
        ios = 0
        namestr = ''
        do while(trim(namestr).ne.trim(varTit))
            bufstr = ''
            call readline(nunitr,bufstr,ios)
            call split(bufstr,'=',namestr)
            if(ios /= 0) then
                write(*,*) 'config error: ',varTit, ' not found'
                stop
            endif
        enddo
        if(len(valstr).lt.len(trim(bufstr))) then
            write(*,*) 'error: ','cant read value, len(valstr) < len(value)'
            stop
        endif
        valstr = trim(bufstr)

    end subroutine txt_findvar



    subroutine parse(str,delims,args,nargs)

        ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
        ! the delimiters contained in the string 'delims'. Preceding a delimiter in
        ! 'str' by a backslash (\) makes this particular instance not a delimiter.
        ! The integer output variable nargs contains the number of arguments found.

        character(len=*) :: str,delims
        character(len=len_trim(str)) :: strsav
        character(len=*),dimension(:) :: args

        strsav=str
        call compact(str)
        na=size(args)
        do i=1,na
            args(i)=' '
        end do
        nargs=0
        lenstr=len_trim(str)
        if(lenstr==0) return
        k=0

        do
            if(len_trim(str) == 0) exit
            nargs=nargs+1
            call split(str,delims,args(nargs))
            call removebksl(args(nargs))
        end do
        str=strsav

    end subroutine parse

    !**********************************************************************

    subroutine compact(str)

        ! Converts multiple spaces and tabs to single spaces; deletes control characters;
        ! removes initial spaces.

        character(len=*):: str
        character(len=1):: ch
        character(len=len_trim(str)):: outstr

        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=' '
        isp=0
        k=0

        do i=1,lenstr
            ch=str(i:i)
            ich=iachar(ch)

            select case(ich)

                case(9,32)     ! space or tab character
                    if(isp==0) then
                        k=k+1
                        outstr(k:k)=' '
                    end if
                    isp=1

                case(33:)      ! not a space, quote, or control character
                    k=k+1
                    outstr(k:k)=ch
                    isp=0

            end select

        end do

        str=adjustl(outstr)

    end subroutine compact

    !**********************************************************************

    subroutine removesp(str)

        ! Removes spaces, tabs, and control characters in string str

        character(len=*):: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr

        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=' '
        k=0

        do i=1,lenstr
            ch=str(i:i)
            ich=iachar(ch)
            select case(ich)
                case(0:32)  ! space, tab, or control character
                    cycle
                case(33:)
                    k=k+1
                    outstr(k:k)=ch
            end select
        end do

        str=adjustl(outstr)

    end subroutine removesp

    !**********************************************************************

    subroutine value_str(str,valstr,ios)
        implicit none
        ! Converts quoted string to a string

        character(len=*) :: str
        character(len=*) :: valstr
        integer          :: ios

        read(str,*,iostat=ios) valstr

        return
    end subroutine value_str

    !**********************************************************************

    subroutine value_dr(str,rnum,ios)

        ! Converts number string to a double precision real number

        character(len=*)::str
        real(kr8)::rnum
        integer :: ios

        ilen=len_trim(str)
        ipos=scan(str,'Ee')
        if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
            ios=3
            return
        end if
        read(str,*,iostat=ios) rnum

    end subroutine value_dr

    !**********************************************************************

    subroutine value_sr(str,rnum,ios)

        ! Converts number string to a single precision real number

        character(len=*)::str
        real(kr4) :: rnum
        real(kr8) :: rnumd

        call value_dr(str,rnumd,ios)
        if( abs(rnumd) > huge(rnum) ) then
            ios=15
            return
        end if
        if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
        rnum=sngl(rnumd)

    end subroutine value_sr

    !**********************************************************************

    subroutine value_di(str,inum,ios)

        ! Converts number string to a double precision integer value

        character(len=*)::str
        integer(ki8) :: inum
        real(kr8)    :: rnum
        integer(ki8) :: hnum     ! huge integer
        real(kr8)    :: rhnum    ! huge integer converted to real

        call value_dr(str,rnum,ios)
        hnum  = huge(inum)
        rhnum = real(hnum,kr8)
        if(abs(rnum)>rhnum) then
            ios=15
            return
        end if
        inum=nint(rnum,ki8)

    end subroutine value_di

    !**********************************************************************

    subroutine value_si(str,inum,ios)

        ! Converts number string to a single precision integer value

        character(len=*)::str
        integer(ki4) :: inum
        real(kr8) :: rnum

        call value_dr(str,rnum,ios)
        if(abs(rnum)>huge(inum)) then
            ios=15
            return
        end if
        inum=nint(rnum,ki4)

    end subroutine value_si

    !**********************************************************************

    subroutine shiftstr(str,n)

        ! Shifts characters in in the string 'str' n positions (positive values
        ! denote a right shift and negative values denote a left shift). Characters
        ! that are shifted off the end are lost. Positions opened up by the shift
        ! are replaced by spaces.

        character(len=*):: str

        lenstr=len(str)
        nabs=iabs(n)
        if(nabs>=lenstr) then
            str=repeat(' ',lenstr)
            return
        end if
        if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
        if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right
        return

    end subroutine shiftstr

    !**********************************************************************

    subroutine insertstr(str,strins,loc)

        ! Inserts the string 'strins' into the string 'str' at position 'loc'.
        ! Characters in 'str' starting at position 'loc' are shifted right to
        ! make room for the inserted string. Trailing spaces of 'strins' are
        ! removed prior to insertion

        character(len=*):: str,strins
        character(len=len(str))::tempstr

        lenstrins=len_trim(strins)
        tempstr=str(loc:)
        call shiftstr(tempstr,lenstrins)
        tempstr(1:lenstrins)=strins(1:lenstrins)
        str(loc:)=tempstr
        return

    end subroutine insertstr

    !**********************************************************************

    subroutine delsubstr(str,substr)

        ! Deletes first occurrence of substring 'substr' from string 'str' and
        ! shifts characters left to fill hole. Trailing spaces or blanks are
        ! not considered part of 'substr'.

        character(len=*):: str,substr

        lensubstr=len_trim(substr)
        ipos=index(str,substr)
        if(ipos==0) return
        if(ipos == 1) then
            str=str(lensubstr+1:)
        else
            str=str(:ipos-1)//str(ipos+lensubstr:)
        end if
        return

    end subroutine delsubstr

    !**********************************************************************

    subroutine delall(str,substr)

        ! Deletes all occurrences of substring 'substr' from string 'str' and
        ! shifts characters left to fill holes.

        character(len=*):: str,substr

        lensubstr=len_trim(substr)
        do
            ipos=index(str,substr)
            if(ipos == 0) exit
            if(ipos == 1) then
                str=str(lensubstr+1:)
            else
                str=str(:ipos-1)//str(ipos+lensubstr:)
            end if
        end do
        return

    end subroutine delall

    !**********************************************************************

    function uppercase(str) result(ucstr)

        ! convert string to upper case

        character (len=*):: str
        character (len=len_trim(str)):: ucstr

        ilen=len_trim(str)
        ioffset=iachar('A')-iachar('a')
        iquote=0
        ucstr=str
        do i=1,ilen
            iav=iachar(str(i:i))
            if(iquote==0 .and. (iav==34 .or.iav==39)) then
                iquote=1
                iqc=iav
                cycle
            end if
            if(iquote==1 .and. iav==iqc) then
                iquote=0
                cycle
            end if
            if (iquote==1) cycle
            if(iav >= iachar('a') .and. iav <= iachar('z')) then
                ucstr(i:i)=achar(iav+ioffset)
            else
                ucstr(i:i)=str(i:i)
            end if
        end do
        return

    end function uppercase

    !**********************************************************************

    function lowercase(str) result(lcstr)

        ! convert string to lower case

        character (len=*):: str
        character (len=len_trim(str)):: lcstr

        ilen=len_trim(str)
        ioffset=iachar('A')-iachar('a')
        iquote=0
        lcstr=str
        do i=1,ilen
            iav=iachar(str(i:i))
            if(iquote==0 .and. (iav==34 .or.iav==39)) then
                iquote=1
                iqc=iav
                cycle
            end if
            if(iquote==1 .and. iav==iqc) then
                iquote=0
                cycle
            end if
            if (iquote==1) cycle
            if(iav >= iachar('A') .and. iav <= iachar('Z')) then
                lcstr(i:i)=achar(iav-ioffset)
            else
                lcstr(i:i)=str(i:i)
            end if
        end do
        return

    end function lowercase

    !**********************************************************************

    subroutine readline(nunitr,line,ios)

        ! Reads line from unit=nunitr, ignoring blank lines
        ! and deleting comments beginning with an exclamation point(!)

        character (len=*):: line
        character (len=N_ch_max) :: sbuf
        integer :: icont,ios

        line = ''
        icont = 0
        do
            read(nunitr,'(a)', iostat=ios) sbuf      ! read input line
            if(ios /= 0) then
                if(icont /= 0) then
                    write(*,*) 'error: end of continuous line not found, end of file reached'
                endif
                return
            endif
            sbuf=adjustl(sbuf)
            ipos=index(sbuf,'!')
            if(ipos == 1) cycle
            if(ipos /= 0) then
                sbuf=sbuf(:ipos-1)
            endif
            icont=index(sbuf,'&')
            if(icont /= 0) then
                sbuf=sbuf(:icont-1)
            endif
            ipos = len_trim(line)+1
            call insertstr(line,sbuf,ipos)
            if(icont /= 0) cycle
            if(len_trim(line) /= 0) exit
        end do
        return

    end subroutine readline

    !**********************************************************************

    subroutine match(str,ipos,imatch)

        ! Sets imatch to the position in string of the delimiter matching the delimiter
        ! in position ipos. Allowable delimiters are (), [], {}, <>.

        character(len=*) :: str
        character :: delim1,delim2,ch

        lenstr=len_trim(str)
        delim1=str(ipos:ipos)
        select case(delim1)
            case('(')
                idelim2=iachar(delim1)+1
                istart=ipos+1
                iend=lenstr
                inc=1
            case(')')
                idelim2=iachar(delim1)-1
                istart=ipos-1
                iend=1
                inc=-1
            case('[','{','<')
                idelim2=iachar(delim1)+2
                istart=ipos+1
                iend=lenstr
                inc=1
            case(']','}','>')
                idelim2=iachar(delim1)-2
                istart=ipos-1
                iend=1
                inc=-1
            case default
                write(*,*) delim1,' is not a valid delimiter'
                return
        end select
        if(istart < 1 .or. istart > lenstr) then
            write(*,*) delim1,' has no matching delimiter'
            return
        end if
        delim2=achar(idelim2) ! matching delimiter

        isum=1
        do i=istart,iend,inc
            ch=str(i:i)
            if(ch /= delim1 .and. ch /= delim2) cycle
            if(ch == delim1) isum=isum+1
            if(ch == delim2) isum=isum-1
            if(isum == 0) exit
        end do
        if(isum /= 0) then
            write(*,*) delim1,' has no matching delimiter'
            return
        end if
        imatch=i

        return

    end subroutine match

    !**********************************************************************

    subroutine write_dr(rnum,str,fmt)

        ! Writes double precision real number rnum to string str using format fmt

        real(kr8) :: rnum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) rnum
        str=adjustl(str)

    end subroutine write_dr

    !***********************************************************************

    subroutine write_sr(rnum,str,fmt)

        ! Writes single precision real number rnum to string str using format fmt

        real(kr4) :: rnum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) rnum
        str=adjustl(str)

    end subroutine write_sr

    !***********************************************************************

    subroutine write_di(inum,str,fmt)

        ! Writes double precision integer inum to string str using format fmt

        integer(ki8) :: inum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) inum
        str=adjustl(str)

    end subroutine write_di

    !***********************************************************************

    subroutine write_si(inum,str,fmt)

        ! Writes single precision integer inum to string str using format fmt

        integer(ki4) :: inum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) inum
        str=adjustl(str)

    end subroutine write_si

    !***********************************************************************

    subroutine trimzero(str)

        ! Deletes nonsignificant trailing zeroes from number string str. If number
        ! string ends in a decimal point, one trailing zero is added.

        character(len=*) :: str
        character :: ch
        character(len=10) :: exp

        ipos=scan(str,'eE')
        if(ipos>0) then
            exp=str(ipos:)
            str=str(1:ipos-1)
        endif
        lstr=len_trim(str)
        do i=lstr,1,-1
            ch=str(i:i)
            if(ch=='0') cycle
            if(ch=='.') then
                str=str(1:i)//'0'
                if(ipos>0) str=trim(str)//trim(exp)
                exit
            endif
            str=str(1:i)
            exit
        end do
        if(ipos>0) str=trim(str)//trim(exp)

    end subroutine trimzero

    !**********************************************************************

    subroutine writeq_dr(unit,namestr,value,fmt)

        ! Writes a string of the form <name> = value to unit

        real(kr8) :: value
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr

        call writenum(value,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_dr

    !**********************************************************************

    subroutine writeq_sr(unit,namestr,value,fmt)

        ! Writes a string of the form <name> = value to unit

        real(kr4) :: value
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr

        call writenum(value,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_sr

    !**********************************************************************

    subroutine writeq_di(unit,namestr,ivalue,fmt)

        ! Writes a string of the form <name> = ivalue to unit

        integer(ki8) :: ivalue
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr
        call writenum(ivalue,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_di

    !**********************************************************************

    subroutine writeq_si(unit,namestr,ivalue,fmt)

        ! Writes a string of the form <name> = ivalue to unit

        integer(ki4) :: ivalue
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr
        call writenum(ivalue,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_si

    !**********************************************************************

    function is_letter(ch) result(res)

        ! Returns .true. if ch is a letter and .false. otherwise

        character :: ch
        logical :: res

        select case(ch)
            case('A':'Z','a':'z')
                res=.true.
            case default
                res=.false.
        end select
        return

    end function is_letter

    !**********************************************************************

    function is_digit(ch) result(res)

        ! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

        character :: ch
        logical :: res

        select case(ch)
            case('0':'9')
                res=.true.
            case default
                res=.false.
        end select
        return

    end function is_digit

    !**********************************************************************

    subroutine split(str,delims,before,sep)

        ! Routine finds the first instance of a character from 'delims' in the
        ! the string 'str'. The characters before the found delimiter are
        ! output in 'before'. The characters after the found delimiter are
        ! output in 'str'. The optional output character 'sep' contains the
        ! found delimiter. A delimiter in 'str' is treated like an ordinary
        ! character if it is preceded by a backslash (\). If the backslash
        ! character is desired in 'str', then precede it with another backslash.

        character(len=*) :: str,delims,before
        character,optional :: sep
        logical :: pres
        character :: ch,cha

        pres=present(sep)
        str=adjustl(str)
        call compact(str)
        lenstr=len_trim(str)
        if(lenstr == 0) return        ! string str is empty
        k=0
        ibsl=0                        ! backslash initially inactive
        before=' '
        do i=1,lenstr
            ch=str(i:i)
            if(ibsl == 1) then          ! backslash active
                k=k+1
                before(k:k)=ch
                ibsl=0
                cycle
            end if
            if(ch == '\') then          ! backslash with backslash inactive
                k=k+1
                before(k:k)=ch
                ibsl=1
                cycle
            end if
            ipos=index(delims,ch)
            if(ipos == 0) then          ! character is not a delimiter
                k=k+1
                before(k:k)=ch
                cycle
            end if
            if(ch /= ' ') then          ! character is a delimiter that is not a space
                str=str(i+1:)
                if(pres) sep=ch
                exit
            end if
            cha=str(i+1:i+1)            ! character is a space delimiter
            iposa=index(delims,cha)
            if(iposa > 0) then          ! next character is a delimiter
                str=str(i+2:)
                if(pres) sep=cha
                exit
            else
                str=str(i+1:)
                if(pres) sep=ch
                exit
            end if
        end do
        if(i >= lenstr) str=''
        str=adjustl(str)              ! remove initial spaces
        return

    end subroutine split

    !**********************************************************************

    subroutine removebksl(str)

        ! Removes backslash (\) characters. Double backslashes (\\) are replaced
        ! by a single backslash.

        character(len=*):: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr

        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=' '
        k=0
        ibsl=0                        ! backslash initially inactive

        do i=1,lenstr
            ch=str(i:i)
            if(ibsl == 1) then          ! backslash active
                k=k+1
                outstr(k:k)=ch
                ibsl=0
                cycle
            end if
            if(ch == '\') then          ! backslash with backslash inactive
                ibsl=1
                cycle
            end if
            k=k+1
            outstr(k:k)=ch              ! non-backslash with backslash inactive
        end do

        str=adjustl(outstr)

    end subroutine removebksl

    !**********************************************************************

    subroutine txt_read_struct_data_string(cfgUnit,strtit,vartit,varstr,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: vartit
        ! output:
        character(N_ch_max)   :: varstr
        integer        :: ios
        ! local:
        character(128) :: acttit
        character(128) :: actstr
        character(128) :: sbuf
        integer        :: ipos,imatch
        character(N_ch_max) :: str
        character(N_str_ch),dimension(N_str_arg) :: args
        integer :: i,narg

        character(32), dimension(20) :: StrLevel       ! strcture levels titles
        integer                      :: Nlevs          ! number of strcture levels

        ! 1. parcing structure title -> get structure levels
        Nlevs = 0
        sbuf = trim(strtit)
        sbuf = adjustl(sbuf)
        i = index(sbuf,"%")
        do while (i.gt.0)
            Nlevs = Nlevs + 1
            StrLevel(Nlevs) = sbuf(1:i-1)
            sbuf = sbuf(i+1:)
            sbuf = adjustl(sbuf)
            i = index(sbuf,"%")
        end do
        Nlevs = Nlevs + 1
        StrLevel(Nlevs) = trim(sbuf)


        ! 2. get base level (from file)
        sbuf = trim(StrLevel(1))
        call txt_findvar(cfgUnit,sbuf,str)
        ipos = index(str,"{")
        if(ipos.eq.0) then
            ios = 1
            write(*,*) "struct (",trim(sbuf),") read error: left bracket { not found"
            stop
        end if

        call match(str,ipos,imatch)
        if(imatch.eq.0) then
            ios = 1
            write(*,*) "struct (",trim(sbuf),") read error: right bracket } not found"
            stop
        end if
        if(ipos+1.ge.imatch-1) then
            ios = 1
            write(*,*) "struct (",trim(sbuf),") read error: empty structure"
            stop
        end if

        str = str(ipos+1:imatch-1)

        ! 3. go to relevant structure level (parcing string):
        do i = 2,Nlevs
            ipos = index(str,trim(StrLevel(i)))
            str = str(ipos:)
            ipos = index(str,"{")
            call match(str,ipos,imatch)
            if(ipos+1.ge.imatch-1) then
                write(*,*) " txtread: ","error in reading strcture:",trim(str)
                stop
            end if
            str = str(ipos+1:imatch-1)
        end do

        ! 4. replace (from relevant level) other internal structures data with string "next level"
        ipos = index(str,"{")
        do while (ipos.gt.0)
            call match(str,ipos,imatch)
            if(ipos+1.ge.imatch-1) then
                write(*,*) " txtread: ","error in reading strcture:",trim(str)
                stop
            end if
            str = str(:ipos-1)//"'next_level'"//str(imatch+1:)
            sbuf = adjustl(sbuf)
            ipos = index(sbuf,"{")
        end do

        ! parce relevant structure level string (get element)
        call parse(str,";",args,narg)
        varstr = ""
        do i = 1,narg
            call removesp(args(i))
            sbuf = trim(args(i))
            call split(sbuf,"=",acttit)
            actstr = trim(sbuf)
            if(trim(vartit).eq.trim(acttit)) then
                varstr = actstr
            end if
        enddo
        if(len(trim(varstr)).eq.0) then
            write(*,*) " config error: '",trim(vartit),"' not found"
            stop
        end if

        return
    end subroutine txt_read_struct_data_string

    !**********************************************************************

    subroutine txt_read_struct_val_str(cfgUnit,strtit,vartit,varval,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: vartit
        ! output:
        character(*)   :: varval
        integer        :: ios
        ! local:
        character(N_ch_max) :: varstr

        call txt_read_struct_data_string(cfgUnit,strtit,vartit,varstr,ios)
        call txt_value(varstr,varval,ios)

        return
    end subroutine txt_read_struct_val_str

    !**********************************************************************

    subroutine txt_read_struct_val_dr(cfgUnit,strtit,vartit,varval,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: vartit
        ! output:
        real(kr8)      :: varval
        integer        :: ios
        ! local:
        character(N_ch_max) :: varstr

        call txt_read_struct_data_string(cfgUnit,strtit,vartit,varstr,ios)
        call txt_value(varstr,varval,ios)

        return
    end subroutine txt_read_struct_val_dr

    !**********************************************************************

    subroutine txt_read_struct_val_sr(cfgUnit,strtit,vartit,varval,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: vartit
        ! output:
        real(kr4)      :: varval
        integer        :: ios
        ! local:
        character(N_ch_max) :: varstr

        call txt_read_struct_data_string(cfgUnit,strtit,vartit,varstr,ios)
        call txt_value(varstr,varval,ios)

        return
    end subroutine txt_read_struct_val_sr

    !**********************************************************************

    subroutine txt_read_struct_val_di(cfgUnit,strtit,vartit,varval,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: vartit
        ! output:
        integer(ki8)   :: varval
        integer        :: ios
        ! local:
        character(N_ch_max) :: varstr

        call txt_read_struct_data_string(cfgUnit,strtit,vartit,varstr,ios)
        call txt_value(varstr,varval,ios)

        return
    end subroutine txt_read_struct_val_di

    !**********************************************************************

    subroutine txt_read_struct_val_si(cfgUnit,strtit,vartit,varval,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: vartit
        ! output:
        integer(ki4)   :: varval
        integer        :: ios
        ! local:
        character(N_ch_max) :: varstr

        call txt_read_struct_data_string(cfgUnit,strtit,vartit,varstr,ios)
        call txt_value(varstr,varval,ios)

        return
    end subroutine txt_read_struct_val_si

    !**********************************************************************

    subroutine txt_read_struct_arr_str(cfgUnit,strtit,arrtit,array,num,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: arrtit
        ! output:
        character(*), dimension(:) :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        ! local:
        character(N_ch_max)                         :: arrstr
        character(N_ch_max), dimension(size(array)) :: strdim
        integer             :: i,ipos,imatch

        call txt_read_struct_data_string(cfgUnit,strtit,arrtit,arrstr,ios)
        ipos = 1
        imatch = 0
        call match(arrstr,ipos,imatch)
        arrstr = arrstr(ipos+1:imatch-1)

        call parse(arrstr,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_struct_arr_str

    !**********************************************************************

    subroutine txt_read_struct_arr_dr(cfgUnit,strtit,arrtit,array,num,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: arrtit
        ! output:
        real(kr8), dimension(:)    :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        ! local:
        character(N_ch_max)                         :: arrstr
        character(N_ch_max), dimension(size(array)) :: strdim
        integer             :: i,ipos,imatch

        call txt_read_struct_data_string(cfgUnit,strtit,arrtit,arrstr,ios)
        ipos = 1
        imatch = 0
        call match(arrstr,ipos,imatch)
        arrstr = arrstr(ipos+1:imatch-1)

        call parse(arrstr,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_struct_arr_dr

    !**********************************************************************

    subroutine txt_read_struct_arr_sr(cfgUnit,strtit,arrtit,array,num,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: arrtit
        ! output:
        real(kr4), dimension(:)    :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        ! local:
        character(N_ch_max)                         :: arrstr
        character(N_ch_max), dimension(size(array)) :: strdim
        integer             :: i,ipos,imatch

        call txt_read_struct_data_string(cfgUnit,strtit,arrtit,arrstr,ios)
        ipos = 1
        imatch = 0
        call match(arrstr,ipos,imatch)
        arrstr = arrstr(ipos+1:imatch-1)

        call parse(arrstr,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_struct_arr_sr

    !**********************************************************************

    subroutine txt_read_struct_arr_di(cfgUnit,strtit,arrtit,array,num,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: arrtit
        ! output:
        integer(ki8), dimension(:) :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        ! local:
        character(N_ch_max)                         :: arrstr
        character(N_ch_max), dimension(size(array)) :: strdim
        integer             :: i,ipos,imatch

        call txt_read_struct_data_string(cfgUnit,strtit,arrtit,arrstr,ios)
        ipos = 1
        imatch = 0
        call match(arrstr,ipos,imatch)
        arrstr = arrstr(ipos+1:imatch-1)

        call parse(arrstr,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_struct_arr_di

    !**********************************************************************

    subroutine txt_read_struct_arr_si(cfgUnit,strtit,arrtit,array,num,ios)
        implicit none
        ! input:
        integer :: cfgUnit
        character(*)   :: strtit
        character(*)   :: arrtit
        ! output:
        integer(ki4), dimension(:) :: array      ! array
        integer                    :: num        ! array size
        integer                    :: ios
        ! local:
        character(N_ch_max)                         :: arrstr
        character(N_ch_max), dimension(size(array)) :: strdim
        integer             :: i,ipos,imatch

        call txt_read_struct_data_string(cfgUnit,strtit,arrtit,arrstr,ios)
        ipos = 1
        imatch = 0
        call match(arrstr,ipos,imatch)
        arrstr = arrstr(ipos+1:imatch-1)

        call parse(arrstr,',;',strdim,num)
        do i=1,num
            call txt_value(strdim(i),array(i),ios)
        enddo

        return
    end subroutine txt_read_struct_arr_si


end module module_TxtRead


