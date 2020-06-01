module module_utils

    interface valpos
        module procedure valpos_int1
        module procedure valpos_int2
        module procedure valpos_int4
        module procedure valpos_int8
        module procedure valpos_str
    end interface

    private valpos_int1,valpos_int2,valpos_int4,valpos_int8,valpos_str

contains

    function GetFreeUnit() result(Free_Unit)
        !--------------------------------------------------------------------------------------------------------------------------------
        !!The GetUnit function is used for getting a free logic unit. The users of \LIBVTKIO does not know which is
        !!the logical unit: \LIBVTKIO handels this information without boring the users. The logical unit used is safe-free: if the
        !!program calling \LIBVTKIO has others logical units used \LIBVTKIO will never use these units, but will choice one that is free.
        !--------------------------------------------------------------------------------------------------------------------------------

        implicit none

        !--------------------------------------------------------------------------------------------------------------------------------
        integer, parameter:: I4P = selected_int_kind(9)        ! range $[-2^{31} ,+2^{31}  -1]$
        integer(I4P):: Free_Unit ! free logic unit
        integer(I4P):: n1        ! counter
        integer(I4P):: ios       ! inquiring flag
        logical(4)::   lopen     ! inquiring flag
        !--------------------------------------------------------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------------------------------------------------------
        !!The following is the code snippet of GetUnit function: the units 0, 5, 6, 9 and all non-free units are discarded.
        !!
        !(\doc)codesnippet
        Free_Unit = -1_I4P                                      ! initializing free logic unit
        n1=1_I4P                                                ! initializing counter
        do
            if ((n1/=5_I4P).AND.(n1/=6_I4P).AND.(n1/=9_I4P)) then
                inquire (unit=n1,opened=lopen,iostat=ios)           ! verify logic units
                if (ios==0_I4P) then
                    if (.NOT.lopen) then
                        Free_Unit = n1                                  ! assignment of free logic
                        return
                    endif
                endif
            endif
            n1=n1+1_I4P                                           ! updating counter
        enddo
        return
        !(doc/)codesnippet
        !!GetUnit function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope to public.
        !--------------------------------------------------------------------------------------------------------------------------------
    endfunction GetFreeUnit

!    ! function for upper case
!    function Upper_Case(string)
!        !--------------------------------------------------------------------------------------------------------------------------------
!        !!The Upper\_Case function converts the lower case characters of a string to upper case one. \LIBVTKIO uses this function in
!        !!order to achieve case-insensitive: all character variables used within \LIBVTKIO functions are pre-processed by
!        !!Uppper\_Case function before these variables are used. So the users can call \LIBVTKIO functions whitout pay attention of the
!        !!case of the kwywords passed to the functions: calling the function VTK\_INI with the string \code{E_IO = VTK_INI('Ascii',...)}
!        !!or with the string  \code{E_IO = VTK_INI('AscII',...)} is equivalent.
!        !--------------------------------------------------------------------------------------------------------------------------------
!
!        implicit none
!
!        !--------------------------------------------------------------------------------------------------------------------------------
!        character(len=*), intent(IN):: string     ! string to be converted
!        character(len=len(string))::   Upper_Case ! converted string
!        integer::                      n1         ! characters counter
!        integer::                      idx        ! characters index in list
!        !--------------------------------------------------------------------------------------------------------------------------------
!
!        !--------------------------------------------------------------------------------------------------------------------------------
!        !!The following is the code snippet of Upper\_Case function.
!        !!
!        !(\doc)codesnippet
!        Upper_Case = string
!        do n1=1,len(string)
!            idx = ichar(string(n1:n1))
!            select case(idx)
!                case(97:122)
!                    Upper_Case(n1:n1)=char(idx-32) ! Upper case conversion
!            endselect
!        enddo
!        return
!        !(doc/)codesnippet
!        !!Upper\_Case function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope
!        !!to public.
!        !--------------------------------------------------------------------------------------------------------------------------------
!    endfunction Upper_Case
!
!    function Lower_Case(string)
!        !--------------------------------------------------------------------------------------------------------------------------------
!        !!The Lower\_Case function converts the lower case characters of a string to upper case one. \LIBVTKIO uses this function in
!        !!order to achieve case-insensitive: all character variables used within \LIBVTKIO functions are pre-processed by
!        !!Uppper\_Case function before these variables are used. So the users can call \LIBVTKIO functions whitout pay attention of the
!        !!case of the kwywords passed to the functions: calling the function VTK\_INI with the string \code{E_IO = VTK_INI('Ascii',...)}
!        !!or with the string  \code{E_IO = VTK_INI('AscII',...)} is equivalent.
!        !--------------------------------------------------------------------------------------------------------------------------------
!
!        implicit none
!
!        !--------------------------------------------------------------------------------------------------------------------------------
!        character(len=*), intent(IN):: string     ! string to be converted
!        character(len=len(string))::   Lower_Case ! converted string
!        integer::                      n1         ! characters counter
!        integer::                      idx        ! characters index in list
!        !--------------------------------------------------------------------------------------------------------------------------------
!
!        !--------------------------------------------------------------------------------------------------------------------------------
!        !!The following is the code snippet of Lower\_Case function.
!        !!
!        !(\doc)codesnippet
!        Lower_Case = string
!        do n1=1,len(string)
!            idx = ichar(string(n1:n1))
!            select case(idx)
!                case(65:90)
!                    Lower_Case(n1:n1)=char(idx+32) ! Lower case conversion
!            endselect
!        enddo
!        return
!        !(doc/)codesnippet
!        !!Lower\_Case function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope
!        !!to public.
!        !--------------------------------------------------------------------------------------------------------------------------------
!    endfunction Lower_Case

! ------------------------------------------------------------------
    function Upper_Case (str) result (string)
!   Changes a string to upper case
    implicit none
    ! input
    character(*), intent(in) :: str
    ! output
    character(len(str))      :: string
    ! local
    integer :: ic, i
    character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, len_trim(str)
        ic = index(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

    end function Upper_Case

! ------------------------------------------------------------------
    function Lower_Case (str) result (string)
!   Changes a string to lower case
    implicit none
    ! input
    character(*), intent(in) :: str
    ! output
    character(len(str))      :: string
    ! local
    integer :: ic, i
    character(26), parameter :: cap = 'abcdefghijklmnopqrstuvwxyz'
    character(26), parameter :: low = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, len_trim(str)
        ic = index(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

    end function Lower_Case


    ! --------- Gauss solver for system of linear equations -----------------
    subroutine gaussj(a,n,b,m,ios)
        implicit none
        ! parameters
        integer, parameter :: NMAX=100
        ! input
        real(8) :: a(n,n)    ! input matrix [np x np]
        integer :: n         ! number of rows-columns
        real(8) :: b(n,m)    ! right hand side vector (RHS) -> solution vector at output
        integer :: m         ! number of RHS vectors (solutions), usual case: m=1
        ! addon
        integer :: np        ! np = n
        integer :: mp        ! mp = np
        ! local
        integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
        real(8) big,dum,pivinv
        !character(1) buf
        integer :: ios

        icol = 1
        irow = 1
        np = n
        mp = m
        ios = 0

        do 11 j=1,n
          ipiv(j)=0
    11  continue
        do 22 i=1,n
          big=0.0d0
          do 13 j=1,n
            if(ipiv(j).ne.1)then
              do 12 k=1,n
                if (ipiv(k).eq.0) then
                  if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                  endif
                else if (ipiv(k).gt.1) then
                    !pause 'singular matrix in gaussj'
                    write(*,*) 'singular matrix in gaussj'
                    ios = 1
!                    write(*,*) 'continue? (y/n)'
!                    read(*,*) buf
!                    if(buf.eq.'n') stop
                endif
    12        continue
            endif
    13    continue
          ipiv(icol)=ipiv(icol)+1
          if (irow.ne.icol) then
           do 14 l=1,n
             dum=a(irow,l)
              a(irow,l)=a(icol,l)
              a(icol,l)=dum
    14      continue
            do 15 l=1,m
              dum=b(irow,l)
              b(irow,l)=b(icol,l)
              b(icol,l)=dum
    15      continue
          endif
          indxr(i)=irow
          indxc(i)=icol
          if (a(icol,icol).eq.0.0d0) then
              write(*,*) 'singular matrix in gaussj'
              ios = 1
!              write(*,*) 'continue? (y/n)'
!              read(*,*) buf
!              if(buf.eq.'n') stop
          endif

          pivinv=1.0d0/a(icol,icol)
          a(icol,icol)=1.0d0
          do 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
    16    continue
          do 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
    17    continue
          do 21 ll=1,n
            if(ll.ne.icol)then
              dum=a(ll,icol)
              a(ll,icol)=0.0d0
              do 18 l=1,n
                a(ll,l)=a(ll,l)-a(icol,l)*dum
    18        continue
              do 19 l=1,m
                b(ll,l)=b(ll,l)-b(icol,l)*dum
    19        continue
            endif
    21    continue
    22  continue
        do 24 l=n,1,-1
          if(indxr(l).ne.indxc(l))then
            do 23 k=1,n
              dum=a(k,indxr(l))
              a(k,indxr(l))=a(k,indxc(l))
              a(k,indxc(l))=dum
    23      continue
          endif
    24  continue
        return
    end subroutine gaussj

    ! - - - - - - - - - - - - - - - - - -
    integer function valpos_int1(array,ix)
        ! get position of integer ix in integer array
        implicit none
        ! input
        integer(1), intent(in) :: array(:)
        integer                :: ix
!        integer(1)             :: ix
        ! local
        integer :: arsize
        integer :: i

        arsize = size(array)
        i = 1
        valpos_int1 = 0
        do while(array(i).ne.ix)
            i = i + 1
            if(i.gt.arsize) then
                return
            endif
        enddo
        valpos_int1 = i
        return
    end function valpos_int1

    ! - - - - - - - - - - - - - - - - - -
    integer function valpos_int2(array,ix)
        ! get position of integer ix in integer array
        implicit none
        ! input
        integer(2), intent(in) :: array(:)
        integer                :: ix
!        integer(2)             :: ix
        ! local
        integer :: arsize
        integer :: i

        arsize = size(array)
        i = 1
        valpos_int2 = 0
        do while(array(i).ne.ix)
            i = i + 1
            if(i.gt.arsize) then
                return
            endif
        enddo
        valpos_int2 = i
        return
    end function valpos_int2

    ! - - - - - - - - - - - - - - - - - -
    integer function valpos_int4(array,ix)
        ! get position of integer ix in integer array
        implicit none
        ! input
        integer(4), intent(in) :: array(:)
        integer                :: ix
!        integer(4)             :: ix
        ! local
        integer :: arsize
        integer :: i

        arsize = size(array)
        i = 1
        valpos_int4 = 0
        do while(array(i).ne.ix)
            i = i + 1
            if(i.gt.arsize) then
                return
            endif
        enddo
        valpos_int4 = i
        return
    end function valpos_int4

    ! - - - - - - - - - - - - - - - - - -
    integer function valpos_int8(array,ix)
        ! get position of integer ix in integer array
        implicit none
        ! input
        integer(8), intent(in) :: array(:)
        integer                :: ix
!        integer(8)             :: ix
        ! local
        integer :: arsize
        integer :: i

        arsize = size(array)
        i = 1
        valpos_int8 = 0
        do while(array(i).ne.ix)
            i = i + 1
            if(i.gt.arsize) then
                return
            endif
        enddo
        valpos_int8 = i
        return
    end function valpos_int8

    ! - - - - - - - - - - - - - - - - - -
    integer function valpos_str(array,sbuf)
        ! get position of integer ix in integer array
        implicit none
        character(*), intent(in) :: array(:)
        character(*) :: sbuf
        integer :: arsize
        integer :: i

        arsize = size(array)
        i = 1
        valpos_str = 0
        do while(trim(array(i)).ne.trim(sbuf))
            i = i + 1
            if(i.gt.arsize) then
                return
            endif
        enddo
        valpos_str = i
        return
    end function valpos_str

end module
