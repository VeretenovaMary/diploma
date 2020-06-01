module module_AtomData

    ! -------------------------------------------
    type TElementData
        character(4)  :: title        ! atom title (symbol)
        real(8)       :: A            ! g/mole, atomic/molar weight
        integer       :: Z            ! charge number
        character(3)  :: QSplit       ! "n","nl","nlj"
        integer       :: Nlev         ! max number of levels (main quantum number) for electrons filling
        integer       :: Zmin         ! minimal ionization degree (for current data base)
        integer       :: Zmax         ! maximal ionization degree (for current data base)
        integer       :: NIons        ! number of ions ( = Zmax-Zmin+1)
        integer       :: NQd          ! total number of ion configurations
        integer       :: NQt_bb        ! total number of ion bound-bound transitions (excitations)
        integer       :: NQt_bf        ! total number of ion bound-bound transitions (excitations)
        integer       :: ncells       ! number of electron cells ((Nlev+1)*Nlev/2)
        type(TIonData),pointer :: ion(:) => NULL()   ! data for ions (including neutrals)
    end type TElementData

    !     Q (or q) - symbol for ion configuration definition
    !     n   - main quantum number
    !     l   - orbital quantum number
    !     Q   - whole ion (average)
    !     Qn  - Q by "n" (main quantum number) splitting
    !     Qnl - Qn by "l" (orbital quantum number) splitting
    ! -------------------------------------------
    type TIonData                 ! only main quantum number splitting
        integer       :: ZQ       ! current ion charge
        real(8)       :: EQ       ! eV, Q-ion energy (is taken as {EQnl} lower energy value)
        real(8)       :: dEion    ! eV, ionization energy
        integer       :: Nqd                            ! = Nqd, num. ion configurations Qn (splitting by main q.num.)
        integer       :: NQt_bb                          ! = Nqt, num. transitions           (splitting by main q.num.)
        integer       :: NQt_bf                          ! = Nqt, num. transitions           (splitting by main q.num.)
        type(TIonQdData),pointer :: Qd(:) => NULL()      ! ion configurations data
    end type TIonData

    ! -------------------------------------------
    type TIonQdData                 ! only main quantum number splitting (data from '<sub>qd<k>.ch0' files)
        real(8)                   :: EQd              ! configuration average energy
        real(8)                   :: gQd              ! configuration statistical weight
        integer(1),allocatable    :: cells(:)         ! electrons setling in orbitals
        integer                   :: iList
    end type TIonQdData

    ! -------------------------------------------
    type TAtomQdList                 ! only main quantum number splitting (data from '<sub>qd<k>.ch0' files)
        integer                   :: iE,iI,iQ
    end type TAtomQdList

    ! -------------------------------------------
    type TAtomQtbbList                 ! only main quantum number splitting (data from '<sub>qd<k>.ch0' files)
        integer :: iLQd1
        integer :: iLQd2
        real(8) :: dE12
        real(8) :: gf12
        real(8) :: gQd1
        real(8) :: gQd2
    end type TAtomQtbbList

    ! -------------------------------------------
    type TAtomQtbfList                 ! only main quantum number splitting (data from '<sub>qd<k>.ch0' files)
        integer :: iLQd1
        integer :: iLQd2
        real(8) :: dE12
        integer :: ilev                ! main quantum number
        real(8) :: gQd1
        real(8) :: gQd2
    end type TAtomQtbfList


    ! -------------------------------------------
    type TAtomData                 ! only main quantum number splitting (data from '<sub>qd<k>.ch0' files)
        integer                    :: NElems              ! number of atoms (chem. elements)
        type(TElementData),pointer :: elem(:) => NULL()   ! NULL option doesn't work (?)
        integer                    :: NQds
        integer                    :: NQt_bb
        integer                    :: NQt_bf
        type(TAtomQdList),pointer    :: QdList(:)
        type(TAtomQtbbList),pointer  :: QtbbList(:)
        type(TAtomQtbfList),pointer  :: QtbfList(:)
    end type TAtomData

    ! -------------------------------------------
    type TConcData
        !sequence
        character(4)         :: title
        integer              :: NIons
        integer              :: MajIonIdx
        real(8)              :: Ne
        real(8)              :: dCi
        real(8), allocatable :: EQ(:),gQ(:)
        real(8), allocatable :: Ni(:)
        real(8), allocatable :: Ci(:)
    end type TConcData

    ! variables
    type(TAtomData) :: AtomData
    integer, parameter :: AtomNumMax = 10
    character(128), dimension(AtomNumMax) :: AtomicDBPaths

    !      private SetIonSize,SetQnSize,SetCellsSize
    private SetCellsSize

    interface AllocAtomData
        module procedure SetIonSize
        module procedure SetQdSize
    end interface

    ! *****************************************************************
contains
    ! *****************************************************************
    subroutine SetIonSize(Elem,N)
        implicit none
        integer N
        type(TElementData) :: Elem
        allocate(Elem%ion(N))
        return
    end subroutine SetIonSize

    ! **********************************************************
    subroutine SetQdSize(Ion,n1,n2)
        implicit none
        type(TIonData) :: Ion
        integer        :: n1,n2
        integer :: i

        allocate(Ion%Qd(n1))
        do i = 1,n1
            call SetCellsSize(Ion%Qd(i),n2)
        enddo
        return
    end subroutine SetQdSize

    subroutine SetCellsSize(IonQd,n)
        implicit none
        type(TIonQdData) :: IonQd
        integer          :: n
        allocate(IonQd%cells(n))
        return
    end subroutine SetCellsSize

    ! ----------------------------------------------------------------
    subroutine AtomDataRead(iE,AtomName,dbPath)
        use module_Settings
        use module_TxtRead
        implicit none
        ! input
        character(2)                     :: AtomName
        character(*)                     :: dbPath
        ! local
        character(3)              :: AtomQSplit
        integer                   :: iE,j,k,num,ios
        character(256)            :: fpath
        character(512)            :: sbuf
        character(2)              :: chtitle
        real(8)                   :: EQd,gQd
        integer(1),dimension(210) :: cells
        integer                   :: ncls,Nqd,Nqt
        integer                   :: atmUnit
        real(8)                   :: eps

        !write(*,*) ' - - - - - - - - - - - '
        !write(*,*) AtomName
        !write(*,*)
        !write(*,*) trim(AtomName)//": "//"reading inf: "
        atmUnit = 100

        ! 1. reading base atom data

        ! . . . . . read element/DB header data: Nlev, Z, A
        if(len(trim(AtomName)).eq.1) then
            chtitle = trim(AtomName)//'_'
        else
            chtitle = trim(AtomName)
        endif
        fpath = trim(dbPath)//"/"//chtitle//"_inf.ch0"
        AtomData%elem(iE)%title   = AtomName
        call slash_correct(fpath,OSType)
        open(atmUnit,file=fpath)
        read(atmUnit,*) AtomData%elem(iE)%Nlev,eps,AtomData%elem(iE)%A
        AtomData%elem(iE)%Z = int(eps)
        ! ..................................................


        ! . . . define data base "Quantum Spliting": "nl" or "n"? . . .
        !             -> find number of cells in *qd.ch0 file:
        fpath = trim(dbPath)//"/"//chtitle//"qd"
        write(fpath,'(a,i2.2,a)') trim(fpath),0,".ch0"
        call slash_correct(fpath,OSType)
        open(atmUnit+1,file=fpath)
            read(atmUnit+1,"(a)",iostat=ios) sbuf
        close(atmUnit+1)
        call compact(sbuf)
        sbuf = adjustl(sbuf)
        do j = 1,4
            k = index(sbuf," ")
            sbuf(1:) = sbuf(k+1:)
        end do
        ncls = 0
        do while(len_trim(sbuf).gt.0)
            ncls = ncls + 1
            k = index(sbuf," ")
            sbuf(1:) = sbuf(k+1:)
        end do
        !write(*,*) "NumCells = ",ncls
        ! ..............................................


        k = (1+AtomData%elem(iE)%Nlev)*AtomData%elem(iE)%Nlev/2      ! number of electron "cells" ('nl' splitting - defolt)
        j = AtomData%elem(iE)%Nlev                                         ! case 'n'  splitting (num. of electron "cells")
        if(ncls.eq.k) then
            AtomQSplit = "nl"
        elseif(ncls.eq.j) then
            AtomQSplit = "n"
        else
            AtomQSplit = "err"
            write(*,*)   " warning: Quantum Spliting is undefined for element: ",AtomName
        endif
        AtomData%elem(iE)%QSplit = AtomQSplit

        ! ------------
        AtomData%elem(iE)%ncells = ncls
        read(atmUnit,*) AtomData%elem(iE)%Zmin,AtomData%elem(iE)%Zmax
        AtomData%elem(iE)%NIons = AtomData%elem(iE)%Zmax - AtomData%elem(iE)%Zmin + 1
        !call AllocAtomData(AtomData%elem(i),AtomData%elem(i)%NIons)
        num = AtomData%elem(iE)%NIons
        call SetIonSize(AtomData%elem(iE),num)

        ! 2. reading ions data: splitted by main quantum number (Qn data)
        AtomData%elem(iE)%NQd = 0
        AtomData%elem(iE)%NQt_bb = 0
        do k=1,AtomData%elem(iE)%NIons
            read(atmUnit,*,iostat=ios) AtomData%elem(iE)%ion(k)%ZQ,  &
                                       AtomData%elem(iE)%ion(k)%Nqd, &
                                       AtomData%elem(iE)%ion(k)%NQt_bb
            !call AllocAtomData(AtomData%elem(i)%ion(k),AtomData%elem(i)%ion(k)%Nqd)
            Nqd  = AtomData%elem(iE)%ion(k)%Nqd
            Nqt  = AtomData%elem(iE)%ion(k)%NQt_bb
            ncls = AtomData%elem(iE)%ncells
            AtomData%elem(iE)%NQd = AtomData%elem(iE)%NQd + Nqd
            AtomData%elem(iE)%NQt_bb = AtomData%elem(iE)%NQt_bb + Nqt
            call SetQdSize(AtomData%elem(iE)%ion(k),nqd,ncls)
        enddo
        close(atmUnit)

        ! 3. reading
        !write(*,*) trim(AtomName)//": "//"reading Qd: "
        do k=1,AtomData%elem(iE)%NIons
            fpath = trim(dbPath)//"/"//chtitle//"qd"
            write(fpath,'(a,i2.2,a)') trim(fpath),k-1,".ch0"
            call slash_correct(fpath,OSType)
            !  . . . . . set Qn data . . . . .
            open(atmUnit,file=fpath)
            do j=1,AtomData%elem(iE)%ion(k)%Nqd
                cells = 0
                read(atmUnit,*,iostat=ios) num,num,EQd,gQd,cells(1:ncls)
                AtomData%elem(iE)%ion(k)%Qd(j)%EQd   = EQd
                AtomData%elem(iE)%ion(k)%Qd(j)%gQd   = gQd
                AtomData%elem(iE)%ion(k)%Qd(j)%cells(1:ncls) = cells(1:ncls)
            enddo
            close(atmUnit)

            !                 set ion energy (EQ)
            !AtomData%elem(iE)%ion(k)%EQ = AtomData%elem(iE)%ion(k)%Qd(1)%EQn            ! ion energy (is taken = EQ1)
            eps = 0.d0
            do j = 1,AtomData%elem(iE)%ion(k)%Nqd
                EQd = AtomData%elem(iE)%ion(k)%Qd(j)%EQd
                eps = min(eps,EQd)                       ! Q-ion energy (is taken as {EQnl} lower energy value)
            enddo
            AtomData%elem(iE)%ion(k)%EQ = eps
        enddo

        ! . . . . . read element ionization energies:
        fpath = trim(dbPath)//"/"//chtitle//".ion"
        call slash_correct(fpath,OSType)
        open(atmUnit+1,file=trim(fpath))
        read(atmUnit+1,*) sbuf
        !write(*,*) trim(sbuf)
        do j = 1,AtomData%elem(iE)%Z
            read(atmUnit+1,*) AtomData%elem(iE)%ion(j)%dEion
        end do
        close(atmUnit+1)

        return
    end subroutine AtomDataRead

end module module_AtomData

! ********************************************************************
subroutine AtomicDataReading(NumAt,AtomNames)
    use module_AtomData
    implicit none
    ! input
    integer                          :: NumAt
    character(2),dimension(NumAt)    :: AtomNames
    ! local
    character(128)                   :: dbPath
    integer                          :: iE

    AtomData%NElems = NumAt
    allocate(AtomData%elem(AtomData%NElems))


    !      write(*,'(a)',advance='no') 'AtomAverageQNL2QN: '
    !      do iE = 1,AtomData%NElems
    !            write(*,'(1x,a)',advance='no') AtomNames(iE)
    !            dbPath = atomicDBPaths(iE)
    !            call AtomAverageQNL2QN(dbPath,AtomNames(iE))
    !      enddo
    !      write(*,*)

    write(*,*) 'AtomDataReading ...'
    do iE = 1,AtomData%NElems
        dbPath = AtomicDBPaths(iE)
        write(*,*) trim(atomicDBPaths(iE))
        call AtomDataRead(iE,AtomNames(iE),dbPath)
    enddo

    !   set configurations list
    write(*,*)
    write(*,*) 'SetAtomsQdList'
    call SetAtomsQdList()

    !   set bound-bound transitions list
    write(*,*) 'SetAtomsQtbbList'
    call SetAtomsQtbbList()

    !   set bound-free transitions list
    write(*,*) 'SetAtomsQtbfList'
    call SetAtomsQtbfList()

    write(*,*) "AtomicDataReading: done"
    write(*,*)

    return
end subroutine AtomicDataReading

! **********************************************************
subroutine SetAtomsQdList
    use module_AtomData
    implicit none

    integer :: iE,iI,iQ
    integer :: idx

    idx = 0
    do iE = 1,AtomData%NElems
        do iI = 1,AtomData%elem(iE)%NIons
            do iQ = 1,AtomData%elem(iE)%ion(iI)%Nqd
                idx = idx + 1
                AtomData%elem(iE)%ion(iI)%Qd(iQ)%iList = idx
            enddo
        enddo
    enddo

    AtomData%NQds = idx
    write(*,*) 'AtomData%NQds = ',AtomData%NQds
    allocate(AtomData%QdList(AtomData%NQds))
    idx = 0
    do iE = 1,AtomData%NElems
        do iI = 1,AtomData%elem(iE)%NIons
            do iQ = 1,AtomData%elem(iE)%ion(iI)%Nqd
                idx = idx + 1
                AtomData%QdList(idx)%iE = iE
                AtomData%QdList(idx)%iI = iI
                AtomData%QdList(idx)%iQ = iQ
            enddo
        enddo
    enddo

    return
end subroutine SetAtomsQdList

! **********************************************************
subroutine SetAtomsQtbbList()
    ! list of bound-bound processes (excitation)
    use module_Settings
    use module_AtomData
    use module_TxtRead
    implicit none

    ! input
    ! local
    integer :: i,j,k
    integer :: i1,i2,iLQd1,iLQd2
    integer :: idx,ios
    character(100) :: fpath
    character(2)   :: chtitle
    real(8)        :: dE12,gf12,gQd1,gQd2

    AtomData%NQt_bb = 0
    do i = 1,AtomData%NElems
        do j = 1,AtomData%elem(i)%NIons
            AtomData%NQt_bb = AtomData%NQt_bb + AtomData%elem(i)%ion(j)%NQt_bb
        enddo
    enddo
    write(*,*) 'AtomData%NQt_bb = ',AtomData%NQt_bb

    allocate(AtomData%QtbbList(AtomData%NQt_bb))
    idx = 0
    do i = 1,AtomData%NElems
        do j = 1,AtomData%elem(i)%NIons
            if(len(trim(AtomData%elem(i)%title)).eq.1) then
                chtitle = trim(AtomData%elem(i)%title)//'_'
            else
                chtitle = trim(AtomData%elem(i)%title)
            endif
            fpath = trim(AtomicDBPaths(i))//"/"//chtitle//"qt"
            write(fpath,'(a,i2.2,a)') trim(fpath),j-1,".ch0"
            call slash_correct(fpath,OSType)

            open(1,file=fpath)
            !write(*,*) "fpath: ",fpath
            do k = 1,AtomData%elem(i)%ion(j)%NQt_bb
                read(1,*,iostat=ios) i1,i2,dE12,gf12,gQd1,gQd2
                idx = idx + 1
                iLQd1 = AtomData%elem(i)%ion(j)%Qd(i1)%iList
                iLQd2 = AtomData%elem(i)%ion(j)%Qd(i2)%iList

                AtomData%QtbbList(idx)%iLQd1 = iLQd1
                AtomData%QtbbList(idx)%iLQd2 = iLQd2
                AtomData%QtbbList(idx)%dE12  = dE12
                AtomData%QtbbList(idx)%gf12  = gf12
                AtomData%QtbbList(idx)%gQd1  = gQd1
                AtomData%QtbbList(idx)%gQd2  = gQd2
            enddo
            close(1)
        enddo
    enddo

    return
end subroutine SetAtomsQtbbList

! **********************************************************
subroutine SetAtomsQtbfList()
    ! list of bound-free processes (ionization)
    use module_AtomData
    implicit none

    ! local
    integer :: iE1,iI1,iI2,iQ1,iQ2
    integer :: iLQd1,iLQd2
    integer :: idx,ipos,ilev
    real(8) :: dE12,EQ1,EQ2
    integer(1),dimension(210) :: cells
    integer :: ncells,dsum

    AtomData%NQt_bf = 0
    do iE1 = 1,AtomData%NElems
        AtomData%elem(iE1)%ion(:)%NQt_bf = 0
        AtomData%elem(iE1)%NQt_bf        = 0
        do iI1 = 1,AtomData%elem(iE1)%NIons-1
            do iQ1 = 1,AtomData%elem(iE1)%ion(iI1)%Nqd
                do iQ2 = 1,AtomData%elem(iE1)%ion(iI1+1)%Nqd
                    ncells = AtomData%elem(iE1)%ncells
                    cells = 0
                    cells(1:ncells) = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%cells(1:ncells)
                    cells(1:ncells) = cells(1:ncells)-AtomData%elem(iE1)%ion(iI1+1)%Qd(iQ2)%cells(1:ncells)
                    cells = abs(cells)
                    dsum  = sum(cells)
                    if(dsum.eq.1) then
                        AtomData%elem(iE1)%ion(iI1)%NQt_bf = AtomData%elem(iE1)%ion(iI1)%NQt_bf + 1
                        AtomData%elem(iE1)%NQt_bf = AtomData%elem(iE1)%NQt_bf + 1
                        AtomData%NQt_bf = AtomData%NQt_bf + 1
                    endif
                enddo
            enddo
        enddo
    enddo

    write(*,*) 'AtomData%NQt_bf = ',AtomData%NQt_bf
!    pause

    allocate(AtomData%QtbfList(AtomData%NQt_bf))
    idx = 0
    do iE1 = 1,AtomData%NElems
        do iI1 = 1,AtomData%elem(iE1)%NIons-1
            do iQ1 = 1,AtomData%elem(iE1)%ion(iI1)%Nqd
                iI2 = iI1+1
                do iQ2 = 1,AtomData%elem(iE1)%ion(iI2)%Nqd
                    ncells = AtomData%elem(iE1)%ncells
                    cells = 0
                    cells(1:ncells) = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%cells(1:ncells)
                    cells(1:ncells) = cells(1:ncells)-AtomData%elem(iE1)%ion(iI2)%Qd(iQ2)%cells(1:ncells)
                    cells = abs(cells)
                    dsum  = sum(cells)
                    if(dsum.eq.1) then
                        idx = idx + 1
                        ipos = 1
                        do while(cells(ipos).eq.0)
                            ipos = ipos + 1
                        enddo
                        ilev = ipos
                        if(trim(AtomData%elem(iE1)%QSplit).eq.'nl') then
                            ilev = floor(0.5*(sqrt(1.0+8*ipos)-1)-0.1e-5)+1
                        endif

!                        write(*,*) 'cells_1:',AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%cells(1:ncells)
!                        write(*,*) 'cells_2:',AtomData%elem(iE1)%ion(iI2)%Qd(iQ2)%cells(1:ncells)
!                        write(*,*) 'dsum =  ',dsum
!                        write(*,*) 'ipos =  ',ipos
!                        write(*,*) 'ilev =  ',ilev
!                        pause

                        iLQd1 = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%iList
                        iLQd2 = AtomData%elem(iE1)%ion(iI2)%Qd(iQ2)%iList

                        AtomData%QtbfList(idx)%iLQd1 = iLQd1
                        AtomData%QtbfList(idx)%iLQd2 = iLQd2
                        AtomData%QtbfList(idx)%ilev  = ilev

                        AtomData%QtbfList(idx)%gQd1 = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%gQd
                        AtomData%QtbfList(idx)%gQd2 = AtomData%elem(iE1)%ion(iI2)%Qd(iQ2)%gQd
                        EQ1  = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%EQd
                        EQ2  = AtomData%elem(iE1)%ion(iI2)%Qd(iQ2)%EQd
                        dE12 = EQ2-EQ1
                        AtomData%QtbfList(idx)%dE12 = dE12
!                        write(*,*) 'dE12 = ',dE12
!                        if(dE12.lt.0.d0) pause

                    endif
                enddo
            enddo
        enddo
    enddo


    return
end subroutine SetAtomsQtbfList


! **********************************************************
subroutine AtomAverageQNL2QN(AtomicDBPath,AtomName)
    use module_Settings
    use module_AtomData
    use module_TxtRead
    implicit none

    ! - - - input
    character(64)  :: AtomicDBPath
    character(2)   :: AtomName
    ! - - - local
    integer                :: Nlev,Z
    real(8)                :: A
    integer                :: Zmin,Zmax
    integer                :: NQd,NQt
    integer                :: NQda,NQta
    integer                :: NIons

    integer                :: j,k,i1,i2
    integer                :: iQ1,iQ2
    integer                :: idx
    integer                :: iI,iQ
    integer                :: num,ios
    character(2)           :: chtitle
    character(100)         :: fpath,fmtstr
    real(8),allocatable    :: EQnl(:),gQnl(:)
    real(8),allocatable    :: EQn(:),gQn(:)
    real(8)                :: dE12,gf12,gQ1,gQ2
    integer, allocatable   :: MtxQd(:,:),ListQnCells(:,:),QnCells(:)
    integer                :: ncls
    integer(1),allocatable :: cells(:)         ! electrons setling in orbitals
    real(8), allocatable   :: MtxEQt(:,:),MtxGQt(:,:),MtxGfQt(:,:)
    integer                :: inf_unit,infa_unit


    inf_unit  = 10
    infa_unit = 100

    if(len(trim(AtomName)).eq.1) then
        chtitle = trim(AtomName)//'_'
    else
        chtitle = trim(AtomName)
    endif
    fpath = trim(AtomicDBPath)//"/"//chtitle//"_inf.ch0"
    call slash_correct(fpath,OSType)
    open(inf_unit,file=fpath)
    read(inf_unit,*) Nlev,Z,A
    read(inf_unit,*) Zmin,Zmax

    fpath = trim(AtomicDBPath)//"/"//chtitle//"_infa.ch0"
    call slash_correct(fpath,OSType)
    open(infa_unit,file=fpath)
    write(infa_unit,'(2(i10),2x,1pE12.4e3)') Nlev,Z,A
    write(infa_unit,'(2(i10))')           Zmin,Zmax


    ncls = ((1+Nlev)*Nlev)/2             ! number of electron "cells"
    allocate(cells(ncls))
    NIons = Zmax-Zmin+1
    do j = 1,NIons

        if(len(trim(AtomName)).eq.1) then
            chtitle = trim(AtomName)//'_'
        else
            chtitle = trim(AtomName)
        endif
        fpath = trim(AtomicDBPath)//"/"//chtitle//"qd"
        call slash_correct(fpath,OSType)
        write(fpath,'(a,i2.2,a)') trim(fpath),j-1,".ch0"
        !            write(*,'(a)') trim(fpath)
        ! get number of configurations Qnl (Nqd)
        read(inf_unit,*) k,NQd,NQt
        !        open(1,file=fpath)
        !        nQd = 0
        !        do
        !            read(1,*,iostat=ios,end=10)
        !            nQd = nQd + 1
        !        enddo
        !10      close(1)
        allocate(MtxQd(NQd,NQd))
        allocate(ListQnCells(NQd,Nlev))
        allocate(QnCells(Nlev))
        allocate(EQnl(NQd),gQnl(NQd))


        MtxQd       = 0
        ListQnCells = 0
        ! read configurations Qnl and join Qnl with same "n"
        open(1,file=fpath)
        do k = 1,NQd
            read(1,*,iostat=ios) idx,iI,EQnl(k),gQnl(k),cells(1:ncls)

            num        = 0
            QnCells(:) = 0
            ! set n-levels populations
            do i1 = 1,Nlev
                do i2 = 1,i1
                    num = num + 1
                    QnCells(i1) = QnCells(i1) + cells(num)
                enddo
            enddo
            ! compare current levels configuration with others and set Qnl to Qn
            do i1 = 1,NQd
                if(all(QnCells(:).eq.ListQnCells(i1,:))) then
                    do i2 = 1,NQd
                        if(MtxQd(i1,i2).eq.0) then
                            MtxQd(i1,i2) = k
                            exit
                        endif
                    enddo
                    exit
                endif
                if(all(ListQnCells(i1,:).eq.0)) then
                    ListQnCells(i1,:) = QnCells(:)
                    MtxQd(i1,1) = k
                    exit
                endif
            enddo
        enddo
        close(1)


        ! create <>qda.ch0 files - averaged <>qd.ch0 files
        allocate(EQn(NQd),gQn(NQd))
        fpath = trim(AtomicDBPath)//"/"//chtitle//"qd"
        write(fpath,'(a,i2.2,a,a)') trim(fpath),j-1,"a",".ch0"
        call slash_correct(fpath,OSType)
        !            write(*,'(a)') trim(fpath)
        ! get number of configurations Qnl
        open(1,file=fpath)
        NQda = 0
        EQn = 0.d0
        gQn = 0.d0
        do k = 1,NQd
            if(any(MtxQd(k,:).ne.0)) then
                i1 = 1
                do i1 = 1,NQd
                    if(MtxQd(k,i1).ne.0) then
                        idx = MtxQd(k,i1)
                        EQn(k) = EQn(k) + EQnl(idx)*gQnl(idx)
                        gQn(k) = gQn(K) + gQnl(idx)
                    endif
                enddo
                EQn(k) = EQn(k)/gQn(k)
                fmtstr = ''
                write(fmtstr,'(a,i3,a)') '(2(i6),2(1pE15.7e2),',Nlev,'(i3))'
                write(1,fmtstr) k,iI,EQn(k),gQn(k),ListQnCells(k,:)
                NQda = k
            endif
        enddo
        close(1)
        deallocate(QnCells)
        deallocate(ListQnCells)
        deallocate(EQnl,gQnl)


        ! create <>qta.ch0 files - averaged <>qt.ch0 files
        allocate(MtxEQt(NQda,NQda),MtxGQt(NQda,NQda),MtxGfQt(NQda,NQda))
        MtxEQt  = 0
        MtxGQt  = 0
        MtxGFQt = 0

        fpath = trim(AtomicDBPath)//"/"//chtitle//"qt"
        write(fpath,'(a,i2.2,a)') trim(fpath),j-1,".ch0"
        call slash_correct(fpath,OSType)
        !            write(*,'(a)') trim(fpath)
        open(1,file=fpath)
        do k = 1,NQt
            read(1,*,iostat=ios) i1,i2,dE12,gf12,gQ1,gQ2
            iQ1 = -1
            iQ2 = -1
            do iQ = 1,NQd
                do idx = 1,NQd
                    if(MtxQd(iQ,idx).eq.i1) then
                        iQ1 = iQ
                        exit
                    endif
                enddo
            enddo
            do iQ = 1,NQd
                do idx = 1,NQd
                    if(MtxQd(iQ,idx).eq.i2) then
                        iQ2 = iQ
                        exit
                    endif
                enddo
            enddo
            if(any([iQ1,iQ2].eq.-1)) then
                write(*,*) 'error: qt average - configuration is not found'
                stop
            endif

            if(iQ1.ne.iQ2) then
                MtxEQt(iQ1,iQ2)  = MtxEQt(iQ1,iQ2)  + dE12*gQ1*gQ2
                MtxGfQt(iQ1,iQ2) = MtxGfQt(iQ1,iQ2) + gf12*gQ1*gQ2
                MtxGQt(iQ1,iQ2)  = MtxGQt(iQ1,iQ2)  + gQ1*gQ2
            endif
        enddo
        close(1)


        do iQ1 = 1,NQda
            do iQ2 = 1,NQda
                if(MtxGQt(iQ1,iQ2).ne.0) then
                    !MtxEQt(iQ1,iQ2)  = MtxEQt(iQ1,iQ2)/MtxGQt(iQ1,iQ2)
                    MtxEQt(iQ1,iQ2)  = EQn(iQ2)-EQn(iQ1)
                    MtxGfQt(iQ1,iQ2) = MtxGfQt(iQ1,iQ2)/MtxGQt(iQ1,iQ2)
                endif
            enddo
        enddo


        fpath = trim(AtomicDBPath)//"/"//chtitle//"qt"
        write(fpath,'(a,i2.2,a,a)') trim(fpath),j-1,"a",".ch0"
        call slash_correct(fpath,OSType)
        !            write(*,'(a)') trim(fpath)
        open(1,file=fpath)
        fmtstr = '(2(i7),2(1pE15.7e2),2(1pE12.4e2))'
        NQta = 0
        do iQ1 = 1,NQda
            do iQ2 = 1,NQda
                dE12 = MtxEQt(iQ1,iQ2)
                gf12 = MtxGfQt(iQ1,iQ2)
                if(dE12.gt.0) then
                    NQta = NQta + 1
                    write(1,fmtstr) iQ1,iQ2,dE12,gf12,gQn(iQ1),gQn(iQ2)
                endif
            enddo
        enddo
        close(1)

        deallocate(MtxEQt,MtxGQt,MtxGfQt)
        deallocate(MtxQd,EQn,gQn)

        write(infa_unit,'(3(i10))') Zmin+j-1,NQda,NQta
    enddo
    deallocate(cells)
    close(inf_unit)
    close(infa_unit)


    !n1 = floor(0.5*(sqrt(1.0+8*n1)-1)-0.1e-5)+1

    return
end subroutine AtomAverageQNL2QN

! **********************************************************
subroutine outputAtomStructure
    use module_AtomData
    use module_TxtRead
    use module_Settings

    implicit none

    integer        :: iE,iI,iQ
    integer        :: ios
    character(128) :: sbuf
    logical :: dirExists

    do iE = 1,AtomData%NElems
        write(*,*) trim(AtomData%elem(iE)%title),AtomData%elem(iE)%NQd,AtomData%elem(iE)%NQt_bb
        sbuf = trim(OutputDir)//'levStruct'
        call slash_correct(sbuf,OSType)
        inquire(file=trim(sbuf)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
        if(.not.dirExists) call system('mkdir '//trim(sbuf),ios)
        sbuf = trim(OutputDir)//'levStruct\levStr_'//trim(AtomData%elem(iE)%title)//'.txt'
        call slash_correct(sbuf,OSType)

        open(10,file=sbuf)
        write(10,'(a,a)')  '# atomTitle: ',trim(AtomData%elem(iE)%title)
        write(10,'(a,i3)') '# numIons  = ',AtomData%elem(iE)%NIons
        write(10,'(a,i3)')
        do iI = 1,AtomData%elem(iE)%NIons
            write(10,'(a,i3)') '# ion:     ',iI-1
            write(10,'(a,i3)') '# numQds = ',AtomData%elem(iE)%ion(iI)%Nqd
            do iQ = 1,AtomData%elem(iE)%ion(iI)%Nqd
                write(10,'(i4,2x,1pE17.9e3)') (iI-1),AtomData%elem(iE)%ion(iI)%Qd(iQ)%EQd
!                write(10,'(f4.1,2x,1pE17.9e3)') (iI-1)+0.2e0,AtomData%elem(iE)%ion(iI)%Qd(iQ)%EQd
!                write(10,*)
            enddo
            write(10,*)
        enddo
        close(10)
    enddo

    return
end subroutine outputAtomStructure



! **********************************************************
subroutine DeallocAtomicData
    use module_AtomData
    integer        :: iE,iI,iQ

    do iE = 1,AtomData%NElems
        do iI = 1,AtomData%elem(iE)%NIons
            do iQ = 1,AtomData%elem(iE)%ion(iI)%Nqd
                deallocate(AtomData%elem(iE)%ion(iI)%Qd(iQ)%cells)
            enddo
            deallocate(AtomData%elem(iE)%ion(iI)%Qd)
        enddo
        deallocate(AtomData%elem(iE)%ion)
    enddo
    deallocate(AtomData%elem)
    deallocate(AtomData%QdList)

    return
end subroutine DeallocAtomicData



! ============================ subroutines =============================
subroutine init_mod_atomdata(cfgName)
    use module_Settings
    use module_AtomData
    use module_Constants
    use module_TxtRead

    implicit none
    character(*) :: cfgName
    integer :: NElems
    character(100) :: strbuf
    character(10),dimension(10) :: args
    character(2),allocatable    :: AtomNames(:)
    integer :: ios
    integer(i4p)  :: cfgUnit,GetUnit
    !logical :: dirExists

    ! ------------------ create output directory ------------------
    strbuf = trim(OutputDir)
    call slash_end_set(strbuf)
    call slash_correct(strbuf,OSType)
!    inquire(file=trim(strbuf)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
!    if(.not.dirExists) call system('mkdir '//trim(strbuf),ios)


    !   ---------------- reading config: program atomic setup -----------------------------
    cfgUnit = GetUnit()
    open(cfgUnit,file=cfgName)
    ! read atom names
    call txt_read_array(cfgUnit,"AtomNames",args,NElems,ios)
    allocate(AtomNames(NElems))
    read(args(1:NElems),*) AtomNames(1:NElems)   ! operation removes quotes '' if exist ('He'->He)

    close(cfgUnit)
    !   ------------------------ config reading done ---------


    !   --------- reading atomic data from data base ----------------
    call GetAtomicDBPaths(cfgUnit,cfgName,NElems,AtomNames)
    call AtomicDataReading(NElems,AtomNames)
    !call outputAtomStructure()
    !pause

    deallocate(AtomNames)

    return
end subroutine init_mod_atomdata


! ---------------------------------------------
subroutine close_mod_atomdata()
    use module_AtomData
    implicit none

    call DeallocAtomicData()

    return
end subroutine close_mod_atomdata


! ---------------------------------------------
subroutine GetAtomicDBPaths(cfgUnit,cfgName,NElems,AtomNames)
    use module_Settings
    use module_TxtRead
    use module_AtomData
    implicit none
    ! input parameters
    integer                        :: cfgUnit
    character(*)                   :: cfgName
    integer                        :: NElems
    character(2),dimension(NElems) :: AtomNames
    ! local parameters
    integer       :: iE
    integer       :: ios

    open(cfgUnit,file=cfgName)
    do iE = 1,NElems
        call txt_read_struct_value(cfgUnit,'AtomicDataBases',AtomNames(iE),AtomicDBPaths(iE),ios)
    enddo
    close(cfgUnit)

    write(*,*)
    do iE = 1,NElems
        call slash_end_del(atomicDBPaths(iE))
        call slash_correct(atomicDBPaths(iE),OSType)
!        write(*,*) trim(atomicDBPaths(iE))
    enddo

    return
end subroutine GetAtomicDBPaths

! ---------------------------------------------
subroutine ScreenOutput_SetupInf(TeV,Ntot)
    use module_AtomData
    use module_Constants


    implicit none
    ! - - - - - input parameters
    real(8),dimension(AtomData%NElems) :: Ntot
    ! - - - - - local parameters
    real(8) :: TK,TeV
    integer :: i,j,num

    ! screen output
    write(*,*)
    write(*,*) ' --------- input parameters ---------- '
    TK = TeV*eV_erg/k_Bol
    write(*,"(a,1pE9.3e2,a)",advance="no")  'TK   = ',TK,' K'
    write(*,"(a,1pE9.3e2,a)",advance="yes") ' (',TeV,' eV)'

    write(*,"(a)",advance="no")  'Ntot = '
    do i=1,AtomData%NElems
        write(*,'(3x,1pE9.3e2)',advance='no') Ntot(i)
    enddo
    write(*,*)

    write(*,*)
    write(*,"(a)",advance='no') 'elements  :   '
    do i=1,AtomData%NElems
        write(*,'(3x,a2,5x)',advance='no') trim(AtomData%elem(i)%title)
    enddo

    write(*,*)
    write(*,"(a)",advance='no') 'N(1/cm**3): '
    do i=1,AtomData%NElems
        write(*,'(1pE10.2e2)',advance='no') Ntot(i)
    enddo
    write(*,*)

    write(*,"(a)",advance='no') 'Number Ions:  '
    do i=1,AtomData%NElems
        write(*,'(i5,5x)',advance='no') AtomData%elem(i)%NIons
    enddo
    write(*,*)
    write(*,"(a)",advance='no') 'Number Qds:   '
    do i=1,AtomData%NElems
        num = 0
        do j=1,AtomData%elem(i)%NIons
            num = num + AtomData%elem(i)%ion(j)%Nqd
        enddo
        write(*,'(i5,5x)',advance='no') num
    enddo
    write(*,*)
    write(*,"(a)",advance='no') 'Number Qts:   '
    do i=1,AtomData%NElems
        num = 0
        do j=1,AtomData%elem(i)%NIons
            num = num + AtomData%elem(i)%ion(j)%NQt_bb
        enddo
        write(*,'(i5,5x)',advance='no') num
    enddo
    write(*,*)
    write(*,*)
    write(*,*)

    return
end subroutine ScreenOutput_SetupInf
