! By WQY
! <!>  - mark places for check
module module_RKModel
    use module_TxtRead
    use module_RayTrace
    use module_GridFunct
    use module_AtomData
    use module_Mesh
    use module_Spectrum
    use module_LevKinetics

!   ... hydrodynamics mesh ...
    type(TMeshData),        target    :: meshH        ! 2D mesh for hydrodynamics
    type(TGridNodalScalar), target    :: densH        ! density
    type(TGridNodalScalar), target    :: tempH        ! temperature
    type(TGridNodalScalar), target    :: IonDegH      ! ioniztion degree
    type(TGridNodalScalar), target    :: cUdenH       ! integral (by photon energy) radiation energy density
    type(TGridNodalVector), target    :: WfluxH       ! integral (by photon energy) radiation energy flux

!   ... rad.transport mesh ...
    type(TMeshData),        target    :: meshR        ! 2D mesh for radiation transport
    type(TGridNodalScalar), target    :: densR        ! density
    type(TGridNodalScalar), target    :: tempR        ! temperature
    type(TGridNodalScalar), target    :: IonDegR      ! ioniztion degree
    type(TGridNodalScalar), target    :: cUden        ! integral (by photon energy) radiation energy density
    type(TGridNodalVector), target    :: Wflux        ! integral (by photon energy) radiation energy flux
    type(TGridCellVector)             :: kappa        ! spectral absorption coefficient
    type(TGridCellVector)             :: etta         ! spectral emissivity

!   ... interpolators for values: hydro.mesh <-> rad.mesh ...
    type(TValuesInterpolator) :: Interp_HR    ! interpolator: hydro.mesh -> rad.mesh
    type(TValuesInterpolator) :: Interp_RH    ! interpolator: rad.mesh -> hydro.mesh

!   ... ray tracer data ...
!    type(TTraceData),       target    :: tracer
    type(TTraceDataCompact),target  :: tracer
    type(TAngularGridData), target  :: AngGrid

contains

! ============================ subroutines =============================
    subroutine read_trace1D_prm(ftrace,Trace1DData,pMesh,Npts)
        use module_Settings
        use module_Mesh
        implicit none
        ! input:
        character(*)                        :: ftrace
        type(TTrace1DData)                  :: Trace1DData
        type(TMeshData),  pointer           :: pMesh         ! mesh to trace
        integer                             :: Npts
        ! local:
        integer        :: i,ipos,ios
        integer        :: ipt,iface
        integer        :: NFaces
        character(128) :: sbuf
        character(1)   :: ch1


        Trace1DData%LenUnit = 1.0d0
        Trace1DData%pMesh   => pMesh
        NFaces = pMesh%NFaces
        Trace1DData%Nfaces  = NFaces
        Trace1DData%Npoints = Npts
        call alloc_Trace1DData(Trace1DData,Npts,NFaces)

        open(3,                                     &
            file       = trim(ftrace),              &
            form       = 'UNFORMATTED',             &
            action     = 'READ',                    &
            convert    = 'BIG_ENDIAN',              &
            access     = 'STREAM',                  &
            iostat     = ios)

        ! 1. skip mesh filename
        ipos = 0
        read(3,iostat=ios) sbuf
        i    = index(sbuf,char(10))
        ipos = ipos + i
        read(3,pos=ipos,iostat=ios) ch1
        ! 2. skip Nth
        read(3,iostat=ios) sbuf
        i = index(sbuf,char(10))
        ipos = ipos + i
        read(3,pos=ipos,iostat=ios) ch1
        ! 3. skip Ndirs
        read(3,iostat=ios) sbuf
        i    = index(sbuf,char(10))
        ipos = ipos + i
        read(3,pos=ipos,iostat=ios) ch1

!        open(21,file="prms-1.dat")
!        open(22,file="prms-2.dat")
        do ipt = 1,Npts
            do iface = 1,Nfaces
                read(3,iostat=ios) Trace1DData%point(ipt)%icell0(iface)
                read(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%a
                read(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%b
                read(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%c
                read(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%k
                read(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%idir
!                write(21,"(i5,i5,4(1x,1pE12.5e2))") ipt,iface,                &
!                                    Trace1DData%point(ipt)%prm0(iface)%a,     &
!                                    Trace1DData%point(ipt)%prm0(iface)%b,     &
!                                    Trace1DData%point(ipt)%prm0(iface)%c,     &
!                                    Trace1DData%point(ipt)%prm0(iface)%k

                read(3,iostat=ios) Trace1DData%point(ipt)%icell1(iface)
                read(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%a
                read(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%b
                read(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%c
                read(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%k
                read(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%idir
!                write(22,"(i5,i5,4(1x,1pE12.5e2))") ipt,iface,                &
!                                    Trace1DData%point(ipt)%prm1(iface)%a,     &
!                                    Trace1DData%point(ipt)%prm1(iface)%b,     &
!                                    Trace1DData%point(ipt)%prm1(iface)%c,     &
!                                    Trace1DData%point(ipt)%prm1(iface)%k
            end do
        end do
        close(3)
!        close(21)
!        close(22)
!        pause

        return
    end subroutine

    ! ==============================================================
    subroutine check_trace1D_file(icase,ftrace,pMesh,pAngGrid)
        implicit none
        ! input:
        character(*)                    :: ftrace
        type(TMeshData), pointer        :: pMesh
        type(TAngularGridData), pointer :: pAngGrid
        ! output:
        integer                         :: icase
        ! local:
        integer        :: Ndirs,Nth
        character(64)  :: meshTitIni,meshTitAct
        character(128) :: sbuf
        integer        :: n1,n2
        integer        :: i,ipos,ios
        logical        :: file_exists
        character(1)   :: ch1


        icase = 0
        inquire(file=ftrace, exist=file_exists)
        if(file_exists) then

            ! set mesh filename: actual
            meshTitAct = trim(pMesh%meshTit)
            n1 = 1
            n2 = 1
            do while(n1.gt.0 .or. n2.gt.0)
                n1 = index(meshTitAct,"\")
                meshTitAct = meshTitAct(n1+1:)
                n2 = index(meshTitAct,"/")
                meshTitAct = meshTitAct(n2+1:)
            end do
            meshTitAct = adjustl(meshTitAct)

            ! open file for reading header information:
            open(3,                                     &
                file       = "data_in\Trace1DData.dat", &
                form       = 'UNFORMATTED',             &
                action     = 'READ',                    &
                convert    = 'BIG_ENDIAN',              &
                access     = 'STREAM',                  &
                iostat     = ios)

            ! 1. ---------------------
            ipos = 0
            read(3,iostat=ios) sbuf
            i    = index(sbuf,char(10))
            sbuf = sbuf(:i-1)
            ipos = ipos + i
            read(3,pos=ipos,iostat=ios) ch1     ! set pos. in file for next data block read (read 1 char)
            ! set mesh filename: initial
            meshTitIni = trim(sbuf(:i-1))
            i    = index(sbuf,":")
            meshTitIni = trim(sbuf(i+1:))
            n1 = 1
            n2 = 1
            do while(n1.gt.0 .or. n2.gt.0)
                n1 = index(meshTitIni,"\")
                meshTitIni = meshTitIni(n1+1:)
                n2 = index(meshTitIni,"/")
                meshTitIni = meshTitIni(n2+1:)
            end do
            meshTitIni = adjustl(meshTitIni)

            ! 2. ---------------------
            read(3,iostat=ios) sbuf
            i = index(sbuf,char(10))
            sbuf = sbuf(:i-1)
            ipos = ipos + i
            read(3,pos=ipos,iostat=ios) ch1     ! set pos. in file for next data block read (read 1 char)
            ! set Nth
            n1 = index(sbuf,"=")
            sbuf = sbuf(n1+1:)
            read(sbuf,*) Nth

            ! 3. ---------------------
            read(3,iostat=ios) sbuf
            i    = index(sbuf,char(10))
            sbuf = sbuf(:i-1)
            ipos = ipos + i
            read(3,pos=ipos,iostat=ios) ch1     ! set pos. in file for next data block read (read 1 char)
            ! set Ndirs
            n1 = index(sbuf,"=")
            sbuf = sbuf(n1+1:)
            read(sbuf,*) Ndirs



            ! check settings:
            icase = 0
            if(trim(meshTitAct).eq.trim(meshTitIni)) then
                icase = icase + 1
            end if
            if(Nth.eq.pAngGrid%NTheta) then
                icase = icase + 1
            end if
            if(Ndirs.eq.pAngGrid%NDirs) then
                icase = icase + 1
            end if

            close(3)
        endif

        return
    end subroutine

end module

    ! ==============================================================
subroutine init_RKModel(CaseDir,TDPDir,NumStepsTotal,NumStepsPart,RunRadEvery)
    use module_RKModel
    use module_Constants
    use module_OutputMesh
    use module_Interpolate
    implicit none
    ! input:
    character(128)                  :: CaseDir
    ! output:
    character(128)                  :: TDPDir
    integer                         :: NumStepsTotal
    integer                         :: NumStepsPart
    integer                         :: RunRadEvery

    ! local:
    character(128)                  :: cfgName            ! config file name
    type(TMeshData),        pointer :: pmeshH   => null() ! hydro.mesh pointer
    type(TMeshData),        pointer :: pmeshR   => null() ! rad.mesh pointer
!    type(TTraceData),       pointer :: ptrace   => null() ! tracer pointer (extended version)
    type(TTraceDataCompact),pointer :: ptrace   => null() ! tracer pointer (compact version)
    type(TAngularGridData), pointer :: pAngGrid => null() ! angular grid pointer

!     ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!    type TScalarContainer
!        type(TGridNodalScalar), pointer :: scl => null()
!    end type TScalarContainer
!    type TVectorContainer
!        type(TGridNodalVector), pointer :: vec => null()
!    end type TVectorContainer
!    type(TScalarContainer), allocatable, dimension(:) :: ContScl     ! container for scalar grid functions
!    type(TVectorContainer), allocatable, dimension(:) :: ContVec     ! container for vector grid functions
!    integer :: NScls,NVecs

    ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    real(r8p), allocatable :: Rint0(:),Rext0(:)      ! line point coordinates
    real(r8p), allocatable :: Rint(:),Rext(:)        ! line point coordinates
    integer                :: NR0,NZ0
    integer                :: NXpt,NXgr
    integer                :: ipt,num,ios
    character(128)          :: FileMesh
    integer                 :: NZ,NR
    integer                 :: StepR,StepZ
    integer                 :: SymType
    real(r8p), dimension(3) :: SymAxis
    character(16)           :: RayType,LocTit,ElemType
    !character(16)           :: ftit
    real(r8p)               :: BaseLen      ! character length
    integer                 :: cfgUnit,GetUnit
    logical                 :: dirExists


    ! ... get OS type (win, lin...) for slash correction
    call SetOSType()       ! value -> mod_Settings/OSType
    call slash_end_set(CaseDir)
    call slash_correct(CaseDir,OSType)

    TDPDir = trim(CaseDir)//"/TDP/"
    call slash_correct(TDPDir,OSType)
    inquire(file=trim(TDPDir)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
    if(.not.dirExists) call system('mkdir '//trim(TDPDir))


    ! ... prepare for open config.txt file
    cfgUnit   = GetUnit()
    OutputDir = trim(CaseDir)
    cfgName   = trim(OutputDir)//'config.txt'

    ! ... set base settings (dimDen, dimVel, dimTemp, ...)
    call set_base_settings(cfgName)

    ! ... read number of steps from config.txt
    open(cfgUnit,file=cfgName)
        call txt_read_value(cfgUnit,'NT',         NumStepsTotal,ios)
        call txt_read_value(cfgUnit,'NSTEP',      NumStepsPart, ios)
        call txt_read_value(cfgUnit,'RunRadEvery',RunRadEvery, ios)
    close(cfgUnit)

    ! ... set atomic data
    call init_mod_atomdata(cfgName)

    ! ... set spectral data
    call init_mod_spectrum(cfgName)
    NXpt = SpectrData%EnGrid%NX
    NXgr = SpectrData%EnGrid%NX - 1
    write(*,*) "NumSpectalGroups = ",NXgr
    write(*,*)

    ! ... set hydro.mesh (reading parameters from "config.txt")
    open(cfgUnit,file=cfgName)
        call txt_read_struct_value(cfgUnit,'meshH','NZ',NZ0,ios)
        call txt_read_struct_value(cfgUnit,'meshH','NR',NR0,ios)
        call txt_read_struct_value(cfgUnit,'meshH','file',FileMesh,ios)
        call txt_read_struct_value(cfgUnit,'meshH','length',BaseLen,ios)
        call txt_read_struct_value(cfgUnit,'meshH','SymType',SymType,ios)
        call txt_read_struct_array(cfgUnit,'meshH','SymAxis',SymAxis,num,ios)
    close(cfgUnit)

    ! read 2D mesh
    allocate(Rint0(NZ0),Rext0(NZ0))
    open(4,file=trim(OutputDir)//FileMesh)
        read(4,*) Rext0(:),Rint0(:)
    close(4)
    meshH%meshTit = 'simp2Dmesh_H'
    call new_mesh_simple2D(NZ0,NR0,Rint0,Rext0,BaseLen,meshH)
    pmeshH => meshH


    ! ... set rad.mesh (reading parameters from "config.txt")
    open(cfgUnit,file=cfgName)
        call txt_read_struct_value(cfgUnit,'meshR','NZ',NZ0,ios)
        call txt_read_struct_value(cfgUnit,'meshR','NR',NR0,ios)
        call txt_read_struct_value(cfgUnit,'meshR','file',FileMesh,ios)
        call txt_read_struct_value(cfgUnit,'meshR','length',BaseLen,ios)
        call txt_read_struct_value(cfgUnit,'meshR','SymType',SymType,ios)
        call txt_read_struct_array(cfgUnit,'meshR','SymAxis',SymAxis,num,ios)
        call txt_read_struct_value(cfgUnit,'meshR','StepZ',StepZ,ios)
        call txt_read_struct_value(cfgUnit,'meshR','StepR',StepR,ios)
    close(cfgUnit)
    NZ = (NZ0-1)/StepZ+1
    NR = (NR0-1)/StepR+1
    allocate(Rint(NZ),Rext(NZ))
    do ipt = 1,NZ
        Rint(ipt) = Rint0((ipt-1)*StepZ+1)
        Rext(ipt) = Rext0((ipt-1)*StepZ+1)
    end do
    meshR%meshTit = 'simp2Dmesh_R'
    call new_mesh_simple2D(NZ,NR,Rint,Rext,BaseLen,meshR)
    pmeshR => meshR

    ! ... output meshH and meshR
    call output_2DMesh_gmsh(pmeshH,trim(OutputDir)//"meshH.msh")
    call output_2DMesh_gmsh(pmeshR,trim(OutputDir)//"meshR.msh")

    ! ... create interpolators:
    call new_interp_ValExtNodesToNodes(Interp_HR,pmeshH,pmeshR)
    call new_interp_ValExtNodesToNodes(Interp_RH,pmeshR,pmeshH)


    write(*,*)
    write(*,*) "....... new grid functions ........"

!   ... hydrodynamics mesh ...
    call new_grid_function(densH,"PLT",pmeshH)
    call new_grid_function(tempH,"TEM",pmeshH)
    call new_grid_function(IonDegH,"IonDeg",pmeshH)
    ! cU - radiation energy density
    call new_grid_function(cUdenH,"cUden",pmeshH)
    ! Wz - radiation energy flux
    call new_grid_function(WfluxH,"Wflux",3,pmeshH)

!   ... rad.transport mesh ...
    call new_grid_function(densR,"PLT",pmeshR)
    call new_grid_function(tempR,"TEM",pmeshR)
    call new_grid_function(IonDegR,"IonDeg",pmeshR)
    ! cU - radiation energy density
    call new_grid_function(cUden,"cUden",pmeshR)
    ! Wz - radiation energy flux
    call new_grid_function(Wflux,"Wflux",3,pmeshR)
    ! kappa - spectral coefficient
    call new_grid_function(kappa,"kappa",NXgr,pmeshR)
    ! etta  - spectral coefficient
    call new_grid_function(etta, "etta", NXgr,pmeshR)


    ! ... set angular grid for tracing
    write(*,*)
    write(*,*) "....... create angular grid ........"
    call new_angular_grid(cfgName,AngGrid)
    pAngGrid => AngGrid
    call output_angular_grid(AngGrid,[0.0d0,0.0d0,0.0d0])

    RayType  = "long"
    LocTit   = "channel"
    ElemType = "node"

    call new_rays_tracer(RayType,pAngGrid,pmeshR,LocTit,ElemType,tracer)
    ptrace => tracer

    return
end subroutine

    ! ==============================================================
subroutine rad_transport(ctime,NZ,NR,PLT,TEM,CU2D,WR2D,WZ2D)
    use module_Constants
    use module_Mesh
    use module_RKModel
    implicit none
    ! input:
    integer                   :: NZ,NR
    real(4), dimension(NZ,NR) :: PLT,TEM
    real(4)                   :: ctime
    ! output:
    real(4), dimension(NZ,NR) :: CU2D,WR2D,WZ2D
    ! local:
!    type(TTraceData), pointer        :: ptrace
    type(TTraceDataCompact), pointer :: ptrace
    ! local:
    real(8), dimension(meshR%NCells,SpectrData%EnGrid%NX-1) :: CellKppa
    real(8), dimension(meshR%NCells,SpectrData%EnGrid%NX-1) :: CellEtta
    real(8), dimension(SpectrData%EnGrid%NX-1) :: cKppa
    real(8), dimension(SpectrData%EnGrid%NX-1) :: cEtta
    real(8), dimension(AtomData%NElems) :: Ntot
    real(8), dimension(AtomData%NQds)   :: cNq
    real(8)                             :: TK,TeV
    real(8)                             :: Ne,IonDeg
    integer :: i,j,ipt
    integer :: inode,icell
    integer :: idx
    integer :: NCells
    real(r8p) :: wt
    real(r8p) :: dtime
    ! . . . . . . . . . data containers for output . . . . . . . . .
    type TScalarContainer
        type(TGridNodalScalar), pointer :: scl => null()
    end type TScalarContainer
    type TVectorContainer
        type(TGridNodalVector), pointer :: vec => null()
    end type TVectorContainer
    type(TScalarContainer), allocatable, dimension(:) :: ContScl     ! container for scalar grid functions
    type(TVectorContainer), allocatable, dimension(:) :: ContVec     ! container for vector grid functions
    integer :: NScls,NVecs                                           ! number of scalar and vector functions
    ! ...............................................................


    ! ... set densH and tempH
    ipt = 0
    do i = 1,NZ
        do j = 1,NR
            ipt = ipt + 1
            densH%fval(ipt) = PLT(i,j)
            tempH%fval(ipt) = TEM(i,j)
        end do
    end do

    ! ... interpolate: densH->densR, tempH->tempR
    do i = 1,meshR%NNodes
        densR%fval(i) = 0.0d0
        tempR%fval(i) = 0.0d0
        do j = 1,Interp_HR%ielem(i)%nPush
            idx = Interp_HR%ielem(i)%iPush(j)
            wt  = Interp_HR%ielem(i)%wPush(j)
            densR%fval(i) = densR%fval(i) + wt*densH%fval(idx)
            tempR%fval(i) = tempR%fval(i) + wt*tempH%fval(idx)
        end do
    end do


    ! ... calculate transport coefficients: kppa,etta
    !write(*,*)
    !write(*,*) ' - - - - - - - - '
    !write(*,*) ' kappa, etta ... '
    write(*,"(a)",advance='no') 'ke-'


    NCells = meshR%NCells

    do icell = 1,NCells
        !write(*,"(a,a,I3,a)",advance="no") achar(13),'done ',icell*100/NCells," %"

        ! set dens, temp in cells
        Ntot(:) = 0.0e0
        TK      = 0.0e0
        wt      = 1.0e0/meshR%Cell(icell)%nnodes
        do idx = 1,meshR%Cell(icell)%nnodes
            inode = meshR%Cell(icell)%inode(idx)
            Ntot(1) = Ntot(1) + wt*densR%fval(inode)
            TK      = TK      + wt*tempR%fval(inode)
        enddo

        TK   = TK   * dimTem
        TeV  = TK   * k_Bol/eV_erg
        Ntot = Ntot * dimDen

        call SahaKin_TeV_Nt(TeV,Ntot,cNq,Ne,IonDeg)

        call kappa_etta(TeV,Ne,cNq,cKppa,cEtta)
        CellKppa(icell,:) = cKppa(:)
        CellEtta(icell,:) = cEtta(:)
    end do
!    pause


    ! ... run radiation transport
    write(*,"(a)",advance='no') 'r'
    ptrace => tracer
!    call rad_transport_long_chars(ptrace,CellKppa,CellEtta,cUden,Wflux)
    call rad_transport_long_chars_compact(ptrace,CellKppa,CellEtta,cUden,Wflux)


    ! ... back interpolation: cUden->cUdenH, Wflux->WfluxH
    do i = 1,meshH%NNodes
        cUdenH%fval(i)     = 0.0d0
        WfluxH%fvec(i)%cmp = 0.0d0
        do j = 1,Interp_RH%ielem(i)%nPush
            idx = Interp_RH%ielem(i)%iPush(j)
            wt  = Interp_RH%ielem(i)%wPush(j)
            cUdenH%fval(i)     = cUdenH%fval(i)     + wt*cUden%fval(idx)
            WfluxH%fvec(i)%cmp = WfluxH%fvec(i)%cmp + wt*Wflux%fvec(idx)%cmp
        end do
    end do


    ! ... set cU2D, Wr2D, Wz2D
    ipt = 0
    do i = 1,NZ
        do j = 1,NR
            ipt = ipt + 1
            CU2D(i,j) = sngl(cUdenH%fval(ipt))
            WZ2D(i,j) = sngl(WfluxH%fvec(ipt)%cmp(1))
            WR2D(i,j) = sngl(WfluxH%fvec(ipt)%cmp(2))
        end do
    end do


!    ! ... output tdp and radiation
!    NScls = 3
!    allocate(ContScl(1:NScls))
!    ContScl(1)%scl => densR
!    ContScl(2)%scl => tempR
!    ContScl(3)%scl => cUden
!    NVecs = 1
!    allocate(ContVec(1:NVecs))
!    ContVec(1)%vec => Wflux
!    !ftit = "radH"
!    dtime = dble(ctime)
!    call VtkWriteStep("meshR",NScls,NVecs,ContScl,ContVec,dtime)
!    deallocate(ContScl)
!    deallocate(ContVec)

    ! ... output tdp and radiation
    NScls = 3
    allocate(ContScl(1:NScls))
    ContScl(1)%scl => densH
    ContScl(2)%scl => tempH
    ContScl(3)%scl => cUdenH
    NVecs = 1
    allocate(ContVec(1:NVecs))
    ContVec(1)%vec => WfluxH
    !ftit = "radH"
    dtime = dble(ctime)
    call VtkWriteStep("time",NScls,NVecs,ContScl,ContVec,dtime)
    deallocate(ContScl)
    deallocate(ContVec)

    return
end subroutine

    ! ==============================================================
subroutine close_mod_RKModel()
    use module_RKModel
    implicit none

    ! ----------------------------

    call close_grid_function(densR)
    call close_grid_function(tempR)
    call close_grid_function(cUden)
    call close_grid_function(Wflux)
    call close_grid_function(densH)
    call close_grid_function(tempH)
    call close_grid_function(cUdenH)
    call close_grid_function(WfluxH)
    call close_interpolator(Interp_HR)
    call close_interpolator(Interp_RH)
    call close_rays_tracer(tracer)
    call close_angular_grid(AngGrid)

    call close_mesh(meshH)        ! deallocate all variables
    call close_mesh(meshR)        ! deallocate all variables
    call close_spectral_data()    ! close module spectral data
    !call close_mod_LevKinetics()
    call close_mod_atomdata()

    return
end subroutine

! -------------------------------------------------------------
subroutine output_mesh_bounds(mesh)
    use module_Settings
    use module_Mesh
    use module_TxtRead
    implicit none
    ! input
    type(TMeshData) :: mesh
    ! local
    integer :: iloc,iface,inode
    integer :: j,k
    character(32) :: fname

    fname = trim(OutputDir)//'\bounds.dat'
    call slash_bck2frw(fname)
    open(10,file = fname)
    !write(*,*) " NLocs = ",mesh%NLocs
    do iloc = 1,mesh%NLocs
        if(trim(mesh%Locations(iloc)%locType).eq."bound") then
            do j = 1,mesh%Locations(iloc)%NElems
                iface = mesh%Locations(iloc)%ielem(j)
                do k = 1,mesh%Face(iface)%nNodes
                    inode = mesh%Face(iface)%inode(k)
                    write(10,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
                end do
                inode = mesh%Face(iface)%inode(1)
                write(10,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
                write(10,"(3(1x,1pE12.4e3))")
                write(10,"(3(1x,1pE12.4e3))")
            end do
        endif
    end do
    close(10)

    return
end subroutine output_mesh_bounds

! -------------------------------------------------------------
subroutine output_angular_grid(AngGrid,pt0)
    use module_Settings
    use module_Mesh
    use module_RayTrace
    use module_TxtRead
    implicit none
    ! input
    type(TAngularGridData)  :: AngGrid
    real(r8p), dimension(3) :: pt0
    ! local
    real(r8p), dimension(3) :: ptc
    real(r8p) :: lray
    real(r8p) :: cosTh,cosPhi,sinTh,sinPhi
    real(r8p) :: theta,phi
    integer   :: iOx,iOy,iOz
    integer   :: idir
    logical   :: dirExists
    character(128) :: fname

!    pt0(:) = 0.d0
    lray = 1.0d0
    iOx = 1        ! Oy
    iOy = 2        ! Oz
    iOz = 3        ! axis of symmetry (if 1 - Ox is axis)

    fname = trim(OutputDir)//'\angular'
    call slash_correct(fname,OSType)
    inquire( file=trim(fname)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
    if(.not.dirExists) call system('mkdir '//trim(fname))
    fname = trim(OutputDir)//'\angular\angular_grid.dat'
    call slash_correct(fname,OSType)

    open(10,file = fname)
!    write(*,"(a,i4)") "Ndirs = ",AngGrid%NDirs
    do idir = 1,AngGrid%NDirs
        theta = real(AngGrid%theta(idir),r8p)
        phi   = real(AngGrid%phi(idir),r8p)
        cosTh  = cos(theta)
        cosPhi = cos(phi)
        sinTh  = sin(theta)
        sinPhi = sin(phi)
!        write(*,"(2(1x,1pe12.5e2))") cosTh,cosPhi

        ptc(iOx) = pt0(iOx) + lray*sinTh*cosPhi
        ptc(iOy) = pt0(iOy) + lray*sinTh*sinPhi
        ptc(iOz) = pt0(iOz) + lray*cosTh

        write(10,"(3(1x,1pE12.4e3))") pt0
        write(10,"(3(1x,1pE12.4e3))") ptc
        write(10,"(3(1x,1pE12.4e3))")
        write(10,"(3(1x,1pE12.4e3))")
    end do
    close(10)

    return
end subroutine output_angular_grid


! -------------------------------------------------------------
subroutine output_beam_lighting(trace,ipt)
    use module_Settings
    use module_Mesh
    use module_RayTrace
    use module_TxtRead
    implicit none
    ! input
    type(TTraceData) :: trace
    integer          :: ipt
    ! local
    type(TMeshData), pointer :: mesh
    integer :: iloc,iface,inode,idir
    integer :: k
    character(32) :: fname

    mesh => trace%pmesh
    fname = trim(OutputDir)//'\beam_faces.dat'
    call slash_bck2frw(fname)
    open(10,file = fname)
    fname = trim(OutputDir)//'\beam_points.dat'
    call slash_bck2frw(fname)
    open(11,file = fname)
    do idir = 1,trace%pAngGrid%NDirs
        iloc = trace%point(ipt)%rays(idir)%NSegm
        if(iloc.gt.0) then
            iface = trace%point(ipt)%rays(idir)%iface(iloc)
            do k = 1,mesh%Face(iface)%nNodes
                inode = mesh%Face(iface)%inode(k)
                write(10,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
            end do
            inode = mesh%Face(iface)%inode(1)
            write(10,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
            write(10,"(3(1x,1pE12.4e3))")
            write(10,"(3(1x,1pE12.4e3))")
            write(11,"(3(1x,1pE12.4e3))") trace%point(ipt)%pt0%crd
            write(11,"(3(1x,1pE12.4e3))") trace%point(ipt)%rays(idir)%ptOut(iloc)%crd
            write(11,"(3(1x,1pE12.4e3))")
            write(11,"(3(1x,1pE12.4e3))")
        endif
    end do
    close(10)
    close(11)

    return
end subroutine output_beam_lighting

! ---------------------------------------------------------------
subroutine output_point(funit,pt)
    use module_Settings
    use module_Mesh
    implicit none
    ! input
    integer            :: funit
    type(T3DPoint) :: pt
    ! local

    write(funit,"(3(1x,1pE12.4e3))") pt%crd

    return
end subroutine output_point

! ---------------------------------------------------------------
subroutine output_face(funit,iface,mesh)
    use module_Settings
    use module_Mesh
    implicit none
    ! input
    integer :: funit
    type(TMeshData) :: mesh
    ! local
    integer :: iface,inode
    integer :: k

    do k = 1,mesh%Face(iface)%nNodes
        inode = mesh%Face(iface)%inode(k)
        write(funit,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
    end do
    inode = mesh%Face(iface)%inode(1)
    write(funit,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
    write(funit,"(3(1x,1pE12.4e3))")
    write(funit,"(3(1x,1pE12.4e3))")

    return
end subroutine output_face

! ---------------------------------------------------------------
subroutine output_cell(funit,icell,mesh)
    use module_Settings
    use module_Mesh
    implicit none
    ! input
    integer :: funit
    integer         :: icell
    type(TMeshData) :: mesh
    ! local
    integer :: iface,inode
    integer :: i,k

    do i = 1,mesh%Cell(icell)%nFaces
        iface = mesh%Cell(icell)%iface(i)
        do k = 1,mesh%Face(iface)%nNodes
            inode = mesh%Face(iface)%inode(k)
            write(funit,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
        end do
        inode = mesh%Face(iface)%inode(1)
        write(funit,"(3(1x,1pE12.4e3))") mesh%Node(inode)%crd
        write(funit,"(3(1x,1pE12.4e3))")
        write(funit,"(3(1x,1pE12.4e3))")
    end do

    return
end subroutine output_cell
