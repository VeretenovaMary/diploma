! ------------------------------------------------
subroutine vtk_title_set(pmesh,vtkTitle)
    use module_Settings
    use module_AtomData
    use module_Mesh
    IMPLICIT NONE
    ! input
    type(TMeshData) :: pmesh    ! mesh pointer
    ! output
    character(256)  :: vtkTitle
    ! local
    character(128)  :: fmesh                ! mesh filename (from mesh titile)
    integer         :: ivtkVers
    integer         :: meshSym

    character(256)  :: s_buf
    integer :: i


    ! vtk file version information:
    ivtkVers = 1
    vtkTitle = 'header vers.:'
    s_buf = ''
    write(s_buf,'(i2,a1)') ivtkVers,';'
    vtkTitle = trim(vtkTitle)//trim(s_buf)

    ! mesh information:
    write(s_buf,'(a4,i8,a1)') ' NN=',pmesh%NNodes,';'
    vtkTitle = trim(vtkTitle)//trim(s_buf)
    write(s_buf,'(a4,i8,a1)') ' NC=',pmesh%NCells,';'
    vtkTitle = trim(vtkTitle)//trim(s_buf)
    fmesh = trim(pmesh%meshTit)
    vtkTitle = trim(vtkTitle)//' mesh:'//trim(fmesh)//';'
    meshSym  = pmesh%SymType
    selectcase(meshSym)
    case(0)
        vtkTitle = trim(vtkTitle)//' symm: native;'
    case(1)
        vtkTitle = trim(vtkTitle)//' symm: planar;'
    case(2)
        vtkTitle = trim(vtkTitle)//' symm: cylindrical;'
    endselect

    ! components information:
    write(s_buf,'(a8,i8,a1)') ' NAtoms=',AtomData%NElems,';'
    vtkTitle = trim(vtkTitle)//trim(s_buf)
    vtkTitle = trim(vtkTitle)//' chems:'
    do i = 1,AtomData%NElems
        if(i.lt.AtomData%NElems) then
            vtkTitle = trim(vtkTitle)//' '//trim(AtomData%elem(i)%title)//','
        else
            vtkTitle = trim(vtkTitle)//' '//trim(AtomData%elem(i)%title)//';'
        endif
    enddo

    write(s_buf,'(a4,i8,a1)') ' NQds=',AtomData%NQds,';'

    return
end subroutine vtk_title_set

! ------------------------------------------------
subroutine vtk_title_parse(vtkTitle,NNodes,NCells,meshSym)
    use module_Settings
    use module_AtomData
    use module_Mesh
    IMPLICIT NONE
    ! input
    character(*)     :: vtkTitle
    ! output
    integer :: NCells,NNodes
    ! local
    character(256)   :: sbuf
    integer :: meshSym
    integer :: ipos


    ! vtk file version information:
!    ipos = index(vtkTitle,"header vers.:")
!    sbuf = vtkTitle(ipos:)
!    ipos = index(sbuf,":")
!    sbuf = sbuf(ipos+1:)
!    ipos = index(sbuf,";")
!    sbuf = sbuf(1:ipos-1)
!    read(sbuf,*) ivtkVers

    ! get number of cells
    ipos = index(vtkTitle,"NC=")
    sbuf = vtkTitle(ipos:)
    ipos = index(sbuf,"=")
    sbuf = sbuf(ipos+1:)
    ipos = index(sbuf,";")
    sbuf = sbuf(1:ipos-1)
    read(sbuf,*) NCells

    ! get number of nodes
    ipos = index(vtkTitle,"NN=")
    sbuf = vtkTitle(ipos:)
    ipos = index(sbuf,"=")
    sbuf = sbuf(ipos+1:)
    ipos = index(sbuf,";")
    sbuf = sbuf(1:ipos-1)
    read(sbuf,*) NNodes

    ! get mesh file
!    ipos = index(vtkTitle,"mesh:")
!    sbuf = vtkTitle(ipos:)
!    ipos = index(sbuf,":")
!    sbuf = sbuf(ipos+1:)
!    ipos = index(sbuf,";")
!    sbuf = sbuf(1:ipos-1)
!    read(sbuf,*) fmesh
!
!    ! get mesh symmetry
!    ipos = index(vtkTitle,"symm:")
!    sbuf = vtkTitle(ipos:)
!    ipos = index(sbuf,":")
!    sbuf = sbuf(ipos+1:)
!    ipos = index(sbuf,";")
!    sbuf = sbuf(1:ipos-1)
!    sbuf = adjustl(sbuf)
!    selectcase(trim(sbuf))
    meshSym = 0
!    case("native")
!        meshSym = 0
!    case("planar")
!        meshSym = 1
!    case("cylindrical")
!        meshSym = 2
!    case("spherical")
!        meshSym = 3
!    endselect

    return
end subroutine vtk_title_parse

! ------------------------------------------------
subroutine vtk_format_get(fname,vtkFmt)
    use module_Settings
    use module_utils
    IMPLICIT NONE
    ! input
    character(*) :: fname                ! vtk filename
    ! output
    character(8) :: vtkFmt               ! vtk filename
    ! local
    integer          :: ivtkUnit
    character(256)   :: sbuf

    ! vtk file version information:
    ivtkUnit = GetFreeUnit()
    open(ivtkUnit,file=fname)
    read(ivtkUnit,*) sbuf
    read(ivtkUnit,*) sbuf
    read(ivtkUnit,*) vtkFmt
    close(ivtkUnit)

    return
end subroutine vtk_format_get

! -------------------------------------------------------------
subroutine WriteRadFieldsVTK(cUden,Wflux,locTit,CurTime,fileTitle,vtkTitle,fileNameOut)
    ! version 4
    USE module_Mesh             ! KV (including additional grid data)
    USE module_GridFunct
    USE LIB_VTK_IO
    use module_Settings
    use module_utils
    use module_TxtRead

    IMPLICIT NONE
    ! ----- input
    type(TGridNodalVector) :: Wflux
    type(TGridNodalScalar) :: cUden
    character(32)          :: locTit
    real(8)                :: CurTime
    character(32)          :: fileTitle
    character(256)         :: vtkTitle
    ! ----- output
    character(128) :: fileNameOut
    ! ----- local
    type(TMeshData), pointer           :: mesh
    REAL(8), ALLOCATABLE, DIMENSION(:) :: var
    REAL(8), ALLOCATABLE, DIMENSION(:) :: varX,varY,varZ
    REAL(4), ALLOCATABLE, DIMENSION(:) :: X,Y,Z
    INTEGER, ALLOCATABLE, DIMENSION(:) :: connect
    INTEGER, ALLOCATABLE, DIMENSION(:) :: cell_type
    character(32)          :: ftit
    character(32)          :: locType
    integer :: E_IO
    integer :: i,j,k,ncon
    integer :: idx
    real(8) :: dtime
    character(256) :: mesh_topology
    integer        :: locPos,locIdx
    integer        :: NNodes,NCells
    integer        :: inode,ielem


    mesh => cUden%pmesh
    locPos  = valpos(mesh%Locations(:)%locTit,locTit)
    locIdx  = mesh%Locations(locPos)%locIdx
    NNodes  = mesh%Locations(locPos)%nnodes
    NCells  = mesh%Locations(locPos)%NElems
    locType = mesh%Locations(locPos)%locType


    ! set filename
    idx = floor(CurTime*1.0e6)
    fileNameOut = trim(OutputDir)//"/vtk/"//trim(fileTitle)//"_"
    !fileNameOut = trim(fileTitle)//"_"
    write(fileNameOut,'(a,I0.9,".vtk")') trim(fileNameOut),idx
    call slash_bck2frw(fileNameOut)


    ! ------------------------- set title and topology ----------------------------

    ! allocation
    allocate(var(NNodes))
    allocate(varX(NNodes),varY(NNodes),varZ(NNodes))
    allocate(X(NNodes),Y(NNodes),Z(NNodes))
    allocate(connect(10*NCells))
    allocate(cell_type(NCells))


    mesh_topology = 'UNSTRUCTURED_GRID'
!    E_IO = VTK_INI(vtkFmtOut,fileNameOut,vtkTitle,mesh_topology)
    E_IO = VTK_INI('ascii',fileNameOut,vtkTitle,mesh_topology)

!    dtime = CurTime*1.0d+6
    dtime = CurTime*1.0d+0
    E_IO = VTK_TimeFieldData(dtime)


    ! ------------------------- SET GRID ------------------------------
    ! NODES
    do i=1,NNodes
        inode = mesh%Locations(locPos)%inode(i)
        X(i) = sngl(mesh%Node(inode)%crd(1))
        Y(i) = sngl(mesh%Node(inode)%crd(2))
        Z(i) = sngl(mesh%Node(inode)%crd(3))
    enddo
    E_IO = VTK_GEO(NNodes,X,Y,Z)
    ! connections
    j = 0
    do i = 1,NCells
        j = j + 1      ! connect(j) = number of nodes
        if(trim(locType).eq."bound") then
            ielem = mesh%Locations(locPos)%ielem(i)
            ncon  = mesh%Face(ielem)%nNodes
            connect(j) = mesh%Face(ielem)%nNodes
            do k = 1,ncon
                j = j + 1
                idx = mesh%Face(ielem)%inode(k)
                inode = valpos(mesh%Locations(locPos)%inode(:),idx)

                connect(j) = inode-1
            enddo
        elseif(trim(locType).eq."volume") then
            ielem = mesh%Locations(locPos)%ielem(i)
            ncon  = mesh%Face(ielem)%nNodes
            connect(j) = mesh%Cell(ielem)%nNodes
            do k = 1,ncon
                j = j + 1
                idx = mesh%Cell(ielem)%inode(k)
                inode = valpos(mesh%Locations(locPos)%inode(:),idx)
                connect(j) = inode-1
            enddo
        end if
    enddo
    ncon = j

    ! cell type
    do i = 1,NCells
            k = 0
        if(trim(locType).eq."bound") then
            ielem = mesh%Locations(locPos)%ielem(i)
            k     = mesh%Face(ielem)%nNodes                ! 3/4-nodes -> tri/quad = (3/4 position in meshElemsTypesVtk)
        elseif(trim(locType).eq."volume") then
            ielem = mesh%Locations(locPos)%ielem(i)
            j     = mesh%Cell(ielem)%itype                 ! get cell type (indexation as in program)
            k     = valpos(meshElemIdxMesh%idx,j)           ! get cell type position in array
        else
            write(*,*) "error: data_io -> WriteRadFieldsVTK: unknown locType"
            stop
        endif
        cell_type(i) = meshElemIdxVtk(k)%idx     ! set cell type (vtk indexation)
    enddo

    E_IO = VTK_CON(NCells,connect(1:ncon),cell_type)


    ! ------------ DATA HEADER (location: cell, node) ----------------------------
    E_IO = VTK_DAT(NNodes,'node')


    ! ------------------------- data arrays ----------------------------
    ! ------- cUden ----------
    do i = 1,NNodes
        inode = mesh%Locations(locPos)%inode(i)
        var(i) = cUden%fval(inode)
    end do
    ftit = trim(cUden%ftit)
    E_IO = VTK_VAR(NNodes,trim(ftit),var)

    ! ------- Wflux ----------
    do i = 1,NNodes
        inode = mesh%Locations(locPos)%inode(i)
        varX(i) = Wflux%fvec(inode)%cmp(1)
        varY(i) = Wflux%fvec(inode)%cmp(2)
        varZ(i) = Wflux%fvec(inode)%cmp(3)
    end do
    ftit = trim(Wflux%ftit)
    E_IO = VTK_VAR('vect',NNodes,trim(ftit),varX,varY,varZ)

    E_IO = VTK_END()

    deallocate(var)
    deallocate(varX,varY,varZ)
    deallocate(X,Y,Z)
    deallocate(connect)
    deallocate(cell_type)

    return
end subroutine WriteRadFieldsVTK


! -------------------------------------------------------------
subroutine WriteStepVTK(gfun,nfun,CurTime,fileTitle,vtkTitle,fileNameOut)
    ! version 4
    USE module_Mesh             ! KV (including additional grid data)
    USE module_GridFunct
    USE LIB_VTK_IO
    use module_utils
    use module_TxtRead

    IMPLICIT NONE
    type(TMeshData), pointer :: mesh
    integer                  :: nfun
    type container
        type(TGridNodalScalar), pointer :: f
    end type container
    type(container), dimension(:) :: gfun(*)
    character(32)   :: ftit

    REAL(8), ALLOCATABLE, DIMENSION(:) :: var
!    REAL(8), ALLOCATABLE, DIMENSION(:) :: varX,varY,varZ
    REAL(4), ALLOCATABLE, DIMENSION(:) :: X,Y,Z
    INTEGER, ALLOCATABLE, DIMENSION(:) :: connect
    INTEGER, ALLOCATABLE, DIMENSION(:) :: cell_type
    integer :: E_IO
    integer :: i,j,k,ncon
    integer :: idx
    real(8) :: CurTime,dtime
    character(128) :: fileNameOut
    character(256) :: mesh_topology
    character(256) :: vtkTitle
    character(32)  :: fileTitle

    idx = floor(CurTime*1.0e6)
    fileNameOut = trim(OutputDir)//"/vtk/"//trim(fileTitle)//"_"
!    fileNameOut = trim(OutputDir)//"/"//trim(fileTitle)//"_"
!    fileNameOut = trim(fileTitle)//"_"
    write(fileNameOut,'(a,I0.9,".vtk")') trim(fileNameOut),idx
    call slash_bck2frw(fileNameOut)
!    write(*,*) trim(fileNameOut)

    ! ------------------------- set title and topology ----------------------------

    !nfun = size(gfun)
!    write(*,*) 'nfun = ',nfun
!    write(*,*) 'NNodes = ',gfun(1)%f%pmesh%NNodes
    mesh => gfun(1)%f%pmesh
    ! allocation
    allocate(var(mesh%NNodes))
    allocate(X(mesh%NNodes),Y(mesh%NNodes),Z(mesh%NNodes))
    allocate(connect(10*mesh%NCells))
    allocate(cell_type(mesh%NCells))


    mesh_topology = 'UNSTRUCTURED_GRID'
    E_IO = VTK_INI(vtkFmtOut,fileNameOut,vtkTitle,mesh_topology)

!    dtime = CurTime*1.0d+6
    dtime = CurTime*1.0d+0
    E_IO = VTK_TimeFieldData(dtime)


    ! ------------------------- SET GRID ------------------------------
    ! NODES
    do i=1,mesh%NNodes
        X(i) = sngl(mesh%Node(i)%crd(1))
        Y(i) = sngl(mesh%Node(i)%crd(2))
        Z(i) = sngl(mesh%Node(i)%crd(3))
    enddo
    E_IO = VTK_GEO(mesh%NNodes,X,Y,Z)
    ! connections
    j = 0
    do i = 1,mesh%NCells
        j = j + 1      ! connect(j) = number of nodes
        connect(j) = mesh%Cell(i)%nNodes
        do k = 1,mesh%Cell(i)%nNodes
            j = j + 1
            connect(j) = mesh%Cell(i)%inode(k)-1
        enddo
    enddo
    ncon = j

    ! cell type
    do i = 1,mesh%NCells
        j = mesh%Cell(i)%itype                 ! get cell type (indexation as in program)
        k = valpos(meshElemIdxMesh%idx,j)    ! get cell type position in array
        cell_type(i) = meshElemIdxVtk(k)%idx     ! set cell type (vtk indexation)
    enddo

    E_IO = VTK_CON(mesh%NCells,connect(1:ncon),cell_type)


    ! ------------ DATA HEADER (location: cell, node) ----------------------------
    E_IO = VTK_DAT(mesh%NNodes,'node')


    ! ------------------------- data arrays ----------------------------
    ! PCI-1,HRFI,FHZ,FHR,PLT,TEM,VIF,VIR,VIZ,SS,AL
    !nfun = size(gfun)
    do i = 1,nfun
        var(:) = gfun(i)%f%fval(:)
        ftit = trim(gfun(i)%f%ftit)
        E_IO = VTK_VAR(mesh%NNodes,trim(ftit),var)
    enddo

    E_IO = VTK_END()

    deallocate(var)
    deallocate(X,Y,Z)
    deallocate(connect)
    deallocate(cell_type)

    return
end subroutine WriteStepVTK



!  ***********         V T K         *******************************
subroutine kspu2dReadStepVTK(gf,mesh,fname,vtkFmt,gfTit,curTime)
    ! version 4
    USE module_Mesh             ! KV (including additional grid data)
    USE module_GridFunct
    USE LIB_VTK_IO
    use module_utils

    IMPLICIT NONE

    ! input variables
    type(TGridNodalScalar) :: gf
    type(TMeshData)        :: mesh
    character(*)           :: fname
    character(*)           :: gfTit
    character(*)           :: vtkFmt
    real(8)                :: curTime
    ! local variables
    REAL(8), ALLOCATABLE :: var(:),varX(:),varY(:),varZ(:)         ! NCells
    REAL(4), ALLOCATABLE :: X(:),Y(:),Z(:)                         ! NNodes
    INTEGER, ALLOCATABLE :: connect(:)                             ! 5*NCells
    INTEGER, ALLOCATABLE :: cell_type(:)                           ! NCells
    integer :: i,k
    character(128)  :: var_location
    character(32)   :: var_name,var_type,var_prec
    character(1024) :: vtkTitle,meshTopo,s_buffer
    character(50)   :: s_buf
    integer         :: E_IO,E_GRID

    integer         :: vtkVers
    integer         :: NNodes,NCells


    NNodes = mesh%NNodes
    NCells = mesh%NCells

    !write(*,*) ' ------------ VTK_INI_READ ----------- '
    call VTK_INI_READ(vtkFmt,fname,vtkTitle,meshTopo,E_IO)
    ! get NNodes and NCells from vtkTitle
    s_buffer = vtkTitle
    !write(*,*) 'filename: ',trim(fname)
    !write(*,*) 'vtkFmt  : ',trim(vtkFmt)
    !write(*,*) 'vtkTitle: ',trim(vtkTitle)

    ! get vtkVers:
    vtkVers = 0
    i=index(s_buffer,'kspu vtk v.')
    if(i.gt.0) then
        s_buf = s_buffer(i+12:)
        k=index(s_buf,';')
        s_buf = s_buf(:k-1)
        read(s_buf,*) vtkVers
    endif
    !write(*,*) 'vtkVers: ',vtkVers

!    ! get NCells and NNodes from vtkTitle
!        i=index(s_buffer,'NN')
!        s_buf = s_buffer(i:)
!        i=index(s_buf,'=')
!        s_buf = s_buf(i+1:)
!        i=index(s_buf,',')
!        if(i.gt.0) then
!            s_buf = s_buf(:i-1)
!        endif
!        i=index(s_buf,';')
!        if(i.gt.0) then
!            s_buf = s_buf(:i-1)
!        endif
!        read(s_buf,*) vtkNNodes
!        ! ..................................
!        s_buffer = vtkTitle
!        i=index(s_buffer,'NC')
!        s_buf = s_buffer(i:)
!        j=index(s_buf,';')
!        s_buf = s_buf(:j-1)
!        i=index(s_buf,'=')
!        s_buf = s_buf(i+1:)
!        read(s_buf,*) vtkNCells
!        if(vtkNNodes.ne.NNodes) then
!            write(*,*) "error: can't read data from vtk: vtkNNodes != meshNNodes"
!            stop
!        endif

    !write(*,*) ' ------------ VTK_TimeFieldData_READ ----------- '
    call VTK_TimeFieldData_READ(CurTime,E_IO)
    !write(*,*) 'vtkTime: ',CurTime


    allocate(X(NNodes),Y(NNodes),Z(NNodes))
    allocate(var(NNodes),varX(NNodes),varY(NNodes),varZ(NNodes))
    allocate(connect(5*NCells),cell_type(NCells))


    !write(*,*) ' ------------ VTK_GEO_READ ----------- '
    X = 0.e0; Y = 0.e0; Z = 0.e0;
    call VTK_GEO_R4_READ(NNodes,X,Y,Z,E_IO,E_GRID)
!    do i = 1,NNodes
!        GridNodes(i)%crd(1) = dble(X(i))
!        GridNodes(i)%crd(2) = dble(Y(i))
!    enddo


    !write(*,*) ' ------------ VTK_CON_READ ----------- '
    call VTK_CON_READ(NCells,connect,cell_type,E_IO,E_GRID)
!    do i = 1,NCells
!        j = (i-1)*5 + 1
!        do k = 1,4
!            GridCells(i)%node(k) = connect(j+k)+1
!        enddo
!    enddo


    !write(*,*) ' ------------ VTK_DAT_READ ----------- '
    call VTK_DAT_READ(NNodes,var_location,E_IO,E_GRID)


!    !write(*,*) ' ------------ scip unused data ----------- '
!    do i = 1,2
!        call VTK_VAR_SCAL_R8_READ(NNodes,var_name,var,E_IO)
!    enddo
!    call VTK_VAR_VECT_R8_READ(NNodes,var_name,varX,varY,varZ,E_IO)
    ! -------------------------- done sciping -----------------------------

    k = 0
    do while(k.eq.0)
        call VTK_VAR_TYPE_IDENT(NNodes,var_type,var_name,var_prec,E_IO)
        var_type = Upper_Case(var_type)
        var_name = Upper_Case(var_name)
        var_prec = Upper_Case(var_prec)
!        write(*,*)
!        write(*,*) "varType: ",trim(var_type)
!        write(*,*) "varName: ",trim(var_name)
!        write(*,*) "varPrec: ",trim(var_prec)
        select case(trim(Upper_Case(var_type)))
            case("SCALARS")
                call VTK_VAR_SCAL_R8_READ(NNodes,var_name,var,E_IO)
            case("VECTORS")
                call VTK_VAR_VECT_R8_READ(NNodes,var_name,varX,varY,varZ,E_IO)
        endselect
        if(trim(gfTit).eq.trim(var_name)) then
            gf%fval = var
            k = 1
        end if
    end do

    !write(*,*) ' ------------ READ TEM ----------- '
    E_IO = VTK_END()

    deallocate(X,Y,Z)
    deallocate(var,varX,varY,varZ)
    deallocate(connect,cell_type)

    return
end subroutine kspu2dReadStepVTK


!  ***********         V T K         *******************************
subroutine VtkReadStep(fname,NScls,NVecs,ContScl,ContVec,cTime)
    ! version 1
    use module_Settings
    use module_Mesh             ! KV (including additional grid data)
    use module_GridFunct
    use LIB_VTK_IO
    use module_utils

    IMPLICIT NONE

    ! input:
    character(*)           :: fname
    integer                :: NScls,NVecs
    ! input/output:
    type TScalarContainer
        type(TGridNodalScalar), pointer :: scl => null()
    end type TScalarContainer
    type TVectorContainer
        type(TGridNodalVector), pointer :: vec => null()
    end type TVectorContainer
    type(TScalarContainer), dimension(NScls) :: ContScl
    type(TVectorContainer), dimension(NVecs) :: ContVec
    ! local:
    character(8)           :: vtkFmt
    real(8)                :: cTime
    ! local variables
    REAL(8), ALLOCATABLE :: var(:),varX(:),varY(:),varZ(:)         ! NCells
    REAL(4), ALLOCATABLE :: X(:),Y(:),Z(:)                         ! NNodes
    INTEGER, ALLOCATABLE :: connect(:)                             ! 5*NCells
    INTEGER, ALLOCATABLE :: cell_type(:)                           ! NCells
    integer :: j,k
    integer :: iScl,iVec
    character(128)  :: var_location
    character(32)   :: var_name,var_type,var_prec
    character(32)   :: gfTitle
    character(1024) :: vtkTitle,meshTopo
    integer         :: E_IO,E_GRID

    integer         :: meshSym
    integer         :: NNodes,NCells

    !write(*,*) ' ------------ VTK_INI_READ ----------- '
    call vtk_format_get(fname,vtkFmt)
    call VTK_INI_READ(vtkFmt,fname,vtkTitle,meshTopo,E_IO)
    !write(*,*) "vtkTitle: ",trim(vtkTitle)

    call vtk_title_parse(vtkTitle,NNodes,NCells,meshSym)

    !write(*,*) ' ------------ VTK_TimeFieldData_READ ----------- '
    call VTK_TimeFieldData_READ(cTime,E_IO)
    !write(*,*) 'vtkTime: ',cTime


    allocate(X(NNodes),Y(NNodes),Z(NNodes))
    allocate(var(NNodes),varX(NNodes),varY(NNodes),varZ(NNodes))
    allocate(connect(5*NCells),cell_type(NCells))


    !write(*,*) ' ------------ VTK_GEO_READ ----------- '
    X = 0.e0; Y = 0.e0; Z = 0.e0;
    call VTK_GEO_R4_READ(NNodes,X,Y,Z,E_IO,E_GRID)

    !write(*,*) ' ------------ VTK_CON_READ ----------- '
    call VTK_CON_READ(NCells,connect,cell_type,E_IO,E_GRID)



    !write(*,*) ' ------------ VTK_DAT_READ ----------- '
    call VTK_DAT_READ(NNodes,var_location,E_IO,E_GRID)


!    !write(*,*) ' ------------ skip unused data ----------- '

    iScl = 0
    iVec = 0
    k    = 0
    do while(k.lt.NScls+NVecs)
        call VTK_VAR_TYPE_IDENT(NNodes,var_type,var_name,var_prec,E_IO)
        var_type = Upper_Case(var_type)
        var_name = Upper_Case(var_name)
        var_prec = Upper_Case(var_prec)
        select case(trim(Upper_Case(var_type)))
            case("SCALARS")
                call VTK_VAR_SCAL_R8_READ(NNodes,var_name,var,E_IO)
                do iScl = 1,NScls
                    gfTitle = ContScl(iScl)%scl%ftit
                    if(trim(gfTitle).eq.trim(var_name)) then
                        ContScl(iScl)%scl%fval = var
                        k = k + 1
                    end if
                end do
            case("VECTORS")
                call VTK_VAR_VECT_R8_READ(NNodes,var_name,varX,varY,varZ,E_IO)
                do iVec = 1,NVecs
                    gfTitle = ContVec(iVec)%vec%ftit
                    if(trim(gfTitle).eq.trim(var_name)) then
                        do j = 1,NNodes
                            ContVec(iVec)%vec%fvec(j)%cmp(1) = varX(j)
                            ContVec(iVec)%vec%fvec(j)%cmp(2) = varY(j)
                            ContVec(iVec)%vec%fvec(j)%cmp(3) = varZ(j)
                        enddo
                        k = k + 1
                    end if
                end do
        endselect
    end do

    !write(*,*) ' ------------ READ TEM ----------- '
    E_IO = VTK_END()

    deallocate(X,Y,Z)
    deallocate(var,varX,varY,varZ)
    deallocate(connect,cell_type)

    return
end subroutine VtkReadStep


! --------------------------------------------------------------
! -------------------------------------------------------------
subroutine VtkWriteStep(fname,NScls,NVecs,ContScl,ContVec,cTime)
    ! version 1
    use module_Mesh             ! KV (including additional grid data)
    use module_GridFunct
    use LIB_VTK_IO
    use module_utils
    use module_TxtRead
    implicit none
    ! input:
    character(*)           :: fname
    integer                :: NScls,NVecs
    ! input/output:
    type TScalarContainer
        type(TGridNodalScalar), pointer :: scl => null()
    end type TScalarContainer
    type TVectorContainer
        type(TGridNodalVector), pointer :: vec => null()
    end type TVectorContainer
    type(TScalarContainer), dimension(NScls) :: ContScl
    type(TVectorContainer), dimension(NVecs) :: ContVec

    ! local
    type(TMeshData),pointer            :: pmesh
    real(8), allocatable, dimension(:) :: var
    real(8), allocatable, dimension(:) :: vecX,vecY,vecZ
    real(4), allocatable, dimension(:) :: X,Y,Z
    integer, allocatable, dimension(:) :: connect
    integer, allocatable, dimension(:) :: cell_type
    integer :: E_IO
    integer :: i,j,k,ncon
    integer :: idx
    real(8) :: cTime,dtime
    character(128) :: fileNameOut
    character(128) :: fmesh
    character(256) :: mesh_topology
    character(256) :: vtkTitle
    character(32)  :: ftit                 ! function title
    logical        :: dirExists

    idx = floor(cTime*1.0e6)
    fileNameOut = trim(OutputDir)//"/vtk"
    call slash_correct(fileNameOut,OSType)
    inquire( file=trim(fileNameOut)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
    if(.not.dirExists) call system('mkdir '//trim(fileNameOut))

    fileNameOut = trim(OutputDir)//"/vtk/"//trim(fname)//"_"
    write(fileNameOut,'(a,I0.9,".vtk")') trim(fileNameOut),idx
    call slash_correct(fileNameOut,OSType)
!    write(*,*) trim(fileNameOut)

    ! ------------------------- set title and topology ----------------------------

    if(NScls.gt.0) then
        pmesh => ContScl(1)%scl%pmesh
    elseif(NVecs.gt.0) then
        pmesh => ContVec(1)%vec%pmesh
    else
        write(*,*) " VtkWriteError: no fnctions to write"
        stop
    end if

    ! allocation
    allocate(var(pmesh%NNodes))
    allocate(X(pmesh%NNodes),Y(pmesh%NNodes),Z(pmesh%NNodes))
    allocate(vecX(pmesh%NNodes),vecY(pmesh%NNodes),vecZ(pmesh%NNodes))
    allocate(connect(10*pmesh%NCells))
    allocate(cell_type(pmesh%NCells))


    mesh_topology = 'UNSTRUCTURED_GRID'
    fmesh = trim(pmesh%meshTit)

    call vtk_title_set(pmesh,vtkTitle)
    E_IO = VTK_INI(vtkFmtOut,fileNameOut,vtkTitle,mesh_topology)

!    dtime = CurTime*1.0d+6
    dtime = cTime*1.0d+0
    E_IO = VTK_TimeFieldData(dtime)

    ! ------------------------- SET GRID ------------------------------
    ! NODES
    do i=1,pmesh%NNodes
        X(i) = sngl(pmesh%Node(i)%crd(1))
        Y(i) = sngl(pmesh%Node(i)%crd(2))
        Z(i) = sngl(pmesh%Node(i)%crd(3))
    enddo
    E_IO = VTK_GEO(pmesh%NNodes,X,Y,Z)
    ! connections
    j = 0
    do i = 1,pmesh%NCells
        j = j + 1      ! connect(j) = number of nodes
        connect(j) = pmesh%Cell(i)%nNodes
        do k = 1,pmesh%Cell(i)%nNodes
            j = j + 1
            connect(j) = pmesh%Cell(i)%inode(k)-1
        enddo
    enddo
    ncon = j

    ! cell type
    do i = 1,pmesh%NCells
        j = pmesh%Cell(i)%itype                  ! get cell type (indexation as in program)
        k = valpos(meshElemIdxMesh%idx,j)        ! get cell type position in array
        cell_type(i) = meshElemIdxVtk(k)%idx     ! set cell type (vtk indexation)
    enddo

    E_IO = VTK_CON(pmesh%NCells,connect(1:ncon),cell_type)

    ! ------------ DATA HEADER (location: cell, node) ----------------------------
    E_IO = VTK_DAT(pmesh%NNodes,'node')

    ! ------------------------- data arrays ----------------------------
    ! PCI-1,HRFI,FHZ,FHR,PLT,TEM,VIF,VIR,VIZ,SS,AL

    do i = 1,NScls
        var(:) = ContScl(i)%scl%fval(:)
        ftit = trim(ContScl(i)%scl%ftit)
        E_IO = VTK_VAR(pmesh%NNodes,trim(ftit),var)
    enddo

    do i = 1,NVecs
        do j = 1,ContVec(i)%vec%Nval
            vecX(j) = ContVec(i)%vec%fvec(j)%cmp(1)
            vecY(j) = ContVec(i)%vec%fvec(j)%cmp(2)
            vecZ(j) = ContVec(i)%vec%fvec(j)%cmp(3)
        enddo
        ftit = trim(ContVec(i)%vec%ftit)
        E_IO = VTK_VAR('vect',pmesh%NNodes,trim(ftit),vecX,vecY,vecZ)
    enddo

    E_IO = VTK_END()

    deallocate(var,vecX,vecY,vecZ)
    deallocate(X,Y,Z)
    deallocate(connect)
    deallocate(cell_type)

    return
end subroutine VtkWriteStep


