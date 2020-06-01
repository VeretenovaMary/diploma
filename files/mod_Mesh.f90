module module_Mesh
    use module_Settings
    use module_Geom

! <!>  - mark places for check
    integer, parameter       :: meshRepUnit = 500
    character(32), parameter :: meshRepName = 'setup_mesh.err'
    integer                  :: meshNumErrors
    integer                  :: meshNumWarnings

    ! parameters for 2D quad mesh
    integer, parameter  :: NumCellNodes = 8  ! number of nodes for mesh cells
    integer, parameter  :: NumCellEdges = 12
    integer, parameter  :: NumCellFaces = 6
    integer, parameter  :: NumCellNeigs = 6
    integer, parameter  :: NumEdgeNodes = 2
    integer, parameter  :: NumFaceNodes = 4
    integer, parameter  :: NumFaceEdges = 4
    integer, parameter  :: NumFaceCells = 2
    integer, parameter  :: NumNodeEdges = 64     ! <!> better 128
    integer, parameter  :: NumNodeFaces = 64     ! <!> better 128
    integer, parameter  :: NumNodeCells = 64     ! <!> better 128
    integer, parameter  :: NumNodeLocs  = 8      ! possible number of node locations
    integer, parameter  :: NumFaceLocs  = 1      ! just boundary location index
    integer, parameter  :: NumCellLocs  = 1      ! just volume   location index
    integer, parameter  :: NumCoords    = 3

    ! elements: point/line/triangle/quad/tetrahedron/hexahedron/prism/pyramid/
    type TMeshElemIndex
        character(32) :: tit
        integer       :: idx
    end type TMeshElemIndex

    ! "mesh module" format elements association
    type(TMeshElemIndex), parameter, dimension(8) :: meshElemIdxMesh = [TMeshElemIndex("points",       1), &
                                                                        TMeshElemIndex("lines",        2), &
                                                                        TMeshElemIndex("triangles",    3), &
                                                                        TMeshElemIndex("quadrangles",  4), &
                                                                        TMeshElemIndex("tetrahedrons", 5), &
                                                                        TMeshElemIndex("hexahedrons",  6), &
                                                                        TMeshElemIndex("prisms",       7), &
                                                                        TMeshElemIndex("pyramids",     8)]
    ! "gmsh" format elements association
    type(TMeshElemIndex), parameter, dimension(8) :: meshElemIdxGmsh = [TMeshElemIndex("points",      15), &
                                                                        TMeshElemIndex("lines",        1), &
                                                                        TMeshElemIndex("triangles",    2), &
                                                                        TMeshElemIndex("quadrangles",  3), &
                                                                        TMeshElemIndex("tetrahedrons", 4), &
                                                                        TMeshElemIndex("hexahedrons",  5), &
                                                                        TMeshElemIndex("prisms",       6), &
                                                                        TMeshElemIndex("pyramids",     7)]
    ! "Vtk" format elements association
    type(TMeshElemIndex), parameter, dimension(8) :: meshElemIdxVtk  = [TMeshElemIndex("points",       1), &
                                                                        TMeshElemIndex("lines",        3), &
                                                                        TMeshElemIndex("triangles",    5), &
                                                                        TMeshElemIndex("quadrangles",  9), &
                                                                        TMeshElemIndex("tetrahedrons",10), &
                                                                        TMeshElemIndex("hexahedrons", 12), &
                                                                        TMeshElemIndex("prisms",      13), &
                                                                        TMeshElemIndex("pyramids",    14)]

    integer, parameter, dimension(8) :: meshElemsNodesNum  = [1,2,3,4,4,8,6,5]
    integer, parameter, dimension(8) :: meshElemsEdgesNum  = [0,1,3,4,6,12,9,8]
    integer, parameter, dimension(8) :: meshElemsFacesNum  = [0,0,1,1,4,6,5,5]
    integer, parameter, dimension(8) :: meshElemsDim       = [0,1,2,2,3,3,3,3]
    integer                          :: NMeshes        = 0                            ! number of loaded meshes


    ! array of association: bndCondTit <-> bndCondIdx
    ! gmsh all physical indexes and names
    character(30), dimension(50) :: physGroupTit   ! gmsh physical names
    integer      , dimension(50) :: physGroupIdx   ! gmsh physical indexes
    integer                      :: physNGroups    ! number of bound  groups

    ! ========= MESH NODES ====================
    type TNodeData
        real(r8p),dimension(NumCoords)   :: crd                 ! node coordinates (x,y,z)
        integer                          :: ncells              ! number of cells adjucent to node
        integer,pointer                  :: icell(:) => null()  ! cells adjucent to node
        !integer,dimension(NumNodeEdges) :: iedge               ! edge adjucent to node
        !integer                         :: nedges              ! number of edges adjucent to node
        integer,dimension(NumNodeLocs)   :: locIds              ! location index (index of physical group in gmsh mesh file)
    end type TNodeData


    ! ========= MESH EDGES ====================
    type TEdgeData
        integer,dimension(NumEdgeNodes)  :: inode       ! list of nodes
        real(r8p)                        :: length      ! edge length
    end type TEdgeData


    ! ========= MESH FACES ====================
    type TFaceData
        integer(i1p)                     :: itype        ! face type (1-point,2-line,3-triangle or 4-quadrangle)
        integer,dimension(NumFaceCells)  :: icell        ! list of cells (for bnd. edge: cell(2) = - Edge_Num = - Bnd_Edge_Num)
!        integer,dimension(NumFaceNodes)  :: inode        ! list of nodes (max nodes = 4 in case of quad)
        integer,pointer                  :: inode(:) => null()  ! list of nodes (max nodes = 4 in case of quad)
        integer                          :: nnodes       ! number of nodes
        integer,pointer                  :: iedge        ! list of edges (max edges = 4 in case of quad)
        integer                          :: nedges       ! number of edges
        integer,dimension(NumFaceEdges)  :: iEdgeDir     ! edge orientation for face (+1 or -1)
        real(r8p)                        :: area         ! face area
        real(r8p),dimension(NumCoords)   :: center       ! face center
        real(r8p),dimension(NumCoords)   :: vnorm        ! normal vector (3d vector)
        integer                          :: locIdx       ! location index (index of physical group in gmsh mesh file)
        integer                          :: bndType      ! boundary conditions type (0 if internal face)
    end type TFaceData

    ! ========= MESH BOUNDARY FACES ====================
    type TBndFaceData
        integer         :: locIdx             ! group index
        integer         :: nfaces             ! number of faces
        integer,pointer :: iface(:) => NULL() ! list of boundary faces for current group
    end type TBndFaceData


    ! ========= MESH CELLS ====================
    type TCellData
        integer                         :: itype               ! cell type
        integer,pointer                 :: iface(:) => null()  ! list of faces
        integer                         :: nfaces              ! number of faces
!        integer,dimension(NumCellNodes) :: inode              ! list of nodes
        integer,pointer                 :: inode(:) => null()  ! list of nodes
        integer                         :: nnodes      ! number of nodes
        integer,dimension(NumCellFaces) :: iFaceDir    ! face normal orientation for cell (+1=out,  -1=in)
        integer,dimension(NumCellNeigs) :: ineighbor   ! list of neighbor cells
        real(r8p),dimension(NumCoords)  :: center      ! coordinates of the cell center
        real(r8p)                       :: volume      ! cell volume
        integer                         :: locIdx      ! location index (index of physical group in gmsh mesh file)
    end type TCellData


    ! ========= MESH locations data ====================
    type TLocationData
        character(32)     :: locTit       ! location title
        integer           :: locIdx       ! location index
        character(8)      :: locType      ! (bound/volume)
        integer           :: nnodes       ! number of nodes for the location
        integer           :: nelems       ! number of elements (cells/faces) for the location
        integer,pointer   :: inode(:) => NULL()     ! list of location nodes for current group, size = nnodes
        integer,pointer   :: ielem(:) => NULL()     ! list of location cells/faces for current group, size = ncells/nfaces
    end type TLocationData


    ! ========= MESH ====================
    type TMeshData
        integer                        :: meshDim             ! = 1,2,3 (1D,2D,3D)
        character(64)                  :: meshTit             ! mesh title (filename)
        integer                        :: meshIdx             ! mesh index
        integer                        :: SymType             ! mesh symmetry (1 - plane, 2 - cylindrical, 3 - spherical)
        real(r8p),dimension(NumCoords) :: SymAxis             ! mesh symmetry (1 - plane, 2 - cylindrical, 3 - spherical)
        integer                        :: NNodes              ! number of mesh nodes
        integer                        :: NEdges              ! number of mesh edges
        integer                        :: NFaces              ! number of mesh faces
        integer                        :: NBFaces             ! number of mesh bound faces
        integer                        :: NCells              ! number of mesh cells
        integer                        :: NLocs               ! number of mesh locations (physical regions)
        type(TNodeData),pointer        :: Node(:) => NULL()  ! NULL doesn't work (?)
        type(TEdgeData),pointer        :: Edge(:) => NULL()  ! NULL doesn't work (?)
        type(TFaceData),pointer        :: Face(:) => NULL()  ! NULL doesn't work (?)
        type(TCellData),pointer        :: Cell(:) => NULL()  ! NULL doesn't work (?)
        type(TBndFaceData),pointer     :: BFaces(:) => NULL()  ! boundary faces information
        ! program physical indexes and names (locations)
        type(TLocationData),pointer    :: Locations(:)  => NULL()  ! NULL doesn't work (?)
        integer                        :: LocNBnds     ! number of bound  groups
        character(32), dimension(50)   :: LocBndTit    ! names of bounds  (groups of bound faces)
        integer, dimension(50)         :: LocBndIdx    ! index of bounds  (groups of bound faces)
        integer                        :: LocNVols     ! number of volume groups
        character(32), dimension(50)   :: LocVolTit    ! names of volumes (groups of cells)
        integer, dimension(50)         :: LocVolIdx    ! index of volumes (groups of cells)
    end type TMeshData

!    ! all finctions in module are private (not usable from outside) except public (see list below)

    ! --------------------------------
    type TFaceAtNode
        integer :: iface                               ! face serial idx
        integer :: imin                                ! index of base (minimal) node index (=inode)
        integer :: imax                                ! index of maximal  node index (for compare)
        integer :: ipre                                ! index of previous node index (for compare)
        integer :: icount                              ! number of incomings
    end type
    type TFaceListAtNode
        type(TFaceAtNode), dimension(NumNodeFaces) :: ListPos
    end type
    ! --------------------------------


    interface valpos
        module procedure valpos_int1
        module procedure valpos_int2
        module procedure valpos_int4
        module procedure valpos_int8
    end interface

    private
    private :: valpos

!    ! list of public subroutines and functions
    public :: new_mesh,close_mesh
    public :: T3DPoint
    public :: meshElemIdxMesh,meshElemIdxGmsh,meshElemIdxVtk
    public :: TMeshData,TNodeData,TEdgeData,TFaceData,TCellData
    public :: TBndFaceData,TLocationData
    public :: NumCellNodes,NumCellEdges,NumCellFaces
    public :: NumFaceNodes,NumFaceEdges
    public :: NMeshes
    public :: mesh_1DMeshCalcul,mesh_2DMeshCalcul,mesh_3DMeshCalcul
    public :: allocate_cell,allocate_face,allocate_node
    public :: allocBFacesList


CONTAINS

    ! - - - - mesh module close - - - - - - - - - - - - - - - - - - - - - - -
    subroutine close_mesh(mesh)
        implicit none
        type(TMeshData) :: mesh
        integer         :: i

        ! deallocate faces
        do i = 1,mesh%NFaces
            deallocate(mesh%Face(i)%inode)
        end do
        deallocate(mesh%Face)

        ! deallocate cells
        do i = 1,mesh%NCells
            deallocate(mesh%Cell(i)%inode)
            deallocate(mesh%Cell(i)%iface)
        end do
        deallocate(mesh%Cell)

        ! deallocate nodes
        do i = 1,mesh%NNodes
           deallocate(mesh%Node(i)%icell)
        enddo
        deallocate(mesh%Node)

        ! deallocate BFaces
        do i = 1,mesh%locNBnds
           deallocate(mesh%BFaces(i)%iface)
        enddo
        deallocate(mesh%BFaces)
        if(associated(mesh%Locations)) then
            do i = 1,mesh%nLocs
               deallocate(mesh%Locations(i)%inode)
               deallocate(mesh%Locations(i)%ielem)
            enddo
            deallocate(mesh%Locations)
        endif
!        deallocate(mesh%Edge)
        return
    end subroutine close_mesh


    ! - - - -  mesh setup (initialization)- - - - - - - - - - - - - - - - - - - - - - -
    subroutine new_mesh(meshFPath,meshSymType,meshSymAxis,mesh)
        ! subroutine do:
        ! 1. open report file
        ! 2. mesh initialization: set mesh data to "zero" (or other initial value)
        ! 3. find mesh file type from name (.msh or .vtk)
        ! 4. read mesh data from file
        ! 5. make report about warnings or errors and close report file

        implicit none
        ! input:
        character(128)          :: meshFPath       ! mesh file path
        integer                 :: meshSymType     ! mesh space symmetry  (1 - plane, 2 - cylindrical, 3 - spherical)
        real(r8p), dimension(3) :: meshSymAxis
        ! output:
        type(TMeshData)         :: mesh            ! mesh object
        ! local:
        integer        :: meshFUnit
        character(4)   :: meshFType
        character(128) :: strRep


        ! 1. - - - open report file
        call mesh_ReportInit(meshFPath)

        ! 2. set variables to initial position (values)
        NMeshes = NMeshes + 1
        mesh%meshTit = mesh_get_ftitle(meshFPath)
        mesh%meshIdx = NMeshes

        meshFUnit = meshRepUnit + 1
        mesh%SymType    = meshSymType
        mesh%SymAxis(:) = meshSymAxis(:)
        mesh%NNodes = 0
        mesh%NEdges = 0
        mesh%NFaces = 0
        mesh%NBFaces = 0
        mesh%NCells  = 0
        mesh%locNBnds = 0
        mesh%locNVols = 0
        mesh%locBndIdx = 0
        mesh%locVolIdx = 0
        mesh%locBndTit(:) = 'none'
        mesh%locVolTit(:) = 'none'

        ! temporal arrays:
        physNGroups     = 0
        physGroupIdx(:) = 0
        physGroupTit(:) = 'none'

        ! report information:
        meshNumErrors   = 0
        meshNumWarnings = 0


        ! 3. find mesh file description (mesh file type: .msh or .vtk?)
        call mesh_getType(meshFPath,meshFType)

        ! 4. reading mesh data from file
        write(*,*) 'load mesh: ',trim(meshFPath)
        if(trim(meshFType).eq.'msh') then
            write(*,*) 'mesh type: gmsh'
            call mesh_ReadGMSH(mesh,meshFPath,meshFUnit)
        else
            write(strRep,'(a)')  'error: unknown mesh format'
            write(*,*) trim(strRep)
            call str2report(strRep)
            stop
        endif

        ! 5. mesh processing (set geometrical data)
        selectcase(mesh%meshDim)
        case(1)
            call mesh_1DMeshCalcul(mesh)
        case(2)
            call mesh_2DMeshCalcul(mesh)
        case(3)
            call mesh_3DMeshCalcul(mesh)
        endselect

        ! 6. finalise report
        write(*,*)
        write(*,'(a,i4,a,i4,a)') 'setup: mesh: ',meshNumErrors,' errors, ',meshNumWarnings,' warnings'
        call mesh_ReportClose()

        return
    end subroutine new_mesh


    ! ---------------------------------------------------------
    function mesh_get_ftitle(fpath) result (ftitle)
        ! input:
        character(*)  :: fpath
        ! output:
        character(32) :: ftitle
        ! local:
        integer :: i
        character(128) :: sbuf

        sbuf = trim(fpath)
        i = 1
        do while(i.ne.0)
            i = index(sbuf,"/")
            sbuf = sbuf(i+1:)
        end do
        i = 1
        do while(i.ne.0)
            i = index(sbuf,"\")
            sbuf = sbuf(i+1:)
        end do
        ftitle = trim(sbuf)
        return

    end function





    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine mesh_1DMeshCalcul(mesh)
        implicit none
        type(TMeshData) :: mesh
        integer :: j,k
        integer :: ipos
        integer :: n1,n2
        integer :: ncells
        integer :: icell,iface,inode
        integer :: locIdx,locPos
        ! prameters:
        real(r8p), parameter :: pi = 3.14159265e0_r8p
        real(8),dimension(NumCoords) :: vec1,vecL
        real(8),dimension(NumCoords) :: pt1
        real(8),dimension(NumCoords) :: fnorm
        real(R8P), parameter :: val0 = 1.e-8_R8P
        real(8) :: eps,elemVol,elemArea
        real(8) :: r1,r2
        ! 1. nodes <- coords correction
        ! 2. nodes <- cells
        ! 3. cells <- centers
        ! 4. faces <- areas
        ! 5. faces <- centers
        ! 6. faces <- normals
        ! 7. cells <- volumes
        ! 8. cells <- normals orientations
        ! 9-... locations data


        ! --- 1. node <- coordinates correction (very small -> 0)
        do k = 1,mesh%NNodes
            do j = 1,3
                eps = mesh%Node(k)%crd(j)
                if(abs(eps).le.val0) then
                    mesh%Node(k)%crd(j) = 0.e0
                end if
            enddo
        enddo


        ! --- 2. nodes <- cells
        ! 2.a - allocate node(inode)%cells
        do inode = 1,mesh%NNodes
            ncells = mesh%Node(inode)%ncells
            call allocate_node(mesh%Node(inode),ncells)
            mesh%Node(inode)%icell(:) = 0
        end do
        ! 2.b - set node(inode)%cells
        do icell = 1,mesh%NCells
            do j = 1,mesh%Cell(icell)%nNodes
                inode = mesh%Cell(icell)%inode(j)
                ipos = valpos(mesh%Node(inode)%icell,0)
                mesh%Node(inode)%icell(ipos) = icell
            end do
        end do


        ! --- 3. cell <- centers
        do icell = 1,mesh%NCells
            pt1 = 0.d0
            do j = 1,mesh%Cell(icell)%nNodes
                inode = mesh%Cell(icell)%inode(j)
                pt1   = pt1 + mesh%Node(inode)%crd
            enddo
            pt1 = pt1/mesh%Cell(icell)%nNodes
            mesh%Cell(icell)%center = pt1
        enddo


        ! --- 4. faces <- areas
        do iface = 1,mesh%NFaces
            elemArea = 0.e0_r8p
            if(mesh%SymType.eq.0) then                ! native symmetry (as is)
                elemArea = 1.e0_r8p
            elseif(mesh%SymType.eq.1) then            ! planar symmetry
                elemArea = 1.e0_r8p
            elseif(mesh%SymType.eq.2) then            ! clindrical symmetry <!> only if axis is along face normal
                elemArea = 1.e0_r8p
            elseif(mesh%SymType.eq.3) then            ! spherical symmetry
                eps = vectLength(mesh%Face(iface)%center)
                elemArea = 4*pi*eps**2
            else                                      ! other symmetry
                write(*,*) " error: unsupported symmetry type for 2D mesh: ",mesh%SymType
                stop
            endif
            mesh%Face(iface)%area = elemArea
        enddo


        ! --- 5. faces <- centers
        do iface = 1,mesh%NFaces
            pt1 = 0.e0
            do j = 1,mesh%Face(iface)%nnodes
                n1 = mesh%Face(iface)%inode(j)
                pt1 = pt1 + mesh%Node(n1)%crd
            enddo
            pt1 = pt1/mesh%Face(iface)%nNodes
            mesh%Face(iface)%center = pt1
        enddo

        ! --- 6. faces <- normals  (<!> only for 1D case)
        do iface = 1,mesh%NFaces
            ! 1st point
            icell = mesh%Face(iface)%icell(1)
            pt1(:) = 0.e0_r8p
            if(icell.gt.0) then
                pt1 = mesh%Cell(icell)%center
            else
                pt1 = mesh%Face(iface)%center
            end if
            ! 2nd point and vector
            icell = mesh%Face(iface)%icell(2)
            vec1(:) = 0.e0_r8p
            if(icell.gt.0) then
                vec1 = mesh%Cell(icell)%center - pt1
            else
                vec1 = mesh%Face(iface)%center - pt1
            end if

            eps = vectLength(vec1)
            fnorm = vec1/eps
            mesh%Face(iface)%vnorm = fnorm
        enddo


        ! --- 7. cell <- volumes       (2D case)
        elemVol = 0.e0_r8p
        do icell = 1,mesh%NCells
            if(mesh%SymType.eq.0) then                ! native symmetry (as is)
                n1 = mesh%Cell(icell)%iface(1)
                n2 = mesh%Cell(icell)%iface(2)
                elemArea = 0.5*(mesh%Face(n1)%area+mesh%Face(n2)%area)
                vecL = 0.5*(mesh%Face(n2)%center-mesh%Face(n1)%center)
                eps  = vectLength(vec1)
                elemVol = eps * elemArea
            elseif(mesh%SymType.eq.1) then            ! planar symmetry
                n1 = mesh%Cell(icell)%iface(1)
                n2 = mesh%Cell(icell)%iface(2)
                elemArea = 0.5*(mesh%Face(n1)%area+mesh%Face(n2)%area)
                vecL = 0.5*(mesh%Face(n2)%center-mesh%Face(n1)%center)
                eps  = vectLength(vec1)
                elemVol = eps * elemArea
            elseif(mesh%SymType.eq.2) then            ! clindrical symmetry <!> only if axis is along face normal
                n1 = mesh%Cell(icell)%iface(1)
                n2 = mesh%Cell(icell)%iface(2)
                elemArea = 0.5*(mesh%Face(n1)%area+mesh%Face(n2)%area)
                vecL = 0.5*(mesh%Face(n2)%center-mesh%Face(n1)%center)
                eps  = vectLength(vec1)
                elemVol = eps * elemArea
            elseif(mesh%SymType.eq.3) then            ! spherical symmetry
                n1 = mesh%Cell(icell)%iface(1)
                n2 = mesh%Cell(icell)%iface(2)
                r1 = vectLength(mesh%Face(n1)%center)
                r2 = vectLength(mesh%Face(n2)%center)
                elemVol = 4.0e0_r8p/3.0e0_r8p*pi*abs(r2**3-r1**3)
            else                                      ! other symmetry
                write(*,*) " error: unsupported symmetry type for 2D mesh: ",mesh%SymType
                stop
            endif
            mesh%Cell(icell)%volume = elemVol
        enddo

        ! --- 8. cell <- face normal orientation
        do icell = 1,mesh%NCells
            do j = 1,mesh%Cell(icell)%nFaces
                iface = mesh%Cell(icell)%iface(j)
                fnorm = mesh%Face(iface)%vnorm
                ! set direction from inside 1st face cell "icells(1)"
                vec1 = mesh%Face(iface)%center - mesh%Cell(icell)%center
                ! scalar production:
                eps = vectScalProd(vec1,fnorm)
                ! check scalar prod. sign and set direction
                if(eps.gt.0.e0) then
                    mesh%Cell(icell)%iFaceDir(j) = +1_i1p
                else
                    mesh%Cell(icell)%iFaceDir(j) = -1_i1p
                endif
            enddo
        enddo


        ! set mesh locations data
        mesh%NLocs = mesh%locNBnds + mesh%locNVols
        write(*,*) " NLocs = ",mesh%NLocs
        allocate(mesh%Locations(mesh%NLocs))
        ! 1. set boundary locations
        mesh%Locations(:)%locIdx = 0
        mesh%Locations(:)%locTit = 'undefined'
        locIdx = 0
        do j = 1,mesh%locNBnds
            locIdx = locIdx + 1
            mesh%Locations(locIdx)%locTit = mesh%locBndTit(j)
            mesh%Locations(locIdx)%locIdx = mesh%locBndIdx(j)
            mesh%Locations(locIdx)%locType = 'bound'
!            write(*,"(a,a32,i3)") "bound: ",trim(mesh%Locations(locIdx)%locTit),mesh%locBndIdx(j)
        enddo
        do j = 1,mesh%locNVols
            locIdx = locIdx + 1
            mesh%Locations(locIdx)%locTit = mesh%locVolTit(j)
            mesh%Locations(locIdx)%locIdx = mesh%locVolIdx(j)
            mesh%Locations(locIdx)%locType = 'volume'
!            write(*,"(a,a)") "volume: ",trim(mesh%Locations(locIdx)%locTit)
        enddo

        ! calculate number of nodes for locations
        mesh%Locations(:)%nnodes = 0
        do j=1,mesh%NNodes
            k = 1
            do while(mesh%Node(j)%locIds(k).gt.0)
                locIdx = mesh%Node(j)%locIds(k)
                locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
                mesh%Locations(locPos)%nnodes = mesh%Locations(locPos)%nnodes + 1
                k = k + 1
            end do
        end do

        ! calculate number of elements for locations
        !   a) bounds
        mesh%Locations(:)%NElems = 0
        do iface = 1,mesh%NFaces
            locIdx = mesh%Face(iface)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
            endif
        enddo
        !   b) volumes
        do icell = 1,mesh%NCells
            locIdx = mesh%Cell(icell)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
            endif
        enddo
        ! allocate arrays of elements
        do j = 1,mesh%NLocs
            call mesh_set_location_size(mesh%Locations(j),mesh%Locations(j)%nnodes,mesh%Locations(j)%NElems)
        enddo

        ! set nodes for locations
        mesh%Locations(:)%nnodes = 0        ! use as node counter
        do j=1,mesh%NNodes
            k = 1
            do while(mesh%Node(j)%locIds(k).gt.0)
                locIdx = mesh%Node(j)%locIds(k)
                locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
                if(locPos.gt.0) then
                    mesh%Locations(locPos)%nnodes = mesh%Locations(locPos)%nnodes + 1
                    n1 = mesh%Locations(locPos)%nnodes
                    mesh%Locations(locPos)%inode(n1) = j
                endif
                k = k + 1
            end do
        end do

        ! set elements for locations
        mesh%Locations(:)%NElems = 0        ! use as element counter
        !   a) bounds (faces)
        do iface = 1,mesh%NFaces
            locIdx = mesh%Face(iface)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
                n1 = mesh%Locations(locPos)%NElems
                mesh%Locations(locPos)%ielem(n1) = iface
            endif
        enddo
        !   b) volumes (cells)
        do icell = 1,mesh%NCells
            locIdx = mesh%Cell(icell)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
                n1 = mesh%Locations(locPos)%NElems
                mesh%Locations(locPos)%ielem(n1) = icell
            endif
        enddo

        return
    end subroutine mesh_1DMeshCalcul



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine mesh_2DMeshCalcul(mesh)
        implicit none
        type(TMeshData) :: mesh
        integer :: j,k
        integer :: ipos
        integer :: n1,n2
        integer :: ncells,nnodes
        integer :: icell,iface,inode
        integer :: itype,locIdx,locPos
        real(8),dimension(NumCoords) :: vec1,vecL
        real(8),dimension(NumCoords) :: pt1,ptL
        real(8),dimension(NumCoords) :: fnorm
        real(R8P), parameter :: val0 = 1.e-8_R8P
        integer, dimension(NumCellNodes) :: elemNodes
        real(8) :: eps,elemVol,elemArea
        real(8) :: r1,r2,z1,z2
        ! 1. nodes <- coords correction
        ! 2. nodes <- cells
        ! 3. cells <- centers
        ! 4. cells <- volumes
        ! 5. faces <- centers
        ! 6. faces <- areas
        ! 7. faces <- normals
        ! 8. cells <- normals orientations
        ! 9-... locations data


        ! --- 1. node <- coordinates correction (very small -> 0)
        do k = 1,mesh%NNodes
            do j = 1,3
                eps = mesh%Node(k)%crd(j)
                if(abs(eps).le.val0) then
                    mesh%Node(k)%crd(j) = 0.e0
                end if
            enddo
        enddo


        ! --- 2. nodes <- cells
        ! 2.a - allocate node(inode)%cells
        mesh%Node(:)%ncells = 0
        do icell = 1,mesh%NCells
            do j = 1,mesh%Cell(icell)%nnodes
                inode = mesh%Cell(icell)%inode(j)
                mesh%Node(inode)%ncells = mesh%Node(inode)%ncells + 1
            end do
        end do

        do inode = 1,mesh%NNodes
            ncells = mesh%Node(inode)%ncells
            call allocate_node(mesh%Node(inode),ncells)
            mesh%Node(inode)%icell(:) = 0
        end do


        ! 2.b - set node(inode)%cells
        do icell = 1,mesh%NCells
            do j = 1,mesh%Cell(icell)%nNodes
                inode = mesh%Cell(icell)%inode(j)
                ipos = valpos(mesh%Node(inode)%icell,0)
                mesh%Node(inode)%icell(ipos) = icell
            end do
        end do


        ! --- 3. cell <- centers
        do icell = 1,mesh%NCells
            pt1 = 0.d0
            do j = 1,mesh%Cell(icell)%nNodes
                inode = mesh%Cell(icell)%inode(j)
                pt1   = pt1 + mesh%Node(inode)%crd
            enddo
            pt1 = pt1/mesh%Cell(icell)%nNodes
            mesh%Cell(icell)%center = pt1
        enddo


        ! --- 4. cell <- volumes       (2D case)
        do icell = 1,mesh%NCells
            itype = mesh%Cell(icell)%itype
            nnodes = mesh%Cell(icell)%nNodes
            elemNodes(1:nnodes) = mesh%Cell(icell)%inode(1:nnodes)
            elemVol = 0.0d0
            if(mesh%SymType.eq.0) then                ! native symmetry (as is)
                call meshGetFaceScalArea(itype,elemNodes,mesh,elemArea)
                elemVol = 1.d0 * elemArea
            elseif(mesh%SymType.eq.1) then            ! planar symmetry
                call meshGetFaceScalArea(itype,elemNodes,mesh,elemArea)
                elemVol = 1.d0 * elemArea
            elseif(mesh%SymType.eq.2) then            ! clindrical symmetry
                ptL(:)  = 0.0d0
                vecL(:) = mesh%SymAxis(:)
                elemVol = 0.0d0
                do j = 1,mesh%Cell(icell)%nnodes-1
                    n1 = mesh%Cell(icell)%inode(j)
                    pt1 = mesh%Node(n1)%crd
                    r1  = point_line_distance(ptL,vecL,pt1)
                    z1  = vectScalProd(vecL,pt1)
                    n2 = mesh%Cell(icell)%inode(j+1)
                    pt1 = mesh%Node(n2)%crd
                    r2  = point_line_distance(ptL,vecL,pt1)
                    z2  = vectScalProd(vecL,pt1)
                    elemVol = elemVol + 1.0d0/6.0d0*(r1**2+r1*r2+r2**2)*(z2-z1)      ! devided by 2*pi
                end do
                elemVol = abs(elemVol)
            else                                      ! other symmetry
                write(*,*) " error: unsupported symmetry type for 2D mesh: ",mesh%SymType
                stop
            endif
            mesh%Cell(icell)%volume = elemVol
        enddo


        ! --- 5. faces <- centers
        do iface = 1,mesh%NFaces
            pt1 = 0.e0
            do j = 1,mesh%Face(iface)%nnodes
                n1 = mesh%Face(iface)%inode(j)
                pt1 = pt1 + mesh%Node(n1)%crd
            enddo
            pt1 = pt1/mesh%Face(iface)%nNodes
            mesh%Face(iface)%center = pt1
        enddo


        ! --- 6. faces <- normals  (<!> only for 2D case, only if mesh plane is (x,y))
        do iface = 1,mesh%NFaces
            n1 = mesh%Face(iface)%inode(1)
            n2 = mesh%Face(iface)%inode(2)
            vec1 = mesh%Node(n2)%crd-mesh%Node(n1)%crd
            fnorm = 0.d0
            fnorm(1) = +vec1(2)/dsqrt(vec1(1)**2+vec1(2)**2)    ! external (right) normal
            fnorm(2) = -vec1(1)/dsqrt(vec1(1)**2+vec1(2)**2)    ! external (right) normal
            mesh%Face(iface)%vnorm = fnorm
        enddo


        ! --- 7. faces <- areas
        do iface = 1,mesh%NFaces
            n1 = mesh%Face(iface)%inode(1)
            n2 = mesh%Face(iface)%inode(2)
            elemArea = 0.0d0
            if(mesh%SymType.eq.0) then                ! native symmetry (as is)
                vec1 = mesh%Node(n2)%crd-mesh%Node(n1)%crd
                elemArea = 1.d0 * vectLength(vec1)
            elseif(mesh%SymType.eq.1) then            ! planar symmetry
                vec1 = mesh%Node(n2)%crd-mesh%Node(n1)%crd
                elemArea = 1.d0 * vectLength(vec1)
            elseif(mesh%SymType.eq.2) then            ! clindrical symmetry
                vec1 = mesh%Node(n2)%crd-mesh%Node(n1)%crd
                ptL(:)  = 0.0d0
                vecL(:) = mesh%SymAxis(:)
                elemVol = 0.0d0
                pt1 = mesh%Node(n1)%crd
                r1  = point_line_distance(ptL,vecL,pt1)
                z1  = vectScalProd(vecL,pt1)
                pt1 = mesh%Node(n2)%crd
                r2  = point_line_distance(ptL,vecL,pt1)
                z2  = vectScalProd(vecL,pt1)
                elemArea = 0.50d0*(r1+r2)*vectLength(vec1)      ! devided by 2*pi
                elemArea = abs(elemArea)
            else                                      ! other symmetry
                write(*,*) " error: unsupported symmetry type for 2D mesh: ",mesh%SymType
                stop
            endif

        enddo


        ! --- set cells <- face normal orintations
        do icell = 1,mesh%NCells
            do j = 1,mesh%Cell(icell)%nFaces
                iface = mesh%Cell(icell)%iface(j)
                fnorm = mesh%Face(iface)%vnorm
                ! set direction from inside 1st face cell "icells(1)"
                vec1 = mesh%Face(iface)%center - mesh%Cell(icell)%center
                ! scalar production:
                eps = vectScalProd(vec1,fnorm)
                ! check scalar prod. sign and set direction
                if(eps.gt.0.e0) then
                    mesh%Cell(icell)%iFaceDir(j) = +1_i1p
                else
                    mesh%Cell(icell)%iFaceDir(j) = -1_i1p
                endif
            enddo
        enddo


        ! set mesh locations data
        mesh%NLocs = mesh%locNBnds + mesh%locNVols
        write(*,*) " NLocs = ",mesh%NLocs
        allocate(mesh%Locations(mesh%NLocs))
        ! 1. set boundary locations
        mesh%Locations(:)%locIdx = 0
        mesh%Locations(:)%locTit = 'undefined'
        locIdx = 0
        do j = 1,mesh%locNBnds
            locIdx = locIdx + 1
            mesh%Locations(locIdx)%locTit = mesh%locBndTit(j)
            mesh%Locations(locIdx)%locIdx = mesh%locBndIdx(j)
            mesh%Locations(locIdx)%locType = 'bound'
!            write(*,"(a,a32,i3)") "bound: ",trim(mesh%Locations(locIdx)%locTit),mesh%locBndIdx(j)
        enddo
        do j = 1,mesh%locNVols
            locIdx = locIdx + 1
            mesh%Locations(locIdx)%locTit = mesh%locVolTit(j)
            mesh%Locations(locIdx)%locIdx = mesh%locVolIdx(j)
            mesh%Locations(locIdx)%locType = 'volume'
!            write(*,"(a,a)") "volume: ",trim(mesh%Locations(locIdx)%locTit)
        enddo


        ! calculate number of nodes for locations
        mesh%Locations(:)%nnodes = 0
        do j=1,mesh%NNodes
            k = 1
            do while(mesh%Node(j)%locIds(k).gt.0)
                locIdx = mesh%Node(j)%locIds(k)
                locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
                mesh%Locations(locPos)%nnodes = mesh%Locations(locPos)%nnodes + 1
                k = k + 1
            end do
        end do


        ! calculate number of elements for locations
        !   a) bounds
        mesh%Locations(:)%NElems = 0
        do iface = 1,mesh%NFaces
            locIdx = mesh%Face(iface)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
            endif
        enddo


        !   b) volumes
        do icell = 1,mesh%NCells
            locIdx = mesh%Cell(icell)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
            endif
        enddo
        ! allocate arrays of elements
        do j = 1,mesh%NLocs
            call mesh_set_location_size(mesh%Locations(j),mesh%Locations(j)%nnodes,mesh%Locations(j)%NElems)
        enddo

        ! set nodes for locations
        mesh%Locations(:)%nnodes = 0        ! use as node counter
        do j=1,mesh%NNodes
            k = 1
            do while(mesh%Node(j)%locIds(k).gt.0)
                locIdx = mesh%Node(j)%locIds(k)
                locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
                if(locPos.gt.0) then
                    mesh%Locations(locPos)%nnodes = mesh%Locations(locPos)%nnodes + 1
                    n1 = mesh%Locations(locPos)%nnodes
                    mesh%Locations(locPos)%inode(n1) = j
                endif
                k = k + 1
            end do
        end do

        ! set elements for locations
        mesh%Locations(:)%NElems = 0        ! use as element counter
        !   a) bounds (faces)
        do iface = 1,mesh%NFaces
            locIdx = mesh%Face(iface)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
                n1 = mesh%Locations(locPos)%NElems
                mesh%Locations(locPos)%ielem(n1) = iface
            endif
        enddo
        !   b) volumes (cells)
        do icell = 1,mesh%NCells
            locIdx = mesh%Cell(icell)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
                n1 = mesh%Locations(locPos)%NElems
                mesh%Locations(locPos)%ielem(n1) = icell
            endif
        enddo

        return
    end subroutine mesh_2DMeshCalcul



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine mesh_3DMeshCalcul(mesh)
        implicit none
        type(TMeshData) :: mesh
        integer :: j,k
        integer :: ipos
        integer :: icell,iface,inode
        integer :: ncells,nnodes,nelems
        integer :: locIdx,locPos
        integer :: itype
        real(8),dimension(NumCoords) :: vec1
        real(8),dimension(NumCoords) :: pt1
        real(8),dimension(NumCoords) :: fnorm
        real(R8P), parameter :: val0 = 1.e-8_r8p
        integer, dimension(NumCellNodes) :: elemNodes
        integer, dimension(NumFaceNodes) :: faceNodes
        real(8) :: eps,elemVol,elemArea
        type(T3DPoint) :: vectArea

        ! 1. nodes <- coords correction
        ! 2. nodes <- cells
        ! 3. cells <- centers
        ! 4. cells <- volumes
        ! 5. faces <- centers
        ! 6. faces <- areas
        ! 7. faces <- normals
        ! 8. cells <- normals orientations
        ! 9-... locations data


        ! --- 1. node <- coordinates correction (very small -> 0)
        do k = 1,mesh%NNodes
            do j = 1,3
                eps = mesh%Node(k)%crd(j)
                if(abs(eps).le.val0) then
                    mesh%Node(k)%crd(j) = 0.e0_r8p
                end if
            enddo
        enddo


        ! --- 2. nodes <- cells
        ! 2.a - allocate node(inode)%cells
        do inode = 1,mesh%NNodes
            ncells = mesh%Node(inode)%ncells
            call allocate_node(mesh%Node(inode),ncells)
            mesh%Node(inode)%icell(:) = 0
        end do
        ! 2.b - set node(inode)%cells
        do icell = 1,mesh%NCells
            do j = 1,mesh%Cell(icell)%nNodes
                inode = mesh%Cell(icell)%inode(j)
                ipos = valpos(mesh%Node(inode)%icell,0)
                mesh%Node(inode)%icell(ipos) = icell
            end do
        end do


        ! --- 3. cell <- centers
        do icell = 1,mesh%NCells
            pt1 = 0.d0
            do j = 1,mesh%Cell(icell)%nNodes
                inode = mesh%Cell(icell)%inode(j)
                pt1   = pt1 + mesh%Node(inode)%crd
            enddo
            pt1 = pt1/mesh%Cell(icell)%nNodes
            mesh%Cell(icell)%center = pt1
        enddo
        ! --- 4. cell <- volumes
        do icell = 1,mesh%NCells
            itype = mesh%Cell(icell)%itype
            nnodes = mesh%Cell(icell)%nNodes
            elemNodes(1:nnodes) = mesh%Cell(icell)%inode(1:nnodes)
            call meshGetCellVolume(itype,elemNodes,mesh,elemVol)
            mesh%Cell(icell)%volume = elemVol
        enddo


        ! --- 5. faces <- centers
        do iface = 1,mesh%NFaces
            ! face <- center
            pt1 = 0.e0
            do j = 1,mesh%Face(iface)%nNodes
                inode = mesh%Face(iface)%inode(j)
                pt1   = pt1 + mesh%Node(inode)%crd
            enddo
            nnodes = mesh%Face(iface)%nNodes
            pt1 = pt1/nnodes
            mesh%Face(iface)%center = pt1
        enddo

        ! --- 6.-7. faces <- normals, areas
        do iface = 1,mesh%NFaces
            ! face <- normal
            itype     = int(mesh%Face(iface)%itype)
            nnodes    = mesh%Face(iface)%nNodes
            faceNodes(1:nnodes) = mesh%Face(iface)%inode(1:nnodes)
            call meshGetFaceVectArea(itype,faceNodes,mesh,vectArea)
            call meshGetFaceScalArea(itype,faceNodes,mesh,elemArea)
            mesh%Face(iface)%area = elemArea
            fnorm = vectArea%crd/vectLength(vectArea)

            ! set direction from inside 1st face cell "icell(1)"
            icell = mesh%Face(iface)%icell(1)
            vec1 = mesh%Face(iface)%center - mesh%Cell(icell)%center
            eps = vectScalProd(vec1,fnorm)

            if(eps.lt.0.e0) then
                nnodes = mesh%Face(iface)%nnodes
                faceNodes(:) = 0
                do j = 1,nnodes
                    faceNodes(nnodes+1-j) = mesh%Face(iface)%inode(j)  ! reverse nodes list
                end do
                mesh%Face(iface)%inode(1:nnodes) = faceNodes(1:nnodes)
                fnorm = -fnorm                                          ! reverse normal
            endif
            mesh%Face(iface)%vnorm = fnorm
        enddo

        ! --- set cells <- face normal orintations
        do icell = 1,mesh%NCells
            do j = 1,mesh%Cell(icell)%nFaces
                iface = mesh%Cell(icell)%iface(j)
                fnorm = mesh%Face(iface)%vnorm
                ! set direction from inside 1st face cell "icells(1)"
                vec1 = mesh%Face(iface)%center - mesh%Cell(icell)%center
                ! scalar production:
                eps = vectScalProd(vec1,fnorm)
                ! check scalar prod. sign and set direction
                if(eps.gt.0.e0) then
                    mesh%Cell(icell)%iFaceDir(j) = +1_i1p
                else
                    mesh%Cell(icell)%iFaceDir(j) = -1_i1p
                endif
            enddo
        enddo


        ! set mesh locations data
        mesh%NLocs = mesh%locNBnds + mesh%locNVols
        write(*,*) " NLocs = ",mesh%NLocs
        allocate(mesh%Locations(mesh%NLocs))
        ! 1. set boundary locations
        mesh%Locations(:)%locIdx = 0
        mesh%Locations(:)%locTit = 'undefined'
        locIdx = 0
        do j = 1,mesh%locNBnds
            locIdx = locIdx + 1
            mesh%Locations(locIdx)%locTit = mesh%locBndTit(j)
            mesh%Locations(locIdx)%locIdx = mesh%locBndIdx(j)
            mesh%Locations(locIdx)%locType = 'bound'
!            write(*,"(a,a32,i3)") "bound: ",trim(mesh%Locations(locIdx)%locTit),mesh%locBndIdx(j)
        enddo
        do j = 1,mesh%locNVols
            locIdx = locIdx + 1
            mesh%Locations(locIdx)%locTit = mesh%locVolTit(j)
            mesh%Locations(locIdx)%locIdx = mesh%locVolIdx(j)
            mesh%Locations(locIdx)%locType = 'volume'
!            write(*,"(a,a)") "volume: ",trim(mesh%Locations(locIdx)%locTit)
        enddo

        ! calculate number of nodes for locations
        mesh%Locations(:)%nnodes = 0
        do j=1,mesh%NNodes
            k = 1
            do while(mesh%Node(j)%locIds(k).gt.0)
                locIdx = mesh%Node(j)%locIds(k)
                locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
                mesh%Locations(locPos)%nnodes = mesh%Locations(locPos)%nnodes + 1
                k = k + 1
            end do
        end do

        ! calculate number of elements for locations
        !   a) bounds
        mesh%Locations(:)%NElems = 0
        do iface = 1,mesh%NFaces
            locIdx = mesh%Face(iface)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
            endif
        enddo
        !   b) volumes
        do icell = 1,mesh%NCells
            locIdx = mesh%Cell(icell)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
            endif
        enddo
        ! allocate arrays of elements
        do j = 1,mesh%NLocs
            call mesh_set_location_size(mesh%Locations(j),mesh%Locations(j)%nnodes,mesh%Locations(j)%NElems)
        enddo

        ! set nodes for locations
        mesh%Locations(:)%nnodes = 0        ! use as node counter
        do j=1,mesh%NNodes
            k = 1
            do while(mesh%Node(j)%locIds(k).gt.0)
                locIdx = mesh%Node(j)%locIds(k)
                locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
                if(locPos.gt.0) then
                    mesh%Locations(locPos)%nnodes = mesh%Locations(locPos)%nnodes + 1
                    nnodes = mesh%Locations(locPos)%nnodes
                    mesh%Locations(locPos)%inode(nnodes) = j
                endif
                k = k + 1
            end do
        end do

        ! set elements for locations
        mesh%Locations(:)%NElems = 0        ! use as element counter
        !   a) bounds (faces)
        do iface = 1,mesh%NFaces
            locIdx = mesh%Face(iface)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
                nelems = mesh%Locations(locPos)%NElems
                mesh%Locations(locPos)%ielem(nelems) = iface
            endif
        enddo
        !   b) volumes (cells)
        do icell = 1,mesh%NCells
            locIdx = mesh%Cell(icell)%locIdx
            locPos = valpos(mesh%Locations(:)%locIdx,locIdx)
            if(locPos.gt.0) then
                mesh%Locations(locPos)%NElems = mesh%Locations(locPos)%NElems + 1
                nelems = mesh%Locations(locPos)%NElems
                mesh%Locations(locPos)%ielem(nelems) = icell
            endif
        enddo

        return
    end subroutine mesh_3DMeshCalcul

    ! allocate node data
    subroutine allocate_node(node,ncells)
        implicit none
        ! input:
        type(TNodeData) :: node
        integer         :: ncells
        ! . . . . . . . . . . . . . .
        allocate(node%icell(ncells))

        return
    end subroutine

    ! allocate face data
    subroutine allocate_face(face,nnodes)
        implicit none
        ! input:
        type(TFaceData) :: face
        integer         :: nnodes
        ! . . . . . . . . . . . . . .
        face%nnodes = nnodes
        allocate(face%inode(nnodes))
!        allocate(face%iedge(nnodes))

        return
    end subroutine

    ! allocate cell data
    subroutine allocate_cell(cell,nnodes,nfaces)
        implicit none
        ! input:
        type(TCellData) :: cell
        integer         :: nnodes
        integer         :: nfaces
        ! . . . . . . . . . . . . . .
        cell%nnodes = nnodes
        cell%nfaces = nfaces
        allocate(cell%inode(nnodes))
        allocate(cell%iface(nfaces))

        return
    end subroutine


    ! set mesh location data sizes
    subroutine mesh_set_location_size(meshLoc,nnodes,NElems)
        implicit none
        type(TLocationData) :: meshLoc
        integer             :: nnodes
        integer             :: NElems

        allocate(meshLoc%inode(nnodes))
        allocate(meshLoc%ielem(NElems))

        return
    end subroutine mesh_set_location_size


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    SUBROUTINE mesh_ReadGMSH(mesh,meshFName,meshFUnit)
        ! subroutine do:
        ! 1. get mesh dimension (read elements from meshFName and check their dimension: 1D/2D/3D) (1x file read)
        ! 2. read mesh physical groups labels (indexes (integer) and names (character(30)))
        ! 3. mesh (1D/2D/3D) data reading                                                          (2x file read)
        ! 4. process mesh data - connections (edges for nodes, orientations, neighbors ...)
        ! 5. mesh geometrical calculations (element lengths, areas, volumess, normals, ...)

        implicit none
        character(128)         :: meshFName
        integer                :: meshFUnit
        integer                :: meshDim
        integer, dimension(33) :: GmshElemsList
        type(TMeshData) :: mesh

        ! - - - - - - - - - - - - - - -


        ! 1. - - - find mesh elements, mesh dimension (1D/2D/3D)
        call gmsh_GetMeshElemsList(meshFName,meshFUnit,GmshElemsList)
        call gmsh_GetMeshDim(GmshElemsList,meshDim)
        mesh%meshDim = meshDim
        write(*,*) 'meshDim = ',meshDim

        ! 2. - - - read physical groups
        call gmsh_GetMeshPhysGroups(meshFUnit,meshFName)

        ! 3. read mesh data (nodes and elements)
        selectcase(meshDim)
        case(1)
            call gmsh_1DMeshRead(mesh,meshFUnit,meshFName,GmshElemsList)
        case(2)
            call gmsh_2DMeshRead(mesh,meshFUnit,meshFName,GmshElemsList)
        case(3)
            call gmsh_3DMeshRead(mesh,meshFUnit,meshFName,GmshElemsList)
        endselect

        return
    end subroutine mesh_ReadGMSH

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine gmsh_GetMeshElemsList(meshFName,meshFUnit,GmshElemsList)
        ! subroutine do:
        ! read elements from meshFName and check their dimension, get max: 1D/2D/3D (1x file read)
        implicit none
        character(128)      :: meshFName
        integer             :: meshFUnit
        integer             :: idx,curType
        integer             :: NElems
        integer             :: i,ios,ipos
        character(128)      :: sbuf
        character(128)      :: strRep
        integer, dimension(33)       :: GmshElemsList

        open(meshFUnit,file=meshFName)
        ! go to elements list
        i   = 0
        ios = 0
        do while(i.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            if(ios.ne.0) then
                write(strRep,'(a)') 'error - mesh: elements list not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif
            i = index(sbuf,'Elements')
        enddo
        read(meshFUnit,*) NElems

        ! find max element type index
        GmshElemsList(:) = 0
        do i = 1,NElems
            read(meshFUnit,*) idx,curType
            ipos = valpos(meshElemIdxGmsh(:)%idx,curType)
            GmshElemsList(ipos) = GmshElemsList(ipos) + 1

        enddo
        close(meshFUnit)

        ! screen output: mesh elements
        write(*,'(1x,a)') 'mesh elements:'
        do i = 1,8
            if(GmshElemsList(i).gt.0) then
                curType = meshElemIdxGmsh(i)%idx
                write(*,'(5x,i7,1x,a32)') GmshElemsList(i),meshElemIdxGmsh(i)%tit
            endif
        enddo

        return
    end subroutine gmsh_GetMeshElemsList

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine gmsh_GetMeshDim(GmshElemsList,meshDim)
        ! subroutine do:
        ! read elements from meshFName and check their dimension, get max: 1D/2D/3D (1x file read)
        implicit none
        integer             :: meshDim
        character(128)      :: strRep
        integer, dimension(33)       :: GmshElemsList

        if(ANY(GmshElemsList(5:8).gt.0)) then
            meshDim = 3
        elseif(ANY(GmshElemsList(3:4).gt.0)) then
            meshDim = 2
        elseif(ANY(GmshElemsList(1:2).gt.0)) then
            meshDim = 1
        else
            write(strRep,'(a,i3,a)') 'error - mesh: found unsupported elements'
            write(*,*) trim(strRep)
            call str2report(strRep)
            meshNumErrors   = meshNumErrors   + 1
        endif

        return
    end subroutine gmsh_GetMeshDim

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine gmsh_GetMeshPhysGroups(meshFUnit,meshFName)
        ! subroutine do:
        ! read mesh physical groups labels
        implicit none

        integer                      :: meshFUnit
        character(128)               :: meshFName
        integer,dimension(128)       :: phGrpIndex
        character(32),dimension(128) :: phGrpTitle
        integer             :: i,j,k,ios
        character(128)      :: sbuf
        character(128)      :: strRep

        open(meshFUnit,file=meshFName)
        ! - - - find physical groups
        i   = 0
        j   = 0
        ios = 0
        do while(i.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            i = index(sbuf,'PhysicalNames')
            j = index(sbuf,'Nodes')
            if(j.gt.0) then
                strRep = 'warning - mesh: physical groups list not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumWarnings = meshNumWarnings + 1
                exit
            endif
        enddo

        physNGroups = 0
        if(i.gt.0) then
            ! read physical groups
            read(meshFUnit,*) physNGroups
            do j = 1,physNGroups
                read(meshFUnit,*) k,phGrpIndex(j),phGrpTitle(j)
                physGroupIdx(j) = phGrpIndex(j)
                physGroupTit(j) = trim(phGrpTitle(j))
            enddo
            j = 0
        endif
        close(meshFUnit)

        return
    end subroutine gmsh_GetMeshPhysGroups

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine gmsh_1DMeshRead(mesh,meshFUnit,meshFName,GmshElemsList)
        ! subroutine do:
        ! 1. run specific (1D/2D/3D) subroutine for meshFName reading (nodes and elements)         (2x file read)

        implicit none

        integer                      :: meshFUnit
        character(128)               :: meshFName
        integer                :: i,j,k,ios,iel
        integer                :: NElems
        character(128)         :: sbuf
        character(128)         :: strRep
        integer                :: idx,iElemType
        integer                :: nElemTags,nElemNodes,nElemFaces,nFaceNodes
        integer, dimension(5)  :: idsTags
        integer                :: nodeIdx,faceIdx,cellIdx
        integer                :: locIdx,locPos
        integer                :: n1,n2,inode
        integer                :: ElemDim
        integer                :: ipos,icount
        integer, dimension(NumCellNodes)              :: idsNodes
        integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
        integer, dimension(NumCellFaces)              :: CellFacesTypes
        integer, dimension(NumFaceNodes)              :: FaceNodesList
        integer, dimension(33) :: GmshElemsList

        type(TFaceListAtNode), pointer :: FaceListAtNode(:) => null()


        type(TMeshData) :: mesh



        open(meshFUnit,file=meshFName)
        ! - - - find nodes list
        j   = 0
        ios = 0
        do while(j.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            j = index(sbuf,'Nodes')
            if(ios.ne.0) then
                strRep = 'error - mesh: nodes list begining not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif
        enddo

        ! ----- Set Mesh Sizes
        read(meshFUnit,*) mesh%NNodes
        mesh%NCells = 0
        do i = 2,2
            mesh%NCells = mesh%NCells + GmshElemsList(i)
        end do
        mesh%NBFaces = 0
        do i = 1,1
            mesh%NBFaces = mesh%NBFaces + GmshElemsList(i)
        end do
        mesh%NFaces = mesh%NBFaces
        do i = 2,2
            mesh%NFaces = mesh%NFaces + GmshElemsList(i)*2
        end do
        mesh%NFaces = mesh%NFaces/2

        allocate(mesh%Node(mesh%NNodes))           ! mesh setup (allocation NNodes)
        allocate(mesh%Cell(mesh%NCells))           ! mesh setup (allocation)
        allocate(mesh%Face(mesh%NFaces))           ! mesh setup (allocation)
!        allocate(mesh%Edge(mesh%NEdges))           ! mesh setup (allocation)


        write(*,*)
        write(*,*) 'mesh%NNodes  = ',mesh%NNodes
        write(*,*) 'mesh%NEdges  = ',mesh%NEdges
        write(*,*) 'mesh%NFaces  = ',mesh%NFaces
        write(*,*) 'mesh%NBFaces = ',mesh%NBFaces
        write(*,*) 'mesh%NCells  = ',mesh%NCells

        ! set values to initial position
        ! nodes
        !mesh%Node(:)%ncells   = 0
        !mesh%Node(:)%nedges   = 0
        do i = 1,mesh%NNodes
            mesh%Node(i)%locIds(:) = 0
        enddo
        ! edges
        do i = 1,mesh%NEdges
            mesh%Edge(i)%inode(:) = 0
        enddo
        ! faces
        do i = 1,mesh%NFaces
!            mesh%Face(i)%iedge(:)  = 0
            mesh%Face(i)%icell(:)  = 0
            mesh%Face(i)%locIdx    = 0
        enddo
        mesh%Face(:)%nEdges   = 0   ! initial value (to count number of edges)
        mesh%Face(:)%nNodes   = 0   ! initial value (to count number of nodes)
        ! cells
        do i = 1,mesh%NCells
            mesh%Cell(i)%locIdx    = 0
        enddo


        ! FaceListAtNode - initialization
        allocate(FaceListAtNode(mesh%NNodes))         ! special ids
        do inode = 1,mesh%NNodes
            FaceListAtNode(inode)%ListPos(:)%icount = 0
            FaceListAtNode(inode)%ListPos(:)%imin  = 0
            FaceListAtNode(inode)%ListPos(:)%imax  = 0
            FaceListAtNode(inode)%ListPos(:)%ipre  = 0
            FaceListAtNode(inode)%ListPos(:)%iface = 0
        enddo

        ! ----- Read Nodes coordinates
        do i = 1,mesh%NNodes
            read(meshFUnit,*) j,mesh%Node(i)%crd
        enddo

        ! - - - find elements list
        i   = 0
        ios = 0
        do while(i.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            i = index(sbuf,'Elements')
            if(ios.ne.0) then
                write(strRep,*) 'error - mesh: elements list begining not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif
        enddo
        read(meshFUnit,*) NElems

        ! first read of elements (get number of elements, separate bndPhysicalIds/volPhysicalIds)
        faceIdx = 0
        cellIdx = 0
        mesh%Node(:)%ncells = 0
        do iel = 1,NElems
            read(meshFUnit,'(a)') sbuf
            read(sbuf,*) idx,iElemType,nElemTags,idsTags(1:nElemTags)

            ! set element attributes (location idx and name)
            nElemNodes = 0
            j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
            nElemNodes = meshElemsNodesNum(j)
            if(meshElemIdxGmsh(1)%idx.eq.iElemType)            then       ! node in GMSH - skip
                ElemDim = 0
                if(ANY(mesh%locBndIdx(:).eq.idsTags(1))) then
                    continue
                else
                    mesh%locNBnds = mesh%locNBnds + 1
                    j = valpos(physGroupIdx,idsTags(1))
                    if(j.eq.0) then
                        write(strRep,'(a,i7,a)') 'warning - mesh: for physical index',idsTags(1),' physical name not found'
                        meshNumWarnings = meshNumWarnings + 1
                        call str2report(strRep)
                    else
                        mesh%locBndTit(mesh%locNBnds) = trim(physGroupTit(j))
                    endif
                    mesh%locBndIdx(mesh%locNBnds) = idsTags(1)
                endif
            elseif(meshElemIdxGmsh(2)%idx.eq.iElemType)        then       ! 2D case: boundary faces in gmsh
                ElemDim = 1
                if(ANY(mesh%locVolIdx(:).eq.idsTags(1))) then
                    continue
                else
                    mesh%locNVols = mesh%locNVols + 1
                    j = valpos(physGroupIdx,idsTags(1))
                    if(j.eq.0) then
                        write(strRep,'(a,i7,a)') 'warning - mesh: for physical index',idsTags(1),' physical name not found'
                        meshNumWarnings = meshNumWarnings + 1
                        call str2report(strRep)
                    else
                        mesh%locVolTit(mesh%locNVols) = trim(physGroupTit(j))
                    endif
                    mesh%locVolIdx(mesh%locNVols)   = idsTags(1)
                endif
            else
                write(strRep,*) 'error - mesh: unknown element type = ',iElemType
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif


            ! set element geometry data
            read(sbuf,*) idx,iElemType,nElemTags,idsTags(1:nElemTags),idsNodes(1:nElemNodes)
            j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
            nElemNodes = meshElemsNodesNum(j)
            ElemDim = meshElemsDim(j)
            if(ElemDim.eq.0) then                ! boundary faces
                ! check element type and set number of nodes for read from sbuf
                faceIdx = faceIdx + 1
                locIdx = idsTags(1)
                locPos = valpos(mesh%locBndIdx,idsTags(1))
                if(locPos.eq.0) then
                    write(*,*) 'error: location ',idsTags(1),' is undefined'
                    stop
                end if
                if(nElemNodes.eq.1) then
                    mesh%Face(faceIdx)%itype = 1    ! point
                else
                    write(*,*) "error: unknown face type for 1D: nnodes = ",nElemNodes
                end if
                ! set faces <- nodes
                mesh%Face(faceIdx)%nNodes = nElemNodes
                call allocate_face(mesh%Face(faceIdx),nElemNodes)
                mesh%Face(faceIdx)%inode(1:nElemNodes) = idsNodes(1:nElemNodes)
                ! create new face in faces list
                FaceNodesList(:) = 0
                FaceNodesList(1:nElemNodes) = mesh%Face(faceIdx)%inode(1:nElemNodes)
                call FaceIdentify(FaceNodesList,FaceListAtNode,mesh%NNodes,inode,ipos,icount)
                if(icount.eq.1) then          ! new face
                    if(icount.eq.2) then
                        write(*,*) "error in mesh: boundary face doublicating"
                        stop
                    end if
                    FaceListAtNode(inode)%ListPos(ipos)%iface = faceIdx
                end if
                ! set faces <- location (bound index)
                mesh%Face(faceIdx)%locIdx = locIdx       ! new (ordered) physical group index (=/= gmsh index)
                ! set nodes <- location index
                do i = 1,nElemNodes
                    nodeIdx = idsNodes(i)
                    k = valpos(mesh%Node(nodeIdx)%locIds,locIdx)
                    if(k.eq.0) then
                        k = valpos(mesh%Node(nodeIdx)%locIds,0)
                        if(k.eq.0) then
                            write(*,*) "error: no place to set node location index"
                            stop
                        else
                            mesh%Node(nodeIdx)%locIds(k) = locIdx
                        end if
                    endif
                end do
            elseif(ElemDim.eq.1) then         ! cell
                ! check element type and set number of nodes for read from sbuf
                cellIdx = cellIdx + 1
                ipos = valpos(meshElemIdxGmsh%idx,iElemType)
                locPos = valpos(mesh%locVolIdx,idsTags(1))
                locIdx = idsTags(1)
                if(locIdx.eq.0) then
                    write(*,*) 'error: location ',idsTags(1),' is undefined'
                    stop
                end if

                nElemFaces = meshElemsNodesNum(ipos)
                mesh%Cell(cellIdx)%nNodes = nElemNodes
                mesh%Cell(cellIdx)%nFaces = nElemFaces
                mesh%Cell(cellIdx)%itype = meshElemIdxMesh(ipos)%idx
                ! set cells <- nodes, locations
                call allocate_cell(mesh%Cell(cellIdx),nElemNodes,nElemFaces)
                mesh%Cell(cellIdx)%inode(1:nElemNodes) = idsNodes(1:nElemNodes)
                mesh%Cell(cellIdx)%locIdx = locIdx     ! new (ordered) physical group index (=/= gmsh index)

                ! set nodes <- cells, locIdx  (first create edge)
                do i = 1,nElemNodes
                    nodeIdx = idsNodes(i)
                    mesh%Node(nodeIdx)%ncells = mesh%Node(nodeIdx)%ncells + 1
                    k = valpos(mesh%Node(nodeIdx)%locIds,locIdx)
                    if(k.eq.0) then
                        k = valpos(mesh%Node(nodeIdx)%locIds,0)
                        if(k.eq.0) then
                            write(*,*) "error: no place to set node location index"
                            stop
                        else
                            mesh%Node(nodeIdx)%locIds(k) = locIdx
                        end if
                    endif
                enddo

                ! set faces <- nodes, cells    (first create face)
                ! set cells <- faces           (first create face)
                j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
                nElemFaces = meshElemsNodesNum(j)
                call meshGet1DElementFaces(j,idsNodes,CellFacesNodes,CellFacesTypes)
                do j = 1,nElemFaces
                    FaceNodesList(:) = 0
                    FaceNodesList(:) = CellFacesNodes(j,:)
                    nFaceNodes = meshElemsNodesNum(CellFacesTypes(j))
!                    write(*,*) iel,'FaceNodesList: ',FaceNodesList(1:nFaceNodes)
!                    pause
                    call FaceIdentify(FaceNodesList,FaceListAtNode,mesh%NNodes,inode,ipos,icount)
                    if(icount.eq.1) then                            ! it is a new face (not boundary)
                        faceIdx = faceIdx + 1
                        FaceListAtNode(inode)%ListPos(ipos)%iface = faceIdx
                        mesh%Face(faceIdx)%icell(1) = cellIdx
                        ! get current face number of nodes by face type
                        n1 = CellFacesTypes(j)
                        n2 = valpos(meshElemsNodesNum,n1)
                        mesh%Face(faceIdx)%nNodes = n2
                        if(n2.eq.1) then
                            mesh%Face(faceIdx)%itype = 1    ! point
                        else
                            write(*,*) "error: unknown face type: nnodes = ",n2
                        end if

                        call allocate_face(mesh%Face(faceIdx),n2)
                        mesh%Face(faceIdx)%inode(1:n2) = FaceNodesList(1:n2)
                        mesh%Face(faceIdx)%locIdx = 0         ! it is new not boundary face
                    else             ! case: face already exists
                        idx = FaceListAtNode(inode)%ListPos(ipos)%iface
                        mesh%Face(idx)%icell(2) = cellIdx
                    endif
                    idx = FaceListAtNode(inode)%ListPos(ipos)%iface
                    mesh%Cell(cellIdx)%iface(j) = idx
                enddo
            endif
        enddo

        call checkFaceCount(mesh,FaceListAtNode)

        ! set neighbours
        do k = 1,mesh%NFaces
            if(k.le.mesh%NBFaces) then
                mesh%Face(k)%icell(1) = mesh%Face(k)%icell(2)
                mesh%Face(k)%icell(2) = -k
                i = mesh%Face(k)%icell(1)
                if(i.eq.0) then
                        write(*,*) ' iface = ',k
                        idx = mesh%Face(k)%inode(1)
                        write(*,*) " inode = ",idx
                        write(*,"(a,3(1x,1pE12.4e3))") " loc: ",mesh%Node(idx)%crd(:)
                        write(*,*) ' cells = ',mesh%Face(k)%icell
                        write(*,*) ' nodes = ',mesh%Face(k)%inode
                end if
                j = valpos(mesh%Cell(i)%iface,k)
                mesh%Cell(i)%ineighbor(j) = -k
            else
                i = mesh%Face(k)%icell(1)
                j = mesh%Face(k)%icell(2)

                if(i.eq.0) then
                        write(*,*) ' iface = ',k
                        idx = mesh%Face(k)%inode(1)
                        write(*,*) " inode = ",idx
                        write(*,"(a,3(1x,1pE12.4e3))") " loc: ",mesh%Node(idx)%crd(:)
                        write(*,*) ' cells = ',mesh%Face(k)%icell
                        write(*,*) ' nodes = ',mesh%Face(k)%inode
                end if

                idx = valpos(mesh%Cell(i)%iface,k)
                mesh%Cell(i)%ineighbor(idx) = j
                idx = valpos(mesh%Cell(j)%iface,k)
                mesh%Cell(j)%ineighbor(idx) = i
            endif
        enddo


!        open(45,file="2dmesh.dat")
!        do i = 1,mesh%NFaces
!            do j = 1,mesh%Face(i)%nNodes
!                n1 = mesh%Face(i)%inode(j)
!                write(45,"(3(1pE12.4e2))") mesh%Node(n1)%crd
!            end do
!            n1 = mesh%Face(i)%inode(1)
!            write(45,"(3(1pE12.4e2))") mesh%Node(n1)%crd
!            write(45,"(3(1pE12.4e2))")
!            write(45,"(3(1pE12.4e2))")
!        end do
!        close(45)


        !pause

        close(meshFUnit)
        deallocate(FaceListAtNode)


        if(physNGroups.ne.mesh%locNBnds+mesh%locNVols) then
            write(strRep,'(a)') 'warning - mesh: physNGroups =/= locNBnds+locNVols '
            meshNumWarnings = meshNumWarnings + 1
            call str2report(strRep)

            write(*,*) ' warning: in locations indexes and names setup'
        end if

        ! - - - set boundary faces - - - - - - - - - - - - - - - - - - -
        ! 1. calculate number of faces for each group
        allocate(mesh%BFaces(mesh%locNBnds))
        mesh%BFaces(:)%nfaces = 0
        do i = 1,mesh%NFaces
            locIdx = mesh%Face(i)%locIdx
            if(locIdx.ne.0) then
                locPos = valpos(mesh%locBndIdx,locIdx)
                mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
            end if
        end do
        ! 2. allocation of boundary group faces list
        do i = 1,mesh%locNBnds
            call allocBFacesList(mesh%BFaces(i),mesh%BFaces(i)%nfaces)
        end do
        ! 3. set faces indexes
        mesh%BFaces(:)%nfaces = 0
        do i = 1,mesh%NFaces
            locIdx = mesh%Face(i)%locIdx
            if(locIdx.ne.0) then
                locPos = valpos(mesh%locBndIdx,locIdx)
                mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
                k = mesh%BFaces(locPos)%nfaces
                mesh%BFaces(locPos)%iface(k) = i
                mesh%BFaces(locPos)%locIdx = locIdx
            end if
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        write(*,*) "1D_reading - done..."

        return
    end subroutine gmsh_1DMeshRead



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine gmsh_2DMeshRead(mesh,meshFUnit,meshFName,GmshElemsList)
        ! subroutine do:
        ! 1. run specific (1D/2D/3D) subroutine for meshFName reading (nodes and elements)         (2x file read)

        implicit none

        integer                      :: meshFUnit
        character(128)               :: meshFName
        integer                :: i,j,k,ios,iel
        integer                :: NElems
        character(128)         :: sbuf
        character(128)         :: strRep
        integer                :: idx,iElemType
        integer                :: nElemTags,nElemNodes,nElemFaces,nFaceNodes
        integer, dimension(5)  :: idsTags
        integer                :: nodeIdx,faceIdx,cellIdx
        integer                :: locIdx,locPos
        integer                :: n1,n2,inode
        integer                :: ElemDim
        integer                :: ipos,icount
        integer, dimension(NumCellNodes)              :: idsNodes
        integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
        integer, dimension(NumCellFaces)              :: CellFacesTypes
        integer, dimension(NumFaceNodes)              :: FaceNodesList
        integer, dimension(33) :: GmshElemsList

        type(TFaceListAtNode), pointer :: FaceListAtNode(:) => null()


        type(TMeshData) :: mesh



        open(meshFUnit,file=meshFName)
        ! - - - find nodes list
        j   = 0
        ios = 0
        do while(j.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            j = index(sbuf,'Nodes')
            if(ios.ne.0) then
                strRep = 'error - mesh: nodes list begining not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif
        enddo

        ! ----- Set Mesh Sizes
        read(meshFUnit,*) mesh%NNodes
        mesh%NCells = 0
        do i = 3,4
            mesh%NCells = mesh%NCells + GmshElemsList(i)
        end do
        mesh%NBFaces = 0
        do i = 2,2
            mesh%NBFaces = mesh%NBFaces + GmshElemsList(i)
        end do
        mesh%NFaces = mesh%NBFaces
        do i = 3,4
            mesh%NFaces = mesh%NFaces + GmshElemsList(i)*meshElemsEdgesNum(i)
        end do
        mesh%NFaces = mesh%NFaces/2

        allocate(mesh%Node(mesh%NNodes))           ! mesh setup (allocation NNodes)
        allocate(mesh%Cell(mesh%NCells))           ! mesh setup (allocation)
        allocate(mesh%Face(mesh%NFaces))           ! mesh setup (allocation)
!        allocate(mesh%Edge(mesh%NEdges))           ! mesh setup (allocation)


        write(*,*)
        write(*,*) 'mesh%NNodes  = ',mesh%NNodes
        write(*,*) 'mesh%NEdges  = ',mesh%NEdges
        write(*,*) 'mesh%NFaces  = ',mesh%NFaces
        write(*,*) 'mesh%NBFaces = ',mesh%NBFaces
        write(*,*) 'mesh%NCells  = ',mesh%NCells

        ! set values to initial position
        ! nodes
        !mesh%Node(:)%ncells   = 0
        !mesh%Node(:)%nedges   = 0
        do i = 1,mesh%NNodes
            mesh%Node(i)%locIds(:) = 0
        enddo
        ! edges
        do i = 1,mesh%NEdges
            mesh%Edge(i)%inode(:) = 0
        enddo
        ! faces
        do i = 1,mesh%NFaces
!            mesh%Face(i)%iedge(:)  = 0
            mesh%Face(i)%icell(:)  = 0
            mesh%Face(i)%locIdx    = 0
        enddo
        mesh%Face(:)%nEdges   = 0   ! initial value (to count number of edges)
        mesh%Face(:)%nNodes   = 0   ! initial value (to count number of nodes)
        ! cells
        do i = 1,mesh%NCells
            mesh%Cell(i)%locIdx    = 0
        enddo


        ! FaceListAtNode - initialization
        allocate(FaceListAtNode(mesh%NNodes))         ! special ids
        do inode = 1,mesh%NNodes
            FaceListAtNode(inode)%ListPos(:)%icount = 0
            FaceListAtNode(inode)%ListPos(:)%imin  = 0
            FaceListAtNode(inode)%ListPos(:)%imax  = 0
            FaceListAtNode(inode)%ListPos(:)%ipre  = 0
            FaceListAtNode(inode)%ListPos(:)%iface = 0
        enddo

        ! ----- Read Nodes coordinates
        do i = 1,mesh%NNodes
            read(meshFUnit,*) j,mesh%Node(i)%crd
        enddo

        ! - - - find elements list
        i   = 0
        ios = 0
        do while(i.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            i = index(sbuf,'Elements')
            if(ios.ne.0) then
                write(strRep,*) 'error - mesh: elements list begining not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif
        enddo
        read(meshFUnit,*) NElems

        ! first read of elements (get number of elements, separate bndPhysicalIds/volPhysicalIds)
        faceIdx = 0
        cellIdx = 0
        mesh%Node(:)%ncells = 0
        do iel = 1,NElems
            read(meshFUnit,'(a)') sbuf
            read(sbuf,*) idx,iElemType,nElemTags,idsTags(1:nElemTags)

            ! set element attributes (location idx and name)
            nElemNodes = 0
            j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
            nElemNodes = meshElemsNodesNum(j)
            if(meshElemIdxGmsh(1)%idx.eq.iElemType)            then       ! node in GMSH - skip
                ElemDim = 0
            elseif(meshElemIdxGmsh(2)%idx.eq.iElemType)        then       ! 2D case: boundary faces in gmsh
                ElemDim = 1
                if(ANY(mesh%locBndIdx(:).eq.idsTags(1))) then
                    continue
                else
                    mesh%locNBnds = mesh%locNBnds + 1
                    j = valpos(physGroupIdx,idsTags(1))
                    if(j.eq.0) then
                        write(strRep,'(a,i7,a)') 'warning - mesh: for physical index',idsTags(1),' physical name not found'
                        meshNumWarnings = meshNumWarnings + 1
                        call str2report(strRep)
                    else
                        mesh%locBndTit(mesh%locNBnds) = trim(physGroupTit(j))
                    endif
                    mesh%locBndIdx(mesh%locNBnds) = idsTags(1)
                endif
            elseif(ANY(meshElemIdxGmsh(3:4)%idx.eq.iElemType)) then     ! 2D case: cells in gmsh
                ElemDim = 2
                if(ANY(mesh%locVolIdx(:).eq.idsTags(1))) then
                    continue
                else
                    mesh%locNVols = mesh%locNVols + 1
                    j = valpos(physGroupIdx,idsTags(1))
                    if(j.eq.0) then
                        write(strRep,'(a,i7,a)') 'warning - mesh: for physical index',idsTags(1),' physical name not found'
                        meshNumWarnings = meshNumWarnings + 1
                        call str2report(strRep)
                    else
                        mesh%locVolTit(mesh%locNVols) = trim(physGroupTit(j))
                    endif
                    mesh%locVolIdx(mesh%locNVols)   = idsTags(1)
                endif
            else
                write(strRep,*) 'error - mesh: unknown element type = ',iElemType
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif


            ! set element geometry data
            read(sbuf,*) idx,iElemType,nElemTags,idsTags(1:nElemTags),idsNodes(1:nElemNodes)
            j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
            nElemNodes = meshElemsNodesNum(j)
            ElemDim = meshElemsDim(j)
            if(ElemDim.eq.1) then                ! boundary faces
                ! check element type and set number of nodes for read from sbuf
                faceIdx = faceIdx + 1
                locIdx = idsTags(1)
                locPos = valpos(mesh%locBndIdx,idsTags(1))
                if(locPos.eq.0) then
                    write(*,*) 'error: location ',idsTags(1),' is undefined'
                    stop
                end if
                if(nElemNodes.eq.2) then
                    mesh%Face(faceIdx)%itype = 2    ! line
                else
                    write(*,*) "error: unknown face type for 2D: nnodes = ",nElemNodes
                end if
                ! set faces <- nodes
                mesh%Face(faceIdx)%nNodes = nElemNodes
                call allocate_face(mesh%Face(faceIdx),nElemNodes)
                mesh%Face(faceIdx)%inode(1:nElemNodes) = idsNodes(1:nElemNodes)
                ! create new face in faces list
                FaceNodesList(:) = 0
                FaceNodesList(1:nElemNodes) = mesh%Face(faceIdx)%inode(1:nElemNodes)
                call FaceIdentify(FaceNodesList,FaceListAtNode,mesh%NNodes,inode,ipos,icount)
                if(icount.eq.1) then          ! new face
                    if(icount.eq.2) then
                        write(*,*) "error in mesh: boundary face doublicating"
                        stop
                    end if
                    FaceListAtNode(inode)%ListPos(ipos)%iface = faceIdx
                end if
                ! set faces <- location (bound index)
                mesh%Face(faceIdx)%locIdx = locIdx       ! new (ordered) physical group index (=/= gmsh index)
                ! set nodes <- location index
                do i = 1,nElemNodes
                    nodeIdx = idsNodes(i)
                    k = valpos(mesh%Node(nodeIdx)%locIds,locIdx)
                    if(k.eq.0) then
                        k = valpos(mesh%Node(nodeIdx)%locIds,0)
                        if(k.eq.0) then
                            write(*,*) "error: no place to set node location index"
                            stop
                        else
                            mesh%Node(nodeIdx)%locIds(k) = locIdx
                        end if
                    endif
                end do
            elseif(ElemDim.eq.2) then         ! cell
                ! check element type and set number of nodes for read from sbuf
                cellIdx = cellIdx + 1
                ipos = valpos(meshElemIdxGmsh%idx,iElemType)
                locPos = valpos(mesh%locVolIdx,idsTags(1))
                locIdx = idsTags(1)
                if(locIdx.eq.0) then
                    write(*,*) 'error: location ',idsTags(1),' is undefined'
                    stop
                end if

                nElemFaces = meshElemsEdgesNum(ipos)
                mesh%Cell(cellIdx)%nNodes = nElemNodes
                mesh%Cell(cellIdx)%nFaces = nElemFaces
                mesh%Cell(cellIdx)%itype = meshElemIdxMesh(ipos)%idx
                ! set cells <- nodes, locations
                call allocate_cell(mesh%Cell(cellIdx),nElemNodes,nElemFaces)
                mesh%Cell(cellIdx)%inode(1:nElemNodes) = idsNodes(1:nElemNodes)
                mesh%Cell(cellIdx)%locIdx = locIdx     ! new (ordered) physical group index (=/= gmsh index)

                ! set nodes <- cells, locIdx  (first create edge)
                do i = 1,nElemNodes
                    nodeIdx = idsNodes(i)
                    mesh%Node(nodeIdx)%ncells = mesh%Node(nodeIdx)%ncells + 1
                    k = valpos(mesh%Node(nodeIdx)%locIds,locIdx)
                    if(k.eq.0) then
                        k = valpos(mesh%Node(nodeIdx)%locIds,0)
                        if(k.eq.0) then
                            write(*,*) "error: no place to set node location index"
                            stop
                        else
                            mesh%Node(nodeIdx)%locIds(k) = locIdx
                        end if
                    endif
                enddo

                ! set faces <- nodes, cells    (first create face)
                ! set cells <- faces           (first create face)
                j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
                nElemFaces = meshElemsEdgesNum(j)
                call meshGet2DElementFaces(j,idsNodes,CellFacesNodes,CellFacesTypes)
                do j = 1,nElemFaces
                    FaceNodesList(:) = 0
                    FaceNodesList(:) = CellFacesNodes(j,:)
                    nFaceNodes = meshElemsNodesNum(CellFacesTypes(j))
!                    write(*,*) iel,'FaceNodesList: ',FaceNodesList(1:nFaceNodes)
!                    pause
                    call FaceIdentify(FaceNodesList,FaceListAtNode,mesh%NNodes,inode,ipos,icount)
                    if(icount.eq.1) then                            ! it is a new face (not boundary)
                        faceIdx = faceIdx + 1
                        FaceListAtNode(inode)%ListPos(ipos)%iface = faceIdx
                        mesh%Face(faceIdx)%icell(1) = cellIdx
                        ! get current face number of nodes by face type
                        n1 = CellFacesTypes(j)
                        n2 = valpos(meshElemsNodesNum,n1)
                        mesh%Face(faceIdx)%nNodes = n2
                        if(n2.eq.2) then
                            mesh%Face(faceIdx)%itype = 2    ! line
                        else
                            write(*,*) "error: unknown face type: nnodes = ",n2
                        end if

                        call allocate_face(mesh%Face(faceIdx),n2)
                        mesh%Face(faceIdx)%inode(1:n2) = FaceNodesList(1:n2)
                        mesh%Face(faceIdx)%locIdx = 0         ! it is new not boundary face
                    else             ! case: face already exists
                        idx = FaceListAtNode(inode)%ListPos(ipos)%iface
                        mesh%Face(idx)%icell(2) = cellIdx
                    endif
                    idx = FaceListAtNode(inode)%ListPos(ipos)%iface
                    mesh%Cell(cellIdx)%iface(j) = idx
                enddo
            endif
        enddo

        call checkFaceCount(mesh,FaceListAtNode)

        ! set neighbours
        do k = 1,mesh%NFaces
            if(k.le.mesh%NBFaces) then
                mesh%Face(k)%icell(1) = mesh%Face(k)%icell(2)
                mesh%Face(k)%icell(2) = -k
                i = mesh%Face(k)%icell(1)
                if(i.eq.0) then
                        write(*,*) ' iface = ',k
                        idx = mesh%Face(k)%inode(1)
                        write(*,*) " inode = ",idx
                        write(*,"(a,3(1x,1pE12.4e3))") " loc: ",mesh%Node(idx)%crd(:)
                        write(*,*) ' cells = ',mesh%Face(k)%icell
                        write(*,*) ' nodes = ',mesh%Face(k)%inode
                end if
                j = valpos(mesh%Cell(i)%iface,k)
                mesh%Cell(i)%ineighbor(j) = -k
            else
                i = mesh%Face(k)%icell(1)
                j = mesh%Face(k)%icell(2)

                if(i.eq.0) then
                        write(*,*) ' iface = ',k
                        idx = mesh%Face(k)%inode(1)
                        write(*,*) " inode = ",idx
                        write(*,"(a,3(1x,1pE12.4e3))") " loc: ",mesh%Node(idx)%crd(:)
                        write(*,*) ' cells = ',mesh%Face(k)%icell
                        write(*,*) ' nodes = ',mesh%Face(k)%inode
                end if

                idx = valpos(mesh%Cell(i)%iface,k)
                mesh%Cell(i)%ineighbor(idx) = j
                idx = valpos(mesh%Cell(j)%iface,k)
                mesh%Cell(j)%ineighbor(idx) = i
            endif
        enddo


!        open(45,file="2dmesh.dat")
!        do i = 1,mesh%NFaces
!            do j = 1,mesh%Face(i)%nNodes
!                n1 = mesh%Face(i)%inode(j)
!                write(45,"(3(1pE12.4e2))") mesh%Node(n1)%crd
!            end do
!            n1 = mesh%Face(i)%inode(1)
!            write(45,"(3(1pE12.4e2))") mesh%Node(n1)%crd
!            write(45,"(3(1pE12.4e2))")
!            write(45,"(3(1pE12.4e2))")
!        end do
!        close(45)


        !pause

        close(meshFUnit)
        deallocate(FaceListAtNode)


        if(physNGroups.ne.mesh%locNBnds+mesh%locNVols) then
            write(strRep,'(a)') 'warning - mesh: physNGroups =/= locNBnds+locNVols '
            meshNumWarnings = meshNumWarnings + 1
            call str2report(strRep)

            write(*,*) ' warning: in locations indexes and names setup'
        end if

        ! - - - set boundary faces - - - - - - - - - - - - - - - - - - -
        ! 1. calculate number of faces for each group
        allocate(mesh%BFaces(mesh%locNBnds))
        mesh%BFaces(:)%nfaces = 0
        do i = 1,mesh%NFaces
            locIdx = mesh%Face(i)%locIdx
            if(locIdx.ne.0) then
                locPos = valpos(mesh%locBndIdx,locIdx)
                mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
            end if
        end do
        ! 2. allocation of boundary group faces list
        do i = 1,mesh%locNBnds
            call allocBFacesList(mesh%BFaces(i),mesh%BFaces(i)%nfaces)
        end do
        ! 3. set faces indexes
        mesh%BFaces(:)%nfaces = 0
        do i = 1,mesh%NFaces
            locIdx = mesh%Face(i)%locIdx
            if(locIdx.ne.0) then
                locPos = valpos(mesh%locBndIdx,locIdx)
                mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
                k = mesh%BFaces(locPos)%nfaces
                mesh%BFaces(locPos)%iface(k) = i
                mesh%BFaces(locPos)%locIdx = locIdx
            end if
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        write(*,*) "2D_reading - done..."

        return
    end subroutine gmsh_2DMeshRead



    ! - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine gmsh_3DMeshRead(mesh,meshFUnit,meshFName,GmshElemsList)
        ! subroutine do:
        ! 1. run specific (1D/2D/3D) subroutine for meshFName reading (nodes and elements)         (2x file read)

        implicit none

        integer                      :: meshFUnit
        character(128)               :: meshFName
        integer                :: i,j,k,ios,iel
        integer                :: NElems
        character(128)         :: sbuf
        character(128)         :: strRep
        integer                :: idx,iElemType
        integer                :: nElemTags,nElemNodes,nElemFaces,nFaceNodes
        integer, dimension(5)  :: idsTags
        integer                :: nodeIdx,faceIdx,cellIdx
        integer                :: locIdx,locPos
        integer                :: n1,n2,inode
        integer                :: ElemDim
        integer                :: ipos,icount
        integer, dimension(NumCellNodes)              :: idsNodes
        integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
        integer, dimension(NumCellFaces)              :: CellFacesTypes
        integer, dimension(NumFaceNodes)              :: FaceNodesList
        integer, dimension(33) :: GmshElemsList

        type(TFaceListAtNode), pointer :: FaceListAtNode(:) => null()


        type(TMeshData) :: mesh



        open(meshFUnit,file=meshFName)
        ! - - - find nodes list
        j   = 0
        ios = 0
        do while(j.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            j = index(sbuf,'Nodes')
            if(ios.ne.0) then
                strRep = 'error - mesh: nodes list begining not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif
        enddo

        ! ----- Set Mesh Sizes
        read(meshFUnit,*) mesh%NNodes
        mesh%NCells = 0
        do i = 5,8
            mesh%NCells = mesh%NCells + GmshElemsList(i)
        end do
        mesh%NBFaces = 0
        do i = 3,4
            mesh%NBFaces = mesh%NBFaces + GmshElemsList(i)
        end do
        mesh%NFaces = mesh%NBFaces
        do i = 5,8
            mesh%NFaces = mesh%NFaces + GmshElemsList(i)*meshElemsFacesNum(i)
        end do
        mesh%NFaces = mesh%NFaces/2

        allocate(mesh%Node(mesh%NNodes))           ! mesh setup (allocation NNodes)
        allocate(mesh%Cell(mesh%NCells))           ! mesh setup (allocation)
        allocate(mesh%Face(mesh%NFaces))           ! mesh setup (allocation)
!        allocate(mesh%Edge(mesh%NEdges))           ! mesh setup (allocation)


        write(*,*)
        write(*,*) 'mesh%NNodes  = ',mesh%NNodes
        write(*,*) 'mesh%NEdges  = ',mesh%NEdges
        write(*,*) 'mesh%NFaces  = ',mesh%NFaces
        write(*,*) 'mesh%NBFaces = ',mesh%NBFaces
        write(*,*) 'mesh%NCells  = ',mesh%NCells

        ! set values to initial position
        ! nodes
        !mesh%Node(:)%ncells   = 0
        !mesh%Node(:)%nedges   = 0
        do i = 1,mesh%NNodes
            mesh%Node(i)%locIds(:) = 0
        enddo
        ! edges
        do i = 1,mesh%NEdges
            mesh%Edge(i)%inode(:) = 0
        enddo
        ! faces
        do i = 1,mesh%NFaces
!            mesh%Face(i)%iedge(:)  = 0
            mesh%Face(i)%icell(:)  = 0
            mesh%Face(i)%locIdx    = 0
        enddo
        mesh%Face(:)%nEdges   = 0   ! initial value (to count number of edges)
        mesh%Face(:)%nNodes   = 0   ! initial value (to count number of nodes)
        ! cells
        do i = 1,mesh%NCells
            mesh%Cell(i)%locIdx    = 0
        enddo


        ! FaceListAtNode - initialization
        allocate(FaceListAtNode(mesh%NNodes))         ! special ids
        do inode = 1,mesh%NNodes
            FaceListAtNode(inode)%ListPos(:)%icount = 0
            FaceListAtNode(inode)%ListPos(:)%imin  = 0
            FaceListAtNode(inode)%ListPos(:)%imax  = 0
            FaceListAtNode(inode)%ListPos(:)%ipre  = 0
            FaceListAtNode(inode)%ListPos(:)%iface = 0
        enddo

        ! ----- Read Nodes coordinates
        do i = 1,mesh%NNodes
            read(meshFUnit,*) j,mesh%Node(i)%crd
        enddo

        ! - - - find elements list
        i   = 0
        ios = 0
        do while(i.eq.0 .and. ios.eq.0)
            read(meshFUnit,*,iostat=ios) sbuf
            i = index(sbuf,'Elements')
            if(ios.ne.0) then
                write(strRep,*) 'error - mesh: elements list begining not found'
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif
        enddo
        read(meshFUnit,*) NElems

        ! first read of elements (get number of elements, separate bndPhysicalIds/volPhysicalIds)
        faceIdx = 0
        cellIdx = 0
        mesh%Node(:)%ncells = 0
        do iel = 1,NElems
            read(meshFUnit,'(a)') sbuf
            read(sbuf,*) idx,iElemType,nElemTags,idsTags(1:nElemTags)

            ! set element attributes (location idx and name)
            nElemNodes = 0
            j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
            nElemNodes = meshElemsNodesNum(j)
            if(meshElemIdxGmsh(1)%idx.eq.iElemType)            then       ! node in GMSH - skip
                ElemDim = 0
            elseif(meshElemIdxGmsh(2)%idx.eq.iElemType)        then       ! 3D case: line in GMSH - skip
                ElemDim = 1
            elseif(ANY(meshElemIdxGmsh(3:4)%idx.eq.iElemType)) then       ! 3D case: boundary faces in gmsh (tri or quad)
                ElemDim = 2
                if(ANY(mesh%locBndIdx(:).eq.idsTags(1))) then
                    continue
                else
                    mesh%locNBnds = mesh%locNBnds + 1
                    j = valpos(physGroupIdx,idsTags(1))
                    if(j.eq.0) then
                        write(strRep,'(a,i7,a)') 'warning - mesh: for physical index',idsTags(1),' physical name not found'
                        meshNumWarnings = meshNumWarnings + 1
                        call str2report(strRep)
                    else
                        mesh%locBndTit(mesh%locNBnds) = trim(physGroupTit(j))
                    endif
                    mesh%locBndIdx(mesh%locNBnds) = idsTags(1)
                endif
            elseif(ANY(meshElemIdxGmsh(5:8)%idx.eq.iElemType)) then   ! 3D case: cells in GMSH
                ElemDim = 3
                if(ANY(mesh%locVolIdx(:).eq.idsTags(1))) then
                    continue
                else
                    mesh%locNVols = mesh%locNVols + 1
                    j = valpos(physGroupIdx,idsTags(1))
                    if(j.eq.0) then
                        write(strRep,'(a,i7,a)') 'warning - mesh: for physical index',idsTags(1),' physical name not found'
                        meshNumWarnings = meshNumWarnings + 1
                        call str2report(strRep)
                    else
                        mesh%locVolTit(mesh%locNVols) = trim(physGroupTit(j))
                    endif
                    mesh%locVolIdx(mesh%locNVols)   = idsTags(1)
                endif
            else
                write(strRep,*) 'error - mesh: unknown element type = ',iElemType
                write(*,*) trim(strRep)
                call str2report(strRep)
                meshNumErrors   = meshNumErrors   + 1
                stop
            endif


            ! set element geometry data
            read(sbuf,*) idx,iElemType,nElemTags,idsTags(1:nElemTags),idsNodes(1:nElemNodes)
            j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
            nElemNodes = meshElemsNodesNum(j)
            ElemDim = meshElemsDim(j)
            if(ElemDim.eq.2) then                ! boundary faces
                ! check element type and set number of nodes for read from sbuf
                faceIdx = faceIdx + 1
                locIdx = idsTags(1)
                locPos = valpos(mesh%locBndIdx,idsTags(1))
                if(locPos.eq.0) then
                    write(*,*) 'error: location ',idsTags(1),' is undefined'
                    stop
                end if
                if(nElemNodes.eq.3) then
                    mesh%Face(faceIdx)%itype = 3    ! triangle
                elseif(nElemNodes.eq.4) then
                    mesh%Face(faceIdx)%itype = 4    ! quadrangle
                else
                    write(*,*) "error: unknown face type: nnodes = ",nElemNodes
                end if
                ! set faces <- nodes
                mesh%Face(faceIdx)%nNodes = nElemNodes
                call allocate_face(mesh%Face(faceIdx),nElemNodes)
                mesh%Face(faceIdx)%inode(1:nElemNodes) = idsNodes(1:nElemNodes)
                ! create new face in faces list
                FaceNodesList(:) = 0
                FaceNodesList(1:nElemNodes) = mesh%Face(faceIdx)%inode(1:nElemNodes)
                call FaceIdentify(FaceNodesList,FaceListAtNode,mesh%NNodes,inode,ipos,icount)
                if(icount.eq.1) then          ! new face
                    if(icount.eq.2) then
                        write(*,*) "error in mesh: boundary face doublicating"
                        stop
                    end if
                    FaceListAtNode(inode)%ListPos(ipos)%iface = faceIdx
                end if
                ! set faces <- location (bound index)
                mesh%Face(faceIdx)%locIdx = locIdx       ! new (ordered) physical group index (=/= gmsh index)
                ! set nodes <- location index
                do i = 1,nElemNodes
                    nodeIdx = idsNodes(i)
                    k = valpos(mesh%Node(nodeIdx)%locIds,locIdx)
                    if(k.eq.0) then
                        k = valpos(mesh%Node(nodeIdx)%locIds,0)
                        if(k.eq.0) then
                            write(*,*) "error: no place to set node location index"
                            stop
                        else
                            mesh%Node(nodeIdx)%locIds(k) = locIdx
                        end if
                    endif
                end do
            elseif(ElemDim.eq.3) then         ! cell
                ! check element type and set number of nodes for read from sbuf
                cellIdx = cellIdx + 1
                ipos = valpos(meshElemIdxGmsh%idx,iElemType)
                locPos = valpos(mesh%locVolIdx,idsTags(1))
                locIdx = idsTags(1)
                if(locIdx.eq.0) then
                    write(*,*) 'error: location ',idsTags(1),' is undefined'
                    stop
                end if

                nElemFaces = meshElemsFacesNum(ipos)
                mesh%Cell(cellIdx)%nNodes = nElemNodes
                mesh%Cell(cellIdx)%nFaces = nElemFaces
                mesh%Cell(cellIdx)%itype = meshElemIdxMesh(ipos)%idx
                ! set cells <- nodes, locations
                call allocate_cell(mesh%Cell(cellIdx),nElemNodes,nElemFaces)
                mesh%Cell(cellIdx)%inode(1:nElemNodes) = idsNodes(1:nElemNodes)
                mesh%Cell(cellIdx)%locIdx = locIdx     ! new (ordered) physical group index (=/= gmsh index)

                ! set nodes <- cells, locIdx  (first create edge)
                do i = 1,nElemNodes
                    nodeIdx = idsNodes(i)
                    mesh%Node(nodeIdx)%ncells = mesh%Node(nodeIdx)%ncells + 1
                    k = valpos(mesh%Node(nodeIdx)%locIds,locIdx)
                    if(k.eq.0) then
                        k = valpos(mesh%Node(nodeIdx)%locIds,0)
                        if(k.eq.0) then
                            write(*,*) "error: no place to set node location index"
                            stop
                        else
                            mesh%Node(nodeIdx)%locIds(k) = locIdx
                        end if
                    endif
                enddo

                ! set faces <- nodes, cells    (first create face)
                ! set cells <- faces           (first create face)
                j = valpos(meshElemIdxGmsh(:)%idx,iElemType)
                nElemFaces = meshElemsFacesNum(j)
                call meshGet3DElementFaces(j,idsNodes,CellFacesNodes,CellFacesTypes)
                do j = 1,nElemFaces
                    FaceNodesList(:) = 0
                    FaceNodesList(:) = CellFacesNodes(j,:)
                    nFaceNodes = meshElemsNodesNum(CellFacesTypes(j))
!                    write(*,*) iel,'FaceNodesList: ',FaceNodesList(1:nFaceNodes)
!                    pause
                    call FaceIdentify(FaceNodesList,FaceListAtNode,mesh%NNodes,inode,ipos,icount)
                    if(icount.eq.1) then                            ! it is a new face (not boundary)
                        faceIdx = faceIdx + 1
                        FaceListAtNode(inode)%ListPos(ipos)%iface = faceIdx
                        mesh%Face(faceIdx)%icell(1) = cellIdx
                        ! get current face number of nodes by face type
                        n1 = CellFacesTypes(j)
                        n2 = valpos(meshElemsNodesNum,n1)
                        mesh%Face(faceIdx)%nNodes = n2
                        if(n2.eq.3) then
                            mesh%Face(faceIdx)%itype = 3    ! triangle
                        elseif(n2.eq.4) then
                            mesh%Face(faceIdx)%itype = 4    ! quadrangle
                        else
                            write(*,*) "error: unknown face type: nnodes = ",n2
                        end if

                        call allocate_face(mesh%Face(faceIdx),n2)
                        mesh%Face(faceIdx)%inode(1:n2) = FaceNodesList(1:n2)
                        mesh%Face(faceIdx)%locIdx = 0         ! it is new not boundary face
                    else             ! case: face already exists
                        idx = FaceListAtNode(inode)%ListPos(ipos)%iface
                        mesh%Face(idx)%icell(2) = cellIdx
                    endif
                    idx = FaceListAtNode(inode)%ListPos(ipos)%iface
                    mesh%Cell(cellIdx)%iface(j) = idx
                enddo
            endif
        enddo

        call checkFaceCount(mesh,FaceListAtNode)

        ! set neighbours
        do k = 1,mesh%NFaces
            if(k.le.mesh%NBFaces) then
                mesh%Face(k)%icell(1) = mesh%Face(k)%icell(2)
                mesh%Face(k)%icell(2) = -k
                i = mesh%Face(k)%icell(1)
                if(i.eq.0) then
                        write(*,*) ' iface = ',k
                        idx = mesh%Face(k)%inode(1)
                        write(*,*) " inode = ",idx
                        write(*,"(a,3(1x,1pE12.4e3))") " loc: ",mesh%Node(idx)%crd(:)
                        write(*,*) ' cells = ',mesh%Face(k)%icell
                        write(*,*) ' nodes = ',mesh%Face(k)%inode
                end if
                j = valpos(mesh%Cell(i)%iface,k)
                mesh%Cell(i)%ineighbor(j) = -k
            else
                i = mesh%Face(k)%icell(1)
                j = mesh%Face(k)%icell(2)

                if(i.eq.0) then
                        write(*,*) ' iface = ',k
                        idx = mesh%Face(k)%inode(1)
                        write(*,*) " inode = ",idx
                        write(*,"(a,3(1x,1pE12.4e3))") " loc: ",mesh%Node(idx)%crd(:)
                        write(*,*) ' cells = ',mesh%Face(k)%icell
                        write(*,*) ' nodes = ',mesh%Face(k)%inode
                end if

                idx = valpos(mesh%Cell(i)%iface,k)
                mesh%Cell(i)%ineighbor(idx) = j
                idx = valpos(mesh%Cell(j)%iface,k)
                mesh%Cell(j)%ineighbor(idx) = i
            endif
        enddo


!        open(45,file="3dmesh.dat")
!        do i = 1,mesh%NFaces
!            do j = 1,mesh%Face(i)%nNodes
!                n1 = mesh%Face(i)%inode(j)
!                write(45,"(3(1pE12.4e2))") mesh%Node(n1)%crd
!            end do
!            n1 = mesh%Face(i)%inode(1)
!            write(45,"(3(1pE12.4e2))") mesh%Node(n1)%crd
!            write(45,"(3(1pE12.4e2))")
!            write(45,"(3(1pE12.4e2))")
!        end do
!        close(45)


        !pause

        close(meshFUnit)
        deallocate(FaceListAtNode)


        if(physNGroups.ne.mesh%locNBnds+mesh%locNVols) then
            write(strRep,'(a)') 'warning - mesh: physNGroups =/= locNBnds+locNVols '
            meshNumWarnings = meshNumWarnings + 1
            call str2report(strRep)

            write(*,*) ' warning: in locations indexes and names setup'
        end if

        ! - - - set boundary faces - - - - - - - - - - - - - - - - - - -
        ! 1. calculate number of faces for each group
        allocate(mesh%BFaces(mesh%locNBnds))
        mesh%BFaces(:)%nfaces = 0
        do i = 1,mesh%NFaces
            locIdx = mesh%Face(i)%locIdx
            if(locIdx.ne.0) then
                locPos = valpos(mesh%locBndIdx,locIdx)
                mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
            end if
        end do
        ! 2. allocation of boundary group faces list
        do i = 1,mesh%locNBnds
            call allocBFacesList(mesh%BFaces(i),mesh%BFaces(i)%nfaces)
        end do
        ! 3. set faces indexes
        mesh%BFaces(:)%nfaces = 0
        do i = 1,mesh%NFaces
            locIdx = mesh%Face(i)%locIdx
            if(locIdx.ne.0) then
                locPos = valpos(mesh%locBndIdx,locIdx)
                mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
                k = mesh%BFaces(locPos)%nfaces
                mesh%BFaces(locPos)%iface(k) = i
                mesh%BFaces(locPos)%locIdx = locIdx
            end if
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        write(*,*) "3D_reading - done..."

        return
    end subroutine gmsh_3DMeshRead


    ! allocation of boundary group faces list
    subroutine allocBFacesList(BndFaceData,nfaces)
        implicit none
        type(TBndFaceData) :: BndFaceData
        integer            :: nfaces

        allocate(BndFaceData%iface(nfaces))

        return
    end subroutine

    ! - - - - - - - - - - - - - - - - - -
    subroutine mesh_getType(meshFName,meshType)
        implicit none
        character(128) :: meshFName
        character(4)   :: meshType
        integer        :: i

        i = len_trim(meshFName)
        do while(meshFName(i:i).ne.'.' .and. i.gt.0)
            i = i-1
        enddo
        meshType = meshFName(i+1:len_trim(meshFName))

        return
    end subroutine mesh_getType

    ! - - - - - - - - - - - - - - - - - -
    subroutine mesh_ReportInit(meshFName)
        implicit none
        character(128) :: meshFName

        ! delete old file if it is exist
        open(meshRepUnit,file = meshRepName)
        close(meshRepUnit,status='delete')

        ! create new file
        open(meshRepUnit,file = meshRepName)
        write(meshRepUnit,*) 'mesh: ',meshFName
        close(meshRepUnit)

        return
    end subroutine mesh_ReportInit

    ! - - - - - - - - - - - - - - - - - -
    subroutine mesh_ReportClose()
        implicit none

        if(meshNumWarnings.gt.0) then
!            stop
            continue
        endif
        if(meshNumErrors.gt.0) then
            stop
        endif
        ! clean file if no errors/warnings occured
        if(meshNumErrors+meshNumWarnings.eq.0) then
            ! delete file with mesh name
            open(meshRepUnit,file = meshRepName)
            close(meshRepUnit,status='delete')
            ! create new clean file
            open(meshRepUnit,file = meshRepName)
            close(meshRepUnit)
        endif

        return
    end subroutine mesh_ReportClose


    ! - - - - - - - - - - - - - - - - - -
    subroutine str2report(sbuf)
        implicit none
        character(128) :: sbuf

        open(meshRepUnit,file=meshRepName,position='append')
        write(meshRepUnit,*) trim(sbuf)
        close(meshRepUnit)

        return
    end subroutine str2report



    ! ! element nodes list -> element faces list
    subroutine meshGet1DElementFaces(iElemType,ElemNodes,CellFacesNodes,iFacesTypes)
        implicit none
        integer                                       :: iElemType
        integer, dimension(NumCellNodes)              :: ElemNodes
        integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
        integer, dimension(NumCellFaces)              :: iFacesTypes

        CellFacesNodes(:,:) = 0
        iFacesTypes(:)  = 0
        if(iElemType.eq.1)      then           ! point
            continue
        elseif(iElemType.eq.2) then            ! line
            CellFacesNodes(1,1) = ElemNodes(1)
            iFacesTypes(1) = 1
            CellFacesNodes(2,1) = ElemNodes(2)
            iFacesTypes(2) = 1
        elseif(iElemType.eq.3) then            ! triangle
            continue
        elseif(iElemType.eq.4) then            ! quadrangle
            continue
        else
            continue
        end if

        return
    end subroutine meshGet1DElementFaces


!    ! element nodes list -> element edges list
!    subroutine meshGet2DElementEdges(iElemType,ElemNodes,EdgeNodes)
!        implicit none
!        integer                            :: iElemType
!        integer, dimension(NumCellNodes)   :: ElemNodes
!        integer, dimension(NumCellEdges,2) :: EdgeNodes
!
!        EdgeNodes(:,:) = 0
!        if(iElemType.eq.1)      then           ! point
!            continue
!        elseif(iElemType.eq.2) then           ! line
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!        elseif(iElemType.eq.3) then           ! triangle
!            EdgeNodes(1,1) = ElemNodes(1)                  !  1
!            EdgeNodes(1,2) = ElemNodes(2)                  !  |\
!            EdgeNodes(2,1) = ElemNodes(2)                  !  | \
!            EdgeNodes(2,2) = ElemNodes(3)                  !  |  \
!            EdgeNodes(3,1) = ElemNodes(3)                  !  2---3
!            EdgeNodes(3,2) = ElemNodes(1)
!        elseif(iElemType.eq.4) then           ! quadrangle
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!            EdgeNodes(2,1) = ElemNodes(2)
!            EdgeNodes(2,2) = ElemNodes(3)
!            EdgeNodes(3,1) = ElemNodes(3)
!            EdgeNodes(3,2) = ElemNodes(4)
!            EdgeNodes(4,1) = ElemNodes(4)
!            EdgeNodes(4,2) = ElemNodes(1)
!        else
!            continue
!        end if
!
!        return
!    end subroutine meshGet2DElementEdges

    ! ! element nodes list -> element faces list
    subroutine meshGet2DElementFaces(iElemType,ElemNodes,CellFacesNodes,iFacesTypes)
        implicit none
        integer                                       :: iElemType
        integer, dimension(NumCellNodes)              :: ElemNodes
        integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
        integer, dimension(NumCellFaces)              :: iFacesTypes

        CellFacesNodes(:,:) = 0
        iFacesTypes(:)  = 0
        if(iElemType.eq.1)      then           ! point
            continue
        elseif(iElemType.eq.2) then           ! line
            continue
        elseif(iElemType.eq.3) then           ! triangle
            CellFacesNodes(1,1) = ElemNodes(1)
            CellFacesNodes(1,2) = ElemNodes(2)
            iFacesTypes(1) = 2
            CellFacesNodes(2,1) = ElemNodes(2)
            CellFacesNodes(2,2) = ElemNodes(3)
            iFacesTypes(2) = 2
            CellFacesNodes(3,1) = ElemNodes(3)
            CellFacesNodes(3,2) = ElemNodes(1)
            iFacesTypes(3) = 2
        elseif(iElemType.eq.4) then           ! quadrangle
            CellFacesNodes(1,1) = ElemNodes(1)
            CellFacesNodes(1,2) = ElemNodes(2)
            iFacesTypes(1) = 2
            CellFacesNodes(2,1) = ElemNodes(2)
            CellFacesNodes(2,2) = ElemNodes(3)
            iFacesTypes(2) = 2
            CellFacesNodes(3,1) = ElemNodes(3)
            CellFacesNodes(3,2) = ElemNodes(4)
            iFacesTypes(3) = 2
            CellFacesNodes(4,1) = ElemNodes(4)
            CellFacesNodes(4,2) = ElemNodes(1)
            iFacesTypes(4) = 2
        else
            continue
        end if

        return
    end subroutine meshGet2DElementFaces

!    ! element nodes list -> element edges list
!    subroutine meshGet3DElementEdges(iElemType,ElemNodes,EdgeNodes)
!        implicit none
!        integer                            :: iElemType
!        integer, dimension(NumCellNodes)   :: ElemNodes
!        integer, dimension(NumCellEdges,2) :: EdgeNodes
!
!        EdgeNodes(:,:) = 0
!        if(iElemType.eq.1)      then           ! point
!            continue
!        elseif(iElemType.eq.2) then           ! line
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!        elseif(iElemType.eq.3) then           ! triangle
!            EdgeNodes(1,1) = ElemNodes(1)                  !  1
!            EdgeNodes(1,2) = ElemNodes(2)                  !  |\
!            EdgeNodes(2,1) = ElemNodes(2)                  !  | \
!            EdgeNodes(2,2) = ElemNodes(3)                  !  |  \
!            EdgeNodes(3,1) = ElemNodes(3)                  !  2---3
!            EdgeNodes(3,2) = ElemNodes(1)
!        elseif(iElemType.eq.4) then           ! quadrangle
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!            EdgeNodes(2,1) = ElemNodes(2)
!            EdgeNodes(2,2) = ElemNodes(3)
!            EdgeNodes(3,1) = ElemNodes(3)
!            EdgeNodes(3,2) = ElemNodes(4)
!            EdgeNodes(4,1) = ElemNodes(4)
!            EdgeNodes(4,2) = ElemNodes(1)
!        elseif(iElemType.eq.5) then           ! tetrahedron
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!            EdgeNodes(2,1) = ElemNodes(2)
!            EdgeNodes(2,2) = ElemNodes(3)
!            EdgeNodes(3,1) = ElemNodes(3)
!            EdgeNodes(3,2) = ElemNodes(1)
!            EdgeNodes(4,1) = ElemNodes(1)
!            EdgeNodes(4,2) = ElemNodes(4)
!            EdgeNodes(5,1) = ElemNodes(2)
!            EdgeNodes(5,2) = ElemNodes(4)
!            EdgeNodes(6,1) = ElemNodes(3)
!            EdgeNodes(6,2) = ElemNodes(4)
!        elseif(iElemType.eq.6) then           ! hexahedron
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!            EdgeNodes(2,1) = ElemNodes(2)
!            EdgeNodes(2,2) = ElemNodes(3)
!            EdgeNodes(3,1) = ElemNodes(3)
!            EdgeNodes(3,2) = ElemNodes(4)
!            EdgeNodes(4,1) = ElemNodes(4)
!            EdgeNodes(4,2) = ElemNodes(1)
!            EdgeNodes(5,1) = ElemNodes(5)
!            EdgeNodes(5,2) = ElemNodes(6)
!            EdgeNodes(6,1) = ElemNodes(6)
!            EdgeNodes(6,2) = ElemNodes(7)
!            EdgeNodes(7,1) = ElemNodes(7)
!            EdgeNodes(7,2) = ElemNodes(8)
!            EdgeNodes(8,1) = ElemNodes(8)
!            EdgeNodes(8,2) = ElemNodes(5)
!            EdgeNodes(9,1)  = ElemNodes(1)
!            EdgeNodes(9,2)  = ElemNodes(5)
!            EdgeNodes(10,1) = ElemNodes(2)
!            EdgeNodes(10,2) = ElemNodes(6)
!            EdgeNodes(11,1) = ElemNodes(3)
!            EdgeNodes(11,2) = ElemNodes(7)
!            EdgeNodes(12,1) = ElemNodes(4)
!            EdgeNodes(12,2) = ElemNodes(8)
!        elseif(iElemType.eq.7) then           ! prism
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!            EdgeNodes(2,1) = ElemNodes(2)
!            EdgeNodes(2,2) = ElemNodes(3)
!            EdgeNodes(3,1) = ElemNodes(3)
!            EdgeNodes(3,2) = ElemNodes(1)
!            EdgeNodes(4,1) = ElemNodes(4)
!            EdgeNodes(4,2) = ElemNodes(5)
!            EdgeNodes(5,1) = ElemNodes(5)
!            EdgeNodes(5,2) = ElemNodes(6)
!            EdgeNodes(6,1) = ElemNodes(6)
!            EdgeNodes(6,2) = ElemNodes(7)
!            EdgeNodes(7,1) = ElemNodes(1)
!            EdgeNodes(7,2) = ElemNodes(4)
!            EdgeNodes(8,1) = ElemNodes(2)
!            EdgeNodes(8,2) = ElemNodes(5)
!            EdgeNodes(9,1) = ElemNodes(3)
!            EdgeNodes(9,2) = ElemNodes(6)
!        elseif(iElemType.eq.8) then           ! piramid
!            EdgeNodes(1,1) = ElemNodes(1)
!            EdgeNodes(1,2) = ElemNodes(2)
!            EdgeNodes(2,1) = ElemNodes(2)
!            EdgeNodes(2,2) = ElemNodes(3)
!            EdgeNodes(3,1) = ElemNodes(3)
!            EdgeNodes(3,2) = ElemNodes(4)
!            EdgeNodes(4,1) = ElemNodes(4)
!            EdgeNodes(4,2) = ElemNodes(1)
!            EdgeNodes(5,1) = ElemNodes(1)
!            EdgeNodes(5,2) = ElemNodes(5)
!            EdgeNodes(6,1) = ElemNodes(2)
!            EdgeNodes(6,2) = ElemNodes(5)
!            EdgeNodes(7,1) = ElemNodes(3)
!            EdgeNodes(7,2) = ElemNodes(5)
!            EdgeNodes(8,1) = ElemNodes(4)
!            EdgeNodes(8,2) = ElemNodes(5)
!        end if
!
!        return
!    end subroutine meshGet3DElementEdges

    ! ! element nodes list -> element faces list
    subroutine meshGet3DElementFaces(iElemType,ElemNodes,CellFacesNodes,iFacesTypes)
        implicit none
        integer                                       :: iElemType
        integer, dimension(NumCellNodes)              :: ElemNodes
        integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
        integer, dimension(NumCellFaces)              :: iFacesTypes

        CellFacesNodes(:,:) = 0
        iFacesTypes(:)  = 0
        if(iElemType.eq.1)      then          ! point
            continue
        elseif(iElemType.eq.2) then           ! line
            continue
        elseif(iElemType.eq.3) then           ! triangle
            CellFacesNodes(1,1) = ElemNodes(1)
            CellFacesNodes(1,2) = ElemNodes(2)
            CellFacesNodes(1,3) = ElemNodes(3)
            iFacesTypes(1) = 3
        elseif(iElemType.eq.4) then           ! quadrangle
            CellFacesNodes(1,1) = ElemNodes(1)
            CellFacesNodes(1,2) = ElemNodes(2)
            CellFacesNodes(1,3) = ElemNodes(3)
            CellFacesNodes(1,4) = ElemNodes(4)
            iFacesTypes(1) = 4
        elseif(iElemType.eq.5) then           ! tetrahedron
            call meshGetTetrahedronFaces(ElemNodes,CellFacesNodes,iFacesTypes)
        elseif(iElemType.eq.6) then           ! hexahedron
            call meshGetHexahedronFaces(ElemNodes,CellFacesNodes,iFacesTypes)
        elseif(iElemType.eq.7) then           ! prism
            call meshGetPrismFaces(ElemNodes,CellFacesNodes,iFacesTypes)
        elseif(iElemType.eq.8) then           ! piramid
            call meshGetPyramidFaces(ElemNodes,CellFacesNodes,iFacesTypes)
        end if

        return
    end subroutine meshGet3DElementFaces


    subroutine meshGetCellVolume(iElemType,elemNodes,mesh,elemVolume)
        ! by spliting pyramide in two tetrahedrons:
        implicit none
        ! input
        integer :: iElemType
        integer, dimension(NumCellNodes) :: elemNodes
        type(TMeshData)                  :: mesh
        ! output
        real(r8p) :: elemVolume
        ! local

        elemVolume = 0.e0_r8p
        if(iElemType.eq.5) then
            call meshGetTetrahedronVolume(elemNodes,mesh,elemVolume)
        elseif(iElemType.eq.6) then
            call meshGetHexahedronVolume(elemNodes,mesh,elemVolume)
        elseif(iElemType.eq.7) then
            call meshGetPrismVolume(elemNodes,mesh,elemVolume)
        elseif(iElemType.eq.8) then
            call meshGetPyramidVolume(elemNodes,mesh,elemVolume)
        else
            write(*,*) " error - meshGetCellVolume: unknown",iElemType
            stop
        end if

        return
    end subroutine meshGetCellVolume

    ! - - - - - - - - - - - - - - - - - -
    subroutine meshGetFaceScalArea(iElemType,elemNodes,mesh,elemArea)
        ! by spliting pyramide in two tetrahedrons:
        implicit none
        ! input
        integer :: iElemType
        integer, dimension(NumFaceNodes) :: elemNodes
        type(TMeshData)                  :: mesh
        ! output
        real(r8p) :: elemArea
        ! local

        elemArea = 0.e0_r8p
        if(iElemType.eq.3) then
            call meshGetTriangleScalArea(elemNodes,mesh,elemArea)
        elseif(iElemType.eq.4) then
            call meshGetQuadrangleScalArea(elemNodes,mesh,elemArea)
        else
            write(*,*) " error - meshGetFaceArea: unknown element type: ",iElemType
            stop
        end if

        return
    end subroutine meshGetFaceScalArea

    ! - - - - - - - - - - - - - - - - - -
    subroutine meshGetFaceVectArea(iElemType,elemNodes,mesh,vectArea)
        ! by spliting pyramide in two tetrahedrons:
        implicit none
        ! input
        integer :: iElemType
        integer, dimension(NumFaceNodes) :: elemNodes
        type(TMeshData)                  :: mesh
        ! output
        type(T3DPoint)               :: vectArea
        ! local

        vectArea%crd = 0.e0_r8p
        if(iElemType.eq.3) then
            call meshGetTriangleVectArea(elemNodes,mesh,vectArea)
        elseif(iElemType.eq.4) then
            call meshGetQuadrangleVectArea(elemNodes,mesh,vectArea)
        else
            write(*,"(a,i2)") " error - meshGetFaceNormVect: unknown element type ",iElemType
            stop
        end if

        return
    end subroutine meshGetFaceVectArea


    ! - - - - - - - - - - - - - - - - - -
    subroutine FaceIdentify(NodesList,FaceListAtNode,NNodes,inode,ipos,icount)
        implicit none
        integer,dimension(NumFaceNodes) :: NodesList
        integer :: i,NNodes
        integer :: imin,imax,ipre
        integer :: imaxList,ipreList
        integer :: inode                ! index of base (lower) node
        integer :: ipos                 ! index of face position in list FaceAtNode(inode)%ipos
        integer :: icount

        type(TFaceListAtNode), dimension(NNodes) :: FaceListAtNode

        imin = NodesList(1)
        imax = 0
        do i = 1,NumFaceNodes
            if(NodesList(i).gt.0) then
                if(imax.lt.NodesList(i)) then
                    imax = NodesList(i)
                endif
                if(imin.gt.NodesList(i)) then
                    imin = NodesList(i)
                endif
            endif
        enddo
        ipre = 0
        do i = 1,NumFaceNodes
            if(ipre.lt.NodesList(i) .and. NodesList(i).lt.imax) then
                ipre = NodesList(i)
            endif
        enddo

        inode  = imin
        ipos   = 0
        icount = 0
        i = 0
!        do i = 1,NumNodeFaces
        do while(i.lt.NumNodeFaces .and. ipos.eq.0)
            i = i + 1
!            write(*,*) " inode = ",inode
!            write(*,*) " i     = ",i
            ipreList = FaceListAtNode(inode)%ListPos(i)%ipre
            imaxList = FaceListAtNode(inode)%ListPos(i)%imax
            icount   = FaceListAtNode(inode)%ListPos(i)%icount
!            write(*,*) " icount = ",icount
            if(imax.eq.imaxList .and. ipre.eq.ipreList) then      ! face exists
                if(ipos.eq.0) then
                    ipos = i
                else
                    write(*,*) "warning: duplicate face"
                endif
            elseif(ipos.eq.0 .and. icount.eq.0) then
                ipos = i
            end if
        enddo

        if(ipos.eq.0) then
            write(*,*) "error: not enougth NumNodeFaces"
            stop
        else
            icount = FaceListAtNode(inode)%ListPos(ipos)%icount
        end if

        ! create new face if it is not exists
        icount = icount + 1
        if(icount.eq.1) then
            FaceListAtNode(inode)%ListPos(ipos)%imin = inode
            FaceListAtNode(inode)%ListPos(ipos)%ipre = ipre
            FaceListAtNode(inode)%ListPos(ipos)%imax = imax
            FaceListAtNode(inode)%ListPos(ipos)%icount = icount
        end if
        if(icount.eq.2) then
            FaceListAtNode(inode)%ListPos(ipos)%icount = icount
        end if

        if(icount.gt.2) then
!            pc(:) = 0.e0
!            nnodes = 0
!            do i=1,NumFaceNodes
!                if(NodesList(i).gt.0) then
!                    j = NodesList(i)
!                    nnodes = nnodes + 1
!                    pc(:) = pc(:) + mesh%Node(j)%crd(:)
!                end if
!            end do
!            pc = pc/nnodes
!            write(*,"(a,3(1pE12.4e2))") "warning: duplicating face at:",pc(:)
            write(*,"(a,3(1pE12.4e2))") "warning: duplicating face"
            !pause
        end if

        return
    end subroutine FaceIdentify

    ! - - - - - - - - - - - - - - - - - -
    subroutine checkFaceCount(mesh,FaceListAtNode)
        implicit none
        ! input:
        type(TMeshData)                               :: mesh
        type(TFaceListAtNode), dimension(mesh%NNodes) :: FaceListAtNode
        ! local:
        integer :: inode,num
        integer :: ierr1,ierr3
        integer :: ipos                 ! index of face position in list FaceAtNode(inode)%ipos
        integer :: icount

        num  = 0
        ierr1 = 0
        ierr3 = 0
        do inode = 1,mesh%NNodes
            do ipos = 1,NumNodeFaces
                icount = FaceListAtNode(inode)%ListPos(ipos)%icount
                if(icount.gt.0) then
                    num = num + 1
                    if(icount.eq.1) then
                        ierr1 = ierr1 + 1
                    end if
                    if(icount.gt.2) then
                        ierr3 = ierr3 + 1
                    end if
                end if
            enddo
        enddo

        write(*,"(a,i7,a)") "error: face icount = 1: ",ierr1," faces"
        write(*,"(a,i7,a)") "error: face icount = 3: ",ierr3," faces"
        write(*,*) "Number of faces: ",num,"/",mesh%NFaces

        return
    end subroutine checkFaceCount

    ! - - - - - - - - - - - - - - - - - -
    integer function valpos_int1(array,ix)
        ! get position of integer ix in integer array
        implicit none
        integer(1), intent(in) :: array(:)
        integer :: ix,arsize
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
        integer(2), intent(in) :: array(:)
        integer :: ix,arsize
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
        integer(4), intent(in) :: array(:)
        integer :: ix,arsize
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
        integer(8), intent(in) :: array(:)
        integer :: ix,arsize
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
    integer function valpos_int(array,ix)
        ! get position of integer ix in integer array
        implicit none
        integer, intent(in) :: array(:)
        integer :: ix,arsize
        integer :: i

        arsize = size(array)
        i = 1
        valpos_int = 0
        do while(array(i).ne.ix)
            i = i + 1
            if(i.gt.arsize) then
                return
            endif
        enddo
        valpos_int = i
        return
    end function valpos_int

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

end module module_Mesh

! -------------------------------------------------------------
subroutine mesh_output_vtk(fileNameOut,mesh)
    ! version 4
    use module_Settings
    use module_Mesh             ! KV (including additional grid data)
    use LIB_VTK_IO
    use module_TxtRead

    IMPLICIT NONE
    type(TMeshData) :: mesh
    REAL(8), DIMENSION(mesh%NCells)   :: var
    REAL(4), DIMENSION(mesh%NNodes)   :: X,Y,Z
    INTEGER, DIMENSION(5*mesh%NCells) :: connect
    INTEGER, DIMENSION(mesh%NCells)   :: cell_type
    integer :: E_IO
    integer :: i,j,k
    real(8) :: CurTime,dtime
    character(*) :: fileNameOut
    character(64) :: pathOut
    character(64) :: mesh_topology
    character(256) :: vtkTitle

    CurTime = 0.d0
    pathOut = trim(OutputDir)//trim(fileNameOut)//'.vtk'
    call slash_bck2frw(pathOut)
    write(*,*) 'pathOut: ',trim(pathOut)
    vtkTitle = ''
    ! ------------------------- set title and topology ----------------------------
    mesh_topology = 'UNSTRUCTURED_GRID'
    E_IO = VTK_INI(vtkFmtOut,pathOut,vtkTitle,mesh_topology)

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
    do i = 1,mesh%NCells
        j = 5*(i-1)+1
        connect(j) = 4
        do k = 1,4
            connect(j+k) = mesh%Cell(i)%inode(k)-1
        enddo
    enddo

    ! cell type
    do i = 1,mesh%NCells
        cell_type(i) = 9   ! => VTK quad type
    enddo

    E_IO = VTK_CON(mesh%NCells,connect,cell_type)


    ! ------------ DATA HEADER (location: cell, node) ----------------------------
    E_IO = VTK_DAT(mesh%NCells,'cell')

    ! ------------------------- PCI ----------------------------
    ! PCI-1,HRFI,FHZ,FHR,PLT,TEM,VIF,VIR,VIZ,SS,AL
    do i = 1,mesh%NCells
        var(i) = i
    enddo
    E_IO = VTK_VAR(mesh%NNodes,'cellId',var)

    ! --------------------------- VTK END ----------------------------
    E_IO = VTK_END()

    return
end subroutine mesh_output_vtk


! --------------------------------------------
function point_in_cell(mesh,pt,icell,rtol)
    ! <!> - check, works not correctly
    use module_Mesh
    use module_Geom
    implicit none
    ! input:
    type(TMeshData)      :: mesh
    real(8),dimension(3) :: pt
    real(8)              :: rtol         ! tolerance near bound
    ! result:
    real(8)              :: point_in_cell
    ! local:
    real(8)              :: val0
    logical              :: bndFace             ! if .true. - external
    real(8),dimension(3) :: fnorm
    real(8),dimension(3) :: ptc
    real(8),dimension(3) :: vec
    integer :: icell
    integer :: i
    integer :: idir
    integer :: iface
    real(8) :: eps,clen

    point_in_cell = -1.0e10
    if(icell.gt.mesh%NCells .or. icell.eq.0) then
!        write(*,*) 'error /point_in_cell/: icell = ',icell,'>',mesh%NCells," = NCells"
        return
    endif

    ptc = mesh%Cell(icell)%center
    vec = pt - ptc
    clen = vectLength(vec)
    point_in_cell = clen

    do i = 1,mesh%Cell(icell)%nFaces
        iface = mesh%Cell(icell)%iface(i)
        if(mesh%Face(iface)%icell(2).gt.0) then
            bndFace = .false.
        else
            bndFace = .true.
        end if
        idir  = mesh%Cell(icell)%iFaceDir(i)
        fnorm = (-1)*mesh%Face(iface)%vnorm*idir       ! internal normal
        ptc = mesh%Face(iface)%center
        vec = pt - ptc
        eps = vectScalProd(vec,fnorm)             !    [    x <-vec ] nor->  : eps > 0 - in cell

        if(bndFace) then
            val0 = -rtol                ! <!> thin place
        else
            val0 = 0.0d0
        end if
        if(eps.lt.val0) then            ! should be "eps.gt.0.e0" <!>
            point_in_cell = -1.0e0*clen
            return
        endif
    enddo

    return
end function point_in_cell

! --------------------------------------------
subroutine get_point_cell(mesh,pt,rtol,icell,rdist)
    ! <!> - check, works not correctly
    use module_Mesh
    use module_Geom
    implicit none
    ! input parameters
    type(TMeshData)      :: mesh
    real(8),dimension(3) :: pt
    real(8)              :: rtol        ! tolerance near bound
    ! output parameters
    integer              :: icell
    ! local parameters
    real(8)              :: rdist       ! distance near bound
    real(8)              :: point_in_cell
    integer              :: tcell       ! temporal cell
    real(8)              :: tdist       ! temporal distance

    icell = 1
    rdist = point_in_cell(mesh,pt,icell,rtol)
    tcell = icell
    tdist = rdist

    do while(rdist.lt.0.d0)
        icell = icell + 1
        if(icell.gt.mesh%NCells) then
            icell = tcell
            rdist = tdist
            return
        endif

        rdist = point_in_cell(mesh,pt,icell,rtol)
        if(rdist.gt.tdist) then
            tdist = rdist
            tcell = icell
        end if
    enddo

    return
end subroutine get_point_cell


! --------------------------------------------
subroutine get_point_node(mesh,pt,inode,dist)
    ! <!> - check, works not correctly
    use module_Settings
    use module_Geom
    use module_Mesh
    implicit none
    ! input parameters
    type(TMeshData)      :: mesh
    real(8),dimension(3) :: pt
    ! output parameters
    integer              :: inode
    ! local parameters
    integer              :: i
    real(8),dimension(3) :: ptn
    real(R8P) :: eps,dist

    dist = 1.e10
    inode = 0
    do i = 1,mesh%NNodes
        ptn = mesh%Node(i)%crd
        eps = point_point_distance(pt,ptn)
        if(eps.lt.dist) then
            dist  = eps
            inode = i
        endif
    enddo

    return
end subroutine get_point_node

! --------------------------------------------
subroutine get_point_node_from_list(mesh,nnodes,nodelist,pt,inode,dist)
    ! <!> - check, works not correctly
    use module_Settings
    use module_Geom
    use module_Mesh
    implicit none
    ! input parameters
    type(TMeshData)           :: mesh
    integer                   :: nnodes
    integer,dimension(nnodes) :: nodelist
    real(8),dimension(3)      :: pt
    ! output parameters
    integer                   :: inode
    ! local parameters
    integer              :: i,jnode
    real(8),dimension(3) :: ptn
    real(R8P) :: eps,dist

    dist  = 1.e10
    inode = 0
    do i = 1,nnodes
        jnode = nodelist(i)
        if(jnode.le.0) then
            continue
        else
            ptn = mesh%Node(jnode)%crd
            eps = point_point_distance(pt,ptn)
            if(eps.lt.dist) then
                dist  = eps
                inode = jnode
            endif
        end if
    enddo

    return
end subroutine get_point_node_from_list


! --------------------------------------------
subroutine get_point_face_from_list(mesh,nfaces,facelist,pt,iface,dist)
    ! <!> - check, works not correctly
    use module_Settings
    use module_Geom
    use module_Mesh
    implicit none
    ! input parameters
    type(TMeshData)           :: mesh
    integer                   :: nfaces
    integer,dimension(nfaces) :: facelist
    real(8),dimension(3)      :: pt
    ! output parameters
    integer                   :: iface
    ! local parameters
    integer              :: j,iflag
    integer              :: jface           ! element index
    real(8),dimension(3) :: fnorm                 ! edge normal vector
    real(8),dimension(3) :: vec                   ! edge normal vector
    real(r8p) :: eps,dist,sprod

    iface = 0
    dist  = 1.e10
    eps   = 0.0e0_r8p

    do j = 1,nfaces
        jface = facelist(j)
        if(jface.le.0) then
            continue
        else
            iflag = 1
            fnorm = mesh%Face(jface)%vnorm
            vec   = pt - mesh%Face(jface)%center
            sprod = vectScalProd(fnorm,vec)
            sprod = abs(sprod)
            if(sprod.lt.dist) then
                iface = jface
                dist  = sprod
            end if
        end if
    enddo

    return
end subroutine get_point_face_from_list


! --------------------------------------------
subroutine get_cell_linsize(mesh,icell,dl)
    ! <!> - check, works not correctly
    use module_Settings
    use module_Geom
    use module_Mesh
    implicit none
    ! input parameters
    type(TMeshData)           :: mesh
    integer                   :: icell
    ! output parameters
    real(r8p)                 :: dl
    ! local parameters
    integer                :: i,iface           ! element index
    real(r8p)              :: dist
    real(r8p),dimension(3) :: vec

    dl = 1.0e10_r8p
    do i = 1,mesh%Cell(icell)%nfaces
        iface = mesh%Cell(icell)%iface(i)
        vec = mesh%Cell(icell)%center - mesh%Face(iface)%center
        dist = vectLength(vec)
        if(dist.lt.dl) then
            dl = dist
        end if
    end do
    dl = 2*dl
    return
end subroutine get_cell_linsize


! --------------------------------------------
function point_in_cell_mapped(mesh,pt,icell,rtol,mapper)
    ! <!> - check, works not correctly
    use module_Mesh
    use module_Geom
    implicit none
    ! input:
    type(TMeshData)      :: mesh
    real(8),dimension(3) :: pt
    real(8)              :: rtol       ! tolerance near bound
    ! result:
    logical :: point_in_cell_mapped
    ! optional:
    external :: mapper
    optional :: mapper
    ! local:
    real(8)              :: val0
    logical              :: bndFace             ! if .true. - external
    real(8),dimension(3) :: fnorm
    real(8),dimension(3) :: ptc
    real(8),dimension(3) :: vec
    integer :: icell
    integer :: i
    integer :: idir
    integer :: iface
    real(8) :: eps

    point_in_cell_mapped = .true.
    if(icell.gt.mesh%NCells .or. icell.eq.0) then
!        write(*,*) 'error /point_in_cell_mapped/: icell = ',icell,'>',mesh%NCells," = NCells"
        return
    endif
    do i = 1,mesh%Cell(icell)%nFaces
        iface = mesh%Cell(icell)%iface(i)
        if(mesh%Face(iface)%icell(2).gt.0) then
            bndFace = .false.
        else
            bndFace = .true.
        end if
        ! normal vector
        idir  = mesh%Cell(icell)%iFaceDir(i)
        fnorm = mesh%Face(iface)%vnorm*idir       ! external normal
        if(present(mapper)) then
            call mapper(fnorm)
        endif
        ! face center
        ptc = mesh%Face(iface)%center

        vec = pt - ptc
        if(present(mapper)) then
            call mapper(vec)
        endif

        eps = vectScalProd(vec,fnorm)             !    [    x <-vec ] nor->  : eps < 0 - in cell
!        write(*,"(a,4(1pE12.4e3))") 'pt      = ',pt
!        write(*,"(a,4(1pE12.4e3))") 'fcenter = ',mesh%Face(iface)%center
!        write(*,"(a,4(1pE12.4e3))") 'fnorm   = ',fnorm*idir
!        write(*,"(a,1(1pE12.4e3))") ' -> eps = ',eps
!        pause
        val0 = 0.d0
        if(bndFace) then
            val0 = rtol
        endif
        if(eps.gt.val0) then            ! should be "eps.gt.0.e0" <!>
            point_in_cell_mapped = .false.
            return
        endif

    enddo

    return
end function point_in_cell_mapped

! --------------------------------------------
subroutine get_point_face_from_list_mapped(mesh,nfaces,facelist,pt,iface,dist,mapper)
    ! <!> - check, works not correctly
    use module_Settings
    use module_Geom
    use module_Mesh
    implicit none
    ! input parameters
    type(TMeshData)           :: mesh
    integer                   :: nfaces
    integer,dimension(nfaces) :: facelist
    real(8),dimension(3)      :: pt
    ! output parameters
    integer                   :: iface
    ! optional:
    external :: mapper
    optional :: mapper
    ! local parameters
    integer              :: j,iflag
    integer              :: jface           ! element index
    real(8),dimension(3) :: fnorm                 ! edge normal vector
    real(8),dimension(3) :: cvec                  ! edge center
    real(8),dimension(3) :: vec                   ! edge normal vector
    real(r8p) :: eps,dist,sprod

    iface = 0
    dist  = 1.e10
    eps   = 0.0e0_r8p

    do j = 1,nfaces
        jface = facelist(j)
        if(jface.le.0) then
            continue
        else
            iflag = 1
            fnorm = mesh%Face(jface)%vnorm
            if(present(mapper)) then
                call mapper(fnorm)
            endif

            cvec   = mesh%Face(jface)%center
            if(present(mapper)) then
                call mapper(cvec)
            endif
            vec   = pt - cvec

            sprod = vectScalProd(fnorm,vec)
            sprod = abs(sprod)
            if(sprod.lt.dist) then
                iface = jface
                dist  = sprod
            end if
        end if
    enddo

    return
end subroutine get_point_face_from_list_mapped
