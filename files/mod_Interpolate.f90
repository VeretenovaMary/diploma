module module_Interpolate
    use module_Settings
    use module_Geom
    use module_Mesh
    use module_utils

    ! --- weights for function value interpolation ---
    type TWeightsForValues                         ! weights for funct. value interpolation
        integer(1)                   :: nPush                ! number of source elements
        integer, pointer             :: iPush(:) => null()   ! indexes of source elements (nodes or cells)
        real(r8p), pointer           :: wPush(:) => null()   ! weights of source elements (nodes or cells)
    end type TWeightsForValues

    ! --- weights for function spatial derivatives interpolation ---
    type TWeightsForDerivs                                 ! weights for funct. spatial derivs interpolation
        integer(1)                   :: nPush                ! number of source elements
        integer, pointer             :: iPush(:) => null()   ! indexes of source elements (nodes or cells)
        type(T3DPoint), pointer  :: wPush(:) => null()   ! weights of source elements (nodes or cells)
    end type TWeightsForDerivs

    type TValuesInterpolator
        type(TMeshData), pointer         :: pPushMesh            ! source mesh
        type(TMeshData), pointer         :: pPullMesh            ! acceptor mesh
        integer(1)                       :: PushElemType         ! source mesh elements type: 0-node,1-edge,2-face,3-cell
        integer(1)                       :: PullElemType         ! accept mesh elements type: 0-node,1-edge,2-face,3-cell
        integer                          :: NElems              ! nmber of interpolator elements: nodes/edges/faces/cells
        type(TWeightsForValues), pointer :: ielem(:) => null()  ! funct. weights for values
    end type TValuesInterpolator

    type TDerivsInterpolator
        type(TMeshData), pointer         :: pPushMesh            ! source mesh
        type(TMeshData), pointer         :: pPullMesh            ! acceptor mesh
        integer(1)                       :: PushElemType         ! source mesh elements type: 0-node,1-edge,2-face,3-cell
        integer(1)                       :: PullElemType         ! accept mesh elements type: 0-node,1-edge,2-face,3-cell
        integer                          :: NElems              ! nmber of interpolator elements: "node"/"cell"
        type(TWeightsForDerivs), pointer :: ielem(:) => null()  ! funct. weights for derivs in diff. directions
    end type TDerivsInterpolator


    interface close_interpolator
        module procedure close_interpolator_for_values
        module procedure close_interpolator_for_derivs
    end interface

    !private new_interpolator_for_derivs
    private close_interpolator_for_values,close_interpolator_for_derivs

    public :: mapper_XYZ_to_RZ_dirX
    public :: mapper_XY_to_X

contains

    ! --- close (deallocate) interpolator for values ---
    subroutine close_interpolator_for_values(valsInterp)
        implicit none
        ! input
        type(TValuesInterpolator) :: valsInterp
        ! local
        integer :: i

        do i = 1,valsInterp%NElems
            deallocate(valsInterp%ielem(i)%wPush)
            deallocate(valsInterp%ielem(i)%iPush)
        end do
        deallocate(valsInterp%ielem)

        return
    end subroutine close_interpolator_for_values

    ! --- close (deallocate) interpolator for derivatives ---
    subroutine close_interpolator_for_derivs(dersInterp)
        implicit none
        ! input
        type(TDerivsInterpolator) :: dersInterp
        ! local
        integer :: i

        do i = 1,dersInterp%NElems
            deallocate(dersInterp%ielem(i)%iPush)
            deallocate(dersInterp%ielem(i)%wPush)
        end do
        deallocate(dersInterp%ielem)

        return
    end subroutine close_interpolator_for_derivs

    ! --- create new interpolator for values ---
    subroutine allocate_valsInterp_weight(valsInterpElem,nPush)
        implicit none
        type(TWeightsForValues) :: valsInterpElem
        integer                 :: nPush

        allocate(valsInterpElem%iPush(nPush))
        allocate(valsInterpElem%wPush(nPush))
        valsInterpElem%nPush = int(nPush,1)

        return
    end subroutine

    ! --- create new interpolator for values ---
    subroutine allocate_dersInterp_weight(dersInterpElem,nPush)
        implicit none
        type(TWeightsForDerivs) :: dersInterpElem
        integer                 :: nPush

        allocate(dersInterpElem%iPush(nPush))
        allocate(dersInterpElem%wPush(nPush))
        dersInterpElem%nPush = int(nPush,1)

        return
    end subroutine

    ! --- create new interpolator for values ---
    subroutine new_interp_ValExtNodesToNodes(Interp,pPushMesh,pPullMesh,mapper)
        implicit none
        ! input-output
        type(TValuesInterpolator) :: Interp
        type(TMeshData), pointer  :: pPushMesh             ! source mesh
        type(TMeshData), pointer  :: pPullMesh             ! accept mesh
        ! optional:
        external :: mapper
        optional :: mapper
        ! local
        integer :: NElems
        real(8) :: objSizeMb
        integer               :: iPushCell,iPushFace,iPushNode
        integer               :: iPullNode
        integer               :: NPushNodes,NPullNodes
        integer               :: i,istop
        type(T3DPoint)        :: nodes(NumCellFaces+1)     ! coords of face nodes: size = 4 (quad) or 3 (tri)
        type(T3DPoint)        :: pt0
        type(T3DPoint)        :: ptc,vec
        real(r8p)             :: Fwt(NumCellFaces+1)       ! function weights to Fpt
        type(T3DPoint)        :: dFwt(NumCellFaces+1)      ! function weights to dFpt
        integer               :: nfaces,nnodes
        integer, dimension(6) :: faceslist
        integer, dimension(6) :: nodeslist
        real(r8p)             :: rdist,dist,dlen,rdistF,rdistN       ! function weights to Fpt

        write(*,*)
        write(*,*) "....... new interpolator for gf values ........"
        write(*,"(a,a,a,a)") " nodes->nodes    : ",trim(pPushMesh%meshTit)," -> ",trim(pPullMesh%meshTit)


        ! set source mesh elements type
        Interp%PushElemType = 0        ! 0 - node
        Interp%PullElemType = 0        ! 0 - node

        Interp%pPushMesh => pPushMesh
        Interp%pPullMesh => pPullMesh

        ! allocate memory for interpolator:
        ! select interpolator type: node->node, cell->node;
        NElems = pPullMesh%NNodes
        Interp%NElems = pPullMesh%NNodes
        allocate(Interp%ielem(Interp%NElems))


        NPullNodes = pPullMesh%NNodes
        do iPullNode = 1,NPullNodes
            write(*,"(a,a,I3,a)",advance="no") achar(13),' done ',iPullNode*100/NPullNodes," %"
            pt0%crd = pPullMesh%Node(iPullNode)%crd
            if(present(mapper)) then
                call mapper(pPullMesh%Node(iPullNode)%crd,pt0%crd)
            endif

            call get_point_cell(pPushMesh,pt0,1.0d-7,iPushCell,rdist)
            call get_cell_linsize(pPushMesh,iPushCell,dlen)

            if(rdist.ge.0.d0) then
                ! find "Push Face" if exists
                nfaces = pPushMesh%Cell(iPushCell)%nfaces
                faceslist(:)        = 0
                faceslist(1:nfaces) = pPushMesh%Cell(iPushCell)%iface(1:nfaces)
                call get_point_face_from_list(pPushMesh,nfaces,faceslist,pt0,iPushFace,dist)
                rdistF = dist/dlen

                ! find "Push Node" if exists
                nnodes = pPushMesh%Face(iPushFace)%nnodes
                nodeslist(:)        = 0
                nodeslist(1:nnodes) = pPushMesh%Face(iPushFace)%inode(1:nnodes)
                call get_point_node_from_list(pPushMesh,nnodes,nodeslist,pt0,iPushNode,dist)
                rdistN = dist/dlen
                if(rdistN.le.1.0e-2) then
                    NPushNodes = 1
                    call allocate_valsInterp_weight(Interp%ielem(iPullNode),NPushNodes)
                    Interp%ielem(iPullNode)%nPush = int(NPushNodes,1)
                    Interp%ielem(iPullNode)%iPush(1) = iPushNode
                    Fwt(:) = 0.e0_r8p
                    Fwt(1) = 1.e0_r8p
                    Interp%ielem(iPullNode)%wPush(1) = Fwt(1)
                elseif(rdistF.le.1.0e-2) then
                    NPushNodes = pPushMesh%Face(iPushFace)%nnodes
                    call allocate_valsInterp_weight(Interp%ielem(iPullNode),NPushNodes)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        nodes(i)%crd = pPushMesh%Node(iPushNode)%crd
                    end do
! improve for 2 points   call interp_bilinear(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt,istop)  ! improve for 2 points
!                    call interp_least_squares(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt)
                    dist = 0.0d0
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        ptc%crd = pPushMesh%Node(iPushNode)%crd
                        vec%crd = pt0%crd - ptc%crd
                        dist = dist + 1.0d0/(vectLength(vec)+1.0d-4*dlen)
                    enddo
                    Fwt(:)  = 1.e0_r8p/NPushNodes                          ! <!> average needed
                    dFwt(1)%crd(:) = 0.e0_r8p                              ! <!> average needed
                    Interp%ielem(iPullNode)%nPush = int(NPushNodes,1)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        ptc%crd = pPushMesh%Node(iPushNode)%crd
                        vec%crd = pt0%crd - ptc%crd
                        Fwt(i) = 1.0d0/(vectLength(vec)+1.0d-4*dlen)/dist
                        Interp%ielem(iPullNode)%iPush(i) = iPushNode
                        Interp%ielem(iPullNode)%wPush(i) = Fwt(i)
                    enddo
                else
                    NPushNodes = pPushMesh%Cell(iPushCell)%nnodes
                    call allocate_valsInterp_weight(Interp%ielem(iPullNode),NPushNodes)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Cell(iPushCell)%inode(i)
                        nodes(i)%crd = pPushMesh%Node(iPushNode)%crd
                    end do
                    call interp_bilinear(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt,istop)
!                    call interp_least_squares(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt)
!                    Fwt(:)  = 1.e0_r8p/NPushNodes                          ! <!> average needed
                    dFwt(1)%crd(:) = 0.e0_r8p                              ! <!> average needed
                    Interp%ielem(iPullNode)%nPush = int(NPushNodes,1)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Cell(iPushCell)%inode(i)
                        Interp%ielem(iPullNode)%iPush(i) = iPushNode
                        Interp%ielem(iPullNode)%wPush(i) = Fwt(i)
                    enddo
                end if
           else
                nfaces = pPushMesh%Cell(iPushCell)%nfaces
                faceslist(:)        = 0
                faceslist(1:nfaces) = pPushMesh%Cell(iPushCell)%iface(1:nfaces)
                call get_point_face_from_list(pPushMesh,nfaces,faceslist,pt0,iPushFace,dist)
                rdistF = dist/dlen

                ! find "Push Node" if exists
                nnodes = pPushMesh%Face(iPushFace)%nnodes
                nodeslist(:)        = 0
                nodeslist(1:nnodes) = pPushMesh%Face(iPushFace)%inode(1:nnodes)
                call get_point_node_from_list(pPushMesh,nnodes,nodeslist,pt0,iPushNode,dist)
                rdistN = dist/dlen

                if(rdistN.le.1.0e-2) then
                    NPushNodes = 1
                    call allocate_valsInterp_weight(Interp%ielem(iPullNode),NPushNodes)
                    Fwt(:) = 0.e0_r8p
                    Fwt(1) = 1.e0_r8p
                    Interp%ielem(iPullNode)%nPush = int(NPushNodes,1)
                    Interp%ielem(iPullNode)%iPush(1) = iPushNode
                    Interp%ielem(iPullNode)%wPush(1) = Fwt(1)
                else
                    NPushNodes = pPushMesh%Face(iPushFace)%nnodes
                    call allocate_valsInterp_weight(Interp%ielem(iPullNode),NPushNodes)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        nodes(i)%crd = pPushMesh%Node(iPushNode)%crd
                    end do
! improve for 2 points   call interp_bilinear(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt,istop)  ! improve for 2 points
!                    call interp_least_squares(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt)
                    dist = 0.0d0
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        ptc%crd = pPushMesh%Node(iPushNode)%crd
                        vec%crd = pt0%crd - ptc%crd
                        dist = dist + 1.0d0/(vectLength(vec)+1.0d-4*dlen)
                    enddo
                    Fwt(:)  = 0.e0_r8p                          ! <!> average needed
                    Interp%ielem(iPullNode)%nPush = int(NPushNodes,1)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        ptc%crd = pPushMesh%Node(iPushNode)%crd
                        vec%crd = pt0%crd - ptc%crd
                        Fwt(i) = 1.0d0/(vectLength(vec)+1.0d-4*dlen)/dist
                        dFwt(i)%crd(:) = 0.e0_r8p                              ! <!> average needed
                        Interp%ielem(iPullNode)%iPush(i) = iPushNode
                        Interp%ielem(iPullNode)%wPush(i) = Fwt(i)
                    enddo
                end if
            end if
!            write(*,"(a,i3)",advance="no") " NPushNodes = ",Interp%ielem(iPullNode)%nPush
!            if(Interp%ielem(iPullNode)%nPush.ne.1) pause
        enddo

        call valsInterp_get_size(Interp,objSizeMb)
        write(*,"(a,1pE12.4e2,a)") " InterpSize = ",objSizeMb," Mb"

        return
    end subroutine new_interp_ValExtNodesToNodes


    ! --- create new interpolator for values ---
    subroutine new_interp_ValExtNodesToCells(Interp,pPushMesh,pPullMesh,mapper)
        implicit none
        ! input-output
        type(TValuesInterpolator) :: Interp
        type(TMeshData), pointer  :: pPushMesh             ! source mesh
        type(TMeshData), pointer  :: pPullMesh             ! accept mesh
        ! optional:
        external :: mapper
        optional :: mapper
        ! local
        integer :: NElems,Nwts
        real(8) :: objSizeMb
        integer                 :: iPullCell
        integer                 :: iPushCell,iPushFace,iPushNode
        integer                 :: NPushNodes,NPullCells
        integer                 :: i
        type(T3DPoint)          :: nodes(NumCellFaces+1)     ! coords of face nodes: size = 4 (quad) or 3 (tri)
        type(T3DPoint)          :: pt0
        real(r8p)               :: Fwt(NumCellFaces+1)       ! function weights to Fpt
        type(T3DPoint)          :: dFwt(NumCellFaces+1)      ! function weights to dFpt
        integer                 :: nfaces,nnodes
        integer, dimension(6)   :: faceslist
        integer, dimension(6)   :: nodeslist
        real(r8p)               :: rdist,dist,dlen,rdistF,rdistN       ! function weights to Fpt

        write(*,*)
        write(*,*) "....... new interpolator for gf values ........"
        write(*,"(a,a,a,a)") " nodes->nodes    : ",trim(pPushMesh%meshTit)," -> ",trim(pPullMesh%meshTit)


        ! set source mesh elements type
        Interp%PushElemType = 0        ! 0 - node
        Interp%PullElemType = 3        ! 3 - cell

        Interp%pPushMesh => pPushMesh
        Interp%pPullMesh => pPullMesh

        ! allocate memory for interpolator:
        ! select interpolator type: node->node, cell->node;
        NElems        = pPullMesh%NCells
        Interp%NElems = pPullMesh%NCells
        allocate(Interp%ielem(NElems))
        do i = 1,pPullMesh%NCells
            Nwts = NumCellNodes
            call allocate_valsInterp_weight(Interp%ielem(i),Nwts)
        end do

        call valsInterp_get_size(Interp,objSizeMb)
        write(*,"(a,1pE12.4e2,a)") " InterpSize = ",objSizeMb," Mb"

        NPullCells = pPullMesh%NCells
        do iPullCell = 1,NPullCells
            !write(*,"(a,a,I3,a)",advance="no") achar(13),' done ',iPullCell*100/NPullCells," %"
            pt0%crd = pPullMesh%Cell(iPullCell)%center
            if(present(mapper)) then
!                call mapper(pt0%crd)
                call mapper(pt0)          ! works that way
            endif

            ! find "Push Cell"
            call get_point_cell(pPushMesh,pt0,1.0d-7,iPushCell,rdist)
            if(rdist.ge.0.d0) then
                call get_cell_linsize(pPushMesh,iPushCell,dlen)
!                write(*,*)
!                write(*,*)                     " cell = ",iPushCell
!                write(*,"(a,3(1x,1pE12.4e3))") " cntr:  ",pPushMesh%Cell(iPushCell)%center
!                write(*,"(a,3(1x,1pE12.4e3))") " pt:    ",pt0%crd(:)

                ! find "Push Face" if exists
                nfaces = pPushMesh%Cell(iPushCell)%nfaces
                faceslist(:)        = 0
                faceslist(1:nfaces) = pPushMesh%Cell(iPushCell)%iface(1:nfaces)
                call get_point_face_from_list(pPushMesh,nfaces,faceslist,pt0,iPushFace,dist)
                rdistF = dist/dlen
!                write(*,"(a,6(i4))")             "FacesList: ",faceslist(1:nfaces)
!                write(*,"(a)")                   "Face: "
!                write(*,"(i3,a,3(1x,1pE12.4e3))") iPushFace," cntr  = ",pPushMesh%Face(iPushFace)%center
!                write(*,"(i3,a,3(1x,1pE12.4e3))") iPushFace," fdist = ",rdistF,dlen

                ! find "Push Node" if exists
                nnodes = pPushMesh%Face(iPushFace)%nnodes
                nodeslist(:)        = 0
                nodeslist(1:nnodes) = pPushMesh%Face(iPushFace)%inode(1:nnodes)
                call get_point_node_from_list(pPushMesh,nnodes,nodeslist,pt0,iPushNode,dist)
                rdistN = dist/dlen
!                write(*,"(a)")                   "Node: "
!                write(*,"(i3,a,3(1x,1pE12.4e3))") iPushNode," cntr  = ",pPushMesh%Node(iPushNode)%crd
!                write(*,"(i3,a,3(1x,1pE12.4e3))") iPushNode," ndist = ",rdistN,dlen

                if(rdistN.le.1.0e-2) then
                    NPushNodes = 1
                    Interp%ielem(iPullCell)%nPush = int(NPushNodes,1)
                    Interp%ielem(iPullCell)%iPush(1) = iPushNode
                    Fwt(:) = 0.e0_r8p
                    Fwt(1) = 1.e0_r8p
                    Interp%ielem(iPullCell)%wPush(1) = Fwt(1)
!                    write(*,*) " case: Node",NPushNodes
                elseif(rdistF.le.1.0e-2) then
                    NPushNodes = pPushMesh%Face(iPushFace)%nnodes
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        nodes(i)%crd = pPushMesh%Node(iPushNode)%crd
                    end do
! improve for 2 points   call interp_bilinear(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt,istop)  ! improve for 2 points
!                    call interp_least_squares(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt)
                    Fwt(:)  = 1.e0_r8p/NPushNodes                          ! <!> average needed
                    dFwt(1)%crd(:) = 0.e0_r8p                              ! <!> average needed
                    Interp%ielem(iPullCell)%nPush = int(NPushNodes,1)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Face(iPushFace)%inode(i)
                        Interp%ielem(iPullCell)%iPush(i) = iPushNode
                        Interp%ielem(iPullCell)%wPush(i) = Fwt(i)
                    enddo
!                    write(*,*) " case: Face",NPushNodes
                else
                    NPushNodes = pPushMesh%Cell(iPushCell)%nnodes
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Cell(iPushCell)%inode(i)
                        nodes(i)%crd = pPushMesh%Node(iPushNode)%crd
                    end do
!                    call interp_bilinear(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt,istop)
!                    call interp_least_squares(NPushNodes,nodes(1:NPushNodes),pt0,Fwt,dFwt)
                    Fwt(:)  = 1.e0_r8p/NPushNodes                          ! <!> average needed
                    dFwt(1)%crd(:) = 0.e0_r8p                              ! <!> average needed
                    Interp%ielem(iPullCell)%nPush = int(NPushNodes,1)
                    do i = 1,NPushNodes
                        iPushNode    = pPushMesh%Cell(iPushCell)%inode(i)
                        Interp%ielem(iPullCell)%iPush(i) = iPushNode
                        Interp%ielem(iPullCell)%wPush(i) = Fwt(i)
                    enddo
                end if
            else
                write(*,*) "warning: point is outside mesh"
            end if
        enddo
        write(*,*)

        return
    end subroutine new_interp_ValExtNodesToCells


    ! --- create new interpolator for values ---
    subroutine new_interp_ValExtCellsToNodes(Interp,pPushMesh,pPullMesh,mapper)
        implicit none
        ! input-output
        type(TValuesInterpolator) :: Interp
        type(TMeshData), pointer  :: pPushMesh             ! source mesh
        type(TMeshData), pointer  :: pPullMesh             ! accept mesh
        ! optional:
        external :: mapper
        optional :: mapper
        ! local
        logical :: point_in_cell_mapped
        integer :: NElems,Nwts
        real(8) :: objSizeMb
        integer                 :: iPullNode
        integer                 :: iPushCell,iPushFace
        integer                 :: NPullNodes,NPushCells
        integer                 :: i
        type(T3DPoint)          :: pt0
        integer                 :: ncells
        real(r8p)               :: rtol                      ! function weights to Fpt
        integer                 :: nfaces
        integer, dimension(6)   :: faceslist
        real(r8p)               :: rdistF,dist,dlen
        logical                 :: event

        write(*,*)
        write(*,*) "....... new interpolator for gf values ........"
        write(*,"(a,a,a,a)") " nodes->nodes    : ",trim(pPushMesh%meshTit)," -> ",trim(pPullMesh%meshTit)


        ! set source mesh elements type
        Interp%PushElemType = 3        ! 3 - cell
        Interp%PullElemType = 0        ! 0 - node

        Interp%pPushMesh => pPushMesh
        Interp%pPullMesh => pPullMesh

        ! allocate memory for interpolator:
        ! select interpolator type: node->node, cell->node;
        NElems        = pPullMesh%NNodes
        Interp%NElems = pPullMesh%NNodes
        allocate(Interp%ielem(NElems))
        do i = 1,pPullMesh%NNodes
            Nwts = 1        ! 1 cell for each node
            call allocate_valsInterp_weight(Interp%ielem(i),Nwts)
        end do

        call valsInterp_get_size(Interp,objSizeMb)
        write(*,"(a,1pE12.4e2,a)") " InterpSize = ",objSizeMb," Mb"

        NPullNodes = pPullMesh%NNodes
        NPushCells = pPushMesh%NCells
        do iPullNode = 1,NPullNodes
            !write(*,"(a,a,I3,a)",advance="no") achar(13),' done ',iPullNode*100/NPullNodes," %"

            pt0%crd = pPullMesh%Node(iPullNode)%crd
            rtol = 1.0e-2_r8p
            ! get number of cells (to allocate memory)
            ncells = 0
            iPushCell = 1
            event = point_in_cell_mapped(pPushMesh,pt0,iPushCell,rtol,mapper)
            do while(.not.event .and. iPushCell.le.pPushMesh%NCells)
                iPushCell = iPushCell + 1
                event = point_in_cell_mapped(pPushMesh,pt0,iPushCell,rtol,mapper)
            enddo

            if(point_in_cell_mapped(pPushMesh,pt0,iPushCell,rtol,mapper)) then
                call get_cell_linsize(pPushMesh,iPushCell,dlen)

                ! find "Push Face" if exists
                nfaces = pPushMesh%Cell(iPushCell)%nfaces
                faceslist(:)        = 0
                faceslist(1:nfaces) = pPushMesh%Cell(iPushCell)%iface(1:nfaces)
                call get_point_face_from_list_mapped(pPushMesh,nfaces,faceslist,pt0,iPushFace,dist,mapper)
                rdistF = dist/dlen

                if(rdistF.le.1.0e-1) then

                    ncells = 2
                    Nwts = ncells        ! number of weights
                    call allocate_valsInterp_weight(Interp%ielem(iPullNode),Nwts)

                    ! <!> better to make least squares interpolation weights
                    ! set 1st cell and weight
                    Interp%ielem(iPullNode)%iPush(1) = iPushCell
                    Interp%ielem(iPullNode)%wPush(1) = 1.0e0_r8p/Nwts
                    ! set 2nd cell and weight
                    if(pPushMesh%Face(iPushFace)%icell(2).eq.iPushCell) then
                        iPushCell = pPushMesh%Face(iPushFace)%icell(1)
                        Interp%ielem(iPullNode)%iPush(2) = iPushCell
                        Interp%ielem(iPullNode)%wPush(2) = 1.0e0_r8p/Nwts
                    elseif(pPushMesh%Face(iPushFace)%icell(2).lt.0) then   ! boundary case
                        Interp%ielem(iPullNode)%iPush(2) = iPushCell
                        Interp%ielem(iPullNode)%wPush(2) = 1.0e0_r8p/Nwts
                    else
                        iPushCell = pPushMesh%Face(iPushFace)%icell(2)
                        Interp%ielem(iPullNode)%iPush(2) = iPushCell
                        Interp%ielem(iPullNode)%wPush(2) = 1.0e0_r8p/Nwts
                    end if
!                    write(*,*) " case: Face",Nwts
                else
                    ncells = 1
                    Nwts = ncells        ! number of weights
                    call allocate_valsInterp_weight(Interp%ielem(iPullNode),Nwts)
                    ! set 1st cell and weight
                    Interp%ielem(iPullNode)%iPush(1) = iPushCell
                    Interp%ielem(iPullNode)%wPush(1) = 1.0e0_r8p/Nwts
!                    write(*,*) " case: Cell",Nwts
                end if
            else
                write(*,*) "warning: point is outside mesh"
            end if
        enddo

        write(*,*)

        return
    end subroutine new_interp_ValExtCellsToNodes

!    ! --- create new interpolator for derivatives ---
!    subroutine new_interpolator_for_derivs(dersInterp,pPushMesh,pPullMesh,sPushElem,sPullElem)
!        implicit none
!        ! input-output
!        type(TDerivsInterpolator) :: dersInterp
!        type(TMeshData), pointer  :: pPushMesh             ! source mesh
!        type(TMeshData), pointer  :: pPullMesh             ! accept mesh
!        character(*)              :: sPushElem             ! source mesh element type: "node"/"cell"
!        character(*)              :: sPullElem             ! accept mesh element type: "node"/"cell"
!        ! local
!!        integer :: NElems
!!        real(8) :: objSizeMb
!!        integer               :: i,icell,inei,iface,nnodes
!!        type(T3DPoint)    :: nodes(NumCellFaces+1)     ! coords of face nodes: size = 4 (quad) or 3 (tri)
!!        type(T3DPoint)    :: pt0
!!        real(r8p)             :: Fwt(NumCellFaces+1)       ! function weights to Fpt
!!        type(T3DPoint)    :: dFwt(NumCellFaces+1)      ! function weights to dFpt
!
!        write(*,*)
!        write(*,*) "....... new interpolator for gf derivatives ........"
!
!        ! set source mesh elements type
!        select case(trim(Lower_Case(sPushElem)))
!        case("node")
!            dersInterp%PushElemType = 0
!        case("edge")
!            dersInterp%PushElemType = 1
!        case("face")
!            dersInterp%PushElemType = 2
!        case("cell")
!            dersInterp%PushElemType = 3
!        end select
!        ! set acceptor mesh elements type
!        select case(trim(Lower_Case(sPullElem)))
!        case("node")
!            dersInterp%PullElemType = 0
!        case("edge")
!            dersInterp%PullElemType = 1
!        case("face")
!            dersInterp%PullElemType = 2
!        case("cell")
!            dersInterp%PullElemType = 3
!        end select
!
!
!        dersInterp%pPushMesh => pPushMesh
!        dersInterp%pPullMesh => pPullMesh
!
!!        allocate(dersInterp%ielem(pmesh%NCells))
!!
!!        call dersInterp_get_size(dersInterp,objSizeMb)
!!        write(*,"(a,1pE12.4e2,a)") " dersInterpSize = ",objSizeMb," Mb"
!!
!!        NElems = pmesh%NCells
!!!        open(54,file="InterpWeights.dat")
!!        do icell = 1,NElems
!!            write(*,"(a,a,I3,a)",advance="no") achar(13),' done ',icell*100/NElems," %"
!!            pt0%crd = pmesh%Cell(icell)%center
!!            nnodes = 1
!!            nodes(nnodes)%crd = pmesh%Cell(icell)%center
!!!            if(pt0%crd(3).gt.40.5 .and. pt0%crd(3).lt.50.5) then
!!!                write(54,*) " --------------- "
!!!                write(54,*) "icell = ",icell
!!!            endif
!!            !do i = 1,NumCellFaces
!!            do i = 1,pmesh%Cell(icell)%nFaces
!!                inei = pmesh%Cell(icell)%ineighbor(i)
!!                if(inei.gt.0) then
!!                    nnodes = nnodes + 1
!!                    nodes(nnodes)%crd = pmesh%Cell(inei)%center
!!                elseif(inei.lt.0) then   ! external face -> virtual (boundary) cell
!!                    nnodes = nnodes + 1
!!                    iface = pmesh%Cell(icell)%iface(i)
!!                    nodes(nnodes)%crd = pmesh%Face(iface)%center
!!                endif
!!            enddo
!!            call interp_least_squares(nnodes,nodes(1:nnodes),pt0,Fwt,dFwt)
!!
!!            nnodes = 1
!!            dersInterp%ielem(icell)%wt(0)%crd = dFwt(nnodes)%crd
!!            do i = 1,NumCellFaces
!!!            do i = 1,pmesh%Cell(icell)%nFaces
!!                dersInterp%ielem(icell)%wt(i)%crd(:) = 0.e0_r8p
!!                inei = pmesh%Cell(icell)%ineighbor(i)
!!                if(inei.ne.0) then
!!                    nnodes = nnodes + 1
!!                    dersInterp%ielem(icell)%wt(i)%crd = dFwt(nnodes)%crd
!!                endif
!!!                if(pt0%crd(3).gt.40.5 .and. pt0%crd(3).lt.50.5) then
!!!                    write(54,*) "dFwt: ",dFwt(nnodes)%crd
!!!                endif
!!            enddo
!!        enddo
!!!        close(54)
!!!        pause
!!        write(*,*)
!
!        return
!    end subroutine new_interpolator_for_derivs

! -------------------------------------------
    subroutine valsInterp_get_size(valsInterp,objSizeMb)
        implicit none
        ! input
        type(TValuesInterpolator) :: valsInterp
        ! output
        real(8)    :: objSizeBytes,objSizeMb
        ! local
        integer    :: i

        objSizeBytes = 0
        ! . . . . . . . size of data in interpolator root variables . . . . . . . .
        objSizeBytes = objSizeBytes + real(sizeof(valsInterp),8)
        ! . . . . . . . size of data in gfun arrays: . . . . . . . . .
        do i = 1,valsInterp%NElems
            objSizeBytes = objSizeBytes + real(sizeof(valsInterp%ielem(i)%iPush),8)
            objSizeBytes = objSizeBytes + real(sizeof(valsInterp%ielem(i)%wPush),8)
        end do
        ! .................................................................
        objSizeMb = objSizeBytes * 1.e-6
        return
    end subroutine

! -------------------------------------------
    subroutine dersInterp_get_size(dersInterp,objSizeMb)
        implicit none
        ! input
        type(TDerivsInterpolator) :: dersInterp
        ! output
        real(8)    :: objSizeBytes,objSizeMb
        ! local
        integer    :: i

        objSizeBytes = 0
        ! . . . . . . . size of data in interpolator root variables . . . . . . . .
        objSizeBytes = objSizeBytes + real(sizeof(dersInterp),8)
        ! . . . . . . . size of data in gfun arrays: . . . . . . . . .
        do i = 1,dersInterp%NElems
            objSizeBytes = objSizeBytes + real(sizeof(dersInterp%ielem(i)%iPush),8)
            objSizeBytes = objSizeBytes + real(sizeof(dersInterp%ielem(i)%wPush),8)
        end do
        ! .................................................................
        objSizeMb = objSizeBytes * 1.e-6
        return
    end subroutine


!    function interpolate_scalar_to_gradient(gfScl,dersInterp) result(gradScl)
!        use module_GridFunct
!        implicit none
!        ! input
!        type(TGridCellScalar)     :: gfScl
!        type(TDerivsInterpolator) :: dersInterp
!        ! output
!        type(TGridCellVector) :: gradScl
!        ! local
!        integer :: i,j
!        integer :: icell,inei
!        real(rPgfn) :: eps
!
!        write(*,*) "size(gfScl) = ",size(gfScl%fval)
!        do icell = 1,gfScl%pmesh%NCells
!            do i = 1,3
!                eps = gfScl%fval(icell)*dersInterp%ielem(icell)%wt(0)%crd(i)
!                do j = 1,NumCellFaces
!                    inei = gfScl%pmesh%Cell(icell)%ineighbor(j)
!                    if(inei.gt.0) then
!                        eps = eps + gfScl%fval(inei)*dersInterp%ielem(icell)%wt(j)%crd(i)
!                    endif
!                end do
!                gradScl%fvec(icell)%cmp(i) = eps
!            end do
!        end do
!
!        return
!    end function

! -------------------------------------------
    subroutine interp_arithm_mean(nnodes,Fwt)
        use module_Settings
        implicit none
        ! input parameters
        integer               :: nnodes
        ! output variables
        real(r8p)           :: Fwt(nnodes)

        Fwt(:) = 1.e0_r8p/nnodes

        return
    end subroutine interp_arithm_mean

! -------------------------------------------
!    subroutine interp_bilinear(nnodes,nodes,pt0,Fwt,dFwt)
! --------------------------------------------------------------------------
    subroutine interp_bilinear(nnodes,nodes,pt,Fwt,dFwt,istop)
        use module_Settings
        use module_Geom
        implicit none
        ! input parameters
        integer               :: nnodes
        type(T3DPoint)    :: nodes(nnodes)     ! coords of face nodes: size = 4 (quad) or 3 (tri)
        type(T3DPoint)    :: pt
        ! output variables
        real(r8p)             :: Fwt(nnodes)       ! function weights to Fpt
        type(T3DPoint)    :: dFwt(nnodes)      ! function weights to dFpt
        ! local variables
        real(r8p), parameter :: val0 = 1.e-6_r8p
        type(T3DPoint), dimension(4) :: quad      ! variable for quad (or tri: if quad(4)=quad(1))
!        real(r8p),          dimension(4) :: fnval     ! variable for quad (or tri: if quad(4)=quad(1))
        type(T3DPoint) :: pt1,pt2,pt3,pt4,ptc
        type(T3DPoint) :: e12,e23,e43,e14,e1p
        type(T3DPoint) :: e_nf
        type(T3DPoint) :: e_p,e_c
        type(T3DPoint) :: e_t,e_s,e_st,e_nst
        real(r8p)          :: elen,eps
        real(r8p)          :: a_p,a_s,a_t,a_st
        real(r8p)          :: b_p,b_s,b_t,b_st
        real(r8p)          :: a,b,c,D
        real(r8p)          :: s,t
!        real(r8p)          :: f12,f43,f14,f23
        real(r8p)            :: rcase
        integer              :: i,istop


        quad(1:nnodes) = nodes(1:nnodes)
        if(nnodes.eq.2) then              ! case of line
            quad(3) = quad(2)
            quad(4) = quad(1)
        elseif(nnodes.eq.3) then              ! case of triangle
            quad(4) = quad(1)
        elseif(nnodes.ne.4) then
            write(*,"(a,i3)") 'error - interp_bilinear: unsupported element type with nnodes = ',nnodes
            stop
        end if
        istop = 0

!        pt1 = nodes(1)
!        pt2 = nodes(2)
!        pt3 = nodes(3)
!        pt4 = nodes(4)
        pt1 = quad(1)
        pt2 = quad(2)
        pt3 = quad(3)
        pt4 = quad(4)
        ptc%crd = 0.25_r8p*(pt1%crd+pt2%crd+pt3%crd+pt4%crd)

        e12%crd = pt2%crd - pt1%crd
        e23%crd = pt3%crd - pt2%crd
        e43%crd = pt3%crd - pt4%crd
        e14%crd = pt4%crd - pt1%crd
        e1p%crd = pt%crd  - pt1%crd
        e_c%crd = ptc%crd - pt1%crd
        e_nf = vectVectProd(e12,e23)
        e_nf%crd = e_nf%crd/vectLength(e_nf)

        e_p%crd  = e1p%crd
        e_t%crd  = e12%crd
        e_s%crd  = e14%crd
        e_st%crd = e43%crd-e12%crd

        elen = vectLength(e_st)/vectLength(e12)
        if(elen.gt.val0) then
            a_p = vectScalProd(e_p,e_t)
            a_t = vectScalProd(e_t,e_t)
            a_s = vectScalProd(e_s,e_t)
            a_st = vectScalProd(e_st,e_t)

            e_nst = vectVectProd(e_st,e_nf)
            e_nst%crd = e_nst%crd/vectLength(e_nst)
            b_p = vectScalProd(e_p,e_nst)
            b_t = vectScalProd(e_t,e_nst)
            b_s = vectScalProd(e_s,e_nst)
            b_st = vectScalProd(e_st,e_nst)       ! by definition = 0

            eps = vectLength(e12)
            if(abs(b_s).gt.val0*eps) then
                a = a_st*b_t
                b = a_s*b_t-a_t*b_s-a_st*b_p
                c = a_p*b_s-a_s*b_p
                D = b**2 - 4.0_r8p*a*c
                if(abs(a).gt.val0) then
                    rcase = 1.11
                    t = 0.5_r8p*(-b+dsqrt(D))/a
                    if(t.gt.1.e0_r8p .or. t.lt.0.e0_r8p) then
                        rcase = 1.12
                        t = 0.5_r8p*(-b-dsqrt(D))/a
                    end if
                else
                    rcase = 1.2
                    t = -c/b
                end if
                s = (b_p-b_t*t)/b_s
            elseif(abs(b_t).gt.val0*eps) then
                a = a_st*b_s
                b = a_t*b_s-a_s*b_t-a_st*b_p
                c = a_p*b_t-a_t*b_p
                D = b**2 - 4.0_r8p*a*c
                if(abs(a).gt.val0) then
                    rcase = 2.11
                    s = 0.5_r8p*(-b+dsqrt(D))/a
                    if(s.gt.1.e0_r8p .or. s.lt.0.e0_r8p) then
                        rcase = 2.12
                        s = 0.5_r8p*(-b-dsqrt(D))/a
                    end if
                else
                    rcase = 2.2
                    s = -c/b
                end if
                t = (b_p-b_s*s)/b_t
            else
                rcase = 3.0
                t = b_p/b_t
                s = (a_p-a_t*t)/(a_s+a_st*t)
            endif
        else
            a_p = vectScalProd(e_p,e_t)
            a_t = vectScalProd(e_t,e_t)
            a_s = vectScalProd(e_s,e_t)

            b_p = vectScalProd(e_p,e_s)
            b_t = vectScalProd(e_t,e_s)
            b_s = vectScalProd(e_s,e_s)

            s = (a_t*b_p-a_p*b_t)/(a_t*b_s-a_s*b_t)
            t = (a_p-a_s*s)/a_t
        endif

        ! function weights
        Fwt(1) = (1.0_r8p - t)*(1.0_r8p - s)
        Fwt(2) = t*(1.0_r8p - s)
        Fwt(3) = t*s
        Fwt(4) = s*(1.0_r8p - t)

        ! remove small values (max >= 0.25)
        do i = 1,4
            if(abs(Fwt(i)).lt.val0) then
                Fwt(i) = 0.e0_r8p
            end if
        end do

        if(ANY(Fwt(:).lt.0)) then
!            ... check for small negative values (exclude) ...
            write(*,"(a,4(1pE12.4e2))")
            write(*,"(a,f4.2)")         "case:  ",rcase
            write(*,"(a,4(1pE12.4e3))") "a    = ",a
            write(*,"(a,4(1pE12.4e3))") "b    = ",b
            write(*,"(a,4(1pE12.4e3))") "c    = ",c
            write(*,"(a,4(1pE12.4e3))") "D    = ",D
            write(*,"(a,4(1pE12.4e3))") "s    = ",s
            write(*,"(a,4(1pE12.4e3))") "t    = ",t
            write(*,"(a,4(1pE12.4e3))") "Fwt:   ",Fwt(:)
            write(*,"(a,4(1pE12.4e3))") "elen = ",elen
            write(*,"(a,4(1pE12.4e3))") "len(e_st) = ",vectLength(e_st)
            write(*,"(a,4(1pE12.4e3))") "len(e12)  = ",vectLength(e12)
            write(*,"(a,4(1pE12.4e3))") "b_s  = ",b_s
            write(*,"(a,4(1pE12.4e3))") "b_t  = ",b_t
            stop

            istop = 1
        end if

        ! gradient calculation
        e_s%crd = e12%crd + (e23%crd-e14%crd)*s
        e_s%crd = e_s%crd/vectLength(e_s)
        e_t%crd = e14%crd + (e43%crd-e12%crd)*t
        e_t%crd = e_t%crd/vectLength(e_t)
        !dfpt%crd = (f43-f12)*e_s%crd + (f23-f14)*e_t%crd

        do i = 1,nnodes
            dFwt(i)%crd = 0.0_r8p
        enddo
!        write(*,"(a)") " - - - - - - - - - - - - "
!        write(*,"(2(1pE12.4e2),1pE18.10e2)") 0.e0,0.e0,fnval(1)
!        write(*,"(2(1pE12.4e2),1pE18.10e2)") 1.e0,0.e0,fnval(2)
!        write(*,"(2(1pE12.4e2),1pE18.10e2)") 1.e0,1.e0,fnval(3)
!        write(*,"(2(1pE12.4e2),1pE18.10e2)") 0.e0,1.e0,fnval(4)
!        write(*,"(2(1pE12.4e2),1pE18.10e2)")
!        write(*,"(2(1pE12.4e2),1pE18.10e2)") t,s,fpt

!        open(24,file="bilin_point.dat")
!        write(24,"(4(1pE12.4e2))") 0.e0,0.e0,fnval(1)
!        write(24,"(4(1pE12.4e2))") 1.e0,0.e0,fnval(2)
!        write(24,"(4(1pE12.4e2))") 1.e0,1.e0,fnval(3)
!        write(24,"(4(1pE12.4e2))") 0.e0,1.e0,fnval(4)
!        write(24,"(4(1pE12.4e2))")
!        write(24,"(4(1pE12.4e2))") t,s,fpt
!        close(24)
!        b = t
!        c = s
!
!        open(24,file="bilin_surf.dat")
!        s = 0.e0_r8p
!        do while(s.le.1.0e0_r8p+val0)
!            t = 0.e0_r8p
!            do while(t.le.1.0e0_r8p+val0)
!                f12 = (fnval(2)-fnval(1))*t+fnval(1)
!                f43 = (fnval(3)-fnval(4))*t+fnval(4)
!                a = f12 + (f43-f12)*s
!                write(24,"(4(1pE12.4e2))") t,s,a
!                t = t + 0.1e0_r8p
!            end do
!            s = s + 0.1e0_r8p
!            write(24,"(4(1pE12.4e2))")
!        end do
!
!        t = b
!        s = c
!!        write(24,"(4(1pE12.4e2))") t,s,fpt
!        f12 = (fnval(2)-fnval(1))*t+fnval(1)
!        f43 = (fnval(3)-fnval(4))*t+fnval(4)
!        f14 = (fnval(4)-fnval(1))*s+fnval(1)
!        f23 = (fnval(3)-fnval(2))*s+fnval(2)
!        !a = f12 + (f43-f12)*s
!        !a = f14 + (f23-f14)*t
!        b = (f23-f14)     ! df/dt
!        c = (f43-f12)     ! df/ds
!        D = dsqrt(b**2+c**2)
!        b = 0.5*b/D
!        c = 0.5*c/D
!        a = fpt + (f23-f14)*b + (f43-f12)*c
!
!        write(24,"(4(1pE12.4e2))") t+b,s+c,a
!        close(24)
!        pause

        return
    end subroutine interp_bilinear

! ---------------------------------------------------------------
    subroutine interp_least_squares(nnodes,nodes,pt0,Fwt,dFwt)
        use module_Settings
        use module_Geom
        implicit none
        ! input parameters
        integer              :: nnodes
        type(T3DPoint)   :: nodes(nnodes)     ! coords of face nodes: size = 4 (quad) or 3 (tri)
        type(T3DPoint)   :: pt0
        ! output variables
        real(r8p)            :: Fwt(nnodes)       ! function weights to Fpt
        type(T3DPoint)   :: dFwt(nnodes)      ! function weights to dFpt
        ! local variables
        real(r8p), parameter   :: val0 = 1.e-8_r8p
        type(T3DPoint)     :: pti,dri
        real(r8p)              :: Dr1(3)
        real(r8p)              :: Dr2(3,3)
        real(8)                :: Amtx(4,4),Rmtx(4,4)
        real(8)                :: bvec(4)
        real(r8p)              :: Smtx(4,4)
        integer                :: inode,i,j
        integer                :: ios

        Dr1(:)   = 0.e0_r8p
        Dr2(:,:) = 0.e0_r8p

        do inode = 1,nnodes
            pti = nodes(inode)
            dri%crd = pti%crd-pt0%crd
            do i = 1,3
                Dr1(i) = Dr1(i) + dri%crd(i)
                do j = 1,3
                    Dr2(i,j) = Dr2(i,j) + dri%crd(i)*dri%crd(j)
                enddo
            enddo
        enddo

        ! set full matrix
        do i = 1,3
            Amtx(i,4) = real(Dr1(i),8)
            Amtx(4,i) = real(Dr1(i),8)
            do j = 1,3
                Amtx(i,j) = real(Dr2(i,j),8)
            enddo
        enddo
        Amtx(4,4) = real(nnodes,8)
        ! set RHS vector
        bvec(:) = 1.0d0

!        do i = 1,3
!            if(.not.any(Amtx(i,:).ne.0.e0)) then
!                write(*,"(a,i1,a,10(1x,1pE12.4e3))") " Amtx (",i,",:) = ",Amtx(i,:)
!            end if
!        enddo

        Rmtx = Amtx                    ! reverse matrix initial values
        call gaussj(Rmtx,4,bvec,1,ios)
        Smtx = real(Rmtx,r8p)

        if(ios.gt.0) then
            write(*,"(a,10(1x,1pE12.4e3))")
            do i = 1,nnodes
                write(*,"(a,10(1x,1pE12.4e3))") " nodes: ",nodes(i)%crd
            enddo
            write(*,"(a,10(1x,1pE12.4e3))")
            write(*,"(a,10(1x,1pE12.4e3))") "    |",Amtx(1,:)
            write(*,"(a,10(1x,1pE12.4e3))") "    |",Amtx(2,:)
            write(*,"(a,10(1x,1pE12.4e3))") " A= |",Amtx(3,:)
            write(*,"(a,10(1x,1pE12.4e3))") "    |",Amtx(4,:)
            stop
        end if

!        open(15,file="Amtx-Rmtx.dat")
!        do i = 1,4
!            write(15,"(10(1x,1pE12.4e3))") Amtx(i,:)
!        end do
!        write(15,"(10(1x,1pE12.4e3))")
!        write(15,"(10(1x,1pE12.4e3))")
!        do i = 1,4
!            write(15,"(10(1x,1pE12.4e3))") Rmtx(i,:)
!        end do
!        close(15)


        do inode = 1,nnodes
            pti = nodes(inode)
            dri%crd = pti%crd-pt0%crd
            Fwt(inode)      = 0.e0_r8p
            dFwt(inode)%crd = 0.e0_r8p
            do i = 1,3
                Fwt(inode) = Fwt(inode) + Smtx(4,i)*dri%crd(i)
                do j = 1,3
                    dFwt(inode)%crd(j) = dFwt(inode)%crd(j) + Smtx(j,i)*dri%crd(i)
                enddo
            enddo
            Fwt(inode) = Fwt(inode) + Smtx(4,4)
            do j = 1,3
                dFwt(inode)%crd(j) = dFwt(inode)%crd(j) + Smtx(j,4)
            enddo
        end do

        return
    end subroutine

    ! mapper
    subroutine mapper_XYZ_to_RZ_dirX(pt)
        use module_Settings
        implicit none
        ! input/output:
        real(r8p), dimension(3) :: pt      ! point coordinates to convert
        ! local:
        real(r8p), dimension(3) :: pt0

        pt0 = pt
        pt(1) = pt0(1)
        pt(2) = sqrt(pt0(2)**2+pt0(3)**2)
        pt(3) = 0.0e0_r8p
        return
    end subroutine

    ! mapper
    subroutine mapper_XY_to_X(pt)
        use module_Settings
        implicit none
        ! input/output:
        real(r8p), dimension(3) :: pt      ! point coordinates to convert
        ! local:
        real(r8p), dimension(3) :: pt0

        pt0 = pt
        pt(1) = pt0(1)
        pt(2) = 0.0e0_r8p
        pt(3) = 0.0e0_r8p
        return
    end subroutine

! -------------------------------------------------------
end module module_interpolate


