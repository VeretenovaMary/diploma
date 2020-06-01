subroutine new_mesh_simple1D(NX,X0,X1,mesh)
    use module_Settings
    use module_Mesh
    use module_utils
    implicit none
    ! input:
    integer         :: NX      ! number of mesh nodes
    real(r8p)       :: X0      ! left bound coordinate
    real(r8p)       :: X1      ! right bound coordinate
    ! output:
    type(TMeshData) :: mesh
    ! local:
    integer      :: i,k
    integer      :: locIdx,locPos
    integer      :: ncells,nfaces,nnodes
    real(r8p)    :: dX      ! coordinate step



    ! set number of elements (nodes, faces, cells):
    NMeshes      = NMeshes + 1                       ! add new mesh (counter)
    mesh%meshIdx = NMeshes                           ! set mesh index
    mesh%meshTit = "simp1D"                          ! mesh title
    mesh%meshDim = 1                                 ! 1D mesh
    mesh%SymType = 0                                 ! native (default)
    mesh%SymAxis = [1.0e0_r8p,0.0e0_r8p,0.0e0_r8p]   ! native (default)
    mesh%NNodes  = NX
    mesh%NCells  = NX-1
    mesh%NBFaces = 2
    mesh%NFaces  = NX
    allocate(mesh%Node(mesh%NNodes))           ! mesh setup (allocation NNodes)
    allocate(mesh%Cell(mesh%NCells))           ! mesh setup (allocation)
    allocate(mesh%Face(mesh%NFaces))           ! mesh setup (allocation)

    ! set nodes:
    write(*,*) "set nodes:"
    dX = (X1-X0)/(NX-1)
    do i = 1,NX
        mesh%Node(i)%crd(:) = [dX*(i-1),0.0e0_r8p,0.0e0_r8p]      ! node coordinate
        ncells = 2                                                ! number of cells adjacent to current node
        mesh%Node(i)%ncells = ncells                              ! number of cells adjacent to current node
        call allocate_node(mesh%Node(i),ncells)               ! allocate memory for adjcent cells information
        mesh%Node(i)%icell(1) = i-1
        mesh%Node(i)%icell(2) = i
        if(i.eq.1) then                      ! if node is left boundary node
            mesh%Node(i)%icell(1) =  i
            mesh%Node(i)%icell(2) = -i
        end if
        if(i.eq.NX) then                     ! if node is right boundary node
            mesh%Node(i)%icell(1) =   i-1
            mesh%Node(i)%icell(2) = -(i-1)
        end if
    end do


    ! set faces:
    write(*,*) "set faces:"
    do i = 1,NX
        nnodes = 1                               ! number of nodes that defines the face
        mesh%Face(i)%nnodes = nnodes             ! number of nodes that defines the face
        call allocate_face(mesh%Face(i),nnodes)  ! allocate memory for face nodes information
        mesh%Face(i)%inode(1) = i
        mesh%Face(i)%itype    = 1                ! point defines the face
        mesh%Face(i)%nedges   = 0                ! number of edges in face
        mesh%Face(i)%icell(:) = 0
        mesh%Face(i)%icell(1) = i-1
        mesh%Face(i)%icell(2) = i
        if(i.eq.1) then                          ! if node is left boundary node
            mesh%Face(i)%icell(1) =  i
            mesh%Face(i)%icell(2) = -i
        end if
        if(i.eq.NX) then                     ! if node is right boundary node
            mesh%Face(i)%icell(1) =   i-1
            mesh%Face(i)%icell(2) = -(i-1)
        end if
    end do


    ! set cells:
    write(*,*) "set cells:"
    do i = 1,NX-1
        nnodes = 2                               ! number of nodes that defines the cell
        mesh%Cell(i)%nnodes = nnodes             ! number of nodes that defines the cell
        nfaces = 2                               ! number of faces that defines the cell
        mesh%Cell(i)%nfaces = nfaces             ! number of faces that defines the cell
        call allocate_cell(mesh%Cell(i),nnodes,nfaces)  ! allocate memory for cell faces and nodes information
        mesh%Cell(i)%inode(1) = i
        mesh%Cell(i)%inode(2) = i+1
        mesh%Cell(i)%iface(1) = i
        mesh%Cell(i)%iface(2) = i+1
    end do

    write(*,*) "set locs:"
    ! set locations data:
    mesh%NLocs = 3
    mesh%LocNBnds = 2
    mesh%LocBndIdx(1) = 1
    mesh%LocBndTit(1) = "inlet"
    mesh%LocBndIdx(2) = 2
    mesh%LocBndTit(2) = "outlet"
    mesh%LocNVols = 1
    mesh%LocVolIdx(1) = 3
    mesh%LocVolTit(1) = "channel"
    ! set node locations:
    do i = 1,NX
        mesh%Node(i)%locIds(:) = 0
        mesh%Node(i)%locIds(1) = 3
        if(i.eq.1) then
            mesh%Node(i)%locIds(2) = 1   ! "inlet" location index
        end if
        if(i.eq.NX) then
            mesh%Node(i)%locIds(2) = 2   ! "outlet" location index
        end if
    end do
    ! set face locations:
    do i = 1,NX
        mesh%Face(i)%locIdx = 0
        if(i.eq.1) then
            mesh%Face(i)%locIdx = 1   ! "inlet" location index
        end if
        if(i.eq.NX) then
            mesh%Face(i)%locIdx = 2   ! "outlet" location index
        end if
    end do
    ! set cell locations:
    do i = 1,NX-1
        mesh%Cell(i)%locIdx = 3       ! channel
    end do

    write(*,*) "mesh_1DMeshCalcul:"
    call mesh_1DMeshCalcul(mesh)

    ! set boundary faces:
    !write(*,*) "set boundary faces:"
    allocate(mesh%BFaces(mesh%LocNBnds))        ! 2 locations: inlet, outlet
    mesh%BFaces(1)%nfaces = 1       ! inlet
    mesh%BFaces(2)%nfaces = 1       ! outlet
    do i = 1,mesh%LocNBnds
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

    return
end subroutine

