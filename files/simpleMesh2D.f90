! ----------------------------------------------------------------------------
subroutine new_mesh_simple2D(NZ,NR,Rint,Rext,DL0,mesh)
    use module_Settings
    use module_Mesh
    implicit none
    ! input:
    integer                          :: NZ,NR
    real(r8p),      dimension(NZ)    :: Rint,Rext      ! line point coordinates
    real(r8p)                        :: DL0
    ! output:
    type(TMeshData)                  :: mesh
    ! local:
    type(T3DPoint), dimension(NZ)    :: LineBot,LineTop      ! line point coordinates
    type(T3DPoint), dimension(NR)    :: LineLef,LineRig      ! line point coordinates
    integer :: i

    do i = 1,NZ
        LineBot(i)%crd(1) = DL0*(i-1)/(NZ-1)
        LineBot(i)%crd(2) = Rint(i)
        LineBot(i)%crd(3) = 0.0d0
        LineTop(i)%crd(1) = DL0*(i-1)/(NZ-1)
        LineTop(i)%crd(2) = Rext(i)
        LineTop(i)%crd(3) = 0.0d0
    end do

    do i = 1,NR
        LineLef(i)%crd(1) = 0.d0
        LineLef(i)%crd(2) = (Rext(1)-Rint(1))*(i-1)/(NR-1) + Rint(1)
        LineLef(i)%crd(3) = 0.0d0
        LineRig(i)%crd(1) = DL0
        LineRig(i)%crd(2) = (Rext(NZ)-Rint(NZ))*(i-1)/(NR-1) + Rint(NZ)
        LineRig(i)%crd(3) = 0.0d0
    end do

    call build_mesh_simple2D(NZ,NR,LineBot,LineTop,LineLef,LineRig,mesh)

    return
end subroutine new_mesh_simple2D


! ----------------------------------------------------------------------------
subroutine build_mesh_simple2D(Nx,Ny,LineBot,LineTop,LineLef,LineRig,mesh)
    use module_Settings
    use module_Mesh
    use module_OutputMesh
    use module_utils
    implicit none
    ! input:
    integer                          :: Nx,Ny
    type(T3DPoint), dimension(Nx)    :: LineBot,LineTop      ! line point coordinates
    type(T3DPoint), dimension(Ny)    :: LineLef,LineRig      ! line point coordinates
    type(T3DPoint), dimension(Nx,Ny) :: Point                ! points coordinates
    ! output:
    type(TMeshData)                  :: mesh
    ! local:
    integer        :: i,j,locIdx,locPos
    integer        :: nn,nc,ne
    ! x - attributes
    real(r8p)        :: dx0,dx1,dxc
    real(r8p)        :: Lx0,Lx1,Lxc
    real(r8p)        :: kx0,kx1,kxc
    ! x - attributes
    real(r8p)        :: dy0,dy1,dyc
    real(r8p)        :: Ly0,Ly1,Lyc
    real(r8p)        :: ky0,ky1,kyc
    character(256) :: vtkTit
    character(16)  :: vtkFmt
    integer        :: NNodes,NFaces,NBFaces,NCells

    NNodes  = Nx*Ny
    NCells  = (Nx-1)*(Ny-1)
    NFaces  = 2*(Nx-1)*(Ny-1) + (Nx-1) + (Ny-1)
    NBFaces = 2*(Nx-1) + 2*(Ny-1)

    ! set number of elements (nodes, faces, cells):
    NMeshes      = NMeshes + 1                       ! add new mesh (counter)
    mesh%meshIdx = NMeshes                           ! set mesh index
    if(len_trim(mesh%meshTit).eq.0) then
        mesh%meshTit = "simp2D"                      ! mesh title
    endif
    mesh%meshDim = 2                                 ! 2D mesh
    mesh%SymType = 2                                 ! native (default)
    mesh%SymAxis = [1.0e0_r8p,0.0e0_r8p,0.0e0_r8p]   ! native (default)
    mesh%NNodes  = NNodes
    mesh%NFaces  = NFaces
    mesh%NCells  = NCells
    mesh%NBFaces = NBFaces
    allocate(mesh%Node(mesh%NNodes))           ! mesh setup (allocation NNodes)
    allocate(mesh%Cell(mesh%NCells))           ! mesh setup (allocation)
    allocate(mesh%Face(mesh%NFaces))           ! mesh setup (allocation)
    do ne = 1,mesh%NFaces
        mesh%Face(ne)%nnodes = 2
        call allocate_face(mesh%Face(ne),2)
    enddo
    do nc = 1,mesh%NCells
        mesh%Cell(nc)%nnodes = 4
        mesh%Cell(nc)%nfaces = 4
        call allocate_cell(mesh%Cell(nc),4,4)
    enddo

    ! set grid points
    nn = 0
    Lx0 = LineBot(Nx)%crd(1) - LineBot(1)%crd(1)
    Lx1 = LineTop(Nx)%crd(1) - LineTop(1)%crd(1)
    Ly0 = LineLef(Ny)%crd(2) - LineLef(1)%crd(2)
    Ly1 = LineRig(Ny)%crd(2) - LineRig(1)%crd(2)
    ! set bottom boundary points -> points array
    do i = 1,Nx
        Point(i,1) = LineBot(i)
    end do
    ! set left boundary points -> points array
    do j = 1,Ny
        Point(1,j) = LineLef(j)
    end do

    do i = 2,Nx
        dx0 = LineBot(i)%crd(1) - LineBot(i-1)%crd(1)
        kx0 = dx0/Lx0
        dx1 = LineTop(i)%crd(1) - LineTop(i-1)%crd(1)
        kx1 = dx1/Lx1

        do j = 2,Ny
            kxc = kx0 + (kx1-kx0)/(Ny-1)*(j-1)
            Lxc = LineRig(j)%crd(1) - LineLef(j)%crd(1)
            dxc = Lxc*kxc

            dy0 = LineLef(j)%crd(2) - LineLef(j-1)%crd(2)
            ky0 = dy0/Ly0
            dy1 = LineRig(j)%crd(2) - LineRig(j-1)%crd(2)
            ky1 = dy1/Ly1

            kyc = ky0 + (ky1-ky0)/(Nx-1)*(i-1)
            Lyc = LineTop(i)%crd(2) - LineBot(i)%crd(2)
            dyc = Lyc*kyc

            Point(i,j)%crd(1) = Point(i-1,j)%crd(1) + dxc
            Point(i,j)%crd(2) = Point(i,j-1)%crd(2) + dyc
            Point(i,j)%crd(3) = 0.0e0
        enddo
    enddo

    ! set grid nodes (from points array)
    nn = 0
    do i = 1,Nx
        do j = 1,Ny
            nn = nn + 1
            mesh%Node(nn)%crd(1) = Point(i,j)%crd(1)
            mesh%Node(nn)%crd(2) = Point(i,j)%crd(2)
            mesh%Node(nn)%crd(3) = 0.0e0
            mesh%Node(nn)%locIds(:) = 0
            mesh%Node(nn)%locIds(1) = 5
        enddo
    enddo
!    open(1,file="points.txt")
!    do i = 1,NNodes
!        write(1,"(2(1x,1pE12.5e2))") mesh%Node(i)%crd(:)
!    enddo
!    close(1)


    ! set bound faces:
    NBFaces = 2*(Nx-1) + 2*(Ny-1)
    ne = 0
    do i = 1,Nx-1     ! bottom bound
        ne = ne + 1
        mesh%Face(ne)%locIdx = 1               ! "bottom"
        mesh%Face(ne)%itype  = 2               ! line
        mesh%Face(ne)%inode(1) = Ny*(i-1) + 1
        mesh%Face(ne)%inode(2) = Ny*(i-0) + 1
        mesh%Face(ne)%icell(1) = (Ny-1)*(i-1) + 1
        mesh%Face(ne)%icell(2) = - ne
        mesh%Face(ne)%bndType = 1

        nc = (Ny-1)*(i-1) + 1
        mesh%Cell(nc)%iface(1)     = ne
        mesh%Cell(nc)%iFaceDir     = 1
        mesh%Cell(nc)%inode(1)     = Ny*(i-1) + 1
        mesh%Cell(nc)%inode(2)     = Ny*(i-0) + 1
        mesh%Cell(nc)%ineighbor(1) = - ne
        mesh%Cell(nc)%locIdx       = 5            ! "channel"
        mesh%Cell(nc)%itype        = 4            ! quad
    enddo
    do j = 1,Ny-1     ! right bound
        ne = ne + 1
        mesh%Face(ne)%locIdx = 2               ! "right"
        mesh%Face(ne)%itype  = 2               ! line
        mesh%Face(ne)%inode(1)  = Ny*(Nx-1) + (j+0)
        mesh%Face(ne)%inode(2)  = Ny*(Nx-1) + (j+1)
        mesh%Face(ne)%icell(1)  = (Nx-2)*(Ny-1) + j
        mesh%Face(ne)%icell(2)  = - ne

        nc = (Nx-2)*(Ny-1) + j
        mesh%Cell(nc)%iface(2)     = ne
        mesh%Cell(nc)%iFaceDir     = 1
        mesh%Cell(nc)%inode(2)     = Ny*(Nx-1) + (j+0)
        mesh%Cell(nc)%inode(3)     = Ny*(Nx-1) + (j+1)
        mesh%Cell(nc)%ineighbor(2) = - ne
        mesh%Cell(nc)%locIdx       = 5            ! "channel"
        mesh%Cell(nc)%itype        = 4            ! quad
    enddo
    do i = 1,Nx-1     ! top bound
        ne = ne + 1
        mesh%Face(ne)%locIdx = 3                  ! "top"
        mesh%Face(ne)%inode(1)  = Ny*(Nx-i+1)
        mesh%Face(ne)%inode(2)  = Ny*(Nx-i+0)
        mesh%Face(ne)%icell(1)  = (Ny-1)*(Nx-1-(i-1))
        mesh%Face(ne)%icell(2)  = - ne

        nc = (Ny-1)*(Nx-1-(i-1))
        mesh%Cell(nc)%iface(3)     = ne
        mesh%Cell(nc)%iFaceDir     = 1
        mesh%Cell(nc)%inode(3)     = Ny*(Nx-i+1)
        mesh%Cell(nc)%inode(4)     = Ny*(Nx-i+0)
        mesh%Cell(nc)%ineighbor(3) = - ne
        mesh%Cell(nc)%locIdx       = 5            ! "channel"
        mesh%Cell(nc)%itype        = 4            ! quad
    enddo
    do j = 1,Ny-1     ! left bound
        ne = ne + 1
        mesh%Face(ne)%locIdx   = 4                ! "left"
        mesh%Face(ne)%inode(1)  = Ny - (j-1)
        mesh%Face(ne)%inode(2)  = Ny - (j-0)
        mesh%Face(ne)%icell(1)  = Ny - j
        mesh%Face(ne)%icell(2)  = - ne

        nc = Ny - j
        mesh%Cell(nc)%iface(4)     = ne
        mesh%Cell(nc)%iFaceDir     = 1
        mesh%Cell(nc)%inode(4)     = Ny - (j-1)
        mesh%Cell(nc)%inode(1)     = Ny - (j-0)
        mesh%Cell(nc)%ineighbor(4) = - ne
        mesh%Cell(nc)%locIdx       = 5            ! "channel"
        mesh%Cell(nc)%itype        = 4            ! quad
    enddo
    do i = 1,Nx-2     ! internal edges
        do j = 1,Ny-2
            ne = ne + 1        ! vertical edge
            mesh%Face(ne)%locIdx   = 0                ! internal
            mesh%Face(ne)%inode(1) = Ny*(i-0) + (j+0)
            mesh%Face(ne)%inode(2) = Ny*(i-0) + (j+1)
            mesh%Face(ne)%icell(1) = (Ny-1)*(i-1) + j
            mesh%Face(ne)%icell(2) = (Ny-1)*(i-0) + j
            !   define left cell
            nc = (Ny-1)*(i-1) + j  ! left cell (cell-1) for vertical edge
            mesh%Cell(nc)%iface(2)     = ne
            mesh%Cell(nc)%iFaceDir     = 1
            mesh%Cell(nc)%inode(2)     = Ny*(i-0) + (j+0)
            mesh%Cell(nc)%inode(3)     = Ny*(i-0) + (j+1)
            mesh%Cell(nc)%ineighbor(2) = (Ny-1)*(i-0) + j
            mesh%Cell(nc)%locIdx       = 5                 ! "channel"
            mesh%Cell(nc)%itype        = 4            ! quad
            !   define right cell
            nc = (Ny-1)*(i-0) + j  ! right cell (cell-2) for vertical edge
            mesh%Cell(nc)%iface(4)     = ne
            mesh%Cell(nc)%iFaceDir     = -1
            mesh%Cell(nc)%inode(4)     = Ny*(i-0) + (j+1)
            mesh%Cell(nc)%inode(1)     = Ny*(i-0) + (j+0)
            mesh%Cell(nc)%ineighbor(4) = (Ny-1)*(i-1) + j
            mesh%Cell(nc)%locIdx       = 5                 ! "channel"
            mesh%Cell(nc)%itype        = 4                 ! quad

            ne = ne + 1        ! horizontal edge
            mesh%Face(ne)%locIdx   = 0                ! internal
            mesh%Face(ne)%inode(1) = Ny*(i-0) + (j+1)
            mesh%Face(ne)%inode(2) = Ny*(i-1) + (j+1)
            mesh%Face(ne)%icell(1) = (Ny-1)*(i-1) + j
            mesh%Face(ne)%icell(2) = (Ny-1)*(i-1) + j+1
            !   define bottom cell
            nc = (Ny-1)*(i-1) + j    ! bottom cell (cell-1) for horizontal edge
            mesh%Cell(nc)%iface(3)     = ne
            mesh%Cell(nc)%iFaceDir     = 1
            mesh%Cell(nc)%inode(3)     = Ny*(i-0) + (j+1)
            mesh%Cell(nc)%inode(4)     = Ny*(i-1) + (j+1)
            mesh%Cell(nc)%ineighbor(3) = (Ny-1)*(i-1) + j+1
            mesh%Cell(nc)%locIdx       = 5            ! "channel"
            mesh%Cell(nc)%itype        = 4            ! quad
            !   define top cell
            nc = (Ny-1)*(i-1) + j+1  ! right cell (cell-2) for horizontal edge
            mesh%Cell(nc)%iface(1)     = ne
            mesh%Cell(nc)%iFaceDir     = -1
            mesh%Cell(nc)%inode(1)     = Ny*(i-1) + (j+1)
            mesh%Cell(nc)%inode(2)     = Ny*(i-0) + (j+1)
            mesh%Cell(nc)%ineighbor(1) = (Ny-1)*(i-1) + j
            mesh%Cell(nc)%locIdx       = 5                 ! "channel"
            mesh%Cell(nc)%itype        = 4            ! quad
        enddo
            ! for top row
            ne = ne + 1        ! vertical edge
            mesh%Face(ne)%locIdx   = 0                ! internal
            mesh%Face(ne)%inode(1) = Ny*(i-0) + (j+0)
            mesh%Face(ne)%inode(2) = Ny*(i-0) + (j+1)
            mesh%Face(ne)%icell(1) = (Ny-1)*(i+0)
            mesh%Face(ne)%icell(2) = (Ny-1)*(i+1)
            !   define left cell
            nc = (Ny-1)*(i+0)       ! left cell (cell-1) for vertical edge
            mesh%Cell(nc)%iface(2)     = ne
            mesh%Cell(nc)%iFaceDir     = 1
            mesh%Cell(nc)%inode(2)     = Ny*(i-0) + (j+0)
            mesh%Cell(nc)%inode(3)     = Ny*(i-0) + (j+1)
            mesh%Cell(nc)%ineighbor(2) = (Ny-1)*(i+1)
            mesh%Cell(nc)%locIdx       = 5                  ! "channel"
            mesh%Cell(nc)%itype        = 4            ! quad
            !   define right cell
            nc = (Ny-1)*(i+1)       ! right cell (cell-2) for vertical edge
            mesh%Cell(nc)%iface(4)     = ne
            mesh%Cell(nc)%iFaceDir     = -1
            mesh%Cell(nc)%inode(4)     = Ny*(i-0) + (j+1)
            mesh%Cell(nc)%inode(1)     = Ny*(i-0) + (j+0)
            mesh%Cell(nc)%ineighbor(4) = (Ny-1)*(i+0)
            mesh%Cell(nc)%locIdx       = 5                ! "channel"
            mesh%Cell(nc)%itype        = 4            ! quad

    enddo
    ! for right column
    i = Nx-1
    do j = 1,Ny-2
        ne = ne + 1        ! horizontal edge
        mesh%Face(ne)%locIdx   = 0                 ! internal
        mesh%Face(ne)%inode(1) = Ny*(i-0) + (j+1)
        mesh%Face(ne)%inode(2) = Ny*(i-1) + (j+1)
        mesh%Face(ne)%icell(1)  = (Ny-1)*(i-1) + j
        mesh%Face(ne)%icell(2)  = (Ny-1)*(i-1) + j+1

        !   define bottom cell
        nc = (Ny-1)*(i-1) + j  ! bottom cell (cell-1) for horizontal edge
        mesh%Cell(nc)%iface(3)     = ne
        mesh%Cell(nc)%iFaceDir     = 1
        mesh%Cell(nc)%inode(3)     = Ny*(i-0) + (j+1)
        mesh%Cell(nc)%inode(4)     = Ny*(i-1) + (j+1)
        mesh%Cell(nc)%ineighbor(3) = (Ny-1)*(i-1) + j+1
        mesh%Cell(nc)%locIdx       = 5                 ! "channel"
        mesh%Cell(nc)%itype        = 4            ! quad
        !   define top cell
        nc = (Ny-1)*(i-1) + j+1 ! right cell (cell-2) for horizontal edge
        mesh%Cell(nc)%iface(1)     = ne
        mesh%Cell(nc)%iFaceDir     = -1
        mesh%Cell(nc)%inode(1)     = Ny*(i-1) + (j+1)
        mesh%Cell(nc)%inode(2)     = Ny*(i-0) + (j+1)
        mesh%Cell(nc)%ineighbor(1) = (Ny-1)*(i-1) + j
        mesh%Cell(nc)%locIdx       = 5            ! "channel"
        mesh%Cell(nc)%itype        = 4            ! quad
    enddo

!    call output_2DMesh_gmsh(mesh,"2dmesh.msh")

    vtkTit = "grid for kspu"
    vtkFmt = "ascii"
!    call WriteGridVTK(vtkTit,vtkFmt)

    !write(*,*) "set locs:"
    ! set locations data:
    mesh%NLocs    = 5
    mesh%LocNBnds = 4
    mesh%LocBndIdx(1) = 1
    mesh%LocBndTit(1) = "int_electrode"
    mesh%LocBndIdx(2) = 2
    mesh%LocBndTit(2) = "outlet"
    mesh%LocBndIdx(3) = 3
    mesh%LocBndTit(3) = "ext_electrode"
    mesh%LocBndIdx(4) = 4
    mesh%LocBndTit(4) = "inlet"
    mesh%LocNVols = 1
    mesh%LocVolIdx(1) = 5
    mesh%LocVolTit(1) = "channel"

    call mesh_2DMeshCalcul(mesh)

    ! set boundary faces:
    allocate(mesh%BFaces(mesh%LocNBnds))        ! 2 locations: inlet, outlet
    ! 3. set faces indexes
    mesh%BFaces(:)%nfaces = 0
    do i = 1,mesh%NFaces
        locIdx = mesh%Face(i)%locIdx
        if(locIdx.ne.0) then
            locPos = valpos(mesh%locBndIdx,locIdx)
            mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
        end if
    end do
    do i = 1,mesh%LocNBnds
        call allocBFacesList(mesh%BFaces(i),mesh%BFaces(i)%nfaces)
    end do

    mesh%BFaces(:)%nfaces = 0
    do i = 1,mesh%NFaces
        locIdx = mesh%Face(i)%locIdx
        if(locIdx.ne.0) then
            locPos = valpos(mesh%locBndIdx,locIdx)
            mesh%BFaces(locPos)%nfaces = mesh%BFaces(locPos)%nfaces + 1
            j = mesh%BFaces(locPos)%nfaces
            mesh%BFaces(locPos)%iface(j) = i
            mesh%BFaces(locPos)%locIdx = locIdx
        end if
    end do


    return
end subroutine build_mesh_simple2D

