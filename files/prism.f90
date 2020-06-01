! --------- gmsh prism -------------
!
!                   w
!                   |
!                   |
!                   4
!                  /|\
!                 / | \
!                /  |  \
!               /   |   \
!              /    |    \
!             5-----+-----6
!             |     |     |
!             |    /|\    |
!             |   / | \   |
!             |  /  |  \  |
!             | /   1   \ |
!             |/   / \   \|
!             |   /   \   |
!            /|  /     \  |\
!           u | /       \ | v
!             |/         \|
!             2-----------3

! element nodes list -> element edges list

subroutine meshGetPrismEdges(ElemNodes,EdgeNodes)
    use module_mesh
    implicit none
    integer, dimension(NumCellNodes)   :: ElemNodes
    integer, dimension(NumCellEdges,2) :: EdgeNodes

    EdgeNodes(:,:) = 0
    EdgeNodes(1,1) = ElemNodes(1)
    EdgeNodes(1,2) = ElemNodes(2)

    EdgeNodes(2,1) = ElemNodes(2)
    EdgeNodes(2,2) = ElemNodes(3)

    EdgeNodes(3,1) = ElemNodes(3)
    EdgeNodes(3,2) = ElemNodes(1)

    EdgeNodes(4,1) = ElemNodes(4)
    EdgeNodes(4,2) = ElemNodes(5)

    EdgeNodes(5,1) = ElemNodes(5)
    EdgeNodes(5,2) = ElemNodes(6)

    EdgeNodes(6,1) = ElemNodes(6)
    EdgeNodes(6,2) = ElemNodes(4)

    EdgeNodes(7,1) = ElemNodes(1)
    EdgeNodes(7,2) = ElemNodes(4)

    EdgeNodes(8,1) = ElemNodes(2)
    EdgeNodes(8,2) = ElemNodes(5)

    EdgeNodes(9,1) = ElemNodes(3)
    EdgeNodes(9,2) = ElemNodes(6)

    return
end subroutine meshGetPrismEdges

! ! element nodes list -> element faces list
subroutine meshGetPrismFaces(ElemNodes,CellFacesNodes,iFacesTypes)
    use module_mesh
    implicit none
    integer, dimension(NumCellNodes)              :: ElemNodes
    integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
    integer, dimension(NumCellFaces)              :: iFacesTypes

    CellFacesNodes(1,1) = ElemNodes(1)
    CellFacesNodes(1,2) = ElemNodes(3)
    CellFacesNodes(1,3) = ElemNodes(2)
    iFacesTypes(1) = 3

    CellFacesNodes(2,1) = ElemNodes(4)
    CellFacesNodes(2,2) = ElemNodes(5)
    CellFacesNodes(2,3) = ElemNodes(6)
    iFacesTypes(2) = 3

    CellFacesNodes(3,1) = ElemNodes(1)
    CellFacesNodes(3,2) = ElemNodes(2)
    CellFacesNodes(3,3) = ElemNodes(5)
    CellFacesNodes(3,4) = ElemNodes(4)
    iFacesTypes(3) = 4
    CellFacesNodes(4,1) = ElemNodes(1)
    CellFacesNodes(4,2) = ElemNodes(4)
    CellFacesNodes(4,3) = ElemNodes(6)
    CellFacesNodes(4,4) = ElemNodes(3)
    iFacesTypes(4) = 4
    CellFacesNodes(5,1) = ElemNodes(2)
    CellFacesNodes(5,2) = ElemNodes(3)
    CellFacesNodes(5,3) = ElemNodes(5)
    CellFacesNodes(5,4) = ElemNodes(6)
    iFacesTypes(5) = 4

    return
end subroutine meshGetPrismFaces

! ----------------------------------------------
subroutine meshGetPrismVolume(elemNodes,mesh,elemVolume)
    ! by spliting pyramide in two tetrahedrons:
    use module_mesh
    use module_Geom
    implicit none
    ! input
    integer, dimension(NumCellNodes) :: elemNodes
    type(TMeshData)                  :: mesh

    ! output
    real(r8p) :: elemVolume
    ! local
    integer, dimension(NumCellNodes) :: pyrNodes,tetNodes
    real(r8p) :: pyrVol,tetVol

    elemVolume = 0.e0_r8p

    ! create tetrahedron from prism
    tetNodes(:) = 0
    tetNodes(1) = elemNodes(1)
    tetNodes(2) = elemNodes(2)
    tetNodes(3) = elemNodes(3)
    tetNodes(4) = elemNodes(4)
    call meshGetTetrahedronVolume(tetNodes,mesh,tetVol)
    elemVolume = elemVolume + tetVol

    ! create pyramid from prism
    pyrNodes(:) = 0
    pyrNodes(1) = elemNodes(2)
    pyrNodes(2) = elemNodes(5)
    pyrNodes(3) = elemNodes(6)
    pyrNodes(4) = elemNodes(3)
    pyrNodes(5) = elemNodes(4)
    call meshGetPyramidVolume(pyrNodes,mesh,pyrVol)
    elemVolume = elemVolume + pyrVol

    return
end subroutine meshGetPrismVolume

