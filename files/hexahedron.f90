! --------- gmsh hexahedron -------------
!
!
!                v
!         4----------3
!         |\     ^   |\
!         | \    |   | \
!         |  \   |   |  \
!         |   8------+---7
!         |   |  +---|-- |-> u
!         1---+---\--2   |
!          \  |    \  \  |
!           \ |     \  \ |
!            \|      w  \|
!             5----------6

! element nodes list -> element edges list

subroutine meshGetHexahedronEdges(ElemNodes,EdgeNodes)
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
    EdgeNodes(3,2) = ElemNodes(4)
    EdgeNodes(4,1) = ElemNodes(4)
    EdgeNodes(4,2) = ElemNodes(1)
    EdgeNodes(5,1) = ElemNodes(5)
    EdgeNodes(5,2) = ElemNodes(6)
    EdgeNodes(6,1) = ElemNodes(6)
    EdgeNodes(6,2) = ElemNodes(7)
    EdgeNodes(7,1) = ElemNodes(7)
    EdgeNodes(7,2) = ElemNodes(8)
    EdgeNodes(8,1) = ElemNodes(8)
    EdgeNodes(8,2) = ElemNodes(5)
    EdgeNodes(9,1)  = ElemNodes(1)
    EdgeNodes(9,2)  = ElemNodes(5)
    EdgeNodes(10,1) = ElemNodes(2)
    EdgeNodes(10,2) = ElemNodes(6)
    EdgeNodes(11,1) = ElemNodes(3)
    EdgeNodes(11,2) = ElemNodes(7)
    EdgeNodes(12,1) = ElemNodes(4)
    EdgeNodes(12,2) = ElemNodes(8)

    return
end subroutine meshGetHexahedronEdges

! ! element nodes list -> element faces list
subroutine meshGetHexahedronFaces(ElemNodes,CellFacesNodes,iFacesTypes)
    use module_mesh
    implicit none
    integer, dimension(NumCellNodes)              :: ElemNodes
    integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
    integer, dimension(NumCellFaces)              :: iFacesTypes

    CellFacesNodes(1,1) = ElemNodes(1)
    CellFacesNodes(1,2) = ElemNodes(4)
    CellFacesNodes(1,3) = ElemNodes(3)
    CellFacesNodes(1,4) = ElemNodes(2)
    iFacesTypes(1) = 4
    CellFacesNodes(2,1) = ElemNodes(5)
    CellFacesNodes(2,2) = ElemNodes(6)
    CellFacesNodes(2,3) = ElemNodes(7)
    CellFacesNodes(2,4) = ElemNodes(8)
    iFacesTypes(2) = 4
    CellFacesNodes(3,1) = ElemNodes(1)
    CellFacesNodes(3,2) = ElemNodes(5)
    CellFacesNodes(3,3) = ElemNodes(8)
    CellFacesNodes(3,4) = ElemNodes(4)
    iFacesTypes(3) = 4
    CellFacesNodes(4,1) = ElemNodes(1)
    CellFacesNodes(4,2) = ElemNodes(2)
    CellFacesNodes(4,3) = ElemNodes(6)
    CellFacesNodes(4,4) = ElemNodes(5)
    iFacesTypes(4) = 4
    CellFacesNodes(5,1) = ElemNodes(2)
    CellFacesNodes(5,2) = ElemNodes(3)
    CellFacesNodes(5,3) = ElemNodes(7)
    CellFacesNodes(5,4) = ElemNodes(6)
    iFacesTypes(5) = 4
    CellFacesNodes(6,1) = ElemNodes(3)
    CellFacesNodes(6,2) = ElemNodes(4)
    CellFacesNodes(6,3) = ElemNodes(8)
    CellFacesNodes(6,4) = ElemNodes(7)
    iFacesTypes(6) = 4

    return
end subroutine meshGetHexahedronFaces

! ----------------------------------------------
subroutine meshGetHexahedronVolume(elemNodes,mesh,elemVolume)
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
    integer, dimension(NumCellNodes) :: pyrNodes
    real(r8p) :: pyrVol

    elemVolume = 0.e0_r8p

    ! create 1st pyramid from hexahedron
    pyrNodes(:) = 0
    pyrNodes(1) = elemNodes(1)
    pyrNodes(2) = elemNodes(5)
    pyrNodes(3) = elemNodes(6)
    pyrNodes(4) = elemNodes(2)
    pyrNodes(5) = elemNodes(4)
    call meshGetPyramidVolume(pyrNodes,mesh,pyrVol)
    elemVolume = elemVolume + pyrVol

    ! create 2nd pyramid from hexahedron
    pyrNodes(:) = 0
    pyrNodes(1) = elemNodes(5)
    pyrNodes(2) = elemNodes(8)
    pyrNodes(3) = elemNodes(7)
    pyrNodes(4) = elemNodes(6)
    pyrNodes(5) = elemNodes(4)
    call meshGetPyramidVolume(pyrNodes,mesh,pyrVol)
    elemVolume = elemVolume + pyrVol

    ! create 3d pyramid from hexahedron
    pyrNodes(:) = 0
    pyrNodes(1) = elemNodes(6)
    pyrNodes(2) = elemNodes(7)
    pyrNodes(3) = elemNodes(3)
    pyrNodes(4) = elemNodes(2)
    pyrNodes(5) = elemNodes(4)
    call meshGetPyramidVolume(pyrNodes,mesh,pyrVol)
    elemVolume = elemVolume + pyrVol

    return
end subroutine meshGetHexahedronVolume

