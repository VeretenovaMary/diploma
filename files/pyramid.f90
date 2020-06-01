! -------------- gmsh pyramid -------------
!
!
!                    5
!                  /..\
!                /  / \‘\
!              /   .   ‘.‘\
!            /     / w   \ ‘\
!          /      .  ^    ‘. ‘\
!         1-------/--|------\--4
!          \     .   |       ‘. \
!           \    /   +---->v   \ \
!            \  .     \         ‘.\
!             \ /      \          \\
!              2--------\-----------3
!                        \
!                         u
!
! element nodes list -> element edges list

subroutine meshGetPyramidEdges(ElemNodes,EdgeNodes)
    use module_mesh
    implicit none
    integer, dimension(NumCellNodes)   :: ElemNodes
    integer, dimension(NumCellEdges,2) :: EdgeNodes

    EdgeNodes(1,1) = ElemNodes(1)
    EdgeNodes(1,2) = ElemNodes(2)
    EdgeNodes(2,1) = ElemNodes(2)
    EdgeNodes(2,2) = ElemNodes(3)
    EdgeNodes(3,1) = ElemNodes(3)
    EdgeNodes(3,2) = ElemNodes(4)
    EdgeNodes(4,1) = ElemNodes(4)
    EdgeNodes(4,2) = ElemNodes(1)
    EdgeNodes(5,1) = ElemNodes(1)
    EdgeNodes(5,2) = ElemNodes(5)
    EdgeNodes(6,1) = ElemNodes(2)
    EdgeNodes(6,2) = ElemNodes(5)
    EdgeNodes(7,1) = ElemNodes(3)
    EdgeNodes(7,2) = ElemNodes(5)
    EdgeNodes(8,1) = ElemNodes(4)
    EdgeNodes(8,2) = ElemNodes(5)

    return
end subroutine meshGetPyramidEdges

! ! element nodes list -> element faces list
subroutine meshGetPyramidFaces(ElemNodes,CellFacesNodes,iFacesTypes)
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
    CellFacesNodes(2,1) = ElemNodes(1)
    CellFacesNodes(2,2) = ElemNodes(2)
    CellFacesNodes(2,3) = ElemNodes(5)
    iFacesTypes(2) = 3
    CellFacesNodes(3,1) = ElemNodes(2)
    CellFacesNodes(3,2) = ElemNodes(3)
    CellFacesNodes(3,3) = ElemNodes(5)
    iFacesTypes(3) = 3
    CellFacesNodes(4,1) = ElemNodes(3)
    CellFacesNodes(4,2) = ElemNodes(4)
    CellFacesNodes(4,3) = ElemNodes(5)
    iFacesTypes(4) = 3
    CellFacesNodes(5,1) = ElemNodes(4)
    CellFacesNodes(5,2) = ElemNodes(1)
    CellFacesNodes(5,3) = ElemNodes(5)
    iFacesTypes(5) = 3

    return
end subroutine meshGetPyramidFaces

subroutine meshGetPyramidVolume(elemNodes,mesh,elemVolume)
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
    integer, dimension(NumCellNodes) :: tetNodes
    real(r8p) :: tetVol

    elemVolume = 0.e0_r8p

    ! create 1st tetrahedron from pyramid
    tetNodes(:) = 0
    tetNodes(1) = elemNodes(1)
    tetNodes(2) = elemNodes(4)
    tetNodes(3) = elemNodes(5)
    tetNodes(4) = elemNodes(3)
    call meshGetTetrahedronVolume(tetNodes,mesh,tetVol)
    elemVolume = elemVolume + tetVol

    ! create 2nd tetrahedron from pyramid
    tetNodes(:) = 0
    tetNodes(1) = elemNodes(1)
    tetNodes(2) = elemNodes(3)
    tetNodes(3) = elemNodes(5)
    tetNodes(4) = elemNodes(2)
    call meshGetTetrahedronVolume(tetNodes,mesh,tetVol)
    elemVolume = elemVolume + tetVol
    return
end subroutine meshGetPyramidVolume

