! ----------- gmsh tetrahedron -------------
!
!                          v
!                        /
!                      /
!                    3
!                  / |\
!                /   | ‘\
!              /     |   ‘\
!            /       |     ‘\
!          /         |       ‘\
!         1----------’.--------2 --> u
!          ‘\.        |       /
!             ‘\.     |     /
!               ‘\.   |   /
!                 ‘\. | /
!                     4
!                       \
!                         \
!                           w
!
! element nodes list -> element edges list

subroutine meshGetTetrahedronEdges(ElemNodes,EdgeNodes)
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
    EdgeNodes(4,1) = ElemNodes(1)
    EdgeNodes(4,2) = ElemNodes(4)
    EdgeNodes(5,1) = ElemNodes(2)
    EdgeNodes(5,2) = ElemNodes(4)
    EdgeNodes(6,1) = ElemNodes(3)
    EdgeNodes(6,2) = ElemNodes(4)

    return
end subroutine meshGetTetrahedronEdges

! ! element nodes list -> element faces list
subroutine meshGetTetrahedronFaces(ElemNodes,CellFacesNodes,iFacesTypes)
    use module_mesh
    implicit none
    integer, dimension(NumCellNodes)              :: ElemNodes
    integer, dimension(NumCellFaces,NumFaceNodes) :: CellFacesNodes
    integer, dimension(NumCellFaces)              :: iFacesTypes

    CellFacesNodes(1,1) = ElemNodes(1)
    CellFacesNodes(1,2) = ElemNodes(3)
    CellFacesNodes(1,3) = ElemNodes(2)
    iFacesTypes(1) = 3
    CellFacesNodes(2,1) = ElemNodes(1)
    CellFacesNodes(2,2) = ElemNodes(2)
    CellFacesNodes(2,3) = ElemNodes(4)
    iFacesTypes(2) = 3
    CellFacesNodes(3,1) = ElemNodes(1)
    CellFacesNodes(3,2) = ElemNodes(4)
    CellFacesNodes(3,3) = ElemNodes(3)
    iFacesTypes(3) = 3
    CellFacesNodes(4,1) = ElemNodes(2)
    CellFacesNodes(4,2) = ElemNodes(3)
    CellFacesNodes(4,3) = ElemNodes(4)
    iFacesTypes(4) = 3

    return
end subroutine meshGetTetrahedronFaces

subroutine meshGetTetrahedronVolume(elemNodes,mesh,elemVolume)
    use module_mesh
    use module_Geom
    implicit none
    ! input
    integer, dimension(NumCellNodes) :: elemNodes
    type(TMeshData)                  :: mesh

    ! output
    real(r8p) :: elemVolume
    ! local
    type(T3DPoint) :: vec1,vec2,vec3,vecArea
    integer :: n1,n2

    ! 1st edge
    n1 = elemNodes(1)
    n2 = elemNodes(2)
    vec1%crd = mesh%Node(n2)%crd - mesh%Node(n1)%crd
    ! 2nd edge
    n1 = elemNodes(1)
    n2 = elemNodes(3)
    vec2%crd = mesh%Node(n2)%crd - mesh%Node(n1)%crd
    ! 3d edge
    n1 = elemNodes(1)
    n2 = elemNodes(4)
    vec3%crd = mesh%Node(n2)%crd - mesh%Node(n1)%crd
    ! area vector
    vecArea = vectVectProd(vec1,vec2)
    ! volume scalar
    elemVolume = 1.0e0_r8p/6.0e0_r8p*vectScalProd(vec3,vecArea)

    return
end subroutine meshGetTetrahedronVolume

