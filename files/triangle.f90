! ----------- gmsh triangle -------------
!
!                          v
!                        /
!                      /
!                    3
!                  / \
!                /    ‘\
!              /        ‘\
!            /            ‘\
!          /                ‘\
!         1-------------------2 --> u

!
! element nodes list -> element edges list

! ------------------------------------------------------
subroutine meshGetTriangleEdges(ElemNodes,EdgeNodes)
    use module_mesh
    implicit none
    integer, dimension(NumFaceNodes)   :: ElemNodes
    integer, dimension(NumFaceEdges,2) :: EdgeNodes

    EdgeNodes(:,:) = 0
    EdgeNodes(1,1) = ElemNodes(1)
    EdgeNodes(1,2) = ElemNodes(2)
    EdgeNodes(2,1) = ElemNodes(2)
    EdgeNodes(2,2) = ElemNodes(3)
    EdgeNodes(3,1) = ElemNodes(3)
    EdgeNodes(3,2) = ElemNodes(1)

    return
end subroutine meshGetTriangleEdges

! -------------------------------------------------------------
subroutine meshGetTriangleVectArea(elemNodes,mesh,vectArea)
    use module_mesh
    use module_Geom
    implicit none
    ! input
    integer, dimension(NumFaceNodes) :: elemNodes
    type(TMeshData)                  :: mesh

    ! output
    type(T3DPoint) :: vectArea
    ! local
    type(T3DPoint) :: vec1,vec2
    integer :: n1,n2

    ! 1st edge
    n1 = elemNodes(1)
    n2 = elemNodes(2)
    vec1%crd = mesh%Node(n2)%crd - mesh%Node(n1)%crd
    ! 2nd edge
    n1 = elemNodes(1)
    n2 = elemNodes(3)
    vec2%crd = mesh%Node(n2)%crd - mesh%Node(n1)%crd

    ! triangle area
    vectArea     = vectVectProd(vec1,vec2)
    vectArea%crd = 0.5e0_r8p*vectArea%crd

    return
end subroutine meshGetTriangleVectArea

! ----------------------------------------------------------
subroutine meshGetTriangleScalArea(elemNodes,mesh,elemArea)
    use module_mesh
    use module_Geom
    implicit none
    ! input
    integer, dimension(NumFaceNodes) :: elemNodes
    type(TMeshData)                  :: mesh
    ! output
    real(r8p) :: elemArea
    ! local
    type(T3DPoint) :: vectArea

    call meshGetTriangleVectArea(elemNodes,mesh,vectArea)
    elemArea = vectLength(vectArea)


    return
end subroutine meshGetTriangleScalArea

