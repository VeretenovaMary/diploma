! ----------- gmsh triangle -------------
!
!
!         v
!         ^
!         |
!         4-------------------3
!         |                   |
!         |                   |
!         |                   |
!         |                   |
!         |                   |
!         |                   |
!         1-------------------2 --> u

!
! element nodes list -> element edges list

subroutine meshGetQuadrangleEdges(ElemNodes,EdgeNodes)
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
    EdgeNodes(3,2) = ElemNodes(4)
    EdgeNodes(4,1) = ElemNodes(4)
    EdgeNodes(4,2) = ElemNodes(1)

    return
end subroutine meshGetQuadrangleEdges

! ----------------------------------------------------------------
subroutine meshGetQuadrangleVectArea(elemNodes,mesh,vectArea)
    use module_mesh
    use module_Geom
    implicit none
    ! input
    integer, dimension(NumFaceNodes) :: elemNodes
    type(TMeshData)                  :: mesh
    ! output
    type(T3DPoint) :: vectArea
    ! local
    type(T3DPoint) :: vec
    integer, dimension(NumFaceNodes) :: triNodes

    vectArea%crd = 0.e0_r8p
    ! 1st triangle
    triNodes(:) = 0
    triNodes(1) = elemNodes(1)
    triNodes(2) = elemNodes(2)
    triNodes(3) = elemNodes(3)
    call meshGetTriangleVectArea(triNodes,mesh,vec)
    vectArea%crd = vectArea%crd + vec%crd
    ! 2nd triangle
    triNodes(:) = 0
    triNodes(1) = elemNodes(1)
    triNodes(2) = elemNodes(3)
    triNodes(3) = elemNodes(4)
    call meshGetTriangleVectArea(triNodes,mesh,vec)
    vectArea%crd = vectArea%crd + vec%crd

    return
end subroutine meshGetQuadrangleVectArea

! --------------------------------------------------------------
subroutine meshGetQuadrangleScalArea(elemNodes,mesh,elemArea)
    use module_mesh
    use module_Geom
    implicit none
    ! input
    integer, dimension(NumFaceNodes) :: elemNodes
    type(TMeshData)                  :: mesh
    ! output
    real(r8p)                        :: elemArea
    ! local
    type(T3DPoint)               :: vectArea
    integer, dimension(NumFaceNodes) :: triNodes

    elemArea = 0.e0_r8p
    ! 1st triangle
    triNodes(:) = 0
    triNodes(1) = elemNodes(1)
    triNodes(2) = elemNodes(2)
    triNodes(3) = elemNodes(3)
    call meshGetTriangleVectArea(triNodes,mesh,vectArea)
    elemArea = elemArea + vectLength(vectArea)
    ! 2nd triangle
    triNodes(:) = 0
    triNodes(1) = elemNodes(1)
    triNodes(2) = elemNodes(3)
    triNodes(3) = elemNodes(4)
    call meshGetTriangleVectArea(triNodes,mesh,vectArea)
    elemArea = elemArea + vectLength(vectArea)

    return
end subroutine meshGetQuadrangleScalArea

