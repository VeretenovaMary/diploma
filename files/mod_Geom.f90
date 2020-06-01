module module_Geom
! <!>  - mark places for check
    use module_Settings

    ! --- 3d spatial point type
    type T3DPoint
        real(R8P), dimension(3) :: crd         ! point coordinates
    end type T3DPoint

    interface vectLength
        module procedure vectLength_dim
        module procedure vectLength_spt
    end interface

    interface vectScalProd
        module procedure vectScalProd_dim
        module procedure vectScalProd_spt
    end interface

    interface vectVectProd
        module procedure vectVectProd_dim
        module procedure vectVectProd_spt
    end interface

    private :: vectLength_dim,vectLength_spt
    private :: vectScalProd_dim,vectScalProd_spt
    private :: vectVectProd_dim,vectVectProd_spt

contains

    ! vector length calculation
    ! 1. version for dimension type
    function vectLength_dim(vec) result(vectLength)
        use module_Settings
        implicit none
        real(R8P), dimension(3) :: vec
        real(R8P) :: eps,vectLength
        eps = vec(1)**2+vec(2)**2+vec(3)**2
        vectLength = sqrt(eps)
        return
    end function vectLength_dim
    ! 2. version for T3DPoint type
    function vectLength_spt(vec) result(vectLength)
        use module_Settings
        implicit none
        type(T3DPoint) :: vec
        real(R8P) :: eps,vectLength
        eps = vec%crd(1)**2+vec%crd(2)**2+vec%crd(3)**2
        vectLength = sqrt(eps)
        return
    end function vectLength_spt


    ! vec1.*.vec2 scalar multiplication
    ! 1. version for dimension type
    function vectScalProd_dim(vec1,vec2) result(vectScalProd)
        use module_Settings
        implicit none
        real(R8P), dimension(3) :: vec1,vec2
        real(R8P) :: vectScalProd
        vectScalProd = vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
        return
    end function vectScalProd_dim
    ! 2. version for T3DPoint type
    function vectScalProd_spt(vec1,vec2) result(vectScalProd)
        use module_Settings
        implicit none
        type(T3DPoint) :: vec1,vec2
        real(R8P) :: vectScalProd
        vectScalProd = vec1%crd(1)*vec2%crd(1)+vec1%crd(2)*vec2%crd(2)+vec1%crd(3)*vec2%crd(3)
        return
    end function vectScalProd_spt


    ! vec1.x.vec2: vector multiplication
    ! 1. version for dimension type
    function vectVectProd_dim(vec1,vec2) result(vectVectProd)
        use module_Settings
        implicit none
        real(R8P), dimension(3) :: vec1,vec2
        real(R8P), dimension(3) :: vectVectProd
        vectVectProd(1) = vec1(2)*vec2(3)-vec1(3)*vec2(2)
        vectVectProd(2) = vec1(3)*vec2(1)-vec1(1)*vec2(3)
        vectVectProd(3) = vec1(1)*vec2(2)-vec1(2)*vec2(1)
        return
    end function vectVectProd_dim
    ! 2. version for T3DPoint type
    function vectVectProd_spt(vec1,vec2) result(vectVectProd)
        use module_Settings
        implicit none
        type(T3DPoint) :: vec1,vec2
        type(T3DPoint) :: vectVectProd
        vectVectProd%crd(1) = vec1%crd(2)*vec2%crd(3)-vec1%crd(3)*vec2%crd(2)
        vectVectProd%crd(2) = vec1%crd(3)*vec2%crd(1)-vec1%crd(1)*vec2%crd(3)
        vectVectProd%crd(3) = vec1%crd(1)*vec2%crd(2)-vec1%crd(2)*vec2%crd(1)
        return
    end function vectVectProd_spt

    ! calculating distance between {line} = {ptL}+{vecL}*t and point {pt}
    function point_line_distance(ptL,vecL,pt)
        use module_Settings
        implicit none
        ! input values
        real(R8P), dimension(3) :: ptL,vecL     ! line parameters: point on line and direction vector
        real(R8P), dimension(3) :: pt           ! point coordinates
        ! output values
        real(R8P) :: point_line_distance
        ! local variables
        real(R8P)               :: t           ! line parameter
        real(R8P)               :: eps
        real(R8P), dimension(3) :: vecD

        eps  = vectScalProd(pt-ptL,vecL)
        t    = eps/vectLength(vecL)**2
        vecD = (pt-ptL) - vecL*t
        point_line_distance = vectLength(vecD)

        return
    end function point_line_distance

    ! --------------------------------------------
    function point_point_distance(pt1,pt2)
        ! <!> - check
        implicit none
        real(R8P)    :: point_point_distance
        real(R8P),dimension(3) :: pt1,pt2
        real(R8P)    :: eps
        integer(I4P) :: i

        eps = 0.e0
        do i = 1,3
            eps = eps + (pt1(i)-pt2(i))**2
        enddo
        point_point_distance = dsqrt(eps)

        return
    end function point_point_distance


end module module_Geom
