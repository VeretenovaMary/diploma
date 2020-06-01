module module_GridFunct
! <!>  - mark places for check
    use module_Settings
    use module_Mesh


    ! --- Grid Nodal Scalar ---
    type TGridNodalScalar
        type(TMeshData), pointer  :: pmesh
        character(32)             :: ftit         ! title of the variable
        character(8)              :: units        ! variable units
        integer                   :: Nval         ! number of values
        real(R8P), pointer        :: fval(:)      ! index = index of mesh node
    end type TGridNodalScalar

    ! --- Grid Cell Scalar ---
    type TGridCellScalar
        type(TMeshData), pointer  :: pmesh
        character(32)             :: ftit         ! title of the variable
        character(8)              :: units        ! variable units
        integer                   :: Nval         ! number of values
        real(R8P), pointer        :: fval(:)      ! index = index of mesh cell
    end type TGridCellScalar

    ! --- Grid Points Scalar ---
    type TGridPointsScalar
        type(TMeshData), pointer    :: pmesh
        character(32)               :: ftit         ! title of the variable
        character(8)                :: units        ! variable units
        integer                     :: Nval         ! number of values
        real(R8P), pointer          :: fval(:)      ! index = index of mesh cell
        type(T3DPoint), pointer     :: pt(:)        ! points
        integer, pointer            :: BindId(:)    ! points bind index
        integer(2), pointer         :: BindTp(:)    ! points bind type (class): 0-node,1-edge,2-face,3-cell,4-free
    end type TGridPointsScalar


    ! --- vector in free space (with free dimension) ---
    type TFreeFormVector
        real(R8P), pointer        :: cmp(:)      ! array of components values
    end type TFreeFormVector


    ! --- Grid Nodal Vector ---
    type TGridNodalVector
        type(TMeshData), pointer       :: pmesh => null()
        character(32)                  :: ftit              ! title of the variable
        character(8)                   :: units             ! variable units
        integer                        :: Nval              ! number of values
        integer                        :: Ncmp              ! vector space dimension
        type(TFreeFormVector), pointer :: fvec(:) => null() ! index = index of mesh node
    end type TGridNodalVector


    ! --- Grid Cell Vector ---
    type TGridCellVector
        type(TMeshData), pointer       :: pmesh => null()
        character(32)                  :: ftit              ! title of the variable
        character(8)                   :: units             ! variable units
        integer                        :: Nval              ! number of values
        integer                        :: Ncmp              ! vector space dimension
        type(TFreeFormVector), pointer :: fvec(:) => null() ! index = index of mesh cell
    end type TGridCellVector


    ! --- Grid Points Vector ---
    type TGridPointsVector
        type(TMeshData), pointer       :: pmesh => null()
        character(32)                  :: ftit              ! title of the variable
        character(8)                   :: units             ! variable units
        integer                        :: Nval              ! number of values
        integer                        :: Ncmp              ! vector space dimension
        type(TFreeFormVector), pointer :: fvec(:) => null() ! index = index of mesh cell
        type(T3DPoint), pointer        :: pt(:)             ! points coordinates
        integer, pointer               :: BindId(:)         ! points bind index
        integer(2), pointer            :: BindTp(:)   ! points bind type (class): 0-node,1-edge,2-face,3-cell,4-free
    end type TGridPointsVector


    integer :: gridfNum = 0

    interface new_grid_function
        module procedure new_nodal_scalar
        module procedure new_cell_scalar
        module procedure new_points_scalar
        module procedure new_nodal_vector
        module procedure new_cell_vector
        module procedure new_points_vector
    end interface

    interface close_grid_function
        module procedure close_nodal_scalar
        module procedure close_cell_scalar
        module procedure close_points_scalar
        module procedure close_nodal_vector
        module procedure close_cell_vector
        module procedure close_points_vector
    end interface

    private new_nodal_scalar,new_nodal_vector
    private new_cell_scalar,new_cell_vector
    private new_points_scalar,new_points_vector
    private close_nodal_scalar,close_nodal_vector
    private close_cell_scalar,close_cell_vector
    private close_points_scalar,close_points_vector

    public new_grid_function
    public close_grid_function

contains

    ! create new grid nodal scalar function
    subroutine new_nodal_scalar(gf,ftit,pmesh)
        implicit none
        type(TGridNodalScalar)   :: gf
        character(*)             :: ftit
        type(TMeshData), pointer :: pmesh

        gf%ftit  = trim(ftit)
        gf%pmesh => pmesh
        gf%Nval = pmesh%NNodes
        allocate(gf%fval(pmesh%NNodes))
        gf%fval(:) = 0.d0

        return
    end subroutine new_nodal_scalar

    ! create new grid cell scalar function
    subroutine new_cell_scalar(gf,ftit,pmesh)
        implicit none
        type(TGridCellScalar)    :: gf
        character(*)             :: ftit
        type(TMeshData),pointer  :: pmesh

        gf%ftit  = trim(ftit)
        gf%pmesh => pmesh
        gf%Nval = pmesh%NCells
        allocate(gf%fval(pmesh%NCells))
        gf%fval(:) = 0.d0

        return
    end subroutine new_cell_scalar


    ! create new points scalar function
    subroutine new_points_scalar(gf,ftit,Npts,pmesh)
        implicit none
        type(TGridPointsScalar)  :: gf
        character(*)             :: ftit
        type(TMeshData),pointer  :: pmesh
        integer :: Npts                   ! number of points

        gf%ftit  = trim(ftit)
        gf%pmesh => pmesh
        gf%Nval  = Npts
!        write(*,*) trim(gf%ftit)," -> ",Npts
        allocate(gf%fval(Npts))
        allocate(gf%pt(Npts))
        allocate(gf%BindId(Npts))
        allocate(gf%BindTp(Npts))
        gf%fval(:) = 0.d0
        gf%BindId(:) = 0
        gf%BindTp(:) = 4        ! 4 - indefined (free)

        return
    end subroutine new_points_scalar


    ! create new grid nodal vector function
    subroutine new_nodal_vector(gf,ftit,Ncmp,pmesh)
        implicit none
        type(TGridNodalVector)   :: gf
        character(*)             :: ftit
        type(TMeshData), pointer :: pmesh
        integer :: Ncmp
        integer :: i

        gf%ftit  = trim(ftit)
        gf%pmesh => pmesh
        gf%Nval = pmesh%NNodes
        gf%Ncmp = Ncmp
        allocate(gf%fvec(pmesh%NNodes))
        do i = 1,gf%Nval
            call alloc_free_form_vector(gf%fvec(i),Ncmp)
            gf%fvec(i)%cmp(:) = 0.d0
        end do

        return
    end subroutine new_nodal_vector

    ! create new grid cell vector function
    subroutine new_cell_vector(gf,ftit,Ncmp,pmesh)
        implicit none
        type(TGridCellVector)    :: gf
        character(*)             :: ftit
        type(TMeshData),pointer  :: pmesh
        integer :: Ncmp
        integer :: i

        gf%ftit  = trim(ftit)
        gf%pmesh => pmesh
        gf%Nval = pmesh%NCells
        gf%Ncmp = Ncmp
        allocate(gf%fvec(pmesh%NCells))
        do i = 1,gf%Nval
            call alloc_free_form_vector(gf%fvec(i),Ncmp)
            gf%fvec(i)%cmp(:) = 0.d0
        end do

        return
    end subroutine new_cell_vector


    ! create new grid cell vector function
    subroutine new_points_vector(gf,ftit,Npts,Ncmp,pmesh)
        implicit none
        type(TGridPointsVector)  :: gf
        character(*)             :: ftit
        type(TMeshData),pointer  :: pmesh
        integer :: Npts                   ! number of points
        integer :: Ncmp                   ! number of components
        integer :: i

        gf%ftit  = trim(ftit)
        gf%pmesh => pmesh
        gf%Nval = Npts
        gf%Ncmp = Ncmp
!        write(*,*) trim(gf%ftit)," -> ",Npts
        allocate(gf%fvec(Npts))
        do i = 1,gf%Nval
            call alloc_free_form_vector(gf%fvec(i),Ncmp)
            gf%fvec(i)%cmp(:) = 0.d0
        end do

        allocate(gf%pt(Npts))
        allocate(gf%BindId(Npts))
        allocate(gf%BindTp(Npts))
        gf%BindId(:) = 0
        gf%BindTp(:) = 4              ! 4 - indefined (free)

        return
    end subroutine new_points_vector


    ! allocate free form vector
    subroutine alloc_free_form_vector(fvec,Ncmp)
        implicit none
        type(TFreeFormVector) :: fvec
        integer :: Ncmp

        allocate(fvec%cmp(Ncmp))

        return
    end subroutine alloc_free_form_vector

    ! close nodal scalar function
    subroutine close_nodal_scalar(x)
        implicit none
        type(TGridNodalScalar) :: x

        deallocate(x%fval)

        return
    end subroutine close_nodal_scalar

    ! close cell scalar function
    subroutine close_cell_scalar(x)
        implicit none
        type(TGridCellScalar) :: x

        deallocate(x%fval)
        return
    end subroutine close_cell_scalar

    ! close points scalar function
    subroutine close_points_scalar(x)
        implicit none
        type(TGridPointsScalar) :: x

        deallocate(x%fval)
        deallocate(x%pt)
        deallocate(x%BindId)
        deallocate(x%BindTp)
        return
    end subroutine close_points_scalar

    ! close nodal vector function
    subroutine close_nodal_vector(x)
        implicit none
        type(TGridNodalVector) :: x

        integer :: i
        do i = 1,x%Nval
                deallocate(x%fvec(i)%cmp)
        end do
        deallocate(x%fvec)

        return
    end subroutine close_nodal_vector

    ! close cell vector function
    subroutine close_cell_vector(x)
        implicit none
        type(TGridCellVector) :: x
        integer :: i
        do i = 1,x%Nval
                deallocate(x%fvec(i)%cmp)
        end do

        deallocate(x%fvec)

        return
    end subroutine close_cell_vector

    ! close points vector function
    subroutine close_points_vector(x)
        implicit none
        type(TGridPointsVector) :: x
        integer :: i
        do i = 1,x%Nval
                deallocate(x%fvec(i)%cmp)
        end do

        deallocate(x%fvec)
        deallocate(x%pt)
        deallocate(x%BindId)
        deallocate(x%BindTp)

        return
    end subroutine close_points_vector



end module module_GridFunct
