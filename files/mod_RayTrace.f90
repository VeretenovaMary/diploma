module module_RayTrace
! <!>  - mark places for check
    use module_Settings
    use module_Mesh
    use module_Geom
    use module_interpolate
    use module_TxtRead
    use module_OutputMesh

    integer, parameter       :: raysRepUnit = 500
    character(32), parameter :: raysRepName = 'setup_rays.err'
    integer :: raysNumErrors
    integer :: raysNumWarnings

    ! --- ANGULAR GRID --- : type with angular grid information
    type TAngularGridData
        integer(i4p)             :: NDirs               ! sum number of directions
        integer(i4p)             :: NTheta              ! number of theta directions
        integer(i4p)             :: NPhi                ! number of theta directions
        real(r8p)                :: dOmega              ! spatial angle (case of uniform angle grid)
        real(r8p), dimension(3)  :: oxdir               ! Ox axis direction (local) in global coordinates
        real(r8p), dimension(3)  :: oydir               ! Oy axis direction (local) in global coordinates
        real(r8p), dimension(3)  :: ozdir               ! Oz axis direction (local) in global coordinates
        real(r8p), pointer       :: theta(:)  => NULL() ! size = number of directions
        real(r8p), pointer       :: phi(:)    => NULL() ! size = number of directions
    end type TAngularGridData


    ! --- RAY --- : type with single ray information
    type TRayData
        integer                     :: NSegm               ! number of segments
        real(r4p), pointer          :: segLen(:) => NULL() ! size = number of segments
        integer(i4p), pointer       :: icell(:)  => NULL() ! size = number of segments
        integer(i4p), pointer       :: iface(:)  => NULL() ! size = number of segments
        type(T3DPoint), pointer     :: ptOut(:)  => NULL() ! size = number of segments
    end type TRayData

    ! --- RAY --- : type with single ray information
    type TRayDataCompact
        integer                     :: NSegm               ! number of segments
        real(r4p), pointer          :: segLen(:) => NULL() ! size = number of segments
        integer(i4p), pointer       :: icell(:)  => NULL() ! size = number of segments
        integer(i4p), pointer       :: iface(:)  => NULL() ! size = number of segments
    end type TRayDataCompact

    ! --- spatial point tracing data --- : type with single point tracing data
    type TPointTraceData
        type(T3DPoint)      :: pt0         ! initial (tracing start point)
        integer(i4p)            :: bindId      ! start point bind (node/cell/face) index
        type(TRayData), pointer :: rays(:) => NULL() ! size = number of directions
    end type TPointTraceData

    ! --- spatial point tracing data --- : type with single point tracing data
    type TPointTraceDataCompact
        type(T3DPoint)                 :: pt0         ! initial (tracing start point)
        integer(i4p)                   :: bindId      ! start point bind (node/cell/face) index
        type(TRayDataCompact), pointer :: rays(:) => NULL() ! size = number of directions
    end type TPointTraceDataCompact

    ! --- ray tree data ---
    type TTreeSegm
        type(T3DPoint)      :: pt0               !
        type(T3DPoint)      :: pt1               !
        integer(i4p)            :: icell  = 0        ! cell index ( =0 for outer ray)
        integer(i4p)            :: iface  = 0        ! face index ( =0 for outer ray)
        integer(i4p)            :: inode  = 0        ! node index ( =0 for outer ray)
        real(r8p)               :: segLen = 0.e0_r8p ! ray segment length
        real(r8p), dimension(4) :: weight = 0.e0_r8p ! size = NumFaceNodes
    end type TTreeSegm

    type TRayNodalTree
        type(TTreeSegm), pointer :: segm(:)  => null()  ! segment data
        integer(i4p)             :: nsegm               ! number of segments (nodes)
    end type TRayNodalTree


    ! --- mesh tracing data --- : type with all mesh points tracing data
    type TTraceData
        character(8)                    :: rayTp                ! ray type: "long"/"short"
        type(TAngularGridData), pointer :: pAngGrid => null()   ! pointer to angular grid
        type(TMeshData), pointer        :: pMesh => null()      ! pointer to mesh where tracing is
        integer(i4p)                    :: regIdx               ! tracing region Idx
        character(32)                   :: elemTp               ! start point bind type (node/cell/face)
        integer(i4p)                    :: NPoints              ! number of points in tracing data
        type(TPointTraceData),  pointer :: point(:) => NULL()   ! size = number of directions
        type(TRayNodalTree), pointer    :: rayTree(:) => NULL() ! size = number of directions
    end type TTraceData

    ! --- mesh tracing data --- : type with all mesh points tracing data
    type TTraceDataCompact
        character(8)                    :: rayTp                ! ray type: "long"/"short"
        type(TAngularGridData), pointer :: pAngGrid => null()   ! pointer to angular grid
        type(TMeshData), pointer        :: pMesh => null()      ! pointer to mesh where tracing is
        integer(i4p)                    :: regIdx               ! tracing region Idx
        character(32)                   :: elemTp               ! start point bind type (node/cell/face)
        integer(i4p)                    :: NPoints              ! number of points in tracing data
        type(TPointTraceDataCompact), pointer :: point(:) => NULL()   ! size = number of directions
        type(TRayNodalTree), pointer    :: rayTree(:) => NULL() ! size = number of directions
    end type TTraceDataCompact


    ! --- mesh 1D tracing data --- : type with all mesh points tracing data
    type TTrace1DPRM       ! -1.d0/(a*kappa+b) + c   ! don't forget kappa->kappa*dimLen due to scale
        real(r8p)    :: a
        real(r8p)    :: b
        real(r8p)    :: c
        real(r8p)    :: k
        integer(2)   :: idir
        !real(r8p)    :: dOm
        !real(r8p)    :: dOmz
    end type

    type TTrace1DPoint
        type(TTrace1DPRM), pointer :: prm0(:)   => null()  ! 0th mju momentum: Sum(exp(-tau)*dOm) -> cU
        type(TTrace1DPRM), pointer :: prm1(:)   => null()  ! 1st mju momentum: Sum(exp(-tau)*mju*dOm) -> Wz
        integer(i4p), pointer      :: icell0(:) => null()  ! end of the ray
        integer(i4p), pointer      :: icell1(:) => null()  ! beginning of the ray
    end type

    type TTrace1DData
        type(TMeshData),     pointer :: pMesh    => null()   ! pointer to mesh where tracing is
        type(TTrace1DPoint), pointer :: point(:) => null()   ! size = number of directions
        real(r8p)                    :: LenUnit              ! cm, length unit
        integer(i4p)                 :: Npoints
        integer(i4p)                 :: Nfaces               ! number of mesh faces = Nprms
    endtype
    ! -----------------------------------------------------------------------

    interface new_rays_tracer
        module procedure new_rays_tracer_extended
        module procedure new_rays_tracer_compact
    end interface

    interface close_rays_tracer
        module procedure close_rays_tracer_compact
        module procedure close_rays_tracer_extended
    end interface

    interface allocTracePointRays
        module procedure allocTracePointRaysExtended
        module procedure allocTracePointRaysCompact
    end interface


    interface allocTracePointRayData
        module procedure allocTracePointRayDataExtended
        module procedure allocTracePointRayDataCompact
    end interface




    real(r8p), parameter :: pi = 3.1415926500000000e+00_r8p
    private :: pi

contains


!   --------------------------------------------
    subroutine alloc_Trace1DData(Trace1DData,Npoints,Nfaces)
        implicit none
        ! input variables
        type(TTrace1DData) :: Trace1DData
        integer            :: Npoints
        integer            :: Nfaces
        integer            :: ipt

        Trace1DData%Npoints = int(Npoints,i4p)
        Trace1DData%Nfaces  = int(Nfaces,i4p)
        allocate(Trace1DData%point(Npoints))
        do ipt = 1,Npoints
            call alloc_Trace1DPoint(Trace1DData%point(ipt),Nfaces)
        enddo
        return
    end subroutine alloc_Trace1DData

!   . . . . . . . . . . . . . . . . . . . . . . . .
    subroutine alloc_Trace1DPoint(Trace1DPoint,Nfaces)
        implicit none
        ! input variables
        type(TTrace1DPoint) :: Trace1DPoint
        integer             :: Nfaces

        allocate(Trace1DPoint%prm0(Nfaces))
        allocate(Trace1DPoint%prm1(Nfaces))
        allocate(Trace1DPoint%icell0(Nfaces))
        allocate(Trace1DPoint%icell1(Nfaces))

        return
    end subroutine alloc_Trace1DPoint

!   . . . . . . . . . . . . . . . . . . . . . . . .
    subroutine close_Trace1DData(Trace1DData)
        implicit none
        ! input variables
        type(TTrace1DData) :: Trace1DData
        integer            :: ipt

        do ipt = 1,Trace1DData%Npoints
            call dealloc_Trace1DPoint(Trace1DData%point(ipt))
        enddo
        deallocate(Trace1DData%point)

        return
    end subroutine close_Trace1DData

!   . . . . . . . . . . . . . . . . . . . . . . . .
    subroutine dealloc_Trace1DPoint(Trace1DPoint)
        implicit none
        ! input variables
        type(TTrace1DPoint) :: Trace1DPoint

        deallocate(Trace1DPoint%prm0)
        deallocate(Trace1DPoint%prm1)
        deallocate(Trace1DPoint%icell0)
        deallocate(Trace1DPoint%icell1)

        return
    end subroutine dealloc_Trace1DPoint


!   --------------------------------------------
    subroutine allocRayTree(rayTree,nsegm)
        implicit none
        ! input variables
        type(TRayNodalTree) :: rayTree
        integer :: nsegm

        rayTree%nsegm = nsegm
        allocate(rayTree%segm(nsegm))

        return
    end subroutine allocRayTree

!   --------------------------------------------
    subroutine set_rays_nodal_tree(trace)
        implicit none
        ! input
        type(TTraceData) :: trace           ! tracing data
        ! local
        integer :: idir,ndirs

        allocate(trace%rayTree(trace%pAngGrid%NDirs))
        do idir = 1,trace%pAngGrid%NDirs
            call allocRayTree(trace%rayTree(idir),trace%pMesh%NNodes)
        end do

        write(*,*)
        write(*,*) ' - - - - - - - - '
        write(*,*) " rays traversing ... "
        ndirs = trace%pAngGrid%NDirs
        do idir = 1,ndirs
            call calc_ray_nodal_tree(trace,idir)
        end do
        write(*,"(a,2(a,I3))")

        return
    end subroutine set_rays_nodal_tree

! ---------------------------------------------------------------------
    subroutine calc_ray_nodal_tree(trace,idir)
        implicit none
        ! input
        type(TTraceData) :: trace         ! tree pointer
        ! local
        type(TMeshData),        pointer :: pMesh         ! mesh to trace
        integer        :: i,j,k,iface,inode,icell,jnode,jcell
        integer        :: ipt,idir,isegm,idx,num,iTop,ios
        integer        :: istate
        type TTraverseElems
            integer(i4p), pointer    :: stateCell(:) => NULL()  ! sum number of directions
            integer(i4p), pointer    :: stateNode(:) => NULL()  ! sum number of directions
            integer(i4p), pointer    :: actCell(:) => NULL()  ! sum number of directions
            integer(i4p), pointer    :: treeNode(:) => NULL()  ! sum number of directions
            integer(i4p)             :: iTop = 0
            integer(i4p)             :: nFin = 0
            integer(i4p)             :: iTree = 0
        end type TTraverseElems
        type(TTraverseElems) :: stack
        ! bilinear interpolation:
        integer               :: nnodes,ndirs
        type(T3DPoint)    :: nodes(4)     ! coords of face nodes: size = 4 (quad) or 3 (tri)
        type(T3DPoint)    :: pt
        real(r8p)             :: Fwt(4)       ! function weights to Fpt
        type(T3DPoint)    :: dFwt(4)      ! function weights to dFpt
        integer               :: istop


        pMesh => trace%pMesh


!        write(*,*)
!        write(*,*) "  creating nodes tree ... "
        allocate(stack%stateCell(trace%pMesh%NCells))
        allocate(stack%actCell(trace%pMesh%NCells))
        allocate(stack%stateNode(trace%pMesh%NNodes))
        allocate(stack%treeNode(trace%pMesh%NNodes))

        stack%stateCell(:) = 0
        stack%stateNode(:) = 0
        stack%actCell(:)   = 0
        stack%iTree = 0
        ! traverse through external bounds

        do iface = 1,pMesh%NBFaces
            icell = pmesh%Face(iface)%icell(1)
            do i = 1,pmesh%Face(iface)%nNodes
                inode = pmesh%Face(iface)%inode(i)
                if(trace%point(inode)%rays(idir)%NSegm.eq.0) then
                    if(stack%stateNode(inode).eq.0) then
                        stack%iTree = stack%iTree + 1
                        stack%treeNode(stack%iTree) = inode
                    endif
                    stack%stateNode(inode) = 1
                endif
            end do
        enddo

        do inode = 1,pMesh%NNodes
            if(stack%stateNode(inode).eq.1) then
                do i = 1,pmesh%Node(inode)%ncells
                    icell = pmesh%Node(inode)%icell(i)
                    stack%stateCell(icell) = stack%stateCell(icell) + 1
                end do
            end if
        end do

        stack%iTop = 0
        stack%nFin = 0
        do icell = 1,pMesh%NCells
            if(stack%stateCell(icell).eq.pmesh%Cell(icell)%nNodes) then
                stack%stateCell(icell) = 1
                stack%nFin = stack%nFin + 1
            elseif(stack%stateCell(icell).gt.0) then
                stack%iTop = stack%iTop + 1
                stack%actCell(stack%iTop) = icell
!                write(*,*) icell,' - > ',stack%stateCell(icell)
            endif
        enddo

        ! traverse other nodes
        num = stack%iTop
!        write(*,"(a,10(i5))") "stack ini: ",stack%actCell(1:num)
        isegm = 0
        ndirs = trace%pAngGrid%NDirs
        do while(stack%nFin.lt.pMesh%NCells)
            num = 100*stack%nFin/pMesh%NCells
            write(*,"(a,3(a,I3),a)",advance="no") achar(13),'done ',idir,"/",ndirs,"  stack: ",num,"%"
            isegm = isegm + 1
!            write(*,*) 'iteration: ',isegm
            iTop = stack%iTop

            i = 0
            istate = 0
            do while( i .lt. stack%iTop )
                i = i + 1
                icell = stack%actCell(i)
                do j = 1,pmesh%Cell(icell)%nNodes
                    inode = pmesh%Cell(icell)%inode(j)
                    iface = trace%point(inode)%rays(idir)%iface(1)
                    num = trace%point(inode)%rays(idir)%NSegm
                    idx = 0
                    if(num.gt.0) then
                        idx = 1         ! initial state: face is full
                        do k = 1,pmesh%Face(iface)%nNodes
                            jnode = pmesh%Face(iface)%inode(k)
                            if(stack%stateNode(jnode).eq.0) then
                                idx = 0   ! -> face not full
                            end if
                        end do
                    endif

                    if(idx.eq.1 .and. stack%stateNode(inode).eq.0) then
                        istate = istate + 1         ! check for infinit loop
                        stack%stateNode(inode) = 1
                        stack%iTree = stack%iTree + 1
                        stack%treeNode(stack%iTree) = inode

                        do k = 1,pmesh%Node(inode)%ncells
                            jcell = pmesh%Node(inode)%icell(k)
                            stack%stateCell(jcell) = stack%stateCell(jcell) + 1
                            num = stack%iTop
                            if(.not.any(stack%actCell(1:num).eq.jcell)) then
!                                write(*,"(a,i5,a,10(i5))") "add:",jcell,"//",stack%actCell(1:num)
                                stack%iTop = stack%iTop + 1
                                num = stack%iTop
                                stack%actCell(num) = jcell
!                                write(*,"(a,i5,a,10(i5))") "  ->",jcell,"//",stack%actCell(1:num)
!                                write(*,"(a,i5,a,10(i5))")
                            end if

                            if(stack%stateCell(jcell).eq.pmesh%Cell(jcell)%nNodes) then
                                num = stack%iTop
!                                write(*,"(a,i5,a,10(i5))") "del:",jcell,"//",stack%actCell(1:num)
                                ipt = valpos(stack%actCell(1:num),jcell)
                                do ios = ipt,num-1
                                    stack%actCell(ios) = stack%actCell(ios+1)
                                enddo
                                stack%actCell(num) = 0
                                num = num-1
!                                write(*,"(a,i5,a,10(i5))") "  ->",jcell,"//",stack%actCell(1:num)
!                                write(*,"(a,i5,a,10(i5))")
                                stack%nFin = stack%nFin + 1
                                stack%iTop = stack%iTop - 1
                            end if
                        end do
                    end if
                end do
            end do
            if(istate.eq.0) then
                write(*,*) " infinite loop"
                stop
            end if
!            stack%iTop = iTop
!            pause
        end do

        do i = 1,trace%pMesh%NNodes
            ipt = stack%treeNode(i)
            trace%rayTree(idir)%segm(i)%icell = trace%point(ipt)%rays(idir)%icell(1)
            trace%rayTree(idir)%segm(i)%iface = trace%point(ipt)%rays(idir)%iface(1)
            trace%rayTree(idir)%segm(i)%inode = trace%point(ipt)%bindId
            trace%rayTree(idir)%segm(i)%segLen = trace%point(ipt)%rays(idir)%segLen(1)
            trace%rayTree(idir)%segm(i)%pt0 = trace%point(ipt)%pt0
            trace%rayTree(idir)%segm(i)%pt1 = trace%point(ipt)%rays(idir)%ptOut(1)
        end do

        deallocate(stack%stateCell)
        deallocate(stack%actCell)
        deallocate(stack%stateNode)
        deallocate(stack%treeNode)

        ! nodes interpolation coefficients:
        do i = 1,trace%rayTree(idir)%nsegm
            iface = trace%rayTree(idir)%segm(i)%iface
            if(iface.gt.0) then
                nnodes = trace%pmesh%Face(iface)%nNodes
                do j = 1,4                                  ! <!> 4 -> num nodes
                    nodes(j)%crd = 0.0_r8p
                enddo
                do j = 1,nnodes
                    inode = pmesh%Face(iface)%inode(j)
                    nodes(j)%crd = pmesh%Node(inode)%crd
                end do
                pt = trace%rayTree(idir)%segm(i)%pt1
                istop = 0
                call interp_bilinear(nnodes,nodes,pt,Fwt,dFwt,istop)
                if(istop.eq.1) then
                    call output_point(37,trace%rayTree(idir)%segm(i)%pt0)
                    call output_point(37,trace%rayTree(idir)%segm(i)%pt1)
                    write(37,*) ""
                    write(37,*) ""
                    call output_face(38,iface,trace%pMesh)
                    call output_cell(39,trace%rayTree(idir)%segm(i)%icell,trace%pMesh)
                    flush(37)
                    flush(38)
                    flush(39)
                    stop
                end if
                trace%rayTree(idir)%segm(i)%weight(:) = 0.0_r8p
                trace%rayTree(idir)%segm(i)%weight(1:nnodes) = Fwt(1:nnodes)
                do j = 1,nnodes
                    if(abs(trace%rayTree(idir)%segm(i)%weight(j)).lt.1.d-80) then
                        trace%rayTree(idir)%segm(i)%weight(j) = 0.d0
                    end if
                enddo
            else
                trace%rayTree(idir)%segm(i)%weight(:) = 0.0_r8p
            end if
        end do

        return
    end subroutine calc_ray_nodal_tree

! ---------------------------------------------------------------------
    subroutine new_angular_grid(cfgName,AngGrid)  !
        implicit none
        ! input data
        character(*) :: cfgName
        ! output data
        type(TAngularGridData) :: AngGrid
        ! local variables
        integer(i4p)  :: NTheta
        real(4), allocatable :: dimCosTh(:)     ! theta: [0,pi/2]
        real(4), allocatable :: dimCosPhi(:)    ! phi: [0,pi]
        real(4)      :: dOm
        real(r8p)    :: theta,phi
        real(r8p)    :: cosTh,cosPhi
        real(r8p)    :: sinTh,sinPhi
        real(r8p)    :: cosTh1,cosPhi1
        real(r8p)    :: sinTh1,sinPhi1
        real(r8p)    :: cosTh0,cosPhi0
        real(r8p)    :: sinTh0,sinPhi0
        integer(i4p) :: NQuarts,NDirs,Npt
        integer      :: i,j,k,idx,iQuart
        integer      :: ios
        ! angular grid options
        character(32)             :: AngMethod
        integer, dimension(4)     :: AngQuarts
        real(r8p), dimension(3,3) :: AngMtx,TMtx        ! matrix for basis translation, transposed matrix
        real(r8p), dimension(3)   :: vec0,vec1
        integer :: cfgUnit,GetUnit

        cfgUnit = GetUnit()
        open(cfgUnit,file=cfgName)
        call txt_read_struct_value(cfgUnit,"AngGrid","N_theta",NTheta,ios)
        call txt_read_struct_value(cfgUnit,"AngGrid","method",AngMethod,ios)
        AngQuarts(:) = 0
        call txt_read_struct_array(cfgUnit,"AngGrid","quarters",AngQuarts,NQuarts,ios)

        call txt_read_struct_array(cfgUnit,"AngGrid","oxdir",vec0,i,ios)
        AngGrid%oxdir(:) = vec0(:)
        AngMtx(1,:)      = vec0(:)
        call txt_read_struct_array(cfgUnit,"AngGrid","oydir",vec0,i,ios)
        AngGrid%oydir(:) = vec0(:)
        AngMtx(2,:)      = vec0(:)
        call txt_read_struct_array(cfgUnit,"AngGrid","ozdir",vec0,i,ios)
        AngGrid%ozdir(:) = vec0(:)
        AngMtx(3,:)      = vec0(:)

        TMtx = transpose(AngMtx)

        allocate(dimCosTh(NTheta),dimCosPhi(NTheta*(NTheta+1)))


        if(trim(AngMethod).eq."Carlson") then
            AngGrid%NTheta = 1_i4p*NQuarts*NTheta           !
            AngGrid%NPhi   = 2_i4p*NQuarts*NTheta           ! 2x - because NPhi=2*NTheta in quarter sphere (inactive value)
            NDirs          = NQuarts*NTheta*(NTheta+1_i4p)  ! 2x - for half spatial angle [0,2*pi]
            AngGrid%NDirs  = NDirs
            allocate(AngGrid%theta(NDirs))
            allocate(AngGrid%phi(NDirs))

            ! set angular grid by Carlson quadrature
            call Carlsonquadr(int(NTheta,8),dimCosTh,dimCosPhi,dOm)

            Npt = NTheta*(NTheta+1)
            idx = 0
            write(*,"(a,i4)") " Ndirs = ",NDirs
            do i = 1,NTheta
                do j = 1,2*(NTheta-i+1)
                    idx = idx + 1
                    ! . . . . . turning basis (change coordinates): . . . . .
                    ! auxiliary basis coordinates
                    cosTh  = real(dimCosTh(i),r8p)
                    sinTh  = sqrt(1.e0_r8p-cosTh**2)
                    cosPhi = real(dimCosPhi(idx),r8p)
                    sinPhi = sqrt(1.e0_r8p-cosPhi**2)
                    ! extension for other quarters (quadrants)
                    do k = 1,NQuarts
                        iQuart = AngQuarts(k)
                        theta = acos(cosTh)
                        phi   = acos(cosPhi)
                        select case(iQuart)
                        case(1)                           ! 1st quarter
                            theta = acos(cosTh)
                            phi   = acos(cosPhi)
                        case(2)                           ! 2nd quarter
                            theta = acos(cosTh)
                            !phi   = 2*pi - acos(cosPhi)
                            phi   = pi + acos(cosPhi)
                        case(3)                           ! 3rd quarter
                            theta = pi + acos(cosTh)
                            phi   = acos(cosPhi)
                        case(4)                           ! 4th quarter
                            theta = pi + acos(cosTh)
                            !phi   = 2*pi - acos(cosPhi)
                            phi   = pi + acos(cosPhi)
                        end select

                        cosTh1  = cos(theta)
                        sinTh1  = sin(theta)
                        cosPhi1 = cos(phi)
                        sinPhi1 = sin(phi)
                        !write(*,*) iQuart,"sin(Phi) = ",sinPhi1

                        vec1(1) = sinTh1*cosPhi1
                        vec1(2) = sinTh1*sinPhi1
                        vec1(3) = cosTh1
                        ! native basis coordinates
                        vec0 = matmul(TMtx,vec1)
                        cosTh0 = vec0(3)                     ! >0 for 1/2 quarters, <0 for 3/4 quarters
                        sinTh0 = sqrt(1.e0_r8p-cosTh0**2)    ! always >= 0
                        if(sinTh0.ne.0.e0_r8p) then
                            cosPhi0 = vec0(1)/sinTh0
                            sinPhi0 = vec0(2)/sinTh0
                        else
                            cosPhi0 = 1.e0_r8p
                            sinPhi0 = 0.e0_r8p
                        end if
                        ! ..................................................
                        ! theta and phi in native basis:
                        theta = acos(cosTh0)
                        phi   = acos(cosPhi0)
                        if(sinPhi0.ge.0) then      ! 1st quarter
                            phi   = acos(cosPhi0)
                        elseif(sinPhi0.lt.0) then  ! 2nd quarter
                            phi   = 2*pi - acos(cosPhi0)
                        endif
                        ! put directions for selected sphere quadrants:
                        AngGrid%theta((k-1)*Npt + idx) = theta
                        AngGrid%phi  ((k-1)*Npt + idx) = phi
                    enddo
                enddo
            enddo
            AngGrid%dOmega = dOm
        else
            write(*,*) "error: angular grid unknown method: ",trim(AngMethod)     ! <!> make file for errors output
            stop
        end if
!        pause

        deallocate(dimCosTh,dimCosPhi)


        return
    end subroutine new_angular_grid

    ! - - - - - - - create single ray - - - - - -
    subroutine new_single_ray(theta,phi,AngGrid)  !
        implicit none
        ! input data
        real(r8p)                 :: theta,phi
        ! output data
        type(TAngularGridData) :: AngGrid
        ! local variables
        real(r8p)                 :: dOm
        integer(i4p)              :: NDirs
        real(r8p), dimension(3,3) :: AngMtx        ! matrix for basis translation, transposed matrix
        real(r8p), dimension(3,3) :: TMtx          ! matrix for basis translation, transposed matrix
        real(r8p)    :: cosTh1,cosPhi1
        real(r8p)    :: sinTh1,sinPhi1
        real(r8p)    :: cosTh0,cosPhi0
        real(r8p)    :: sinTh0,sinPhi0
        real(r8p)    :: theta0,phi0
        real(r8p), dimension(3)   :: vec0,vec1


        AngGrid%NTheta = 1
        AngGrid%NPhi   = 1
        NDirs          = 1
        dOm            = 0.e0
        AngGrid%NDirs  = NDirs
        allocate(AngGrid%theta(NDirs))
        allocate(AngGrid%phi(NDirs))

        ! matrix for basis rotation
        AngMtx(1,:)   = AngGrid%oxdir(:)
        AngMtx(2,:)   = AngGrid%oydir(:)
        AngMtx(3,:)   = AngGrid%ozdir(:)
        TMtx          = transpose(AngMtx)


        cosTh1  = cos(theta)
        sinTh1  = sin(theta)
        cosPhi1 = cos(phi)
        sinPhi1 = sin(phi)
        !write(*,*) iQuart,"sin(Phi) = ",sinPhi1

        vec1(1) = sinTh1*cosPhi1
        vec1(2) = sinTh1*sinPhi1
        vec1(3) = cosTh1
        ! native basis coordinates
        vec0 = matmul(TMtx,vec1)
        cosTh0 = vec0(3)                     ! >0 for 1/2 quarters, <0 for 3/4 quarters
        sinTh0 = sqrt(1.e0_r8p-cosTh0**2)    ! always >= 0
        if(sinTh0.ne.0.e0_r8p) then
            cosPhi0 = vec0(1)/sinTh0
            sinPhi0 = vec0(2)/sinTh0
        else
            cosPhi0 = 1.e0_r8p
            sinPhi0 = 0.e0_r8p
        end if
        ! ..................................................
        ! theta and phi in native basis:
        theta0 = acos(cosTh0)
        phi0   = acos(cosPhi0)
        if(sinPhi0.ge.0) then      ! 1st quarter
            phi0   = acos(cosPhi0)
        elseif(sinPhi0.lt.0) then  ! 2nd quarter
            phi0   = 2*pi - acos(cosPhi0)
        endif

        ! put directions for selected sphere quadrants:
        AngGrid%theta(1) = theta0
        AngGrid%phi(1)   = phi0
        AngGrid%dOmega   = 2*pi

        return
    end subroutine new_single_ray

    ! - - - - angular_grid close - - - - - - -
    subroutine close_angular_grid(AngGrid)
        implicit none
        type(TAngularGridData) :: AngGrid

        deallocate(AngGrid%phi,AngGrid%theta)

        return
    end subroutine close_angular_grid

    ! - - - - tracing close - - - - - - -
    subroutine close_rays_tracer_extended(trace)
        implicit none
        type(TTraceData) :: trace
        integer :: idir,ipt
        do ipt = 1,trace%NPoints
           do idir = 1,trace%pAngGrid%NDirs
               deallocate(trace%point(ipt)%rays(idir)%icell)
               deallocate(trace%point(ipt)%rays(idir)%iface)
               deallocate(trace%point(ipt)%rays(idir)%ptOut)
               deallocate(trace%point(ipt)%rays(idir)%segLen)
           enddo
           deallocate(trace%point(ipt)%rays)
        enddo
        deallocate(trace%point)

        if(trim(trace%rayTp).eq."short") then
            do idir = 1,trace%pAngGrid%NDirs
                deallocate(trace%rayTree(idir)%segm)
            enddo
            deallocate(trace%rayTree)
        endif

        return
    end subroutine close_rays_tracer_extended

    ! - - - - tracing close - - - - - - -
    subroutine close_rays_tracer_compact(trace)
        implicit none
        type(TTraceDataCompact) :: trace
        integer :: idir,ipt
        do ipt = 1,trace%NPoints
           do idir = 1,trace%pAngGrid%NDirs
               deallocate(trace%point(ipt)%rays(idir)%icell)
               deallocate(trace%point(ipt)%rays(idir)%iface)
               deallocate(trace%point(ipt)%rays(idir)%segLen)
           enddo
           deallocate(trace%point(ipt)%rays)
        enddo
        deallocate(trace%point)

        if(trim(trace%rayTp).eq."short") then
            do idir = 1,trace%pAngGrid%NDirs
                deallocate(trace%rayTree(idir)%segm)
            enddo
            deallocate(trace%rayTree)
        endif

        return
    end subroutine close_rays_tracer_compact

!   --------------------------------------------
    subroutine allocTracePointRaysExtended(tracePoint,NDirs)
        implicit none
        ! input variables
        type(TPointTraceData) :: tracePoint
        integer :: NDirs

        allocate(tracePoint%rays(NDirs))

        return
    end subroutine allocTracePointRaysExtended

!   --------------------------------------------
    subroutine allocTracePointRaysCompact(tracePoint,NDirs)
        implicit none
        ! input variables
        type(TPointTraceDataCompact) :: tracePoint
        integer :: NDirs

        allocate(tracePoint%rays(NDirs))

        return
    end subroutine allocTracePointRaysCompact

!   --------------------------------------------
    subroutine allocTracePointRayDataExtended(tracePointRayData,NSegm)
        implicit none
        ! input variables
        type(TRayData) :: tracePointRayData
        integer :: NSegm

        allocate(tracePointRayData%icell(NSegm))
        allocate(tracePointRayData%iface(NSegm))
        allocate(tracePointRayData%ptOut(NSegm))
        allocate(tracePointRayData%segLen(NSegm))

        return
    end subroutine allocTracePointRayDataExtended

!   --------------------------------------------
    subroutine allocTracePointRayDataCompact(tracePointRayData,NSegm)
        implicit none
        ! input variables
        type(TRayDataCompact) :: tracePointRayData
        integer :: NSegm

        allocate(tracePointRayData%icell(NSegm))
        allocate(tracePointRayData%iface(NSegm))
        allocate(tracePointRayData%segLen(NSegm))

        return
    end subroutine allocTracePointRayDataCompact

!   --------------------------------------------
    subroutine new_rays_tracer_extended(rayType,pAngGrid,pMesh,locTit,elemType,trace)  !
        use module_OutputMesh
        implicit none
        ! input variables
        type(TAngularGridData),pointer :: pAngGrid        ! angular grid
        type(TMeshData),       pointer :: pMesh           ! mesh to trace
        character(*)                   :: locTit          ! mesh region name where tracing nodes exists
        character(*)                   :: elemType        ! mesh elements from where ray comes (Nodes/Faces/Cells)
        character(*)                   :: rayType         ! mesh elements from where ray comes (Nodes/Faces/Cells)
        ! setup parameters
        integer, parameter   :: maxNumSegm = 10000
        real(r8p), parameter :: val0 = 1.e-8_r8p

        ! output variables
        type(TTraceData)        :: trace           ! tracing data

        ! local variables
        integer                 :: NPoints,NDirs
        type(TRayData)          :: rayData
        real(r8p)               :: segLen

        real(r8p)               :: cosTh,cosPhi,sinTh,sinPhi
        real(r8p)               :: theta,phi
        real(r8p), dimension(3) :: pt0,pt1,ptc,ptt
        real(r8p), dimension(3) :: rayDir,fnorm
        integer(I4P)            :: icell0,icell1
        integer(I4P)            :: iFace1
        integer                 :: icase
        integer :: iOz,iOx,iOy
        integer :: i,j
        integer :: ipt,iel,idir,isegm,Nsegm
        integer       :: locPos,locIdx
        character(32) :: locType
        character(64) :: sfmt
        !logical :: point_in_cell
        real(8) :: point_in_cell
        real(r8p) :: SegmSizeBytes,RaySizeBytes,TracerSizeBytes,TracerSizeMb

        !write(*,*) 'pNDirs = ',pAngGrid%NDirs

        trace%pmesh    => pMesh
        trace%pAngGrid => pAngGrid

        locPos  = 0
        locPos  = valpos(trace%pmesh%Locations%locTit,trim(locTit))   ! region position in list
        locIdx  = trace%pmesh%Locations(locPos)%locIdx
        trace%regIdx = locIdx
        if(locPos.eq.0) then
            write(*,*) 'error: tracing region not found '
            stop
        end if
        locType = trace%pmesh%Locations(locPos)%locType
!        write(*,*) 'locPos   = ',locPos
!        write(*,*) 'locType  = ',locType
!        write(*,*) 'locTitle = ',locTit
!        write(*,*) 'elemType = ',elemType
        if(index(elemType,'node').gt.0) then
            NPoints = trace%pmesh%Locations(locPos)%nnodes
        elseif(index(elemType,'face').gt.0) then
            NPoints = trace%pmesh%Locations(locPos)%NElems
        elseif(index(elemType,'cell').gt.0) then
            NPoints = trace%pmesh%Locations(locPos)%NElems
        else
            write(*,*) 'error: please set in config the elements to trace (node/face/cell) '
            stop
        endif

        trace%NPoints = NPoints
        trace%elemTp  = trim(elemType)
        trace%rayTp   = trim(rayType)
        allocate(trace%point(NPoints))
        NDirs = pAngGrid%NDirs
        do i = 1,trace%NPoints
            call allocTracePointRays(trace%point(i),trace%pAngGrid%NDirs)
        enddo

        iOx = 1        ! Oy
        iOy = 2        ! Oz
        iOz = 3        ! axis of symmetry (if 1 - Ox is axis)
        allocate(rayData%icell(maxNumSegm))
        allocate(rayData%iface(maxNumSegm))
        allocate(rayData%ptOut(maxNumSegm))
        allocate(rayData%segLen(maxNumSegm))


        !sfmt = '("Point(",i5,") = {",3(1pE12.4e2,","),1pE12.4e2,"};")'
        write(*,*)
        write(*,*) ' - - - - - - - - '
        write(*,*) ' tracing ... '

        Nsegm = 0
        TracerSizeBytes = 0.0e0
        TracerSizeMb    = 0.0e0
        SegmSizeBytes = 0.0e0
        SegmSizeBytes = SegmSizeBytes +   real(sizeof(rayData%icell(1)),8)
        SegmSizeBytes = SegmSizeBytes +   real(sizeof(rayData%iface(1)),8)
        SegmSizeBytes = SegmSizeBytes + 3*real(sizeof(rayData%ptOut(1)%crd(1)),8)
        SegmSizeBytes = SegmSizeBytes +   real(sizeof(rayData%segLen(1)),8)
        !write(*,*) "SegmSizeBytes = ",SegmSizeBytes
        !pause

        sfmt = '('//'a,'//'a,I3,a,'//'4x,a,i7,a,i7,'//'4x,a,1pE10.3e2,a'//')'

        do ipt = 1,NPoints
        !do ipt = 6529,6529
            !write(*,"(a,a,I3,a)",advance="no") achar(13),'done ',ipt*100/NPoints," %"
            !write(*,"(a,a,I3,a,i7,a,i7)",advance="no") achar(13),'done ',ipt*100/NPoints," %",ipt," Nsegm = ",Nsegm
            write(*,sfmt,advance="no") achar(13), &
                    'done ',ipt*100/NPoints," %", &
                    'point ',ipt,'/',NPoints,     &
                    'TracerSize = ',TracerSizeMb,' Mb'
            iel = 0
            ptc = 0.e0_r8p
            if(index(elemType,'node').gt.0) then           ! tracing nodes
                iel = pMesh%Locations(locPos)%inode(ipt)   ! node index
                ptc = pmesh%Node(iel)%crd
            elseif(index(elemType,'face').gt.0) then       ! tracing face centers
                iel = pMesh%Locations(locPos)%ielem(ipt)   ! face index
                ptc = pmesh%Face(iel)%center
            elseif(index(elemType,'cell').gt.0) then       ! tracing cell centers
                iel = pMesh%Locations(locPos)%ielem(ipt)   ! cell index
                ptc = pmesh%Cell(iel)%center
            endif

            trace%point(ipt)%pt0%crd = ptc
            trace%point(ipt)%bindId  = iel

            do idir = 1,trace%pAngGrid%NDirs
            !do idir = 42,42
            !do idir = 195,195
                !write(*,"(a,a,I3,a,I7,2x,a,I7)",advance="no") achar(13),'done ',ipt*100/NPoints," %",ipt,'idir = ',idir
                theta = trace%pAngGrid%theta(idir)
                phi   = trace%pAngGrid%phi(idir)
                cosTh  = cos(theta)
                cosPhi = cos(phi)
                sinTh  = sin(theta)
                sinPhi = sin(phi)

                ! ray direction vector
                rayDir(iOx) = sinTh*cosPhi
                rayDir(iOy) = sinTh*sinPhi
                rayDir(iOz) = cosTh

                pt0(iOx) = ptc(iOx) + 1.e2*val0*rayDir(iOx)
                pt0(iOy) = ptc(iOy) + 1.e2*val0*rayDir(iOy)
                pt0(iOz) = ptc(iOz) + 1.e2*val0*rayDir(iOz)

                icase = -1               ! 0/1/2 - ray is (not inside)/(inside)/ trace region
                if(index(elemType,'node').gt.0) then        ! tracing nodes
                    icell0 = -1
                    do i = 1,pmesh%Node(iel)%ncells
                        icell1 = pmesh%Node(iel)%icell(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !ptt(1) = pt0(1)                    ! only 2D-RZ case !!!!!!!!!!!!!!!!!!!! <!>
                        !ptt(2) = sqrt(pt0(2)**2+pt0(3)**2) ! only 2D-RZ case !!!!!!!!!!!!!!!!!!!! <!>
                        !ptt(3) = 0.0e0                     ! only 2D-RZ case !!!!!!!!!!!!!!!!!!!! <!>
                        ptt = pt0
                        call mapper_XYZ_to_RZ_dirX(ptt)
                        if(point_in_cell(pMesh,ptt,icell1,0.d0).gt.0.d0) then
                            !write(*,"(a,3(1x,1pE12.4e3))") 'ptt: ',ptt
                            icell0 = icell1
                            icase = 1
                            !write(*,*) " point in cell: ",icell0
                        end if                                            ! <!> reflection condition
                    end do
                elseif(index(elemType,'face').gt.0) then    ! tracing face centers
                    icell0 = -1
                    icell1 = pmesh%Face(iel)%icell(1)
                    if(point_in_cell(pMesh,pt0,icell1,0.d0).gt.0.d0) then
                        icell0 = icell1
                        icase = 1
                    end if
                elseif(index(elemType,'cell').gt.0) then    ! tracing cell centers
                    iCell0 = iel
                    icase = 1
                endif

!                if(icase.le.0) then
!                    write(*,"(a)") " tracing point is outside mesh: "
!                    write(*,"(3(1x,1pE10.3e2),a,3(1x,1pE10.3e2))") ptc," -> ",pt0
!                end if


                iFace1 = 0
                isegm  = 0
                rayData%NSegm = 0
                rayData%icell(:)     = 0
                rayData%iface(:)     = 0
                rayData%segLen(:)    = 0.e0

!                if(icase.eq.1) then
!                    open(31,file="cell.geo")
!                    call output_point_gmsh(31,pt0)
!                    call output_cell_gmsh(31,iCell0,pMesh)
!                    ptt = ptc
!                    call mapper_XYZ_to_RZ_dirX(ptt)
!                    pt1 = ptt+0.1*rayDir
!                    call mapper_XYZ_to_RZ_dirX(pt1)
!                    call output_line_gmsh(31,ptt,pt1,'initial ray direction')
!                endif

                do while(icase.ge.1)
                    isegm = isegm + 1
                    Nsegm = Nsegm + 1
                    if(pMesh%meshDim.eq.3) then
                        call tracer_3D(pMesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)
                    elseif(pMesh%meshDim.eq.2) then
                        call tracer_RZ(pMesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)
                    endif
                    if(iCell0.eq.iCell1) then
                        write(*,*)
                        write(*,*) 'ray is in the same cell ...'
                        stop
                    endif
                    !write(*,*) " isegm = ",isegm
                    !if(isegm.eq.8) pause
                    !if(isegm.eq.57) pause

                    if(trim(rayType).eq."short") then
                        icase = 0
                    endif

                    rayData%NSegm = isegm
                    rayData%icell(isegm) = iCell0
                    rayData%iface(isegm) = iFace1
                    !rayData%ptOut(isegm)%crd = pt1
                    ptt = pt1-pt0
                    segLen = vectLength(ptt)
                    rayData%segLen(isegm) = real(segLen,r4p)
                    iCell0 = iCell1
                    pt0    = pt1

                enddo

!                if(icase.gt.-1) then
!                    close(31)
!                    pause
!                    write(*,"(a)",advance="no") repeat(char(8),16)
!                    write(*,"(a,i7)",advance="no") " Nsegm = ",Nsegm
!                endif

                rayData%NSegm = isegm
!                call single_ray_output(ipt,idir,ptc,rayData,pMesh,pAngGrid)
                trace%point(ipt)%rays(idir)%NSegm = isegm
                if(isegm.gt.0) then
                    call allocTracePointRayData(trace%point(ipt)%rays(idir),isegm)
                    do j = 1,isegm
                        trace%point(ipt)%rays(idir)%icell(j)  = rayData%icell(j)
                        trace%point(ipt)%rays(idir)%iface(j)  = rayData%iface(j)
                        trace%point(ipt)%rays(idir)%ptOut(j)  = rayData%ptOut(j)
                        trace%point(ipt)%rays(idir)%segLen(j) = rayData%segLen(j)
                    enddo
                else
                    call allocTracePointRayData(trace%point(ipt)%rays(idir),1)
                    trace%point(ipt)%rays(idir)%icell(1)  = 0
                    trace%point(ipt)%rays(idir)%iface(1)  = 0
                    trace%point(ipt)%rays(idir)%ptOut(1)%crd  = ptc
                    trace%point(ipt)%rays(idir)%segLen(1) = 0.e0_r8p
                endif


                ! calculate tracer size:
                RaySizeBytes = isegm * SegmSizeBytes + real(sizeof(rayData),8)
                TracerSizeBytes = TracerSizeBytes + RaySizeBytes
                TracerSizeMb    = TracerSizeBytes * 1.e-6
            enddo

            write(*,"(a)",advance="no") repeat(char(8),27)
            write(*,"(a,1pE10.3e2,a)",advance="no") " TracerSize = ",TracerSizeMb," Mb"

!            write(*,*)
!            pause
        enddo

        write(*,*)
        write(*,*) "tracing done"
!        pause

        deallocate(rayData%icell)
        deallocate(rayData%iface)
        deallocate(rayData%ptOut)
        deallocate(rayData%segLen)

        return
    end subroutine new_rays_tracer_extended

!   --------------------------------------------
    subroutine new_rays_tracer_compact(rayType,pAngGrid,pMesh,locTit,elemType,trace)  !
        use module_OutputMesh
        implicit none
        ! input variables
        type(TAngularGridData),pointer :: pAngGrid        ! angular grid
        type(TMeshData),       pointer :: pMesh           ! mesh to trace
        character(*)                   :: locTit          ! mesh region name where tracing nodes exists
        character(*)                   :: elemType        ! mesh elements from where ray comes (Nodes/Faces/Cells)
        character(*)                   :: rayType         ! mesh elements from where ray comes (Nodes/Faces/Cells)
        ! setup parameters
        integer, parameter   :: maxNumSegm = 10000
        real(r8p), parameter :: val0 = 1.e-8_r8p

        ! output variables
        type(TTraceDataCompact) :: trace           ! tracing data

        ! local variables
        integer                 :: NPoints,NDirs
        type(TRayData)          :: rayData
        real(r8p)               :: segLen

        real(r8p)               :: cosTh,cosPhi,sinTh,sinPhi
        real(r8p)               :: theta,phi
        real(r8p), dimension(3) :: pt0,pt1,ptc,ptt
        real(r8p), dimension(3) :: rayDir,fnorm
        integer(I4P)            :: icell0,icell1
        integer(I4P)            :: iFace1
        integer                 :: icase
        integer :: iOz,iOx,iOy
        integer :: i,j
        integer :: ipt,iel,idir,isegm,Nsegm
        integer       :: locPos,locIdx
        character(32) :: locType
        character(64) :: sfmt
        !logical :: point_in_cell
        real(8) :: point_in_cell
        real(r8p) :: SegmSizeBytes,RaySizeBytes,TracerSizeBytes,TracerSizeMb

        !write(*,*) 'pNDirs = ',pAngGrid%NDirs

        trace%pmesh    => pMesh
        trace%pAngGrid => pAngGrid

        locPos  = 0
        locPos  = valpos(trace%pmesh%Locations%locTit,trim(locTit))   ! region position in list
        locIdx  = trace%pmesh%Locations(locPos)%locIdx
        trace%regIdx = locIdx
        if(locPos.eq.0) then
            write(*,*) 'error: tracing region not found '
            stop
        end if
        locType = trace%pmesh%Locations(locPos)%locType
!        write(*,*) 'locPos   = ',locPos
!        write(*,*) 'locType  = ',locType
!        write(*,*) 'locTitle = ',locTit
!        write(*,*) 'elemType = ',elemType
        if(index(elemType,'node').gt.0) then
            NPoints = trace%pmesh%Locations(locPos)%nnodes
        elseif(index(elemType,'face').gt.0) then
            NPoints = trace%pmesh%Locations(locPos)%NElems
        elseif(index(elemType,'cell').gt.0) then
            NPoints = trace%pmesh%Locations(locPos)%NElems
        else
            write(*,*) 'error: please set in config the elements to trace (node/face/cell) '
            stop
        endif

        trace%NPoints = NPoints
        trace%elemTp  = trim(elemType)
        trace%rayTp   = trim(rayType)
        allocate(trace%point(NPoints))
        NDirs = pAngGrid%NDirs
        do i = 1,trace%NPoints
            call allocTracePointRays(trace%point(i),trace%pAngGrid%NDirs)
        enddo

        iOx = 1        ! Oy
        iOy = 2        ! Oz
        iOz = 3        ! axis of symmetry (if 1 - Ox is axis)
        allocate(rayData%icell(maxNumSegm))
        allocate(rayData%iface(maxNumSegm))
        allocate(rayData%ptOut(maxNumSegm))
        allocate(rayData%segLen(maxNumSegm))


        !sfmt = '("Point(",i5,") = {",3(1pE12.4e2,","),1pE12.4e2,"};")'
        write(*,*)
        write(*,*) ' - - - - - - - - '
        write(*,*) ' tracing ... '

        Nsegm = 0
        TracerSizeBytes = 0.0e0
        TracerSizeMb    = 0.0e0
        SegmSizeBytes = 0.0e0
        SegmSizeBytes = SegmSizeBytes +   real(sizeof(rayData%icell(1)),8)
        SegmSizeBytes = SegmSizeBytes +   real(sizeof(rayData%iface(1)),8)
        SegmSizeBytes = SegmSizeBytes +   real(sizeof(rayData%segLen(1)),8)
        !write(*,*) "SegmSizeBytes = ",SegmSizeBytes
        !pause

        sfmt = '('//'a,'//'a,I3,a,'//'4x,a,i7,a,i7,'//'4x,a,1pE10.3e2,a'//')'

        do ipt = 1,NPoints
        !do ipt = 525,525
            !write(*,"(a,a,I3,a)",advance="no") achar(13),'done ',ipt*100/NPoints," %"
            !write(*,"(a,a,I3,a,i7,a,i7)",advance="no") achar(13),'done ',ipt*100/NPoints," %",ipt," Nsegm = ",Nsegm
            write(*,sfmt,advance="no") achar(13), &
                    'done ',ipt*100/NPoints," %", &
                    'point ',ipt,'/',NPoints,     &
                    'TracerSize = ',TracerSizeMb,' Mb'
            iel = 0
            ptc = 0.e0_r8p
            if(index(elemType,'node').gt.0) then           ! tracing nodes
                iel = pMesh%Locations(locPos)%inode(ipt)   ! node index
                ptc = pmesh%Node(iel)%crd
            elseif(index(elemType,'face').gt.0) then       ! tracing face centers
                iel = pMesh%Locations(locPos)%ielem(ipt)   ! face index
                ptc = pmesh%Face(iel)%center
            elseif(index(elemType,'cell').gt.0) then       ! tracing cell centers
                iel = pMesh%Locations(locPos)%ielem(ipt)   ! cell index
                ptc = pmesh%Cell(iel)%center
            endif

            trace%point(ipt)%pt0%crd = ptc
            trace%point(ipt)%bindId  = iel

            do idir = 1,trace%pAngGrid%NDirs
            !do idir = 64,64
            !do idir = 195,195
                !write(*,"(a,a,I3,a,I7,2x,a,I7)",advance="no") achar(13),'done ',ipt*100/NPoints," %",ipt,'idir = ',idir
                theta = trace%pAngGrid%theta(idir)
                phi   = trace%pAngGrid%phi(idir)
                cosTh  = cos(theta)
                cosPhi = cos(phi)
                sinTh  = sin(theta)
                sinPhi = sin(phi)

                ! ray direction vector
                rayDir(iOx) = sinTh*cosPhi
                rayDir(iOy) = sinTh*sinPhi
                rayDir(iOz) = cosTh

                pt0(iOx) = ptc(iOx) + 1.e2*val0*rayDir(iOx)
                pt0(iOy) = ptc(iOy) + 1.e2*val0*rayDir(iOy)
                pt0(iOz) = ptc(iOz) + 1.e2*val0*rayDir(iOz)

                icase = -1               ! 0/1/2 - ray is (not inside)/(inside)/ trace region
                if(index(elemType,'node').gt.0) then        ! tracing nodes
                    icell0 = -1
                    do i = 1,pmesh%Node(iel)%ncells
                        icell1 = pmesh%Node(iel)%icell(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !ptt(1) = pt0(1)                    ! only 2D-RZ case !!!!!!!!!!!!!!!!!!!! <!>
                        !ptt(2) = sqrt(pt0(2)**2+pt0(3)**2) ! only 2D-RZ case !!!!!!!!!!!!!!!!!!!! <!>
                        !ptt(3) = 0.0e0                     ! only 2D-RZ case !!!!!!!!!!!!!!!!!!!! <!>
                        ptt = pt0
                        call mapper_XYZ_to_RZ_dirX(ptt)
                        if(point_in_cell(pMesh,ptt,icell1,0.d0).gt.0.d0) then
                            !write(*,"(a,3(1x,1pE12.4e3))") 'ptt: ',ptt
                            icell0 = icell1
                            icase = 1
                            !write(*,*) " point in cell: ",icell0
                        end if                                            ! <!> reflection condition
                    end do
                elseif(index(elemType,'face').gt.0) then    ! tracing face centers
                    icell0 = -1
                    icell1 = pmesh%Face(iel)%icell(1)
                    if(point_in_cell(pMesh,pt0,icell1,0.d0).gt.0.d0) then
                        icell0 = icell1
                        icase = 1
                    end if
                elseif(index(elemType,'cell').gt.0) then    ! tracing cell centers
                    iCell0 = iel
                    icase = 1
                endif

!                if(icase.le.0) then
!                    write(*,"(a)") " tracing point is outside mesh: "
!                    write(*,"(3(1x,1pE10.3e2),a,3(1x,1pE10.3e2))") ptc," -> ",pt0
!                end if


                iFace1 = 0
                isegm  = 0
                rayData%NSegm = 0
                rayData%icell(:)     = 0
                rayData%iface(:)     = 0
                rayData%segLen(:)    = 0.e0

!                if(icase.eq.1) then
!                    open(31,file="cell.geo")
!                    call output_point_gmsh(31,pt0)
!                    call output_cell_gmsh(31,iCell0,pMesh)
!                    ptt = ptc
!                    call mapper_XYZ_to_RZ_dirX(ptt)
!                    pt1 = ptt+0.1*rayDir
!                    call mapper_XYZ_to_RZ_dirX(pt1)
!                    call output_line_gmsh(31,ptt,pt1,'initial ray direction')
!                endif

                do while(icase.ge.1)
                    isegm = isegm + 1
                    Nsegm = Nsegm + 1
                    if(pMesh%meshDim.eq.3) then
                        call tracer_3D(pMesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)
                    elseif(pMesh%meshDim.eq.2) then
                        call tracer_RZ(pMesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)
                    endif
                    if(iCell0.eq.iCell1) then
                        write(*,*)
                        write(*,*) 'ray is in the same cell ...'
                        stop
                    endif
                    !write(*,*) " isegm = ",isegm
                    !if(isegm.eq.8) pause
                    !if(isegm.eq.57) pause

                    if(trim(rayType).eq."short") then
                        icase = 0
                    endif

                    rayData%NSegm = isegm
                    rayData%icell(isegm) = iCell0
                    rayData%iface(isegm) = iFace1
                    rayData%ptOut(isegm)%crd = pt1
                    ptt = pt1-pt0
                    segLen = vectLength(ptt)
                    rayData%segLen(isegm) = real(segLen,r4p)
                    iCell0 = iCell1
                    pt0    = pt1

                enddo

!                if(icase.gt.-1) then
!                    close(31)
!                    pause
!                    write(*,"(a)",advance="no") repeat(char(8),16)
!                    write(*,"(a,i7)",advance="no") " Nsegm = ",Nsegm
!                endif

                rayData%NSegm = isegm
!                call single_ray_output(ipt,idir,ptc,rayData,pMesh,pAngGrid)
                trace%point(ipt)%rays(idir)%NSegm = isegm
                if(isegm.gt.0) then
                    call allocTracePointRayData(trace%point(ipt)%rays(idir),isegm)
                    do j = 1,isegm
                        trace%point(ipt)%rays(idir)%icell(j)  = rayData%icell(j)
                        trace%point(ipt)%rays(idir)%iface(j)  = rayData%iface(j)
                        !trace%point(ipt)%rays(idir)%ptOut(j)  = rayData%ptOut(j)
                        trace%point(ipt)%rays(idir)%segLen(j) = rayData%segLen(j)
                    enddo
                else
                    call allocTracePointRayData(trace%point(ipt)%rays(idir),1)
                    trace%point(ipt)%rays(idir)%icell(1)  = 0
                    trace%point(ipt)%rays(idir)%iface(1)  = 0
                    !trace%point(ipt)%rays(idir)%ptOut(1)%crd  = ptc
                    trace%point(ipt)%rays(idir)%segLen(1) = 0.e0_r8p
                endif


                ! calculate tracer size:
                RaySizeBytes = isegm * SegmSizeBytes + real(sizeof(rayData),8)
                TracerSizeBytes = TracerSizeBytes + RaySizeBytes
                TracerSizeMb    = TracerSizeBytes * 1.e-6
            enddo

            write(*,"(a)",advance="no") repeat(char(8),27)
            write(*,"(a,1pE10.3e2,a)",advance="no") " TracerSize = ",TracerSizeMb," Mb"

!            write(*,*)
!            pause
        enddo

        write(*,*)
        write(*,*) "tracing done"
!        pause

        deallocate(rayData%icell)
        deallocate(rayData%iface)
        deallocate(rayData%ptOut)
        deallocate(rayData%segLen)

        return
    end subroutine new_rays_tracer_compact



!   --------------------------------------------
!   --------------------------------------------
    subroutine new_rays_tracer_points(rayType,pAngGrid,pMesh,Npts,pts,trace)  !
        implicit none
        ! input variables
        type(TAngularGridData),pointer      :: pAngGrid  ! angular grid
        type(TMeshData),       pointer      :: pMesh     ! mesh to trace
        integer                             :: Npts
        type(T3DPoint), dimension(Npts) :: pts
        character(32)                       :: rayType   ! mesh elements from where ray comes (Nodes/Faces/Cells)
        ! setup parameters
        integer, parameter   :: maxNumSegm = 10000
        real(r8p), parameter :: val0 = 1.e-8_r8p

        ! output variables
        type(TTraceData)        :: trace           ! tracing data

        ! local variables
        integer                 :: NPoints
        integer                 :: NDirs
        type(TRayData), target  :: rayData
        real(r8p)               :: segLen

        real(r8p)               :: cosTh,cosPhi,sinTh,sinPhi
        real(r8p)               :: theta,phi
        real(r8p), dimension(3) :: pt0,pt1,ptc,ptt
        real(r8p), dimension(3) :: rayDir,fnorm
        integer(I4P)            :: icell0,icell1
        integer(I4P)            :: iFace1
        integer                 :: icase
        integer :: iOz,iOx,iOy
        integer :: i,j
        integer :: ipt,iel,idir,isegm
        character(32)  :: elemType
        character(128) :: sfmt
        real(8) :: point_in_cell


        NPoints  = Npts
        elemType = 'point'

        !write(*,*) 'pNDirs = ',pAngGrid%NDirs

        trace%pmesh    => pMesh
        trace%pAngGrid => pAngGrid


        trace%NPoints = NPoints
        trace%elemTp  = 'point'
        trace%rayTp   = trim(rayType)
        allocate(trace%point(NPoints))
        NDirs = pAngGrid%NDirs
        do i = 1,trace%NPoints
            call allocTracePointRays(trace%point(i),trace%pAngGrid%NDirs)
        enddo

        iOx = 1        ! Oy
        iOy = 2        ! Oz
        iOz = 3        ! axis of symmetry (if 1 - Ox is axis)
        allocate(rayData%icell(maxNumSegm))
        allocate(rayData%iface(maxNumSegm))
        allocate(rayData%ptOut(maxNumSegm))
        allocate(rayData%segLen(maxNumSegm))



        !sfmt = '("Point(",i5,") = {",3(1pE12.4e2,","),1pE12.4e2,"};")'
        write(*,*)
        write(*,*) ' - - - - - - - - '
        write(*,*) ' tracing ... '

        sfmt = "(a,'done ',I3,' %','    ipt = ',I7,2x,'idir = ',I7)"
        do ipt = 1,NPoints
!        do ipt = 1,1
!            write(*,"(a,a,I3,a)",advance="no") achar(13),'done ',ipt*100/NPoints," %"
            iel = 0
            ptc(:) = pts(ipt)%crd(:)


            trace%point(ipt)%pt0%crd = ptc
            trace%point(ipt)%bindId  = iel

            do idir = 1,trace%pAngGrid%NDirs
                write(*,sfmt,advance="no") achar(13),ipt*100/NPoints,ipt,idir
                theta = trace%pAngGrid%theta(idir)
                phi   = trace%pAngGrid%phi(idir)
                cosTh  = cos(theta)
                cosPhi = cos(phi)
                sinTh  = sin(theta)
                sinPhi = sin(phi)

                ! ray direction vector
                rayDir(iOx) = sinTh*cosPhi
                rayDir(iOy) = sinTh*sinPhi
                rayDir(iOz) = cosTh


                pt0(iOx) = ptc(iOx) + 1.e2*val0*rayDir(iOx)
                pt0(iOy) = ptc(iOy) + 1.e2*val0*rayDir(iOy)
                pt0(iOz) = ptc(iOz) + 1.e2*val0*rayDir(iOz)

                icase = 0               ! 0/1/2 - ray is (not inside)/(inside)/ trace region
                icell0 = -1
                do icell1 = 1,pMesh%NCells
                    if(point_in_cell(pMesh,pt0,icell1,0.d0).gt.0.d0) then
                        icell0 = icell1
                        icase = 1
                    end if
                end do


!                if(icase.le.0) then
!                    write(*,*)
!                    write(*,*) " tracing point is outside mesh"
!                    write(*,"(a,3(1x,1pE12.5e2))") " pt0: ",pt0
!                    pause
!                end if





                iFace1 = 0
                isegm  = 0
                rayData%NSegm = 0
                rayData%icell(:)     = 0
                rayData%iface(:)     = 0
                rayData%segLen(:)    = 0.e0

                do while(icase.ge.1)
                    isegm = isegm + 1
                    if(pMesh%meshDim.eq.3) then
                        call tracer_3D(pMesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)
                    elseif(pMesh%meshDim.eq.2) then
                        call tracer_RZ(pMesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)
                    endif
                    if(iCell0.eq.iCell1) then
                        write(*,*)
                        write(*,*) 'ray is in the same cell ...'
                        stop
                    endif

                    if(trim(rayType).eq."short") then
                        icase = 0
                    endif

                    rayData%NSegm = isegm
                    rayData%icell(isegm) = iCell0
                    rayData%iface(isegm) = iFace1
                    rayData%ptOut(isegm)%crd = pt1
                    ptt = pt1-pt0
                    segLen = vectLength(ptt)
                    rayData%segLen(isegm) = real(segLen,r4p)
                    iCell0 = iCell1
                    pt0    = pt1

!                    if(icase.eq.0) then
!                        if(pMesh%Face(iFace1)%locIdx.ne.12) then
!                            iRefl = iRefl + 1
!                        endif
!                        eps = vectScalProd(rayDir,fnorm)
!                        rayNu  = eps*fnorm
!                        rayTau = rayDir - rayNu
!                        rayDir = rayTau - rayNu
!                        iCell1 = iCell0
!                        icase = 1
!                    end if
!
!                    if(iRefl.lt.7 .and. icase.eq.0) then
!                    rayLen = rayLen + segLen*dimLen
!                    if(rayLen.ge.4.d5) then
!                        icase = 0
!                    end if

                enddo
                rayData%NSegm = isegm
!                call single_ray_output(ipt,idir,ptc,rayData,pMesh,pAngGrid)
                trace%point(ipt)%rays(idir)%NSegm = isegm
                if(isegm.gt.0) then
                    call allocTracePointRayData(trace%point(ipt)%rays(idir),isegm)
                    do j = 1,isegm
                        trace%point(ipt)%rays(idir)%icell(j)  = rayData%icell(j)
                        trace%point(ipt)%rays(idir)%iface(j)  = rayData%iface(j)
                        trace%point(ipt)%rays(idir)%ptOut(j)  = rayData%ptOut(j)
                        trace%point(ipt)%rays(idir)%segLen(j) = rayData%segLen(j)
                    enddo
                else
                    call allocTracePointRayData(trace%point(ipt)%rays(idir),1)
                    trace%point(ipt)%rays(idir)%icell(1)  = 0
                    trace%point(ipt)%rays(idir)%iface(1)  = 0
                    trace%point(ipt)%rays(idir)%ptOut(1)%crd  = ptc
                    trace%point(ipt)%rays(idir)%segLen(1) = 0.e0_r8p
                endif
            enddo
        enddo

        write(*,*)
        write(*,"(a)") " tracing done"
!        pause

        deallocate(rayData%icell)
        deallocate(rayData%iface)
        deallocate(rayData%ptOut)
        deallocate(rayData%segLen)

        return
    end subroutine new_rays_tracer_points

!   --------------------------------------------
    subroutine ray_intersect_face(pt0,RayDir,FaceNodes,NNodes,pt1,FaceNorm,RelDist)  !
        implicit none
        ! input:
        real(r8p), dimension(3)              :: pt0          ! ray start point
        real(r8p), dimension(3)              :: RayDir       ! ray direction
        real(r8p), dimension(NumFaceNodes,3) :: FaceNodes    ! face points coordinates
        integer                              :: NNodes       ! number of face nodes
        ! output:
        real(r8p), dimension(3)        :: pt1          ! ray/face cross point
        real(r8p), dimension(3)        :: FaceNorm     ! face positive normal (clockwise)
        real(r8p)                      :: RelDist      ! relative minimal bound distance (+ if inside, - if otside)
        ! parameters:
        real(r8p), parameter :: val0 = 1.e-8_r8p
        ! local:
        logical                 :: cross
        integer                 :: itri,inode
        real(r8p), dimension(3) :: nd1,nd2,nd3         ! trinlge nodal points
        real(r8p), dimension(3) :: vec1,vec2,vec3      ! trinlge edge vectors
        real(r8p)               :: len1,len2,len3      ! trinlge edge vectors length
        real(r8p), dimension(3) :: cvec,qvec           ! service vectors
        real(r8p), dimension(3) :: dvec01              ! delta vector pt0 -> nd1
        real(r8p)               :: dlen01              ! delta vector pt0 -> nd1 length
        real(r8p)               :: SclPr1,SclPr2       ! scalar products
        real(r8p)               :: FaceSize,diam,eps
        real(r8p), dimension(3,3) :: TriNodes          ! triangle points coordinates
        real(r8p)                 :: t                 ! parameter for ray line: ray(t) = RayDir*t + p0
        real(r8p), dimension(3) :: dist,proj

        ! get face character dimension:
        FaceSize = 0.d0
        do itri = 1,NNodes-2          ! loop through tringles (from face spliting)
            TriNodes(1,:) = FaceNodes(1,:)
            TriNodes(2,:) = FaceNodes(itri+1,:)
            TriNodes(3,:) = FaceNodes(itri+2,:)
            diam = 2*circle_inside_triangle(TriNodes)
            FaceSize = FaceSize + diam
        enddo


        cross   = .false.
        t       = 0.d0
        ! find intersection:
        itri = 0

        do while((.not.cross) .and. (itri.lt.NNodes-2))
            itri   = itri + 1
            ! set triangle:
            nd1(:) = FaceNodes(1,:)
            nd2(:) = FaceNodes(itri+1,:)
            nd3(:) = FaceNodes(itri+2,:)
            ! set triangle normal and edges vectors
            vec1 = nd2-nd1
            len1 = vectLength(vec1)
            vec2 = nd3-nd2
            len2 = vectLength(vec2)
            vec3 = nd1-nd3
            len3 = vectLength(vec3)

            FaceNorm = vectVectProd(vec1,vec2)      ! reversed clockwise (positive, like in mesh module)
            eps      = vectLength(FaceNorm)
            FaceNorm = FaceNorm/eps

            ! set delta vector: pt0 -> nd1:
            dvec01 = nd1 - pt0
            dlen01 = vectLength(dvec01)

            ! get scalar products:
            SclPr1 = vectScalProd(FaceNorm,dvec01)
            SclPr2 = vectScalProd(FaceNorm,RayDir)
            ! check if ray is parallel to triangle:
            t = 0.d0
            if (abs(SclPr2).gt.val0) then
                t = SclPr1/SclPr2
            end if

            if(t.gt.val0) then
                pt1   = pt0 + RayDir*t

                ! check if pt1 is inside face
                qvec = pt1 - nd1
                cvec = vectVectProd(vec1,qvec)/len1
                dist(1) = vectLength(cvec)
                proj(1) = vectScalProd(cvec,FaceNorm)

                qvec = pt1 - nd2
                cvec = vectVectProd(vec2,qvec)/len2
                dist(2) = vectLength(cvec)
                proj(2) = vectScalProd(cvec,FaceNorm)

                qvec = pt1 - nd3
                cvec = vectVectProd(vec3,qvec)/len3
                dist(3) = vectLength(cvec)
                proj(3) = vectScalProd(cvec,FaceNorm)

                RelDist = min(dist(1),dist(2),dist(3))/FaceSize
                eps = 0.d0
                do inode = 1,3
                    if(abs(proj(inode)).gt.val0) then
                        eps = proj(inode)
                    end if
                end do

                ! check if point is inside
                cross = .true.
                do inode = 1,3
                    if(eps*proj(inode).lt.0.d0) then
                        cross = .false.
                    end if
                end do

            end if
        enddo

        ! find minimal relative distance
        RelDist = -1.d0
        if(cross) then
            RelDist = 1.d0           ! 1 = FaceSize/FaceSize (large)
            itri    = 1
            do inode = 1,NNodes
                nd1(:) = FaceNodes(inode,:)
                if(inode.lt.NNodes) then
                    nd2(:) = FaceNodes(inode+1,:)
                else
                    nd2(:) = FaceNodes(1,:)
                endif
                vec1   = nd2-nd1
                eps = point_line_distance(nd1,vec1,pt1)
                eps = eps/FaceSize
                if(eps.lt.RelDist) then
                    RelDist = eps
                    itri = inode/3 + 1
                end if
            end do

            nd1(:) = FaceNodes(1,:)
            nd2(:) = FaceNodes(itri+1,:)
            nd3(:) = FaceNodes(itri+2,:)

            if(RelDist.lt.1.d-4) then
                vec1 = (nd1+nd2+nd3)/3.d0      ! triangle center
                vec2 = vec1 - pt1
                len2 = vectLength(vec2)
!                pt1 = pt1 + 1.d4*val0*vec2
!                RelDist = 1.d4*val0
                pt1 = pt1 + 1.d-4*vec2

                RelDist = 1.d-4
            end if

        end if

        return
    end subroutine

!   --------------------------------------------
    function circle_inside_triangle(TriNodes) result(R)
        implicit none
        ! input:
        real(r8p), dimension(3,3) :: TriNodes    ! triangle points coordinates
        ! output:
        real(r8p)                 :: R
        ! local:
        real(8)                 :: p,len1,len2,len3
        real(r8p), dimension(3) :: vec    ! edge vector

        vec(:) = TriNodes(2,:)-TriNodes(1,:)
        len1   = vectLength(vec)
        vec(:) = TriNodes(3,:)-TriNodes(2,:)
        len2   = vectLength(vec)
        vec(:) = TriNodes(1,:)-TriNodes(3,:)
        len3   = vectLength(vec)
        p = 0.5*(len1 + len2 + len3)
        R = sqrt((p-len1)*(p-len2)*(p-len3)/p)
        return
    end function

!   --------------------------------------------
    subroutine tracer_3D(pmesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)  !
        use module_OutputMesh
        implicit none
        ! input parameters
        type(TMeshData), pointer :: pmesh
        real(r8p), dimension(3)  :: rayDir
        real(r8p), dimension(3)  :: pt0
        integer(I4P)             :: iCell0

        ! output parameters
        real(r8p), dimension(3) :: pt1
        integer(I4P)            :: iCell1
        integer(I4P)            :: iFace1
        real(r8p), dimension(3) :: fnorm         ! face: normal vector
        integer                 :: icase

        ! local parameters
        real(r8p), dimension(NumFaceNodes,3) :: FaceNodes
        real(r8p), dimension(3) :: ptc
        real(r8p) :: RelDist,vproj

        integer :: iFace
        integer :: i,j
        integer :: n1,nnodes
        logical :: cross
        real(r8p), parameter :: val0 = 1.e-8_r8p

        iCell1    = iCell0

        do i = 1,pmesh%Cell(iCell0)%nFaces
            iFace = pmesh%Cell(iCell0)%iface(i)

            ! ----------------- test: ---------------------
            nnodes = pmesh%Face(iFace)%nNodes
            do j = 1,nnodes
                n1 = pmesh%Face(iFace)%inode(j)
                FaceNodes(j,:) = pmesh%Node(n1)%crd(:)
            end do

            RelDist = 0.d0
            call ray_intersect_face(pt0,rayDir,FaceNodes,nnodes,ptc,fnorm,RelDist)
            fnorm = fnorm*pmesh%Cell(iCell0)%iFaceDir(i) ! external for iCell0 normal vector
            vproj = vectScalProd(fnorm,rayDir)

            if(RelDist.gt.0.d0 .and. vproj.gt.0.d0) then
                cross = .true.
                pt1 = ptc
                iCell1 = pmesh%Cell(iCell0)%ineighbor(i)
                iFace1 = iFace
                if(iCell1.lt.0) then
                    icase = 0
                endif
                return
            else
                continue
            end if
        enddo


!        ! correction for small values
!        do i = 1,3
!            if(abs(pt1(i)).le.1.d-8) then
!                pt1(i) = 0.e0
!            end if
!        enddo

        return
    end subroutine tracer_3D


!   --------------------------------------------
    subroutine tracer_RZ(pmesh,rayDir,pt0,iCell0,pt1,iCell1,iFace1,fnorm,icase)  !                             !
        use module_OutputMesh
        use module_Interpolate
        implicit none
        ! input parameters
        type(TMeshData), pointer :: pmesh
        real(r8p), dimension(3)  :: rayDir
        real(r8p), dimension(3)  :: pt0
        integer(I4P)             :: iCell0
        ! output parameters
        real(r8p), dimension(3) :: pt1
        integer(I4P)            :: iCell1
        integer(I4P)            :: iFace1
        real(r8p), dimension(3) :: fnorm           ! face: normal vector
        integer                 :: icase

        ! local parameters
        real(r8p), dimension(3) :: ptc,ptn1,ptn2
        real(r8p), dimension(3) :: vec,ptT
        integer(I4P)            :: iFaceT,iCellT    ! temporal face and cell
        integer(I4P)            :: locIdx,idx
        real(r8p), dimension(3) :: Kr,Br            ! ray line parameters: k*t + b   (in 3D)
        real(r8p), dimension(3) :: tau,ort
        real(r8p), dimension(2) :: Kf,Bf            ! face line parameters: k*t + b  (in 2D)
        real(r8p)               :: Ks,Bs
        real(r8p)               :: a,b,c,D
        real(r8p) :: t1,t2
        real(r8p) :: t,eps
        real(r8p) :: r1
        real(r8p) :: SegLenT,SegLen1            ! segment length

        integer :: iFace,iOrient
        integer :: i,inbr
        integer :: n1,n2
        integer :: iOz,iRx,iRy
        real(r8p), parameter    :: val0 = 1.e-10_r8p
!        logical :: point_in_cell,lres

        ! initial ussumption: free ray/face/Oz orientation
        iRx = 2        ! Oy
        iRy = 3        ! Oz
        iOz = 1        ! axis of symmetry (if 1, then Ox is axis)

        ! coefficients in parametric ray expression
        Kr(1) = rayDir(iRx)         ! ray equation: K - coefficient: x(t) = Kx*t+Bx
        Kr(2) = rayDir(iRy)         ! ray equation: K - coefficient: y(t) = Ky*t+By
        Kr(3) = rayDir(iOz)         ! ray equation: K - coefficient: z(t) = Kz*t+Bz
        Br(1) = pt0(iRx)
        Br(2) = pt0(iRy)
        Br(3) = pt0(iOz)

        ! - - - - -
        iFace1 = 0
        iCell1 = 0
        SegLen1 = 1.0e0/val0
        idx    = 0
        Kf(:) = 0.0e0
        Bf(:) = 0.0e0
        Ks    = 0.0e0
        ort(:) = 0.0e0
        D     = 0.0e0; a = 0.0e0; b = 0.0e0; c = 0.0e0;
        do i = 1,pmesh%Cell(iCell0)%nFaces
            iFaceT = 0
            iCellT = 0
            iFace = pmesh%Cell(iCell0)%iface(i)
            fnorm(1) = pmesh%Face(iFace)%vnorm(iRx)
            fnorm(2) = pmesh%Face(iFace)%vnorm(iRy)
            fnorm(3) = pmesh%Face(iFace)%vnorm(iOz)
            iOrient  = pmesh%Cell(iCell0)%iFaceDir(i)
            fnorm    = fnorm*iOrient
            n1   = pmesh%Face(iFace)%inode(1)
            ptn1 = pmesh%Node(n1)%crd
            n2   = pmesh%Face(iFace)%inode(2)
            ptn2 = pmesh%Node(n2)%crd

            ! face parameters:
            Kf(1) = ptn2(iOz) - ptn1(iOz)                                  ! Kfz
            eps = ptn2(iRx)*ptn2(iRx)+ptn2(iRy)*ptn2(iRy)
            Kf(2) = sqrt(eps)
            eps = ptn1(iRx)*ptn1(iRx)+ptn1(iRy)*ptn1(iRy)
            Kf(2) = Kf(2) - sqrt(eps)                                      ! Kfr
            Bf(1) = ptn1(iOz)                                              ! Bfz
            Bf(2) = sqrt(eps)                                              ! Bfr

            ! ---------------------------------------
            ! 1. case
            if(abs(Kf(1)).le.val0 .and. abs(Kr(3)).gt.val0) then       ! Face |- Oz, Ray |\ Oz
                !write(*,*) ' *** case - 1'
                t = (Bf(1)-Br(3))/Kr(3)
                ptc(iRx) = Kr(1)*t + Br(1)
                ptc(iRy) = Kr(2)*t + Br(2)
                ptc(iOz) = Kr(3)*t + Br(3)
                r1 = sqrt(ptc(iRx)**2+ptc(iRy)**2)
                t1 = sqrt(ptn1(iRx)**2+ptn1(iRy)**2)
                t2 = sqrt(ptn2(iRx)**2+ptn2(iRy)**2)
                eps = (t2-r1)*(t1-r1)
                inbr = pmesh%Cell(iCell0)%ineighbor(i)
                !write(*,"(a,2(1x,1pe12.4e3),a,i7)") 'eps, t = ',eps,t,' NbrCell = ',inbr
                !write(*,"(a,2(1x,1pe12.4e3),a,i7)") 'eps, t = ',eps,t,' NbrCell = ',inbr+pmesh%NBFaces
                if(eps.le.0.e0 .and. t.ge.val0) then
                    ptT = ptc
                    iCellT = pmesh%Cell(iCell0)%ineighbor(i)
                    iFaceT = iFace
                    idx    = i
                    locIdx = pmesh%Face(iFace)%locIdx
                    if(locIdx.gt.0) then
                        icase = 0
                    endif
                endif
            ! 2. case
            elseif(abs(Kf(1)).le.val0 .and. abs(Kr(3)).le.val0) then   ! Face |- Oz, Ray |- Oz
                !write(*,*) ' *** case - 2'
                if(Bf(1).eq.Br(3)) then
                    ptT = ptc
                    iCellT = pmesh%Cell(iCell0)%ineighbor(i)
                    iFaceT = iFace
                    idx    = i
                    inbr = pmesh%Cell(iCell0)%ineighbor(i)
                    !write(*,"(a,1(1x,1pe12.4e3),a,i7)") '-> t_2 = ',t1,' NbrCell = ',inbr
                    !write(*,"(a,1(1x,1pe12.4e3),a,i7)") '-> t_2 = ',t1,' NbrCell = ',inbr+pmesh%NBFaces
                    locIdx = pmesh%Face(iFace)%locIdx
                    !write(*,*) 'bound: ',pmesh%Face(iFace)%locIdx
                    if(locIdx.gt.0) then
                        icase = 0
                    endif
                endif
            ! 3. case
            else                                                       ! Ray |- Oz, Face || Oz
                !write(*,*) ' *** case - 3'
                Ks = Kr(3)/Kf(1)
                Bs = (Br(3)-Bf(1))/Kf(1)
                a = Kr(1)**2 + Kr(2)**2 - (Kf(2)*Ks)**2
                b = 2.e0*(Kr(1)*Br(1)+Kr(2)*Br(2)-Kf(2)**2*Ks*Bs-Kf(2)*Ks*Bf(2))
                c = Br(1)**2 + Br(2)**2 - (Kf(2)*Bs+Bf(2))**2
                D = b**2 - 4.e0*a*c
                !write(*,"(a,4(1x,1pE12.4e3))") 'D,a,b,c: ',D,a,b,c
                if(D.ge.0.e0) then
                    t1 = 0.0e0
                    t2 = 0.0e0
                    if(abs(a).gt.val0) then
                        eps = sqrt(D)
                        t1 = 0.5e0_r8p*(-b-eps)/a
                        t2 = 0.5e0_r8p*(-b+eps)/a
                    elseif(abs(b).gt.val0) then
                        t1 = -c/b
                        t2 = -1.e3
                    endif
                    SegLenT = 1.0e0/val0
                    ptc(iRx) = Kr(1)*t1 + Br(1)
                    ptc(iRy) = Kr(2)*t1 + Br(2)
                    ptc(iOz) = Kr(3)*t1 + Br(3)
                    eps = (ptn1(iOz)-ptc(iOz))*(ptn2(iOz)-ptc(iOz))
                    inbr = pmesh%Cell(iCell0)%ineighbor(i)
                    !write(*,"(a,2(1x,1pe12.4e3),a,i7)") 'eps, t1: ',eps,t1,' NbrCell = ',inbr
                    !write(*,"(a,2(1x,1pe12.4e3),a,i7)") 'eps, t1: ',eps,t1,' NbrCell = ',inbr+pmesh%NBFaces
                    if(eps.le.0.e0 .and. t1.gt.val0) then
                        vec     = ptc-pt0
                        SegLenT = vectLength(vec)
                        ptT     = ptc
                        iCellT = pmesh%Cell(iCell0)%ineighbor(i)
                        iFaceT = iFace
                        idx    = i
                        locIdx = pmesh%Face(iFace)%locIdx
                        if(locIdx.gt.0) then
                            icase = 0
                        endif
                    endif

                    ptc(iRx) = Kr(1)*t2 + Br(1)
                    ptc(iRy) = Kr(2)*t2 + Br(2)
                    ptc(iOz) = Kr(3)*t2 + Br(3)
                    eps = (ptn1(iOz)-ptc(iOz))*(ptn2(iOz)-ptc(iOz))
                    inbr = pmesh%Cell(iCell0)%ineighbor(i)
                    !write(*,"(a,2(1x,1pe12.4e3),a,i7)") 'eps, t2: ',eps,t2,' NbrCell = ',inbr
                    !write(*,"(a,2(1x,1pe12.4e3),a,i7)") 'eps, t2: ',eps,t2,' NbrCell = ',inbr+pmesh%NBFaces
                    if(eps.le.0.e0 .and. t2.gt.val0) then
                        vec     = ptc-pt0
                        if(vectLength(vec).lt.SegLenT) then
                            ptT     = ptc
                            iCellT = pmesh%Cell(iCell0)%ineighbor(i)
                            iFaceT = iFace
                            idx    = i
                            locIdx = pmesh%Face(iFace)%locIdx
                            !write(*,*) " > cross = true"
                            if(locIdx.gt.0) then
                                icase = 0
                            endif
                        endif
                    endif
                endif
            endif

            ! select closest out point
            if(iFaceT.gt.0) then
                vec     = ptT-pt0
                SegLenT = vectLength(vec)
                !write(*,"(a,5(1x,1pe12.4e3))") ' SegLenT = ',SegLenT
                if(SegLenT.lt.SegLen1) then
                    iFace1 = iFaceT
                    iCell1 = iCellT
                    pt1    = ptT
                    SegLen1 = SegLenT
                    !write(*,*) " >>> cross selected"
                endif
            end if
        enddo

!        write(*,*)
!        write(*,*) " - - - - - - - - - - "
!        write(*,*) "cell0 -> cell1: ",iCell0,iCell1
!        write(*,*) "cell0 -> cell1: ",iCell0+pmesh%NBFaces,iCell1+pmesh%NBFaces

        ! calculate distance to shift far from the node
        eps = 1.0e8
        do i = 1,pmesh%Face(iFace1)%nNodes
            n1   = pmesh%Face(iFace1)%inode(i)
            if(i.lt.pmesh%Face(iFace1)%nNodes) then
                n2   = pmesh%Face(iFace1)%inode(i+1)
            else
                n2   = pmesh%Face(iFace1)%inode(1)
            end if
            ptc(1) = pt1(iOz)                         ! pt1 projection to the plane (z,y)
            ptc(2) = dsqrt(pt1(iRx)**2+pt1(iRy)**2)   ! pt1 projection to the plane (z,y)
            ptc(3) = 0.d0                             ! pt1 projection to the plane (z,y)
            ptn1 = pmesh%Node(n1)%crd
            ptn2 = pmesh%Node(n2)%crd
            !tau  = ptn2 - ptn1           ! face tangential vector
            tau  = pmesh%Face(iFace1)%center-ptc           ! face tangential vector to the center of the face
            t1   = vectLength(ptc-ptn1)
            t2   = vectLength(ptn2-ptn1)
            t1   = t1/t2
            eps = min(t1,eps)
            !write(*,"(a,1pE10.3e2)") "t1 = ",t1
            if(t1.lt.1.e-6) then
                !Kf(1) = pt1(iRx)/dsqrt(pt1(iRx)**2+pt1(iRy)**2+val0)
                !Kf(2) = pt1(iRy)/dsqrt(pt1(iRx)**2+pt1(iRy)**2+val0)
                ! tau rotation:
!                tau(iOz) = tau(iOz)
!                tau(iRx) = abs(tau(2)) * pt1(iRx)/dsqrt(pt1(iRx)**2+pt1(iRy)**2+val0)
!                tau(iRy) = abs(tau(2)) * pt1(iRy)/dsqrt(pt1(iRx)**2+pt1(iRy)**2+val0)
!                pt1      = pt1 - 1.0e-2*tau
!                pt1(iRx) = pt1(iRx) + 1.0e-6*Kr(1)
!                pt1(iRy) = pt1(iRy) + 1.0e-6*Kr(2)
!                pt1(iOz) = pt1(iOz) + 1.0e-6*Kr(3)


                ! ptc shift:
                !ptc = ptc + 1.0e-2*tau
                !ptc = pmesh%Node(n1)%crd
                ptc = ptc + 1.0d-4*tau                   ! shift along face
                if(iCell1.gt.0) then
                    ptn2 = pmesh%Cell(iCell1)%center
                    ptc = ptc + 1.0d-4*(ptn2-ptn1)   ! shift across face
                else
                    ort = pmesh%Face(iFace1)%vnorm * pmesh%Cell(iCell0)%iFaceDir(idx)
                    ptc = ptc + 1.0d-7*vectLength(tau)*ort   ! shift across face
                endif
                ! ptc rotate:
                ptn1(iOz) = ptc(1)
                ptn1(iRx) = ptc(2) * pt1(iRx)/dsqrt(pt1(iRx)**2+pt1(iRy)**2+val0)      ! <!>
                ptn1(iRy) = ptc(2) * pt1(iRy)/dsqrt(pt1(iRx)**2+pt1(iRy)**2+val0)      ! <!>

                pt1      = ptn1

                !write(*,*) "done"
            end if
        enddo


!        ! output:
!        if(iFace1.gt.0) then
!            call output_cell_gmsh(31,iCell0,pmesh)
!            ptn1 = pt0
!            call mapper_XYZ_to_RZ_dirX(ptn1)
!            ptn2 = pt1
!            call mapper_XYZ_to_RZ_dirX(ptn2)
!            call output_line_gmsh(31,ptn1,ptn2,'active cross point')
!
!!            ptn1 = pt0
!!            call mapper_XYZ_to_RZ_dirX(ptn1)
!!            ptn2 = pt0+0.1*rayDir
!!            call mapper_XYZ_to_RZ_dirX(ptn2)
!!            call output_line_gmsh(31,ptn1,ptn2,'active cross ray')
!        elseif(iFace1.eq.0) then
!            call output_cell_gmsh(31,iCell0,pmesh)
!            ptn1 = pt0
!            call mapper_XYZ_to_RZ_dirX(ptn1)
!            ptn2 = pt0+0.1*rayDir
!            call mapper_XYZ_to_RZ_dirX(ptn2)
!            call output_line_gmsh(31,ptn1,ptn2,'outsider ray')
!            pause
!        end if
!!        write(*,*)
!        !pause


        return
    end subroutine tracer_RZ


! -------------------------------------------------------------
! ============================ subroutines =============================
    subroutine calc_trace1D_prm(ftrace,Trace1DData,pMesh,pAngGrid,Npts,pts)
        use module_Settings
        use module_Mesh
        implicit none
        ! input:
        character(*)                        :: ftrace
        type(TTrace1DData)                  :: Trace1DData
        type(TMeshData),  pointer           :: pMesh         ! mesh to trace
        type(TAngularGridData),  pointer    :: pAngGrid         ! mesh to trace
        integer                             :: Npts
        type(T3DPoint), dimension(Npts) :: pts
        ! parameters:
        real(8), parameter :: val0 = 0.d-10
        real(8), parameter :: kp0  = 0.d0
        real(8), parameter :: kp1  = 1.d0
        real(8), parameter :: kp2  = 10.d0
        real(8), parameter :: kp3  = 100.d0
        ! local:
        ! --------------- tracer -----------------
        type(TTraceData), target         :: tracer
        type(TTraceData), pointer        :: pTracer
        ! ----------------------------------------
        character(32)  :: rayType
        integer, allocatable :: PtFaceDir(:),PtFaceSeg(:)
        integer :: NDirs,FaceIdx,NCells,NFaces
        integer :: ipt,iface,icell,iseg,idir,idx,ios
        integer :: icell0,icell1
        real(r8p)      :: SegLen,SegLenMin
        real(r8p)      :: theta,phi
        real(r8p),dimension(3) :: RayDir,AxisDir
        real(r8p)      :: f1_0,f1_1,f1_2,f1_3,g1_0,g1_1,g1_2,g1_3
        real(r8p)      :: f2_0,f2_1,f2_2,f2_3,g2_0,g2_1,g2_2,g2_3
        real(r8p)      :: D10,D21,D32,D30,hi,fi,a1,b1,c1,k1,a2,b2,c2,k2
        real(r8p)      :: dOm,dOmz,eps
        character(128) :: meshTit,sbuf
!        real(r8p)    :: kp
!        integer      :: j


        Trace1DData%LenUnit = 1.0d0
!        Trace1DData%pMesh   => pTrace%pMesh
        Trace1DData%pMesh   => pMesh

        NCells  = pMesh%NCells
        NFaces  = pMesh%NFaces
        call alloc_Trace1DData(Trace1DData,Npts,NFaces)

        NDirs = pAngGrid%NDirs
        allocate(PtFaceDir(NDirs),PtFaceSeg(NDirs))

        ! --------------------- ray tracing ----------------------
        write(*,*) ' tracer setup: '
        write(*,*) ' NDirs = ',pAngGrid%NDirs
        rayType  = "long"
        call new_rays_tracer_points(rayType,pAngGrid,pMesh,Npts,pts,tracer)
        pTracer => tracer
        !call output_angular_grid(pAngGrid,pt%crd)

        write(*,*)
        write(*,*) " calculating Trace1DParameters ... "
        do ipt = 1,Npts
!        do ipt = 1,1
!            write(*,"(a,a,I3,a,I5)",advance="no") achar(13),'done ',ipt*100/Npts," %"

            do FaceIdx = 1,pTracer%pMesh%NFaces
!            do FaceIdx = 383,383
                PtFaceDir(:) = 0
                PtFaceSeg(:) = 0
                NDirs        = 0
                SegLenMin    = 1.d10
                do idir = 1,pTracer%pAngGrid%NDirs
                    !write(*,"(a,a,I3,a,I7,2x,a,I7)",advance="no") achar(13),'done ',ipt*100/NPoints," %",ipt,'idir = ',idir
                    do iseg = 1,pTracer%point(ipt)%rays(idir)%NSegm
                        iface = pTracer%point(ipt)%rays(idir)%iface(iseg)
                        if(iface.eq.FaceIdx) then
                            NDirs = NDirs + 1
                            PtFaceDir(NDirs) = idir
                            PtFaceSeg(NDirs) = iseg
                            SegLen = pTracer%point(ipt)%rays(idir)%segLen(iseg)
                            if(SegLen.lt.SegLenMin) then
                                SegLenMin = SegLen
                            end if

                        end if
                    end do
                enddo
                SegLenMin = SegLenMin*0.99d0


!                open(34,file="AproxiExp.dat")
!                do j = 1,1000
!                kp = 0.0e0 + 0.1*(j-1)
!                f1_0 = 0.e0
!                f1_1 = 0.e0
!                f2_0 = 0.e0
!                f2_1 = 0.e0
!                do idx = 1,NDirs
!                    !write(*,"(a,2(1x,a,I5))",advance="no") achar(13),' ipt = ',ipt," FaceIdx = ",FaceIdx,' idir = ',idx
!                    idir = PtFaceDir(idx)
!                    iseg = PtFaceSeg(idx)
!                    SegLen = pTracer%point(ipt)%rays(idir)%segLen(iseg)
!                    icell  = pTracer%point(ipt)%rays(idir)%icell(iseg)
!
!                    icell0 = pTracer%point(ipt)%rays(idir)%icell(1)
!                    icell1 = pTracer%point(ipt)%rays(idir)%icell(iseg)
!
!                    f1_0 = f1_0  + 1.d0
!                    f1_1 = f1_1  + exp(-kp*(SegLen-0.0d0))
!!                    f1_1 = f1_1  + exp(-kp*(SegLen-SegLenMin))
!
!                    ! ray direction vector
!                    theta = pTracer%pAngGrid%theta(idir)
!                    phi   = pTracer%pAngGrid%phi(idir)
!                    RayDir(1) = sin(theta)*cos(phi)
!                    RayDir(2) = sin(theta)*sin(phi)
!                    RayDir(3) = cos(theta)
!                    AxisDir = [1.0d0,0.0d0,0.0d0]
!
!                    eps  = vectScalProd(RayDir,AxisDir)
!                    f2_0 = f2_0 + (-1)*eps
!                    f2_1 = f2_1 + exp(-kp*(SegLen-0.0d0))*(-1)*eps
!!                    f2_1 = f2_1 + exp(-kp*(SegLen-SegLenMin))*(-1)*eps
!                end do
!!                dOm  = f1_0*pTracer%pAngGrid%dOmega
!!                dOmz = f2_0*pTracer%pAngGrid%dOmega
!                dOm  = pTracer%pAngGrid%dOmega
!                dOmz = pTracer%pAngGrid%dOmega
!                write(34,"(5(1x,1pE12.5e2))") kp,dOm*f1_1,abs(dOmz*f2_1)
!                end do
!                close(34)





                ! set initial values:
                Trace1DData%point(ipt)%icell0(FaceIdx) = 0
                Trace1DData%point(ipt)%prm0(FaceIdx)%a = 0.d0    ! idir*exp(-1/(a*x+b)+c) if x=0 or idir=0 - no error
                Trace1DData%point(ipt)%prm0(FaceIdx)%b = 1.d0
                Trace1DData%point(ipt)%prm0(FaceIdx)%c = 1.d0
                Trace1DData%point(ipt)%prm0(FaceIdx)%k = 0.d0
                Trace1DData%point(ipt)%prm0(FaceIdx)%idir = 0
                !Trace1DData%point(ipt)%prm0(FaceIdx)%dOm  = 0.d0
                !Trace1DData%point(ipt)%prm0(FaceIdx)%dOmz = 0.d0
                Trace1DData%point(ipt)%icell1(FaceIdx) = 0
                Trace1DData%point(ipt)%prm1(FaceIdx)%a = 0.d0
                Trace1DData%point(ipt)%prm1(FaceIdx)%b = 1.d0
                Trace1DData%point(ipt)%prm1(FaceIdx)%c = 1.d0
                Trace1DData%point(ipt)%prm1(FaceIdx)%k = 0.d0
                Trace1DData%point(ipt)%prm1(FaceIdx)%idir  = 0
                !Trace1DData%point(ipt)%prm1(FaceIdx)%dOm  = 0.d0
                !Trace1DData%point(ipt)%prm1(FaceIdx)%dOmz = 0.d0

                f1_0 = 0.d0
                f1_1 = 0.d0
                f2_0 = 0.d0
                f2_1 = 0.d0
                f1_2 = 0.d0
                f2_2 = 0.d0
                f1_3 = 0.d0
                f2_3 = 0.d0
                icell0 = 0
                icell1 = 0
                !write(*,*)
                do idx = 1,NDirs
                    !write(*,"(a,2(1x,a,I5))",advance="no") achar(13),' ipt = ',ipt," FaceIdx = ",FaceIdx,' idir = ',idx
                    idir = PtFaceDir(idx)
                    iseg = PtFaceSeg(idx)
                    SegLen = pTracer%point(ipt)%rays(idir)%segLen(iseg)
                    icell  = pTracer%point(ipt)%rays(idir)%icell(iseg)

                    icell0 = pTracer%point(ipt)%rays(idir)%icell(1)
                    icell1 = pTracer%point(ipt)%rays(idir)%icell(iseg)

                    f1_0  = f1_0  + 1.d0
                    f1_1  = f1_1  + exp(-kp1*(SegLen-SegLenMin))
                    f1_2  = f1_2  + exp(-kp2*(SegLen-SegLenMin))
                    f1_3  = f1_3  + exp(-kp3*(SegLen-SegLenMin))

                    ! ray direction vector
                    theta = pTracer%pAngGrid%theta(idir)
                    phi   = pTracer%pAngGrid%phi(idir)
                    RayDir(1) = sin(theta)*cos(phi)
                    RayDir(2) = sin(theta)*sin(phi)
                    RayDir(3) = cos(theta)
                    AxisDir = [1.0d0,0.0d0,0.0d0]

                    eps  = vectScalProd(RayDir,AxisDir)
                    f2_0 = f2_0 + (-1)*eps
                    f2_1 = f2_1 + exp(-kp1*(SegLen-SegLenMin))*(-1)*eps
                    f2_2 = f2_2 + exp(-kp2*(SegLen-SegLenMin))*(-1)*eps
                    f2_3 = f2_3 + exp(-kp3*(SegLen-SegLenMin))*(-1)*eps
                end do

                dOm  = f1_0*pTracer%pAngGrid%dOmega
                dOmz = f2_0*pTracer%pAngGrid%dOmega

                if(NDirs.gt.0) then
                    Trace1DData%point(ipt)%icell0(FaceIdx) = icell0
                    Trace1DData%point(ipt)%icell1(FaceIdx) = icell1

                    ! set sign of the direction
                    Trace1DData%point(ipt)%prm0(FaceIdx)%idir = 1
                    if(f1_1.lt.1.d-80) then
                        f1_1 = 0.d0
                    end if
                    if(f1_2.lt.1.d-80) then
                        f1_2 = 0.d0
                    end if
                    if(f1_3.lt.1.d-80) then
                        f1_3 = 0.d0
                    end if

                    Trace1DData%point(ipt)%prm1(FaceIdx)%idir = 1
                    if(f2_0.lt.-1.d-80) then
                        f2_0 = abs(f2_0)
                        Trace1DData%point(ipt)%prm1(FaceIdx)%idir = -1
                    elseif(f2_0.lt.1.d-80) then
                        f2_0 = 0.d0
                        Trace1DData%point(ipt)%prm1(FaceIdx)%idir =  0
                    else
                        Trace1DData%point(ipt)%prm1(FaceIdx)%idir =  1
                    end if

                    if(f2_1.lt.-1.d-80) then
                        f2_1 = abs(f2_1)
                    elseif(f2_1.lt.1.d-80) then
                        f2_1 = 0.d0
                    end if

                    if(f2_2.lt.-1.d-80) then
                        f2_2 = abs(f2_2)
                    elseif(f2_2.lt.1.d-80) then
                        f2_2 = 0.d0
                    end if

                    if(f2_3.lt.-1.d-80) then
                        f2_3 = abs(f2_3)
                    elseif(f2_3.lt.1.d-80) then
                        f2_3 = 0.d0
                    end if

                    dOm  = pTracer%pAngGrid%dOmega
                    g1_0 =                  log(dOm) + log(f1_0)
                    g1_1 = -kp1*SegLenMin + log(dOm) + log(f1_1)
                    g1_2 = -kp2*SegLenMin + log(dOm) + log(f1_2)
                    g1_3 = -kp3*SegLenMin + log(dOm) + log(f1_3)
                    g2_0 =                  log(dOm) + log(f2_0)
                    g2_1 = -kp1*SegLenMin + log(dOm) + log(f2_1)
                    g2_2 = -kp2*SegLenMin + log(dOm) + log(f2_2)
                    g2_3 = -kp3*SegLenMin + log(dOm) + log(f2_3)

                    D10 = (g1_1-g1_0)/(kp1-kp0)
                    D21 = (g1_2-g1_1)/(kp2-kp1)
                    D32 = (g1_3-g1_2)/(kp3-kp2)

                    ! 1st approximation
                    hi = (D21-D10)/(kp2-kp0)
                    fi = (D32-D21)/(kp3-kp1)
                    if(abs(fi).gt.1.0e-20) then
                        b1 = g1_0
                        a1 = 1.0d0/(kp3-kp0)*(hi/fi-1.0d0)
                        k1 = -hi*(a1*kp1+1.0d0)*(a1*kp2+1.0d0)/a1
                        if(abs(a1*kp1+1.0d0).gt.1.0e-10) then
                            c1 = D10 - k1/(a1*kp1+1.0d0)
                        else
                            a1 = 0.0e0_r8p
                            k1 = 0.0e0_r8p
                            c1 = D10
                        endif
                    else
                        a1 = 0.0d0
                        k1 = 0.0d0
                        b1 = g1_0
                        c1 = (g1_3-g1_0)/(kp3-kp0)
                    end if

                    ! 2nd approximation
                    if(c1.gt.0 .or. a1.lt.0) then
                        D10 = (g1_1-g1_0)/(kp1-kp0)
                        D30 = (g1_3-g1_0)/(kp3-kp0)
                        a1  = (D10-D30)/(g1_3-g1_1)
                        k1  = D10*(a1*kp1+1.0d0)
                        c1  = 0.0d0
                    end if

                    ! 3rd approximation
                    if(a1.lt.0) then
                        D30 = (g1_3-g1_0)/(kp3-kp0)
                        a1  = 0.0d0
                        k1  = 0.0d0
                        c1  = D30
                        b1  = g1_0
                    end if

                    if(isnan(a1)) then
                        write(*,*) ' a1 = ',a1
                        !pause
                    end if

                    if(isnan(b1)) then
                        write(*,*) ' b1 = ',b1
                        !pause
                    end if

                    if(isnan(c1)) then
                        write(*,*) ' a1 = ',a1
                        write(*,*) ' k1 = ',k1
                        write(*,*) ' b1 = ',b1
                        write(*,*) ' c1 = ',c1
                        !pause
                    end if

                    if(isnan(k1)) then
                        write(*,*) ' k1 = ',k1
                        !pause
                    end if

                    Trace1DData%point(ipt)%prm0(FaceIdx)%a = a1
                    Trace1DData%point(ipt)%prm0(FaceIdx)%b = b1
                    Trace1DData%point(ipt)%prm0(FaceIdx)%c = c1
                    Trace1DData%point(ipt)%prm0(FaceIdx)%k = k1

                    D10 = (g2_1-g2_0)/(kp1-kp0)
                    D21 = (g2_2-g2_1)/(kp2-kp1)
                    D32 = (g2_3-g2_2)/(kp3-kp2)

                    ! 1st approximation
                    hi = (D21-D10)/(kp2-kp0)
                    fi = (D32-D21)/(kp3-kp1)
                    if(abs(fi).gt.1.0e-20) then
                        b2 = g2_0
                        a2 = 1.0d0/(kp3-kp0)*(hi/fi-1.0d0)
                        k2 = -hi*(a2*kp1+1.0d0)*(a2*kp2+1.0d0)/a2
                        if(abs(a2*kp1+1.0d0).gt.1.0e-10) then
                            c2 = D10 - k2/(a2*kp1+1.0d0)
                        else
                            a2 = 0.0e0_r8p
                            k2 = 0.0e0_r8p
                            c2 = D10
                        endif
                    else
                        a2 = 0.0d0
                        k2 = 0.0d0
                        b2 = g2_0
                        c2 = (g2_3-g2_0)/(kp3-kp0)
                    end if

                    ! 2nd approximation
                    if(c2.gt.0 .or. a2.lt.0) then
                        D10 = (g2_1-g2_0)/(kp1-kp0)
                        D30 = (g2_3-g2_0)/(kp3-kp0)
                        a2  = (D10-D30)/(g2_3-g2_1)
                        k2  = D10*(a2*kp1+1.0d0)
                        c2  = 0.0d0
                    end if

                    ! 3rd approximation
                    if(a2.lt.0) then
                        D30 = (g2_3-g2_0)/(kp3-kp0)
                        a2  = 0.0d0
                        k2  = 0.0d0
                        c2  = D30
                        b2  = g2_0
                    end if

                    if(isnan(a2)) then
                        write(*,*) ' a2 = ',a2
                        !pause
                    end if

                    if(isnan(b2)) then
                        write(*,*) ' b2 = ',b2
                        !pause
                    end if

                    if(isnan(c2)) then
                        write(*,*) ' a2 = ',a2
                        write(*,*) ' k2 = ',k2
                        write(*,*) ' b2 = ',b2
                        write(*,*) ' c2 = ',c2
                        !pause
                    end if

                    if(isnan(k2)) then
                        write(*,*) ' k2 = ',k2
                        !pause
                    end if

                    Trace1DData%point(ipt)%prm1(FaceIdx)%a = a2
                    Trace1DData%point(ipt)%prm1(FaceIdx)%b = b2
                    Trace1DData%point(ipt)%prm1(FaceIdx)%c = c2
                    Trace1DData%point(ipt)%prm1(FaceIdx)%k = k2

!                    open(34,file="AproxiExp.plt")
!                    write(34,"(a)")           '#set logscale x'
!                    write(34,"(a)")           'set logscale y'
!                    write(34,"(a)")           'set format y "10^{%L}"'
!                    write(34,"(a)")
!                    write(34,"(a)")           '# ------------'
!                    write(34,"(a,1pE12.4e3)") 'a1 = ',a1
!                    write(34,"(a,1pE12.4e3)") 'b1 = ',b1
!                    write(34,"(a,1pE12.4e3)") 'c1 = ',c1
!                    write(34,"(a,1pE12.4e3)") 'k1 = ',k1
!                    write(34,"(a)")           '# ------------'
!                    write(34,"(a,1pE12.4e3)") 'a2 = ',a2
!                    write(34,"(a,1pE12.4e3)") 'b2 = ',b2
!                    write(34,"(a,1pE12.4e3)") 'c2 = ',c2
!                    write(34,"(a,1pE12.4e3)") 'k2 = ',k2
!                    write(34,"(a)")
!                    write(34,"(a)")           'f1(x) = (k1*x)/(a1*x+1)+c1*x+b1'
!                    write(34,"(a)")           'f2(x) = (k2*x)/(a2*x+1)+c2*x+b2'
!                    write(34,"(a)")
!                    write(34,"(a)")           'set xra[0:10000]'
!                    write(34,"(a)")
!                    write(34,"(a)")           'plot "AproxiExp.dat" us 1:2 w l lw 2 ti "ApExp", \'
!                    write(34,"(a)")           '     "AproxiExp.dat" us 1:3 w l lw 2 ti "ApEint", \'
!                    write(34,"(a)")           '                 exp(f1(x)) w l ti "ApExp2", \'
!                    write(34,"(a)")           '                 exp(f2(x)) w l ti "ApEint2"'
!                    write(34,"(a)")
!                    write(34,"(a)")           'pause -1'
!                    close(34)
!
!                    pause

            endif

            end do
        end do
        write(*,*)

        deallocate(PtFaceDir,PtFaceSeg)

        ! output
        ! set mesh filename: actual
        meshTit = trim(pMesh%meshTit)
        idx = 1
        ipt = 1
        do while(idx.gt.0 .or. ipt.gt.0)
            idx = index(meshTit,"\")
            meshTit = meshTit(idx+1:)
            ipt = index(meshTit,"/")
            meshTit = meshTit(ipt+1:)
        end do

        open(3,                           &
            file       = ftrace,          &
            form       = 'UNFORMATTED',   &
            action     = 'WRITE',         &
            convert    = 'BIG_ENDIAN',    &
            access     = 'STREAM',        &
            iostat     = ios)

            sbuf = "meshTit: "//trim(meshTit)
            write(3,iostat=ios) trim(sbuf),char(10)
            write(sbuf,"(a,i5)") " Nth   = ",pTracer%pAngGrid%NTheta
            write(*,"(a,i5)") " Nth   = ",pTracer%pAngGrid%NTheta
            write(3,iostat=ios) trim(sbuf),char(10)
            write(sbuf,"(a,i5)") " Ndirs = ",pTracer%pAngGrid%NDirs
            write(*,"(a,i5)") " Ndirs = ",pTracer%pAngGrid%NDirs
            write(3,iostat=ios) trim(sbuf),char(10)

        !write(*,*) "output ..."
        do ipt = 1,Npts
            do iface = 1,NFaces
                write(3,iostat=ios) Trace1DData%point(ipt)%icell0(iface)
                write(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%a
                write(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%b
                write(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%c
                write(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%k
                write(3,iostat=ios) Trace1DData%point(ipt)%prm0(iface)%idir

                write(3,iostat=ios) Trace1DData%point(ipt)%icell1(iface)
                write(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%a
                write(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%b
                write(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%c
                write(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%k
                write(3,iostat=ios) Trace1DData%point(ipt)%prm1(iface)%idir
            end do
        end do
        close(3)

        pTracer => null()              ! deallocate all variables
        call close_rays_tracer(tracer)   ! deallocate all variables

        write(*,*) " calculating Trace1DParameters: done "

    end subroutine


! -------------------------------------------------------------
    subroutine tracer_rays_output(ifile,idir,trace)
!        use module_Settings
!        use module_Mesh
!        use module_RayTrace
        use module_TxtRead
        implicit none
        ! input parameters
        type(TTraceData), pointer :: trace
        integer :: ifile,idir,ipt
        ! local parameters
        real(r8p), dimension(3)  :: pt0,pt1
        integer :: icell,iface,inode
        integer :: isegm,idx,n1
        integer :: iOz,iOx,iOy
        real(r8p) :: theta,phi
        real(r8p) :: cosTh,cosPhi,sinTh,sinPhi
        character(32) :: sbuf

        iOx = 1        ! Oy
        iOy = 2        ! Oz
        iOz = 3        ! axis of symmetry (if 1 - Ox is axis)

        theta = trace%pAngGrid%theta(idir)
        phi   = trace%pAngGrid%phi(idir)
        cosTh  = cos(theta)
        cosPhi = cos(phi)
        sinTh  = sin(theta)
        sinPhi = sin(phi)

        pt0 = trace%point(1)%pt0%crd
        pt1(iOx) = pt0(iOx) + 1.d+1*sinTh*cosPhi
        pt1(iOy) = pt0(iOy) + 1.d+1*sinTh*sinPhi
        pt1(iOz) = pt0(iOz) + 1.d+1*cosTh


        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_direction-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(34,file=sbuf)
        write(34,"(3(1pE12.4e2))") pt0
        write(34,"(3(1pE12.4e2))") pt1
        close(34)

        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_points-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(35,file=sbuf)
        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_cells-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(36,file=sbuf)
        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_faces-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(37,file=sbuf)

        do ipt = 1,trace%NPoints
            ! ray points
            write(35,"(3(1pE12.4e2))") trace%point(ipt)%pt0%crd
            do isegm = 1,trace%point(ipt)%rays(idir)%NSegm
                write(35,"(3(1pE12.4e2))") trace%point(ipt)%rays(idir)%ptOut(isegm)%crd
            enddo
            write(35,"(3(1pE12.4e2))")
            write(35,"(3(1pE12.4e2))")

            do isegm = 1,trace%point(ipt)%rays(idir)%NSegm
                icell = trace%point(ipt)%rays(idir)%icell(isegm)
                do idx = 1,trace%pmesh%Cell(icell)%nFaces
                    iface = trace%pmesh%Cell(icell)%iface(idx)
                    do n1 = 1,trace%pmesh%Face(iface)%nNodes
                        inode = trace%pmesh%Face(iface)%inode(n1)
                        write(36,"(3(1pE12.4e2))") trace%pmesh%Node(inode)%crd
                    enddo
                    inode = trace%pmesh%Face(iface)%inode(1)
                    write(36,"(3(1pE12.4e2))") trace%pmesh%Node(inode)%crd
                    write(36,"(3(1pE12.4e2))")
                    write(36,"(3(1pE12.4e2))")
                enddo

                iface = trace%point(ipt)%rays(idir)%iface(isegm)
                do n1 = 1,trace%pmesh%Face(iface)%nNodes
                    inode = trace%pmesh%Face(iface)%inode(n1)
                    write(37,"(3(1pE12.4e2))") trace%pmesh%Node(inode)%crd
                enddo
                inode = trace%pmesh%Face(iface)%inode(1)
                write(37,"(3(1pE12.4e2))") trace%pmesh%Node(inode)%crd
                write(37,"(3(1pE12.4e2))")
                write(37,"(3(1pE12.4e2))")
            enddo
        enddo

        close(35)
        close(36)
        close(37)
        return
    end subroutine

! -------------------------------------------------------------
    subroutine single_ray_output(ifile,idir,pt0,ray,pMesh,pAngGrid)
!        use module_Settings
!        use module_Mesh
!        use module_RayTrace
        use module_TxtRead
        implicit none
        ! input parameters
        integer                  :: ifile,idir
        real(r8p), dimension(3)  :: pt0,pt1
        type(TRayData)           :: ray
        type(TMeshData), pointer :: pMesh         ! mesh to trace
        type(TAngularGridData),pointer :: pAngGrid      ! angular grid
        ! local parameters
        integer :: icell,iface,inode
        integer :: isegm,idx,n1
        integer :: iOz,iOx,iOy
        character(32) :: sbuf
        real(r8p) :: theta,phi
        real(r8p) :: cosTh,cosPhi,sinTh,sinPhi

        iOx = 1        ! Oy
        iOy = 2        ! Oz
        iOz = 3        ! axis of symmetry (if 1 - Ox is axis)

        theta = pAngGrid%theta(idir)
        phi   = pAngGrid%phi(idir)
        cosTh  = cos(theta)
        cosPhi = cos(phi)
        sinTh  = sin(theta)
        sinPhi = sin(phi)

        pt1(iOx) = pt0(iOx) + 1.d+1*sinTh*cosPhi
        pt1(iOy) = pt0(iOy) + 1.d+1*sinTh*sinPhi
        pt1(iOz) = pt0(iOz) + 1.d+1*cosTh


        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_direction-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(34,file=sbuf)
        write(34,"(3(1pE12.4e2))") pt0
        write(34,"(3(1pE12.4e2))") pt1
        close(34)


        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_points-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(35,file=sbuf)
        write(35,"(3(1pE12.4e2))") pt0
        do isegm = 1,ray%NSegm
            write(35,"(3(1pE12.4e2))") ray%ptOut(isegm)%crd
        enddo
        close(35)

        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_cells-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(36,file=sbuf)
        write(sbuf,"(a,a,i3.3,a)") trim(OutputDir),'\ray_faces-',ifile,'.dat'
        call slash_bck2frw(sbuf)
        open(37,file=sbuf)
        do isegm = 1,ray%NSegm
            icell = ray%icell(isegm)
            !write(*,*) 'icell = ',icell
            do idx = 1,pmesh%Cell(icell)%nFaces
                iface = pmesh%Cell(icell)%iface(idx)
                do n1 = 1,pmesh%Face(iface)%nNodes
                    inode = pmesh%Face(iface)%inode(n1)
                    write(36,"(3(1pE12.4e2))") pmesh%Node(inode)%crd
                enddo
                inode = pmesh%Face(iface)%inode(1)
                write(36,"(3(1pE12.4e2))") pmesh%Node(inode)%crd
                write(36,"(3(1pE12.4e2))")
                write(36,"(3(1pE12.4e2))")
            enddo

            iface = ray%iface(isegm)
            do n1 = 1,pmesh%Face(iface)%nNodes
                inode = pmesh%Face(iface)%inode(n1)
                write(37,"(3(1pE12.4e2))") pmesh%Node(inode)%crd
            enddo
            inode = pmesh%Face(iface)%inode(1)
            write(37,"(3(1pE12.4e2))") pmesh%Node(inode)%crd
            write(37,"(3(1pE12.4e2))")
            write(37,"(3(1pE12.4e2))")
        enddo
        flush(36)
        close(36)
        flush(37)
        close(37)
        return
    end subroutine

end module module_RayTrace
