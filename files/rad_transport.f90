subroutine rad_transport_short_chars(pTrace,CellKppa,CellEtta,Wflux,cUden)
    use module_Settings
    use module_Constants
    use module_Geom
    use module_RayTrace
    use module_GridFunct
    use module_AtomData
    use module_Spectrum
    use module_utils
    use module_TxtRead
    implicit none

    ! input values
    !type(TTraceData), pointer       :: pTrace
    !type(TGridNodalScalar), pointer :: pDens,pTemp
    type(TTraceData)       :: pTrace
    ! output variables
    type(TGridNodalVector) :: Wflux
    type(TGridNodalScalar) :: cUden
    ! local variables
    real(r8p)               :: I_ptX    ! radiation intencity      at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: dI_ptX   ! add radiation intencity  at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: cU_ptX   ! radiation energy density at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p),dimension(3)  :: W_ptX    ! radiation energy flux    at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p),dimension(pTrace%pAngGrid%NDirs)  :: I_Om ! radiation intencity at (r,:,e)  (r-coords,Om-direction, e-photon energy)

    ! - - - - - - - - - - - - - - - -
    ! local variables
    real(r8p)      :: cosTh,cosPhi
    real(r8p)      :: sinTh,sinPhi
    real(r8p)      :: theta,phi
    real(r8p)      :: dOm
    real(r8p)      :: I_bnd         ! radiation intencity boundary value
    real(r8p)      :: kappa,etta   ! radiation transport coefficients
    real(r8p)      :: tau,ctau     ! optical thickness
    integer(i4p)   :: i,idir,isegm,idx
    integer(i4p)   :: inode,iface,icell
    integer(i4p)   :: NPoints,NDirs
    integer(i4p)   :: iEn,NEnGroups
    real(r8p)      :: X1,X2,dX
    real(r8p)      :: wt
    ! Saha variables
    real(8), dimension(pTrace%pmesh%NCells,SpectrData%EnGrid%NX-1) :: CellKppa
    real(8), dimension(pTrace%pmesh%NCells,SpectrData%EnGrid%NX-1) :: CellEtta
    real(8) :: Nt0,TK0,L0
    real(8), dimension(pTrace%NPoints)   :: I_dir
    !logical :: fileExists

!    write(*,*) "ray%NSegm = ", ray%NSegm

    I_bnd = 0.e0_r8p
    I_ptX = 0.e0_r8p
    tau   = 0.e0_r8p
    ctau  = 0.e0_r8p
    kappa = 1.e0_r8p
    etta  = 1.e0_r8p

    TK0 = dimTem            ! K
    Nt0 = dimDen            ! cm-3
    L0  = dimLen            ! cm

    write(*,*)
    write(*,*) ' - - - - - - - - '
    write(*,*) ' radiation transport ... '

    NDirs   = pTrace%pAngGrid%NDirs
    NPoints = pTrace%NPoints
    NEnGroups = SpectrData%EnGrid%NX-1

    do idir = 1,NDirs                ! cycle by rays (directions)
        write(*,"(a,a,I3,a)",advance="no") achar(13),'done ',idir*100/NDirs," %"
        theta = pTrace%pAngGrid%theta(idir)
        phi   = pTrace%pAngGrid%phi(idir)
        cosTh  = cos(theta)
        cosPhi = cos(phi)
        sinTh  = sin(theta)
        sinPhi = sin(phi)
        dOm    = pTrace%pAngGrid%dOmega


        do iEn = 1,NEnGroups            ! cycle by spectral groups
            I_Om = 0.d0
            X1 = SpectrData%EnGrid%X(iEn)       ! current spectral group left bound
            X2 = SpectrData%EnGrid%X(iEn+1)     ! current spectral group right bound
            dX = X2 - X1                        ! current spectral group width

            I_dir(:) = 0.e0_r8p
            do isegm = 1,NPoints                  ! cycle by mesh points
                iface = pTrace%rayTree(idir)%segm(isegm)%iface
                icell = pTrace%rayTree(idir)%segm(isegm)%icell
                inode = pTrace%rayTree(idir)%segm(isegm)%inode

                I_ptX  = 0.e0_r8p
                dI_ptX = 0.e0_r8p
                if(iface.gt.0) then
                    kappa = CellKppa(icell,iEn)
                    etta  = CellEtta(icell,iEn)
                    ctau = kappa*pTrace%rayTree(idir)%segm(isegm)%segLen*L0

                    I_ptX = 0.e0_r8p
                    do i = 1,pTrace%pmesh%Face(iface)%nNodes
                        wt = pTrace%rayTree(idir)%segm(isegm)%weight(i)
                        idx = pTrace%pmesh%Face(iface)%inode(i)
                        I_ptX = I_ptX + wt*I_dir(idx)
                    end do
                    if(ctau.le.80.d0) then
                        I_ptX = I_ptX*exp(-ctau)
                    else
                        I_ptX = 0.e0_r8p
                    end if

                    dI_ptX = 0.d0
                    if(kappa.gt.0.d0) then
                        dI_ptX = etta/kappa
                    end if
                    if(ctau.le.80.d0) then
                        dI_ptX = dI_ptX*(1.e0_r8p-exp(-ctau))
                    end if
                endif
                I_dir(inode) = I_ptX + dI_ptX

                if(I_dir(inode).lt.1.d-200) then
                    I_dir(inode) = 0.d0
                end if

                ! integration by angle and energy
                cU_ptX   = I_dir(inode)*dOm * dX
                W_ptX(1) = I_dir(inode)*(-1)*sinTh*cosPhi * dOm * dX
                W_ptX(2) = I_dir(inode)*(-1)*sinTh*sinPhi * dOm * dX
                W_ptX(3) = 0.d0

                cUden%fval(inode) = cUden%fval(inode) + cU_ptX
                Wflux%fvec(inode)%cmp(:) = Wflux%fvec(inode)%cmp(:) + W_ptX(:)
            end do
        end do
    enddo
    write(*,*)

    return
end subroutine rad_transport_short_chars

! ----------------------------------------------------------------------
subroutine rad_transport_long_chars(pTrace,CellKppa,CellEtta,cUden,Wflux)
    use module_Settings
    use module_Constants
    use module_Geom
    use module_RayTrace
    use module_GridFunct
    use module_AtomData
    use module_Spectrum
    use module_utils
    use module_TxtRead
    implicit none

    ! input values
    type(TTraceData)       :: pTrace
    real(8), dimension(pTrace%pmesh%NCells,SpectrData%EnGrid%NX-1) :: CellKppa
    real(8), dimension(pTrace%pmesh%NCells,SpectrData%EnGrid%NX-1) :: CellEtta
    ! output variables
    type(TGridNodalScalar) :: cUden
    type(TGridNodalVector) :: Wflux
!    type(TGridNodalVector), optional :: cUsp
    ! local variables
    real(r8p)               :: I_ptX    ! radiation intencity      at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: dI_ptX   ! add radiation intencity  at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: cU_ptX   ! radiation energy density at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p),dimension(3)  :: W_ptX    ! radiation energy flux    at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p),dimension(pTrace%pAngGrid%NDirs)  :: I_Om ! radiation intencity at (r,:,e)  (r-coords,Om-direction, e-photon energy)

    ! - - - - - - - - - - - - - - - -
    ! local variables
    type(T3DPoint) :: pt0
    real(r8p)      :: cosTh,cosPhi
    real(r8p)      :: sinTh,sinPhi
    real(r8p)      :: theta,phi
    real(r8p)      :: dOm
    real(r8p)      :: I_bnd         ! radiation intencity boundary value
    real(r8p)      :: kappa,etta   ! radiation transport coefficients
    real(r8p)      :: tau,ctau     ! optical thickness
    integer(i4p)   :: ipt,idir,isegm
    integer(i4p)   :: inode,icell
    integer(i4p)   :: NPoints,NDirs
    integer(i4p)   :: iEn,NEnGroups
    real(r8p)      :: X1,X2,dX
    real(8)        :: L0
    type(TRayData), pointer              :: pRay => null()          ! ray data
    ! OMP parameters:
    integer     :: threads_count
    integer     :: thread_id
    integer     :: start_idx,stop_idx,dcount
!$  integer     :: omp_get_thread_num
    !integer     :: omp_get_num_threads
    ! time counters:
    real(4) :: deltaTime,iniTime

    iniTime = SECNDS(0.0E0);

    L0  = dimLen            ! cm

    write(*,*)
    write(*,*) ' - - - - - - - - '
    write(*,"(a)",advance='no') ' radiation transport ... '

    NPoints   = pTrace%NPoints
    NEnGroups = SpectrData%EnGrid%NX-1
    NDirs     = pTrace%pAngGrid%NDirs

    threads_count = 1
!$  threads_count = NumThreads
!!$  threads_count = omp_get_num_threads()
    dcount        = NPoints/threads_count
!    write(*,"(a,i2)",advance='no') ' Nth = ',threads_count
    write(*,"(a,i2)") ' NumThreads = ',threads_count

!$OMP   PARALLEL DEFAULT(PRIVATE)                                 &
!$OMP&  SHARED(dcount,threads_count),                             &
!$OMP&  SHARED(NPoints,NEnGroups,NDirs,L0),                       &
!$OMP&  SHARED(cUden,Wflux),                                      &
!$OMP&  SHARED(pTrace,SpectrData,CellKppa,CellEtta)               &
!$OMP&  NUM_THREADS(threads_count)
    thread_id = 0
!$  thread_id = omp_get_thread_num()
    start_idx = thread_id * dcount + 1
    stop_idx  = (thread_id + 1) * dcount
    if(thread_id + 1 .ge. threads_count) then
        stop_idx = NPoints
    endif
    do ipt = start_idx,stop_idx           ! cycle by mesh points, <good for parallel>
!$OMP MASTER
        write(*,"(a,a,I3,a)",advance="no") achar(13),' done ',ipt*100/dcount," %"
!$OMP END MASTER
        pt0 = pTrace%point(ipt)%pt0
        inode = pTrace%point(ipt)%bindId

        do iEn = 1,NEnGroups     ! cycle by spectral groups, <good for parallel>
            I_Om = 0.d0
            X1 = SpectrData%EnGrid%X(iEn)       ! current spectral group left bound
            X2 = SpectrData%EnGrid%X(iEn+1)     ! current spectral group right bound
            dX = X2 - X1                        ! current spectral group width

            do idir = 1,pTrace%pAngGrid%NDirs      ! cycle by rays (directions), <good for parallel>
                pRay => pTrace%point(ipt)%rays(idir)
                I_bnd = 0.e0_r8p
                I_ptX = 0.e0_r8p
                tau   = 0.e0_r8p
                ctau  = 0.e0_r8p
                kappa = 1.e0_r8p
                etta  = 1.e0_r8p
                do isegm = 1,pRay%NSegm
                    icell = pRay%icell(isegm)
                    kappa = CellKppa(icell,iEn)
                    etta  = CellEtta(icell,iEn)
                    ctau = kappa*pRay%segLen(isegm)*L0
                    ! calculate dI_ptX - addon for intencity
                    dI_ptX = 0.d0
                    if(kappa.gt.0.d0) then
                        dI_ptX = etta/kappa
                    end if
                    if(ctau.le.80.d0) then
                        dI_ptX = dI_ptX*(1.e0_r8p-exp(-ctau))
                    end if
                    if(tau.le.80.d0) then
                        dI_ptX = dI_ptX*exp(-tau)
                    else
                        dI_ptX = 0.d0
                    end if
                    I_ptX  = I_ptX + dI_ptX
                    tau  = tau + ctau
                end do
                ! calculate dI_ptX - addon for intencity (at bound)
                dI_ptX = I_bnd
                if(ctau.le.80.d0) then
                    dI_ptX = dI_ptX*(1.e0_r8p-exp(-ctau))
                end if
                if(tau.le.80.d0) then
                    dI_ptX = dI_ptX*exp(-tau)
                else
                    dI_ptX = 0.d0
                end if
                I_ptX  = I_ptX + dI_ptX

                theta = pTrace%pAngGrid%theta(idir)
                phi   = pTrace%pAngGrid%phi(idir)
                cosTh  = cos(theta)
                cosPhi = cos(phi)
                sinTh  = sin(theta)
                sinPhi = sin(phi)

                I_Om(idir) = I_ptX
            enddo

            ! intencity integration
            cU_ptX   = 0.d0
            W_ptX(:) = 0.d0
            dOm = pTrace%pAngGrid%dOmega
            do idir = 1,pTrace%pAngGrid%NDirs       ! cycle by rays (directions), <good for parallel>

                theta = pTrace%pAngGrid%theta(idir)
                phi   = pTrace%pAngGrid%phi(idir)
                cosTh  = cos(theta)
                cosPhi = cos(phi)
                sinTh  = sin(theta)
                sinPhi = sin(phi)

                cU_ptX   = cU_ptX   + I_Om(idir)*dOm
                W_ptX(1) = W_ptX(1) + I_Om(idir)*(-1)*sinTh*cosPhi * dOm
                W_ptX(2) = W_ptX(2) + I_Om(idir)*(-1)*sinTh*sinPhi * dOm
                W_ptX(3) = W_ptX(3) + 0.d0
            end do
            ! global mesh integrals
            cUden%fval(inode)         = cUden%fval(inode)        + cU_ptX*dX
            Wflux%fvec(inode)%cmp(:)  = Wflux%fvec(inode)%cmp(:) + W_ptX(:)*dX
!            if(present(cUsp)) then
!                cUsp%fvec(inode)%cmp(iEn) = cU_ptX
!            endif

        enddo
    enddo
!$OMP   END PARALLEL
    deltaTime = SECNDS(iniTime);
    write(*,"(a,1pE10.3e2,a)") " -> ",deltaTime,' s.'
    write(*,*)

    return
end subroutine rad_transport_long_chars

! ----------------------------------------------------------------------
subroutine rad_transport_long_chars_compact(pTrace,CellKppa,CellEtta,cUden,Wflux)
    use module_Settings
    use module_Constants
    use module_Geom
    use module_RayTrace
    use module_GridFunct
    use module_AtomData
    use module_Spectrum
    use module_utils
    use module_TxtRead
    implicit none

    ! input values
    type(TTraceDataCompact) :: pTrace
    real(8), dimension(pTrace%pmesh%NCells,SpectrData%EnGrid%NX-1) :: CellKppa
    real(8), dimension(pTrace%pmesh%NCells,SpectrData%EnGrid%NX-1) :: CellEtta
    ! output variables
    type(TGridNodalScalar) :: cUden
    type(TGridNodalVector) :: Wflux
    !type(TGridNodalVector), optional :: cUsp
    ! local variables
    real(r8p)               :: I_ptX    ! radiation intencity      at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: dI_ptX   ! add radiation intencity  at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: cU_ptX   ! radiation energy density at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p),dimension(3)  :: W_ptX    ! radiation energy flux    at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p),dimension(pTrace%pAngGrid%NDirs)  :: I_Om ! radiation intencity at (r,:,e)  (r-coords,Om-direction, e-photon energy)

    ! - - - - - - - - - - - - - - - -
    ! local variables
    type(T3DPoint) :: pt0
    real(r8p)      :: cosTh,cosPhi
    real(r8p)      :: sinTh,sinPhi
    real(r8p)      :: theta,phi
    real(r8p)      :: dOm
    real(r8p)      :: I_bnd         ! radiation intencity boundary value
    real(r8p)      :: kappa,etta   ! radiation transport coefficients
    real(r8p)      :: tau,ctau     ! optical thickness
    integer(i4p)   :: ipt,idir,isegm
    integer(i4p)   :: inode,icell
    integer(i4p)   :: NPoints,NDirs
    integer(i4p)   :: iEn,NEnGroups
    real(r8p)      :: X1,X2,dX
    real(8)        :: L0
    type(TRayDataCompact), pointer              :: pRay => null()          ! ray data
    ! OMP parameters:
    integer     :: threads_count
    integer     :: thread_id
    integer     :: start_idx,stop_idx,dcount
!$  integer     :: omp_get_thread_num
    ! time counters:
    real(4) :: deltaTime,iniTime

    iniTime = SECNDS(0.0E0);

    L0  = dimLen            ! cm

    !write(*,*)
    !write(*,*) ' - - - - - - - - '
    !write(*,"(a)",advance='no') ' radiation transport ... '

    NPoints   = pTrace%NPoints
    NEnGroups = SpectrData%EnGrid%NX-1
    NDirs     = pTrace%pAngGrid%NDirs

    threads_count = 1
!$  threads_count = NumThreads
!!$  threads_count = omp_get_num_threads()
    dcount        = NPoints/threads_count
!    write(*,"(a,i2)",advance='no') ' Nth = ',threads_count
    !write(*,"(a,i2)") ' NumThreads = ',threads_count

    write(*,"(a,I3,a)",advance="no") '-',0,"%"

!$OMP   PARALLEL DEFAULT(PRIVATE)                                 &
!$OMP&  SHARED(dcount,threads_count),                             &
!$OMP&  SHARED(NPoints,NEnGroups,NDirs,L0),                       &
!$OMP&  SHARED(cUden,Wflux),                                      &
!$OMP&  SHARED(pTrace,SpectrData,CellKppa,CellEtta)               &
!$OMP&  NUM_THREADS(threads_count)
    thread_id = 0
!$  thread_id = omp_get_thread_num()
    start_idx = thread_id * dcount + 1
    stop_idx  = (thread_id + 1) * dcount
    if(thread_id + 1 .ge. threads_count) then
        stop_idx = NPoints
    endif
    do ipt = start_idx,stop_idx           ! cycle by mesh points, <good for parallel>
!$OMP MASTER
        !write(*,"(a,a,I3,a)",advance="no") achar(13),' done ',ipt*100/dcount," %"
        write(*,"(a,a,I3,a)",advance="no") repeat(char(8),5),'-',ipt*100/dcount,"%"
!$OMP END MASTER
        pt0 = pTrace%point(ipt)%pt0
        inode = pTrace%point(ipt)%bindId

        do iEn = 1,NEnGroups     ! cycle by spectral groups, <good for parallel>
            I_Om = 0.d0
            X1 = SpectrData%EnGrid%X(iEn)       ! current spectral group left bound
            X2 = SpectrData%EnGrid%X(iEn+1)     ! current spectral group right bound
            dX = X2 - X1                        ! current spectral group width

            do idir = 1,pTrace%pAngGrid%NDirs      ! cycle by rays (directions), <good for parallel>
                pRay => pTrace%point(ipt)%rays(idir)
                I_bnd = 0.e0_r8p
                I_ptX = 0.e0_r8p
                tau   = 0.e0_r8p
                ctau  = 0.e0_r8p
                kappa = 1.e0_r8p
                etta  = 1.e0_r8p
                do isegm = 1,pRay%NSegm
                    icell = pRay%icell(isegm)
                    kappa = CellKppa(icell,iEn)
                    etta  = CellEtta(icell,iEn)
                    ctau = kappa*pRay%segLen(isegm)*L0
                    ! calculate dI_ptX - addon for intencity
                    dI_ptX = 0.d0
                    if(kappa.gt.0.d0) then
                        dI_ptX = etta/kappa
                    end if
                    if(ctau.le.80.d0) then
                        dI_ptX = dI_ptX*(1.e0_r8p-exp(-ctau))
                    end if
                    if(tau.le.80.d0) then
                        dI_ptX = dI_ptX*exp(-tau)
                    else
                        dI_ptX = 0.d0
                    end if
                    I_ptX  = I_ptX + dI_ptX
                    tau  = tau + ctau
                end do
                ! calculate dI_ptX - addon for intencity (at bound)
                dI_ptX = I_bnd
                if(ctau.le.80.d0) then
                    dI_ptX = dI_ptX*(1.e0_r8p-exp(-ctau))
                end if
                if(tau.le.80.d0) then
                    dI_ptX = dI_ptX*exp(-tau)
                else
                    dI_ptX = 0.d0
                end if
                I_ptX  = I_ptX + dI_ptX

                theta = pTrace%pAngGrid%theta(idir)
                phi   = pTrace%pAngGrid%phi(idir)
                cosTh  = cos(theta)
                cosPhi = cos(phi)
                sinTh  = sin(theta)
                sinPhi = sin(phi)

                I_Om(idir) = I_ptX
            enddo

            ! intencity integration
            cU_ptX   = 0.d0
            W_ptX(:) = 0.d0
            dOm = pTrace%pAngGrid%dOmega
            do idir = 1,pTrace%pAngGrid%NDirs       ! cycle by rays (directions), <good for parallel>

                theta = pTrace%pAngGrid%theta(idir)
                phi   = pTrace%pAngGrid%phi(idir)
                cosTh  = cos(theta)
                cosPhi = cos(phi)
                sinTh  = sin(theta)
                sinPhi = sin(phi)

                cU_ptX   = cU_ptX   + I_Om(idir)*dOm
                W_ptX(1) = W_ptX(1) + I_Om(idir)*(-1)*sinTh*cosPhi * dOm
                W_ptX(2) = W_ptX(2) + I_Om(idir)*(-1)*sinTh*sinPhi * dOm
                W_ptX(3) = W_ptX(3) + 0.d0
            end do
            ! global mesh integrals
            cUden%fval(inode)         = cUden%fval(inode)        + cU_ptX*dX
            Wflux%fvec(inode)%cmp(:)  = Wflux%fvec(inode)%cmp(:) + W_ptX(:)*dX
!            if(present(cUsp)) then
!                cUsp%fvec(inode)%cmp(iEn) = cU_ptX
!            endif

        enddo
    enddo
!$OMP   END PARALLEL
    write(*,"(a,a)",advance="no") repeat(char(8),5),'     '
    write(*,"(a)",advance="no")   repeat(char(8),5)
    deltaTime = SECNDS(iniTime);
    !write(*,"(a,1pE10.3e2,a)") " -> ",deltaTime,' s.'
    !write(*,*)

    return
end subroutine rad_transport_long_chars_compact

! ----------------------------------------------------------------------
subroutine rad_transport_1D(Trace1DData,pMesh,CellKppa,CellEtta,cUdenSp,WfluxSp)
    use module_Settings
    use module_constants
    use module_geom
    use module_RayTrace
    use module_GridFunct
    use module_AtomData
    use module_Spectrum
    use module_utils
    use module_TxtRead
    use module_Mesh
    implicit none

    ! input values
    !type(TTraceData), pointer       :: pTrace
    !type(TGridNodalScalar), pointer :: pDens,pTemp
    type(TTrace1DData)             :: Trace1DData
    type(TMeshData)                :: pMesh
    ! output variables
    type(TGridCellVector)   :: cUdenSp
    type(TGridCellVector)   :: WfluxSp
    ! local variables
    real(r8p)               :: I_ptX    ! radiation intencity      at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: dI_ptX   ! add radiation intencity  at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: cU_ptX   ! radiation energy density at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: W_ptX    ! radiation energy flux    at (r,Om,e) (r-coords,Om-direction, e-photon energy)
    real(r8p)               :: I_Om     ! radiation intencity at (r,:,e)  (r-coords,Om-direction, e-photon energy)

    ! - - - - - - - - - - - - - - - -
    ! local variables
    real(r8p)      :: I_bnd                      ! radiation intencity boundary value
    real(r8p)      :: kappa,etta                 ! radiation transport coefficients
    real(r8p)      :: ctau                       ! optical thickness
    real(r8p)      :: tau_01,tau_11,tau_02       ! optical thickness
    real(r8p)      :: a0,b0,c0,k0                ! layer parameters for quasi-1D solution
    real(r8p)      :: a1,b1,c1,k1                ! layer parameters for quasi-1D solution
    integer(i4p)   :: ipt,idir,idx
    integer(i4p)   :: icell,iface
    integer(i4p)   :: icell0,icell1
    integer(i4p)   :: imin,imax
    integer(i4p)   :: NPoints
    integer(i4p)   :: iEn,NXgr
    real(r8p)      :: X1,X2,dX
    real(r8p)      :: eps
    ! Saha variables
    real(8), dimension(pMesh%NCells,SpectrData%EnGrid%NX-1) :: CellKppa
    real(8), dimension(pMesh%NCells,SpectrData%EnGrid%NX-1) :: CellEtta
    real(8) :: L0
    ! OMP parameters:
    integer     :: threads_count
    integer     :: thread_id
    integer     :: start_idx,stop_idx,dcount
    !integer     :: omp_get_num_threads
!$  integer     :: omp_get_thread_num

!    write(*,*) "ray%NSegm = ", ray%NSegm

    L0   = dimLen            ! cm
    NXgr = SpectrData%EnGrid%NX-1


!    write(*,*)
!    write(*,*) ' - - - - - - - - '
!    write(*,*) ' radiation transport ... '

    eps = 1.e10_r8p
    NPoints = Trace1DData%NPoints
    NXgr    = SpectrData%EnGrid%NX-1


    threads_count = 1
!$  threads_count = NumThreads
!!$  threads_count = omp_get_num_threads()
    dcount    = NPoints/threads_count
!$OMP   PARALLEL DEFAULT(PRIVATE)                &
!$OMP&  SHARED(dcount,threads_count),            &
!$OMP&  SHARED(NPoints,NXgr,dimLen),             &
!$OMP&  SHARED(cUdenSp,WfluxSp),                 &
!$OMP&  SHARED(SpectrData,Trace1DData)           &
!$OMP&  SHARED(CellKppa,CellEtta)                &
!$OMP&  NUM_THREADS(threads_count)
    thread_id = 0
!$  thread_id = omp_get_thread_num()
! $OMP   BARRIER
    !threads_count = 1
    !thread_id = 0
    start_idx = thread_id * dcount + 1
    stop_idx  = (thread_id + 1) * dcount
    if(thread_id + 1 .ge. threads_count) then
        stop_idx = NPoints
    endif
    !write(*,*) "idx: ",start_idx,stop_idx

    do ipt = start_idx,stop_idx           ! cycle by mesh points, <good for parallel>
!    do ipt = 1,NPoints           ! cycle by mesh points, <good for parallel>
!        write(*,"(a,a,I3,a)",advance="no") achar(13),'done ',ipt*100/NPoints," %"
        !write(*,"(a,a,I3,a)",advance="no") achar(13),'done ',(ipt-start_idx)*100/(stop_idx-start_idx)," %"

        cUdenSp%fvec(ipt)%cmp(:) = 0.d0
        WfluxSp%fvec(ipt)%cmp(:) = 0.d0

        do iEn = 1,NXgr                         ! cycle by spectral groups, <good for parallel>
            I_Om = 0.d0
            X1 = SpectrData%EnGrid%X(iEn)       ! current spectral group left bound
            X2 = SpectrData%EnGrid%X(iEn+1)     ! current spectral group right bound
            dX = X2 - X1                        ! current spectral group width

            cU_ptX = 0.d0
            W_ptX  = 0.d0
            do iface = 1,Trace1DData%Nfaces     ! cycle by rays (directions), <good for parallel>

                icell0 = Trace1DData%point(ipt)%icell0(iface)
                icell1 = Trace1DData%point(ipt)%icell1(iface)
                if(icell0.gt.0) then
                    imin = min(icell0,icell1)
                    imax = max(icell0,icell1)
                    kappa = CellKppa(icell1,iEn)
                    etta  = CellEtta(icell1,iEn)
                    tau_11 = CellKppa(icell1,iEn)*dimLen
                    tau_01 = 0.d0
                    do icell = imin,imax
                        if(icell.ne.icell1) then
                            tau_01 = tau_01 + CellKppa(icell,iEn)*dimLen
                        endif
                    enddo
                    tau_02 = tau_01 + tau_11

                    ! 0. set boundary value
                    I_bnd = 0.e0_r8p
!                    idx = cUdenSp%BindId(1)       ! left boundary face index
                    idx = 1                        ! left boundary face index
                    if(iface.eq.idx .and. kappa.gt.0) then
                        I_bnd = etta/kappa * 0.0d-6
                    end if
!                    idx = cUdenSp%BindId(NPoints) ! right boundary face index
                    idx = Trace1DData%Nfaces       ! right boundary face index
                    if(iface.eq.idx .and. kappa.gt.0) then
                        I_bnd = etta/kappa * 0.0d-6
                    end if

                    ! 1. calculate I_ptX - addon for cU
                    I_ptX = 0.e0_r8p
                    a0 = Trace1DData%point(ipt)%prm0(iface)%a
                    b0 = Trace1DData%point(ipt)%prm0(iface)%b
                    c0 = Trace1DData%point(ipt)%prm0(iface)%c
                    k0 = Trace1DData%point(ipt)%prm0(iface)%k

                    dI_ptX = 0.d0
                    if(kappa.gt.0.d0) then
                        dI_ptX = etta/kappa
                    end if
                    ctau = (k0*tau_02)/(a0*tau_02+1.e0_r8p) - &
                           (k0*tau_01)/(a0*tau_01+1.e0_r8p) + &
                            c0*(tau_02-tau_01)
                    ctau = -ctau
                    if(ctau.le.80.d0) then
                        eps = exp(-ctau)
                        dI_ptX = I_bnd*eps + dI_ptX*(1.e0_r8p-eps)
                    end if

                    ctau = k0*tau_01/(a0*tau_01+1.e0_r8p)+c0*tau_01+b0
                    ctau = -ctau
                    if(ctau.le.80.d0) then
                        dI_ptX = dI_ptX*exp(-ctau)
                    else
                        dI_ptX = 0.d0
                    end if
                    I_ptX  = I_ptX + dI_ptX

                    cU_ptX = cU_ptX   + I_ptX


                    ! 2. calculate I_ptX - addon for Wz
                    I_ptX = 0.e0_r8p
                    a1 = Trace1DData%point(ipt)%prm1(iface)%a
                    b1 = Trace1DData%point(ipt)%prm1(iface)%b
                    c1 = Trace1DData%point(ipt)%prm1(iface)%c
                    k1 = Trace1DData%point(ipt)%prm1(iface)%k
                    idir = Trace1DData%point(ipt)%prm1(iface)%idir

                    dI_ptX = 0.d0
                    if(kappa.gt.0.d0) then
                        dI_ptX = etta/kappa
                    end if
                    ctau = (k1*tau_02)/(a1*tau_02+1.e0_r8p) - &
                           (k1*tau_01)/(a1*tau_01+1.e0_r8p) + &
                            c1*(tau_02-tau_01)
                    ctau = -ctau
                    if(ctau.le.80.d0) then
                        eps = exp(-ctau)
                        dI_ptX = I_bnd*eps + dI_ptX*(1.e0_r8p-eps)
                    end if

                    ctau = k1*tau_01/(a1*tau_01+1.e0_r8p)+c1*tau_01+b1
                    ctau = -ctau
                    if(ctau.le.80.d0) then
                        dI_ptX = dI_ptX*exp(-ctau)
                    else
                        dI_ptX = 0.d0
                    end if
                    I_ptX  = I_ptX + dI_ptX

                    W_ptX  = W_ptX    + idir*I_ptX
                end if

            enddo

           ! global mesh integrals
            cUdenSp%fvec(ipt)%cmp(iEn) = cU_ptX
            WfluxSp%fvec(ipt)%cmp(iEn) = W_ptX

        enddo
    enddo

!$OMP   END PARALLEL

!    write(*,*)

    return
end subroutine rad_transport_1D


! ----------------------------------------------------------------------
subroutine rad_transport_single_ray(pTrace,pDens,pTemp)
!subroutine rad_transport_single_ray(pTrace,pDens,pTemp,I_sp)
    use module_Settings
    use module_Constants
    use module_Geom
    use module_RayTrace
    use module_GridFunct
    use module_AtomData
    use module_Spectrum
    implicit none

    ! input values
    !type(TTraceData), pointer       :: pTrace
    !type(TGridNodalScalar), pointer :: pDens,pTemp
    type(TTraceData)       :: pTrace
    type(TGridNodalScalar) :: pDens,pTemp
    ! output variables
    real(8), dimension(SpectrData%EnGrid%NX-1) :: I_sp
    ! - - - - - - - - - - - - - - - -
    ! local variables
    type(TRayData), pointer :: pRay => null()          ! ray data
    type(T3DPoint)          :: pt0
    real(r8p)      :: cosTh,cosPhi
    real(r8p)      :: sinTh,sinPhi
    real(r8p)      :: theta,phi
    real(r8p)      :: kappa,etta         ! radiation transport coefficients
    real(r8p)      :: dtau,dI            ! radiation transport coefficients
    integer(i4p)   :: ipt,idir,isegm
    integer(i4p)   :: iEn,NEnGroups
    real(r8p)      :: X1,X2,dX
    integer(i4p)   :: icell,nnodes
    real(r8p)      :: eps
    real(r8p)      :: dens,temp,N_all
    real(r8p), dimension(8) :: fvals
    integer,   dimension(8) :: nodes
    ! Saha variables
    real(8), dimension(AtomData%NQds)   :: SahaNi
    real(8), dimension(AtomData%NElems) :: Ntot

    real(8), dimension(SpectrData%EnGrid%NX-1) :: Ibnd_sp
    real(8), dimension(SpectrData%EnGrid%NX-1) :: Kp_sp
    real(8), dimension(SpectrData%EnGrid%NX-1) :: Et_sp
    real(8), dimension(SpectrData%EnGrid%NX-1) :: tau_sp

    real(8) :: TeV,Ne,IonDeg
    real(8) :: Nt0,TK0,L0
    character(64) :: fname,sbuf

!    write(*,*) "ray%NSegm = ", ray%NSegm

    Ibnd_sp(:) = 0.e0_r8p

    TK0 = dimTem             ! K
    Nt0 = dimDen             ! cm-3
    L0  = dimLen             ! cm

    write(*,*)
    write(*,*) ' - - - - - - - - '
    write(*,*) ' single ray radiation transport ... '

    eps = 1.e10_r8p
    Ntot(:)   = 0.e0_r8p
    ipt       = 1
    NEnGroups = SpectrData%EnGrid%NX-1

    do idir = 1,pTrace%pAngGrid%NDirs      ! cycle by rays (directions)
        pRay => pTrace%point(ipt)%rays(idir)
        Ibnd_sp = 0.e0_r8p
        tau_sp(:) = 0.e0_r8p      ! cpectral tau for current ray
        I_sp(:)   = 0.e0_r8p      ! cpectral tau for current ray

        do isegm = 1,pRay%NSegm
            icell = pRay%icell(isegm)

            nnodes = pTrace%pmesh%Cell(icell)%nNodes
            nodes(1:nnodes) = pTrace%pmesh%Cell(icell)%inode(1:nnodes)
            ! average density in cell calculation
            fvals = pDens%fval(nodes(1:nnodes))
            dens = sum(fvals(1:nnodes))/nnodes
            Ntot(1) = Nt0*dens
            ! average temperature in cell calculation
            fvals = pTemp%fval(nodes(1:nnodes))
            temp = sum(fvals(1:nnodes))/nnodes
            TeV     = TK0*temp * k_Bol/eV_erg
            ! Saha populations calculation
            call SahaKin_TeV_Nt(TeV,Ntot,SahaNi,Ne,IonDeg)
            N_all     = sum(Ntot)
            call kappa_etta(TeV,Ne,SahaNi,Kp_sp,Et_sp)

            ! spectral groups cycle: intencity calculation
            do iEn = 1,NEnGroups     ! cycle by spectral groups
                X1 = SpectrData%EnGrid%X(iEn)       ! current spectral group left bound
                X2 = SpectrData%EnGrid%X(iEn+1)     ! current spectral group right bound
                dX = X2 - X1                        ! current spectral group width

                kappa = Kp_sp(iEn)
                etta  = Et_sp(iEn)
                dtau  = kappa*pRay%segLen(isegm)*L0
                ! calculate dInt_X - addon for intencity
                dI = 0.d0
                if(kappa.gt.0.d0) then
                    dI = etta/kappa
                end if
                if(tau_sp(iEn).le.80.d0) then
                    dI = dI*(1.e0_r8p-exp(-tau_sp(iEn)))
                end if
                if(dtau.le.80.d0) then
                    dI = dI*exp(-dtau)
                else
                    dI = 0.d0
                end if
                tau_sp(iEn) = tau_sp(iEn) + dtau
                I_sp(iEn)   = I_sp(iEn)   + dI
            end do
        enddo

        isegm = pRay%NSegm
        do iEn = 1,NEnGroups     ! cycle by spectral groups
            X1 = SpectrData%EnGrid%X(iEn)       ! current spectral group left bound
            X2 = SpectrData%EnGrid%X(iEn+1)     ! current spectral group right bound
            dX = X2 - X1                        ! current spectral group width

            kappa = Kp_sp(iEn)
            etta  = Et_sp(iEn)
            dtau  = kappa*pRay%segLen(isegm)*L0

            ! calculate dI_ptX - addon for intencity (at bound)
            dI = Ibnd_sp(iEn)
            if(tau_sp(iEn).le.80.d0) then
                dI = dI*(1.e0_r8p-exp(-tau_sp(iEn)))
            end if
            if(dtau.le.80.d0) then
                dI = dI*exp(-dtau)
            else
                dI = 0.d0
            end if
            I_sp(iEn)  = I_sp(iEn) + dI
        enddo

        theta = pTrace%pAngGrid%theta(idir)
        phi   = pTrace%pAngGrid%phi(idir)
        cosTh  = cos(theta)
        cosPhi = cos(phi)
        sinTh  = sin(theta)
        sinPhi = sin(phi)

        ! output
        sbuf = 'output\ray_Isp_'
        write(fname,"(a,i2.2,a)") trim(sbuf),idir,'.dat'
        open(18,file=fname)
        write(18,*) "# pt: = ",pt0%crd(:)
        write(18,*) "# th  = ",acos(cosTh)
        write(18,*) "# phi = ",acos(cosPhi)
        do iEn = 1,NEnGroups
            X1 = SpectrData%EnGrid%X(iEn)       ! current spectral group left bound
            X2 = SpectrData%EnGrid%X(iEn+1)     ! current spectral group right bound
            dX = X2 - X1                        ! current spectral group width

            write(18,"(7(1pE12.4e2))") 0.5*(X1+X2),I_sp(iEn)
            write(18,"(7(1pE12.4e2))")
            write(18,"(7(1pE12.4e2))")
        end do
        close(18)
    enddo


    return
end subroutine rad_transport_single_ray
