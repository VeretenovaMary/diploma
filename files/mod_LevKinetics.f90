module module_LevKinetics
    use module_Settings
    use module_AtomData

    type TKinRatesBB
        real(r8p) :: c_ex
        real(r8p) :: c_dex
        real(r8p) :: r_abs
        real(r8p) :: r_em
    end type TKinRatesBB

    type TKinRatesBF
        real(r8p) :: c_ion
        real(r8p) :: c_rec
        real(r8p) :: r_ion
        real(r8p) :: r_rec
    end type TKinRatesBF

    real(r8p), dimension(10) :: KinROpts
    integer, dimension(10)   :: KinIOpts


contains

    ! ==============================================================
    subroutine init_mod_LevKinetics(cfgName)
        use module_Settings
        use module_TxtRead
        implicit none
        ! input:
        character(*) :: cfgName
        ! local:
        integer :: ios
        integer :: NQtbb
        integer :: NQtbf
        integer :: cfgUnit,GetUnit

        cfgUnit = GetUnit()

        NQtbb = AtomData%NQt_bb
        NQtbf = AtomData%NQt_bf

        open(cfgUnit,file=cfgName)
        call txt_read_struct_value(cfgUnit,"LevKinOpts","RadFlag",KinROpts(1),ios)
        call txt_read_struct_value(cfgUnit,"LevKinOpts","CbbFlag",KinROpts(2),ios)
        call txt_read_struct_value(cfgUnit,"LevKinOpts","CbfFlag",KinROpts(3),ios)
        call txt_read_struct_value(cfgUnit,"LevKinOpts","RbbFlag",KinROpts(4),ios)
        call txt_read_struct_value(cfgUnit,"LevKinOpts","RbfFlag",KinROpts(5),ios)
        call txt_read_struct_value(cfgUnit,"LevKinOpts","RelTol",KinROpts(6),ios)
        close(cfgUnit)

!        allocate(KinRatesBB(NQtbb))
!        allocate(KinRatesBF(NQtbf))

        return

    end subroutine


    ! ==============================================================
    subroutine close_mod_LevKinetics()
        implicit none

!        deallocate(KinRatesBB)
!        deallocate(KinRatesBF)

    end subroutine

end module module_LevKinetics



! ==============================================================================
subroutine LevKinSolution(NH,Ntz,Nqz,Nez,cUsp,TKz,dt)
    use module_Constants
    use module_AtomData
    use module_Spectrum
    use module_LevKinetics
    use module_GridFunct
    implicit none
    ! - - - input/output
    integer                :: NH
    type(TGridNodalVector) :: Ntz
    type(TGridNodalVector) :: Nqz
    type(TGridNodalScalar) :: Nez
    type(TGridNodalVector) :: cUsp          ! spectral function at 1d mesh (selected points)
    type(TGridNodalScalar) :: TKz
    real(8)                :: dt

    ! - - - local
    real(8), dimension(AtomData%NElems)        :: Ntot
    real(8)                                    :: Ne,TK,TeV
    real(8), dimension(SpectrData%EnGrid%NX-1) :: cUx
    integer                                    :: NY
    real(8), dimension(AtomData%NQds)          :: Y0,Y1
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: KinRatesBB
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: KinRatesBF
    real(8)                 :: t1,t2
    integer                 :: i
    integer                 :: nsteps
    real(8)                 :: eps12
    real(8)                 :: calc_Ne
    ! OMP parameters:
    integer     :: threads_count
    integer     :: thread_id
    integer     :: start_idx,stop_idx,dcount
    !integer     :: omp_get_num_threads
!$  integer     :: omp_get_thread_num
    ! ---------------------------------------


    t1 = 0.d0
    t2 = dt
    NY = AtomData%NQds

    threads_count = 1
!$  threads_count = NumThreads
!!$  threads_count = omp_get_num_threads()

!$OMP   PARALLEL DEFAULT(PRIVATE), SHARED(t1,t2,cUsp,dcount,NH,NY,Ntz,Nqz,Nez,TKz,threads_count) NUM_THREADS(threads_count)
    thread_id = 0
!$  thread_id = omp_get_thread_num()
!!$OMP   BARRIER
    dcount    = NH/threads_count
    start_idx = thread_id * dcount + 1
    stop_idx  = (thread_id + 1) * dcount
    if(thread_id + 1 .ge. threads_count) then
        stop_idx = NH
    endif
    if(thread_id .eq. 0) then
        start_idx = 2                 ! 1 - is bound
    endif
!    write(*,*)
!    write(*,*) "thread_id = ",thread_id,start_idx,stop_idx

    do i = start_idx,stop_idx           ! cycle by mesh points, <good for parallel>
        cUx(:)  = cUsp%fvec(i)%cmp(:)
!        cUx(:)  = 0.d0

        Ntot(:) = Ntz%fvec(i)%cmp(:)
        TK      = TKz%fval(i)
        TeV     = TK * k_Bol/eV_erg
        Y0(:)   = Nqz%fvec(i)%cmp(:)
        Y1(:) = Y0(:)
        !call SahaKin_TeV_Nt(TeV,Ntot,Y1,Ne,eps12)
        Ne = calc_Ne(Y0)

        !write(*,"(a,i3.3)",advance="no") "-",i
        !call set_LevKinRates_lim(TK,KinRatesBB,KinRatesBF)
        call set_LevKinRates(cUx,TK,Ne,KinRatesBB,KinRatesBF)
        !call gack_solver (NY,t1,Y1,t2,Ne,KinRatesBB,KinRatesBF,eps12,nsteps)
        call gack2_solver(NY,t1,Y1,t2,Ne,KinRatesBB,KinRatesBF,eps12,nsteps)
        !write(*,"(4a)",advance="no") char(8),char(8),char(8),char(8)
        !write(*,"(a)",advance="no")  "    "
        !write(*,"(4a)",advance="no") char(8),char(8),char(8),char(8)

!        call dvode_solver(NY,t1,Y1,t2,IPAR,RPAR,eps12,nsteps)
        call positive_correction(Ntot,Y1,Ne)

        Nqz%fvec(i)%cmp(:) = Y1(:)
        Nez%fval(i)        = Ne
    enddo
!$OMP   END PARALLEL



    return
end subroutine LevKinSolution



! -------------------------------------------------------
subroutine gack_solver(neq,t0,Y0,t1,Ne,RCbb,RCbf,eps,nsteps)
    use module_Settings
    use module_LevKinetics
    implicit none
    integer :: neq
    real(8), dimension(neq) :: Y0
    real(8), dimension(neq) :: u00,u05,u10
    real(8), dimension(neq) :: ust
    real(8) :: t0,t1,tau
    real(8) :: tc
    real(8) :: Ysum0,Ymin,Ysum1
    integer :: nsteps,niter,iflag
    real(8) :: atol,rtol,eps            ! abs/rel tolerance, "zero" (minimal) value
    integer :: i
    ! -----------------------------------
    real(8)                                       :: Ne
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: RCbb
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: RCbf

    !    real(8), dimension(neq) :: vRHS
    !    integer :: k

    niter = 2
    atol = 1.d-4
    rtol = KinROpts(6)
    nsteps = 0
    Ysum0 = 0.d0
    Ysum0 = sum(Y0)

    ust = Y0/1.0d0
    u00 = Y0/1.0d0
    tc = t0
    eps  = 1.d0
    tau = t1 - tc
    do while (tc.lt.t1)
        eps  = 1.d0
        u00 = ust
        iflag = 0
        do while (eps.gt.rtol)
            call gackTauStep(neq,u00,u05,u10,tau,Ne,RCbb,RCbf,iflag,eps)
            if(eps.gt.rtol) then
                tau = tau*0.5d0
            endif
            iflag = 1
        enddo
        ust = u10

        Ysum1 = sum(ust)
        ust = (Ysum0/Ysum1)*ust

        tau = tau*1.1
        if(tc+tau.gt.t1) then
            tau = t1 - tc
            tc = t1
        else
            tc = tc + tau
        endif
        nsteps = nsteps + 1
    enddo

    ! conservativity equation (solution correction):
!    Ysum1 = sum(ust)
!    ust = (Ysum0/Ysum1)*ust
!
!    eps = abs(Ysum0-Ysum1)/Ysum0
!    if(eps.gt.0.5e0) then
!        write(*,*)
!        write(*,"(a,1pE12.4e3)") 'Ysum0 = ',Ysum0
!        write(*,"(a,1pE12.4e3)") 'Ysum1 = ',Ysum1
!        pause
!    end if


    Y0 = ust
    Ymin = 1.0d-20*Ysum0
    do i = 1,neq
        if(Y0(i).lt.Ymin) then
            Y0(i) = Ymin
        end if
    end do

    return
end subroutine gack_solver


! -------------------------------------------------------
subroutine gack2_solver(NY,t0,Y0,t1,Ne,RCbb,RCbf,eps,nsteps)
    use module_Settings
    use module_LevKinetics
    implicit none
    ! input:
    integer                                       :: NY
    real(8)                                       :: t0,t1
    real(8)                                       :: Ne
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: RCbb
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: RCbf
    ! input/output:
    real(8), dimension(NY) :: Y0
    ! output:
    real(8) :: eps                         ! abs/rel tolerance, "zero" (minimal) value
    integer :: nsteps
    ! local:
    real(8), dimension(0:NY) :: ust
    real(8), dimension(0:NY) :: u00,u05,u10
    real(8), dimension(NY)   :: Y1
    integer :: neq
    real(8) :: tc,dt,TStep
    real(8) :: tau,nju,mju
    real(8) :: dlen
    real(8) :: Ysum0,Ymin,Ysum1
    integer :: iflag
    real(8) :: atol,rtol                   ! abs/rel tolerance, "zero" (minimal) value
    integer :: i

    !    real(8), dimension(neq) :: vRHS
    !    integer :: k

    atol = 1.d-2
    rtol = KinROpts(6)
    nsteps = 0
    Ysum0 = sum(Y0)
    ! character values
    neq = NY + 1
    TStep = t1 - t0         ! time step

    ! character scales:
    tau   = t1-t0
    nju   = Ysum0
    mju   = tau/nju

    ust(0)    = 0.0d0
    ust(1:NY) = Y0(1:NY)
    tc  = ust(0)
    eps = 1.d0

    ! get dlen step estimation:
    call gack2GetDLen(NY,Y0,mju,Ne,RCbb,RCbf,dlen)
    dlen = 0.1*dlen

    tc = ust(0)
    do while (tc.lt.TStep)
        eps = 1.d0
        u00 = ust
        iflag = 0
        do while (eps.gt.rtol)
            call gack2LenStep(NY,u00,u05,u10,dlen,tau,mju,Ne,RCbb,RCbf,iflag,eps)
            if(eps.gt.rtol) then
                dlen = dlen*0.5d0
            endif
            !write(*,"(2(1x,1pE12.5e3))") eps,dlen
            iflag = 1
            nsteps = nsteps + 1
        enddo
        ust = u10

        Y1(1:NY) = ust(1:NY)
        Ysum1 = sum(Y1)
        Y1 = (Ysum0/Ysum1)*Y1
        ust(1:NY) = Y1(1:NY)

        !write(*,"(a,3(1x,1pE12.5e3))") "steps: ",TStep,ust(0)

!        write(*,"(a)")                 " . . . . . . . . . "
!        write(*,"(a,4(1x,1pE12.5e3))") " u00: ",u00(1:4)
!        write(*,"(a,4(1x,1pE12.5e3))") " ust: ",ust(1:4)
!        pause
        dt = ust(0)-tc      ! current time step
        tc = ust(0)
        !write(*,"(a,3(1x,1pE12.5e3))") "steps: ",TStep,tc,dt
        !dlen = dlen*2.0
        !dlen = 1.1*dlen
        if(tc+dt.gt.t1) then
            dlen = dlen*(1.e-8 + (t1-tc)/dt)
        else
            dlen = dlen*1.1
        endif
    enddo
    !pause

    Y0   = Y1
    Ymin = 1.0d-20*Ysum0
    do i = 1,NY
        if(Y0(i).lt.Ymin) then
            Y0(i) = Ymin
        end if
    end do

    return
end subroutine gack2_solver


! -------------------------------------------------------
subroutine gackTauStep(neq,u00,u05,u10,tau,Ne,RCbb,RCbf,iflag,eps)
    use module_AtomData
    use module_LevKinetics
    implicit none
    integer :: neq
    real(8), dimension(neq) :: u00,u05,u10
    real(8), dimension(neq) :: us
    real(8), dimension(neq) :: vFi,vPsi,vSig
    real(8) :: tau,dt
    real(8) :: eps,eps0
    integer :: niter,iflag
    integer :: j
    ! -----------------------------------
    real(8)                                       :: Ne
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: RCbb
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: RCbf

    niter = 2

    ! 1. - - -      1*tau step
    dt = tau
    us  = u00
    if(iflag.eq.0) then         ! because u05 is not defined (iterations first enter)
        do j = 1,niter
            us = 0.5d0*(u00 + us)
            call gackFCN (neq,dt,us,Ne,RCbb,RCbf,vFi,vPsi)
            vSig = 1.d0+0.5d0*dt*vFi
            us = (u00 + dt*vPsi*vSig)/(1.d0+dt*vFi*vSig)
        enddo
        u10 = us
    else                        ! u05 is already defined
        u10 = u05
    endif

    ! 2. - - -      2*(tau/2) step
    dt = 0.5d0*tau
    us  = u00
    do j = 1,niter
        us = 0.5d0*(u00 + us)
        call gackFCN (neq,dt,us,Ne,RCbb,RCbf,vFi,vPsi)
        vSig = 1.d0+0.5d0*dt*vFi
        us = (u00 + dt*vPsi*vSig)/(1.d0+dt*vFi*vSig)
    enddo
    u05 = us

    us  = u05
    do j = 1,niter
        us = 0.5d0*(u05 + us)
        call gackFCN (neq,dt,us,Ne,RCbb,RCbf,vFi,vPsi)
        vSig = 1.d0+0.5d0*dt*vFi
        us = (u05 + dt*vPsi*vSig)/(1.d0+dt*vFi*vSig)

    enddo

    ! calculate epsilon
    eps = 0.d0
    do j = 1,neq
        eps0 = 0.d0
        if((u10(j)+us(j)).gt.1.d-5) then
            eps0 = abs(u10(j) - us(j))/(u10(j) + us(j))
        endif
        eps  = max(eps,eps0)
    enddo
    u10 = us


    return
end subroutine gackTauStep

! -------------------------------------------------------
subroutine gack2LenStep(NY,u00,u05,u10,dlen,tau,mju,Ne,RCbb,Rcbf,iflag,eps)
    use module_AtomData
    use module_LevKinetics
    implicit none
    ! input:
    integer :: NY
    ! input/output:
    real(8), dimension(0:NY) :: u00,u05,u10
    ! input:
    real(8)                                       :: dlen
    real(8)                                       :: tau,mju
    real(8)                                       :: Ne
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: RCbb
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: RCbf
    integer                                       :: iflag
    real(8)                                       :: eps
    ! local:
    real(8), dimension(1:NY) :: vFi,vPsi
    real(8), dimension(1:NY) :: Y,FCN
    real(8), dimension(0:NY) :: vFi2,vPsi2,vSig2
    real(8), dimension(0:NY) :: us
    real(8) :: dl,t
    real(8) :: del0,delj
    integer :: niter
    integer :: j

    niter = 2
    t = 0.d0    ! fictive parameter

    ! 1. - - -      1*tau step
    dl = dlen
    us = u00
    if(iflag.eq.0) then         ! because u05 is not defined (iterations first enter)
        do j = 1,niter
            us = 0.5d0*(u00 + us)
            Y(1:NY) = us(1:NY)
            call gackFCN(NY,t,Y,Ne,RCbb,RCbf,vFi,vPsi)
            FCN = -Y*vFi + vPsi
            eps = sqrt(1.d0 + mju**2*sum(FCN**2))
            ! set functions:
            vFi2(0)     = 0.e0
            vPsi2(0)    = tau/eps
            vFi2(1:NY)  = tau*vFi(1:NY)/eps
            vPsi2(1:NY) = tau*vPsi(1:NY)/eps
            vSig2       = 1.d0+0.5d0*dl*vFi2

            us = (u00 + dl*vPsi2*vSig2)/(1.d0+dl*vFi2*vSig2)
        enddo
        u10 = us
    else                        ! u05 is already defined
        u10 = u05
    endif

    ! 2. - - -      2*tau/2 step
    dl = 0.5d0*dlen
    us  = u00
    do j = 1,niter
        us = 0.5d0*(u00 + us)
        Y(1:NY) = us(1:NY)
        call gackFCN(NY,t,Y,Ne,RCbb,RCbf,vFi,vPsi)
        FCN = -Y*vFi + vPsi
        eps = sqrt(1.d0 + mju**2*sum(FCN**2))
        ! set functions:
        vFi2(0)     = 0.e0
        vPsi2(0)    = tau/eps
        vFi2(1:NY)  = tau*vFi(1:NY)/eps
        vPsi2(1:NY) = tau*vPsi(1:NY)/eps
        vSig2       = 1.d0+0.5d0*dl*vFi2

        us = (u00 + dl*vPsi2*vSig2)/(1.d0+dl*vFi2*vSig2)
    enddo
    u05 = us

    us  = u05
    do j = 1,niter
        us = 0.5d0*(u05 + us)
        Y(1:NY) = us(1:NY)
        call gackFCN(NY,t,Y,Ne,RCbb,RCbf,vFi,vPsi)
        FCN = -Y*vFi + vPsi
        eps = sqrt(1.d0 + mju**2*sum(FCN**2))
        ! set functions:
        vFi2(0)     = 0.e0
        vPsi2(0)    = tau/eps
        vFi2(1:NY)  = tau*vFi(1:NY)/eps
        vPsi2(1:NY) = tau*vPsi(1:NY)/eps
        vSig2       = 1.d0+0.5d0*dl*vFi2
        vSig2 = 1.d0+0.5d0*dl*vFi2

        us = (u05 + dl*vPsi2*vSig2)/(1.d0+dl*vFi2*vSig2)
    enddo

    ! calculate tolerance:
    del0 = 0.d0
    if((u10(0)+us(0)).gt.1.d-5) then
        del0 = 2.0*abs(u10(0) - us(0))/(u10(0) + us(0))   ! ~1/3*(u'-u")**2/u"**2
        del0 = 1.0/3.0*del0**2
    endif

    delj = 0.d0
    Y(1:NY) = us(1:NY)
    call gackFCN(NY,t,Y,Ne,RCbb,RCbf,vFi,vPsi)
    FCN = -Y*vFi + vPsi
    do j = 1,NY
        eps = 0.d0
        if((u10(j)+us(j)).gt.1.d-5) then
            eps = 2.0*abs(u10(j) - us(j))/(u10(j) + us(j))
            eps = 1.0/3.0*eps**2
            eps = eps - FCN(j)*del0
            eps = abs(eps)
        endif
        delj = max(eps,delj)
    enddo
    u10 = us

    return
end subroutine gack2LenStep

! ------------------------------------------------
subroutine gack2GetDLen(NY,Y,mju,Ne,RCbb,RCbf,dlen)
    use module_AtomData
    use module_LevKinetics
    implicit none
    ! - - - input/output
    integer                                       :: NY
    real(8), dimension(NY)                        :: Y
    real(8)                                       :: mju
    real(8)                                       :: Ne
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: RCbb
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: RCbf
    ! output:
    real(8) :: dlen
    ! local:
    real(8), dimension(NY)  :: FCN
    real(8), dimension(NY)  :: vFi,vPsi
    real(8) :: t

    t = 0.d0
    call gackFCN(NY,t,Y,Ne,RCbb,RCbf,vFi,vPsi)
    FCN = -Y*vFi + vPsi
    dlen = 1.d0*sqrt(1.d0 + mju**2*sum(FCN**2))

    return
end subroutine gack2GetDLen

! ------------------------------------------------
function calc_Ne(Nq) result(Ne)
    use module_Constants
    use module_AtomData
    implicit none
    ! - - - input/output
    real(8), dimension(AtomData%NQds)    :: Nq
    ! output:
    real(8) :: Ne
    ! - - - local
    integer :: iE1,iI1,iQ1
    integer :: idx1
    real(8) :: Ymin

    Ymin = 1.0d-20*sum(Nq)
    Ne = 0.d0
    do iE1 = 1,AtomData%NElems
        do iI1 = 1,AtomData%elem(iE1)%NIons
            do iQ1 = 1,AtomData%elem(iE1)%ion(iI1)%Nqd
                idx1 = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%iList
                ! Y positive correction
                if(Nq(idx1).lt.Ymin) then
                    Nq(idx1) = Ymin
                end if
                Ne = Ne + Nq(idx1)*AtomData%elem(iE1)%ion(iI1)%ZQ
            enddo
        enddo
    enddo

    return
end function calc_Ne

! ------------------------------------------------
subroutine gackFCN(N,t,Y,Ne,RCbb,RCbf,vFi,vPsi)
    use module_Constants
    use module_AtomData
    use module_LevKinetics
    implicit none
    ! input:
    integer :: N              ! array dimension (number of configurations)
    real(8) :: t              ! time (for time explicit dependence - if exists)
    real(8), dimension(N)                         :: Y
    real(8)                                       :: Ne
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: RCbb
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: RCbf
    ! output:
    real(8), dimension(N)                         :: vFi,vPsi
    ! - - - local
    integer :: i
    integer :: idx1,idx2
    real(8) :: eps1,eps2
    real(8) :: calc_Ne

    ! system of level kinetics equations
    ! dNi/dt = {Km0}_ij*Nj + Ne*{Km1}_ij*Nj + Ne**2*{Km2}_ij*Nj
    ! where:
    !        Km0 = R_ex + R_dex + R_ion
    !        Km1 = C_ex + C_dex + C_ion + R_rec
    !        Km2 = C_rec


    vFi(:)   = 0.d0
    vPsi(:)  = 0.d0
    t        = t

    ! calculate Ne:
    Ne   = calc_Ne(Y)

    ! ionization/recombination part
    do i = 1,AtomData%NQt_bf
        idx1 = AtomData%QtbfList(i)%iLQd1
        idx2 = AtomData%QtbfList(i)%iLQd2

        ! - - - - - collisional processes
        eps1 = Ne    * RCbf(i)%c_ion
        eps2 = Ne**2 * RCbf(i)%c_rec

        vFi(idx1) = vFi(idx1) + eps1
        vFi(idx2) = vFi(idx2) + eps2
        vPsi(idx1) = vPsi(idx1) + Y(idx2)*eps2
        vPsi(idx2) = vPsi(idx2) + Y(idx1)*eps1

        ! - - - - - radiational processes
        eps1 =      RCbf(i)%r_ion
        eps2 = Ne * RCbf(i)%r_rec

        vFi(idx1) = vFi(idx1) + eps1
        vFi(idx2) = vFi(idx2) + eps2
        vPsi(idx1) = vPsi(idx1) + Y(idx2)*eps2
        vPsi(idx2) = vPsi(idx2) + Y(idx1)*eps1


    enddo

    ! excitation/deexcitation part
    do i = 1,AtomData%NQt_bb
        idx1 = AtomData%QtbbList(i)%iLQd1
        idx2 = AtomData%QtbbList(i)%iLQd2

        ! - - - - - collisional processes
        eps1 = Ne * RCbb(i)%c_ex
        eps2 = Ne * RCbb(i)%c_dex

        vFi(idx1) = vFi(idx1) + eps1
        vFi(idx2) = vFi(idx2) + eps2
        vPsi(idx1) = vPsi(idx1) + Y(idx2)*eps2
        vPsi(idx2) = vPsi(idx2) + Y(idx1)*eps1

        ! - - - - - radiational processes
        eps1 = RCbb(i)%r_abs
        eps2 = RCbb(i)%r_em

        vFi(idx1) = vFi(idx1) + eps1
        vFi(idx2) = vFi(idx2) + eps2

        vPsi(idx1) = vPsi(idx1) + Y(idx2)*eps2
        vPsi(idx2) = vPsi(idx2) + Y(idx1)*eps1

    enddo

    return
end subroutine gackFCN


! ------------------------------------------------
subroutine positive_correction(Ntot,Nq,Ne)
    use module_Constants
    use module_AtomData
    implicit none
    ! input/output
    real(8), dimension(AtomData%NElems)  :: Ntot
    real(8), dimension(AtomData%NQds)    :: Nq
    ! output:
    real(8) :: Ne
    ! local
    integer :: iE1,iI1,iQ1
    integer :: idx1
    real(8) :: Ymin,Nt,eps

    Ymin = 1.0d-20*sum(Nq)
    Ne   = 0.d0
    do iE1 = 1,AtomData%NElems
        do iI1 = 1,AtomData%elem(iE1)%NIons
            do iQ1 = 1,AtomData%elem(iE1)%ion(iI1)%Nqd
                idx1 = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%iList
                ! Y positive correction
                if(Nq(idx1).lt.Ymin) then
                    Nq(idx1) = Ymin
                end if
                Ne = Ne + Nq(idx1)*AtomData%elem(iE1)%ion(iI1)%ZQ
            enddo
        enddo
    enddo


    do iE1 = 1,AtomData%NElems
        Nt = 0.d0
        do iI1 = 1,AtomData%elem(iE1)%NIons
            do iQ1 = 1,AtomData%elem(iE1)%ion(iI1)%Nqd
                idx1 = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%iList
                ! Y positive correction
                if(Nq(idx1).lt.Ymin) then
                    Nq(idx1) = Ymin
                end if
                Nt = Nt + Nq(idx1)
            enddo
        enddo
        eps = abs(Ntot(iE1)-Nt)/Ntot(iE1)
        if(eps.gt.0.5e0) then
            write(*,*)
            write(*,*) 'Ntot(iE) = ',Ntot(iE1)
            write(*,*) 'Nt       = ',Nt
            stop
        end if

    enddo


    return
end subroutine positive_correction

!! ============================================================================================================
!subroutine dvode_solver(NEQ,TIN,Y,TOUT,IPAR,RPAR,eps12,nsteps)
!
!    USE dvode_f90m
!
!    IMPLICIT NONE
!    DOUBLE PRECISION ATOL, RTOL, TIN, TOUT, Y, RSTATS
!    INTEGER NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, I
!
!    INTEGER, DIMENSION(10) :: IPAR               ! user integer parameters
!    REAL(8), DIMENSION(10) :: RPAR               ! user real(8) parameters
!    INTEGER, DIMENSION(10) :: UIPAR               ! user integer parameters
!    REAL(8), DIMENSION(10) :: URPAR               ! user real(8) parameters
!    COMMON /InterfPRM/ UIPAR,URPAR
!
!    REAL(8)  :: eps12
!    INTEGER  :: nsteps
!
!    DIMENSION Y(NEQ), ATOL(NEQ), RSTATS(22), ISTATS(31)
!    EXTERNAL :: FCN_DVODE,JAC_DVODE
!
!
!    TYPE (VODE_OPTS) :: OPTIONS
!
!    IERROR = 0
!    ATOL=1.0d-5
!    RTOL=1.0d-5
!    ITASK = 1
!    ISTATE = 1
!!    OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,RELERR=RTOL, &
!!    USER_SUPPLIED_JACOBIAN=.TRUE.)
!
!    UIPAR = IPAR
!    URPAR = 0.d0
!    do i = 1,10
!        URPAR(i) = RPAR(i)
!    enddo
!
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL,      &
!              RELERR=RTOL,USER_SUPPLIED_JACOBIAN=.TRUE.)
!
!    call DVODE_F90(FCN_DVODE,NEQ,Y,TIN,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC_DVODE)
!
!    nsteps = 1
!
!    eps12 = 0.d0
!
!    do i = 1,10
!        RPAR(i) = URPAR(i)
!    enddo
!
!    return
!end subroutine dvode_solver
!
!! ======================================================================
!SUBROUTINE FCN_DVODE(NEQ, T, Y, YDOT)
!    IMPLICIT NONE
!    INTEGER  NEQ
!    REAL(8), DIMENSION(*) :: Y
!    REAL(8), DIMENSION(*) :: YDOT
!    REAL(8) :: T
!    INTEGER, DIMENSION(10) :: UIPAR               ! user integer parameters
!    REAL(8), DIMENSION(10) :: URPAR               ! user real(8) parameters
!    COMMON /InterfPRM/ UIPAR,URPAR
!    REAL(8), DIMENSION(10) :: RPAR               ! user real(8) parameters
!    INTEGER, DIMENSION(10) :: IPAR               ! user integer parameters
!    integer :: i
!
!    IPAR = UIPAR
!    do i = 1,10
!        RPAR(i) = URPAR(i)
!    enddo
!
!    CALL KinFCN(NEQ,T,Y,YDOT,IPAR,RPAR)
!
!    do i = 1,10
!        URPAR(i) = RPAR(i)
!    enddo
!
!    RETURN
!END SUBROUTINE FCN_DVODE
!
!! ============================================================================================================
!SUBROUTINE JAC_DVODE(NEQ, T, Y, ML, MU, PD, NROWPD)
!    IMPLICIT NONE
!    REAL(8) :: PD(NROWPD,NEQ)                       ! JACOBIAN (PARTIAL DIFFERENCES)
!    REAL(8) :: Y(NEQ)
!    REAL(8) :: T
!    INTEGER :: NEQ, ML, MU, NROWPD
!    INTEGER, DIMENSION(10) :: UIPAR               ! user integer parameters
!    REAL(8), DIMENSION(10) :: URPAR               ! user real(8) parameters
!    COMMON /InterfPRM/ UIPAR,URPAR
!    REAL(8), DIMENSION(10) :: RPAR               ! user real(8) parameters
!    INTEGER, DIMENSION(10) :: IPAR               ! user integer parameters
!    integer :: i
!
!    !    ML = ML
!    !    MU = MU
!    IPAR = UIPAR
!    do i = 1,10
!        RPAR(i) = URPAR(i)
!    enddo
!
!    NROWPD = NEQ
!    CALL KinJAC(NEQ,T,Y,PD,IPAR,RPAR)
!
!    do i = 1,10
!        URPAR(i) = RPAR(i)
!    enddo
!
!    RETURN
!END SUBROUTINE JAC_DVODE


! ==============================================================
subroutine set_LevKinRates(cUx,TK,Ne,KinRatesBB,KinRatesBF)
    use module_Constants
    use module_Settings
    use module_GridFunct
    use module_AtomData
    use module_Spectrum
    use module_LevKinetics
    !use module_Spectrum
    implicit none
    ! - - - input/output
    real(8), dimension(SpectrData%EnGrid%NX-1) :: cUx
    real(8)                                    :: TK
    real(8)                                    :: Ne
    ! - - - local
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: KinRatesBB
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: KinRatesBF
    ! . . . . . . . . . . . . . . . .
    integer :: idx1,idx2
    integer :: i,ilev
    real(8) :: TeV,dEeV,gQ1,gQ2,gf12
    real(8) :: C_ION,C_REC,C_EX,C_DEX
    real(8) :: R_ION_int,R_REC_int,R_EX_int,R_DEX_int
    real(8) :: fRad,fCbb,fCbf,fRbb,fRbf
    !    real(8) :: eint

    ! system of level kinetics equations
    ! dNi/dt = {Km0}_ij*Nj + Ne*{Km1}_ij*Nj + Ne**2*{Km2}_ij*Nj
    ! where:
    !        Km0 = R_ex + R_dex + R_ion
    !        Km1 = C_ex + C_dex + C_ion + R_rec
    !        Km2 = C_rec

    TeV = k_Bol*TK/eV_erg

    fRad = KinROpts(1)
    fCbb = KinROpts(2)
    fCbf = KinROpts(3)
    fRbb = KinROpts(4)
    fRbf = KinROpts(5)

    ! ionization/recombination part
    do i = 1,AtomData%NQt_bf
        idx1 = AtomData%QtbfList(i)%iLQd1
        idx2 = AtomData%QtbfList(i)%iLQd2

        dEeV = AtomData%QtbfList(i)%dE12
        gQ1  = AtomData%QtbfList(i)%gQd1
        gQ2  = AtomData%QtbfList(i)%gQd2
        ilev = AtomData%QtbfList(i)%ilev

        ! - - - - - collisional processes
        KinRatesBF(i)%c_ion = 0.e0
        KinRatesBF(i)%c_rec = 0.e0
        if(fCbf.gt.0) then
            KinRatesBF(i)%c_ion = fCbf * C_ION(TeV,dEeV)
            KinRatesBF(i)%c_rec = fCbf * C_REC(TeV,dEeV,gQ1,gQ2)
        endif

        ! - - - - - radiational processes
        KinRatesBF(i)%r_ion = 0.e0
        KinRatesBF(i)%r_rec = 0.e0
        if(fRbf.gt.0) then
            KinRatesBF(i)%r_ion = fRbf * R_ION_int(cUx,dEeV,ilev)
            KinRatesBF(i)%r_rec = fRbf * R_REC_int(cUx,TeV,dEeV,ilev,gQ1,gQ2)
        endif
    enddo

    ! excitation/deexcitation part
    do i = 1,AtomData%NQt_bb
        idx1 = AtomData%QtbbList(i)%iLQd1
        idx2 = AtomData%QtbbList(i)%iLQd2

        dEeV = AtomData%QtbbList(i)%dE12
        gf12 = AtomData%QtbbList(i)%gf12
        gQ1  = AtomData%QtbbList(i)%gQd1
        gQ2  = AtomData%QtbbList(i)%gQd2

        ! - - - - - collisional processes
        KinRatesBB(i)%c_ex  = 0.e0
        KinRatesBB(i)%c_dex = 0.e0
        if(fCbb.gt.0) then
            KinRatesBB(i)%c_ex  = fCbb * C_EX(TeV,dEeV,gf12)
            KinRatesBB(i)%c_dex = fCbb * C_DEX(TeV,dEeV,gf12,gQ1,gQ2)
        endif

        ! - - - - - radiational processes
        KinRatesBB(i)%r_abs = 0.e0
        KinRatesBB(i)%r_em  = 0.e0
        if(fRbb.gt.0) then
            KinRatesBB(i)%r_abs = fRbb * R_EX_int(cUx,Ne,TeV,i)
            KinRatesBB(i)%r_em  = fRbb * R_DEX_int(cUx,Ne,TeV,i)
        endif
    enddo


end subroutine


! ==============================================================
subroutine set_LevKinRates_lim(TK,KinRatesBB,KinRatesBF)
    use module_Constants
    use module_Settings
    use module_AtomData
    use module_LevKinetics
    !use module_Spectrum
    implicit none
    ! - - - input/output
    real(8)                                    :: TK
    ! - - - local
    type(TKinRatesBB), dimension(AtomData%NQt_bb) :: KinRatesBB
    type(TKinRatesBF), dimension(AtomData%NQt_bf) :: KinRatesBF
    ! . . . . . . . . . . . . . . . .
    integer :: idx1,idx2
    integer :: i,ilev
    real(8) :: TeV,dEeV,gQ1,gQ2,gf12
    real(8) :: C_ION,C_REC,C_EX,C_DEX
    real(8) :: R_ION_lim,R_REC_lim,R_EX_lim,R_DEX_lim
    real(8) :: fRad,fCbb,fCbf,fRbb,fRbf
    !    real(8) :: eint

    ! system of level kinetics equations
    ! dNi/dt = {Km0}_ij*Nj + Ne*{Km1}_ij*Nj + Ne**2*{Km2}_ij*Nj
    ! where:
    !        Km0 = R_ex + R_dex + R_ion
    !        Km1 = C_ex + C_dex + C_ion + R_rec
    !        Km2 = C_rec

    TeV = k_Bol*TK/eV_erg

    fRad = KinROpts(1)
    fCbb = KinROpts(2)
    fCbf = KinROpts(3)
    fRbb = KinROpts(4)
    fRbf = KinROpts(5)

    ! ionization/recombination part
    do i = 1,AtomData%NQt_bf
        idx1 = AtomData%QtbfList(i)%iLQd1
        idx2 = AtomData%QtbfList(i)%iLQd2

        dEeV = AtomData%QtbfList(i)%dE12
        gQ1  = AtomData%QtbfList(i)%gQd1
        gQ2  = AtomData%QtbfList(i)%gQd2
        ilev = AtomData%QtbfList(i)%ilev

        ! - - - - - collisional processes
        KinRatesBF(i)%c_ion = 0.e0
        KinRatesBF(i)%c_rec = 0.e0
        if(fCbf.gt.0) then
            KinRatesBF(i)%c_ion = fCbf * C_ION(TeV,dEeV)
            KinRatesBF(i)%c_rec = fCbf * C_REC(TeV,dEeV,gQ1,gQ2)
        endif

        ! - - - - - radiational processes
        KinRatesBF(i)%r_ion = 0.e0
        KinRatesBF(i)%r_rec = 0.e0
        if(fRbf.gt.0) then
            KinRatesBF(i)%r_ion = fRbf * R_ION_lim(fRad,TeV,dEeV,ilev)
            KinRatesBF(i)%r_rec = fRbf * R_REC_lim(fRad,TeV,dEeV,ilev,gQ1,gQ2)
        endif
    enddo

    ! excitation/deexcitation part
    do i = 1,AtomData%NQt_bb
        idx1 = AtomData%QtbbList(i)%iLQd1
        idx2 = AtomData%QtbbList(i)%iLQd2

        dEeV = AtomData%QtbbList(i)%dE12
        gf12 = AtomData%QtbbList(i)%gf12
        gQ1  = AtomData%QtbbList(i)%gQd1
        gQ2  = AtomData%QtbbList(i)%gQd2

        ! - - - - - collisional processes
        KinRatesBB(i)%c_ex  = 0.e0
        KinRatesBB(i)%c_dex = 0.e0
        if(fCbb.gt.0) then
            KinRatesBB(i)%c_ex  = fCbb * C_EX(TeV,dEeV,gf12)
            KinRatesBB(i)%c_dex = fCbb * C_DEX(TeV,dEeV,gf12,gQ1,gQ2)
        endif

        ! - - - - - radiational processes
        KinRatesBB(i)%r_abs = 0.e0
        KinRatesBB(i)%r_em  = 0.e0
        if(fRbb.gt.0) then
            KinRatesBB(i)%r_abs = fRbb * R_EX_lim(fRad,TeV,dEeV,gf12)
            KinRatesBB(i)%r_em  = fRbb * R_DEX_lim(fRad,TeV,dEeV,gf12,gQ1,gQ2)
        endif
    enddo


end subroutine


! ------------------------------------------------
subroutine KinFCN(N,t,Y,Fv,IPAR,RPAR)
    use module_atomdata
    use module_constants
    use module_LevKinetics
    implicit none
    ! - - - input/output
    integer :: N
    real(8) :: t
    real(8), dimension(N)  :: Y,Fv
    integer, dimension(10) :: IPAR
    real(8), dimension(10) :: RPAR

    ! - - - local
    !    integer :: i,j,k
    integer :: iE1,iI1,iQ1
    integer :: idx1,idx2
    integer :: i,ilev
    real(8) :: Ne
    real(8) :: TK,TeV,dEeV,gQ1,gQ2,gf12
    real(8) :: eps1,eps2
    real(8) :: C_ION,C_REC,C_EX,C_DEX
    real(8) :: R_ION_lim,R_REC_lim,R_EX_lim,R_DEX_lim
    real(8) :: fRad,fCbb,fCbf,fRbb,fRbf
    !    real(8) :: eint

    ! system of level kinetics equations
    ! dNi/dt = {Km0}_ij*Nj + Ne*{Km1}_ij*Nj + Ne**2*{Km2}_ij*Nj
    ! where:
    !        Km0 = R_ex + R_dex + R_ion
    !        Km1 = C_ex + C_dex + C_ion + R_rec
    !        Km2 = C_rec

    IPAR = 1*IPAR

    TK  = RPAR(AtomData%NElems+2)
    TeV = k_Bol*TK/eV_erg

    fRad = KinROpts(1)
    fCbb = KinROpts(2)
    fCbf = KinROpts(3)
    fRbb = KinROpts(4)
    fRbf = KinROpts(5)

    Fv  = 0.d0
    t   = t

    Ne = 0.d0
    do iE1 = 1,AtomData%NElems
        do iI1 = 1,AtomData%elem(iE1)%NIons
            do iQ1 = 1,AtomData%elem(iE1)%ion(iI1)%Nqd
                idx1 = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%iList
                Ne = Ne + Y(idx1)*AtomData%elem(iE1)%ion(iI1)%ZQ
            enddo
        enddo
    enddo

    RPAR(AtomData%NElems+1) = Ne

    ! ionization/recombination part
    do i = 1,AtomData%NQt_bf
        idx1 = AtomData%QtbfList(i)%iLQd1
        idx2 = AtomData%QtbfList(i)%iLQd2

        dEeV = AtomData%QtbfList(i)%dE12
        gQ1  = AtomData%QtbfList(i)%gQd1
        gQ2  = AtomData%QtbfList(i)%gQd2
        ilev = AtomData%QtbfList(i)%ilev

                ! - - - - - collisional processes
                eps1 = 0.d0
                eps2 = 0.d0
                if(fCbf.gt.0) then
                    eps1 = fCbf * Ne    * C_ION(TeV,dEeV)
                    eps2 = fCbf * Ne**2 * C_REC(TeV,dEeV,gQ1,gQ2)

                    Fv(idx1) = Fv(idx1) - Y(idx1)*eps1 + Y(idx2)*eps2
                    Fv(idx2) = Fv(idx2) + Y(idx1)*eps1 - Y(idx2)*eps2
                endif

                !                write(*,'(a,1pE12.4e3)') 'C_ION * Ne   = ',eps1
                !                write(*,'(a,1pE12.4e3)') 'C_REC * Ne^2 = ',eps2

                ! - - - - - radiational processes
                eps1 = 0.d0
                eps2 = 0.d0
                if(fRbf.gt.0) then
                    eps1 = fRbf *      R_ION_lim(fRad,TeV,dEeV,ilev)
                    eps2 = fRbf * Ne * R_REC_lim(fRad,TeV,dEeV,ilev,gQ1,gQ2)

                    Fv(idx1) = Fv(idx1) - Y(idx1)*eps1 + Y(idx2)*eps2
                    Fv(idx2) = Fv(idx2) + Y(idx1)*eps1 - Y(idx2)*eps2
                endif

                !                write(*,'(a,1pE12.4e3)') 'R_ION        = ',eps1
                !                write(*,'(a,1pE12.4e3)') 'R_REC * Ne   = ',eps2
                !                write(*,'(a,1pE12.4e3)') '    dEeV/TeV = ',dEeV/TeV
                !                write(*,'(a,1pE12.4e3)') '    exp(t)   = ',dexp(dEeV/TeV)
                !                write(*,'(a,1pE12.4e3)') '   eint(t)   = ',eint(dEeV/TeV)
                !                write(*,'(a,1pE12.4e3)') '        dEeV = ',dEeV
                !                write(*,'(a,1pE12.4e3)') '(Castor,BVJ)-> ',5.2e-14*(dEeV/TeV)**1.5*dexp(dEeV/TeV)*eint(dEeV/TeV)*Ne
                !                pause
                !
                !                write(*,'(a,i3)')            ' - - - Rbf = ',iRbf
                !                write(*,'(a,i3)')            'iRad = ',iRad
                !                write(*,'(a,10(1x,1pE10.3e2))') 'Ne,EQ1,EQ2,gQ1,gQ2:',Ne,EQ1,EQ2,dEeV,gQ1,gQ2
                !                write(*,'(a,10(1pE10.3e2))') 'Y1,Y2           :',Y(idx1),Y(idx2)
                !                write(*,'(a,10(1pE10.3e2))') 'e1,e2           :',eps1,eps2
                !                write(*,'(a,10(1pE10.3e2))') 'Y1*e1,Y2*e2     :',Y(idx1)*eps1,Y(idx2)*eps2
                !                write(*,'(a,10(1pE10.3e2))') 'Y1*e1-Y2*e2     :',Y(idx1)*eps1 - Y(idx2)*eps2
                !                pause

    enddo

    ! excitation/deexcitation part
    do i = 1,AtomData%NQt_bb
        idx1 = AtomData%QtbbList(i)%iLQd1
        idx2 = AtomData%QtbbList(i)%iLQd2

        dEeV = AtomData%QtbbList(i)%dE12
        gf12 = AtomData%QtbbList(i)%gf12
        gQ1  = AtomData%QtbbList(i)%gQd1
        gQ2  = AtomData%QtbbList(i)%gQd2

        ! - - - - - collisional processes
        eps1 = 0.d0
        eps2 = 0.d0
        if(fCbb.gt.0) then
            eps1 = fCbb * Ne * C_EX(TeV,dEeV,gf12)
            eps2 = fCbb * Ne * C_DEX(TeV,dEeV,gf12,gQ1,gQ2)

            Fv(idx1) = Fv(idx1) - Y(idx1)*eps1 + Y(idx2)*eps2
            Fv(idx2) = Fv(idx2) + Y(idx1)*eps1 - Y(idx2)*eps2
        endif

        !        write(*,'(a,1pE12.4e3)') 'C_EX  * Ne   = ',eps1
        !        write(*,'(a,1pE12.4e3)') 'C_DEX * Ne   = ',eps2

        ! - - - - - radiational processes
        eps1 = 0.d0
        eps2 = 0.d0
        if(fRbb.gt.0) then
            eps1 = fRbb * R_EX_lim(fRad,TeV,dEeV,gf12)
            eps2 = fRbb * R_DEX_lim(fRad,TeV,dEeV,gf12,gQ1,gQ2)

            Fv(idx1) = Fv(idx1) - Y(idx1)*eps1 + Y(idx2)*eps2
            Fv(idx2) = Fv(idx2) + Y(idx1)*eps1 - Y(idx2)*eps2
        endif

        !        write(*,'(a,1pE12.4e3)') 'R_EX         = ',eps1
        !        write(*,'(a,1pE12.4e3)') 'R_DEX        = ',eps2
        !        write(*,'(a,2(1pE12.4e3))') 'TeV,dEeV     = ',TeV,dEeV
        !        pause

        !        write(*,'(a,10(1pE10.3e2))') 'Ne,EQ1,EQ2,gQ1,gQ2:',Ne,EQ1,EQ2,gQ1,gQ2
        !        write(*,'(a,10(1pE10.3e2))') 'Y1,Y2           :',Y(idx1),Y(idx2)
        !        write(*,'(a,10(1pE10.3e2))') 'e1,e2           :',eps1,eps2
        !        write(*,'(a,10(1pE10.3e2))') 'Y1*e1,Y2*e2     :',Y(idx1)*eps1,Y(idx2)*eps2
        !        write(*,'(a,10(1pE10.3e2))') 'Y1*e1-Y2*e2     :',Y(idx1)*eps1 - Y(idx2)*eps2
        !        pause

    enddo

    return
end subroutine KinFCN

! ------------------------------------------------
subroutine KinJAC(N,t,Y,dFY,IPAR,RPAR)
    use module_atomdata
    use module_constants
    use module_LevKinetics
    implicit none
    ! - - - input/output
    integer :: N
    real(8) :: t
    real(8), dimension(N)   :: Y
    real(8), dimension(N,N) :: dFY
    integer, dimension(10)  :: IPAR
    real(8), dimension(10)  :: RPAR

    ! - - - local
    !    integer :: i,j,k
    integer :: iE1,iI1,iQ1
    integer :: iE2,iI2
    integer :: idx1,idx2
    integer :: i,ilev
    real(8) :: Ne,dNeY1,dNeY2
    real(8) :: TK,TeV,dEeV,gQ1,gQ2,gf12
    real(8) :: eps1,eps2
    real(8) :: C_ION,C_REC,C_EX,C_DEX
    real(8) :: R_ION_lim,R_REC_lim,R_EX_lim,R_DEX_lim
    real(8) :: fRad,fCbb,fCbf,fRbb,fRbf
    !    real(8) :: eint

    ! system of level kinetics equations
    ! dNi/dt = {Km0}_ij*Nj + Ne*{Km1}_ij*Nj + Ne**2*{Km2}_ij*Nj
    ! where:
    !        Km0 = R_ex + R_dex + R_ion
    !        Km1 = C_ex + C_dex + C_ion + R_rec
    !        Km2 = C_rec

    IPAR = 1*IPAR

    TK  = RPAR(AtomData%NElems+2)
    TeV = k_Bol*TK/eV_erg

    fRad = KinROpts(1)
    fCbb = KinROpts(2)
    fCbf = KinROpts(3)
    fRbb = KinROpts(4)
    fRbf = KinROpts(5)


    dFY = 0.d0
    t   = t

    Ne = 0.d0
    do iE1 = 1,AtomData%NElems
        do iI1 = 1,AtomData%elem(iE1)%NIons
            do iQ1 = 1,AtomData%elem(iE1)%ion(iI1)%Nqd
                idx1 = AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%iList
                Ne = Ne + Y(idx1)*AtomData%elem(iE1)%ion(iI1)%ZQ
            enddo
        enddo
    enddo

    RPAR(AtomData%NElems+1) = Ne

    ! ionization/recombination part
    do i = 1,AtomData%NQt_bf
        idx1 = AtomData%QtbfList(i)%iLQd1
        idx2 = AtomData%QtbfList(i)%iLQd2

        iE1  = AtomData%QdList(idx1)%iE
        iI1  = AtomData%QdList(idx1)%iI
        iE2  = AtomData%QdList(idx2)%iE
        iI2  = AtomData%QdList(idx2)%iI

        dEeV = AtomData%QtbfList(i)%dE12
        gQ1  = AtomData%QtbfList(i)%gQd1
        gQ2  = AtomData%QtbfList(i)%gQd2
        ilev = AtomData%QtbfList(i)%ilev

                dNeY1 = AtomData%elem(iE1)%ion(iI1)%ZQ
                dNeY2 = AtomData%elem(iE1)%ion(iI2)%ZQ

                ! - - - - - collisional processes
                eps1 = 0.d0
                eps2 = 0.d0
                if(fCbf.gt.0) then
                    eps1  = fCbf * C_ION(TeV,dEeV)
                    eps2  = fCbf * C_REC(TeV,dEeV,gQ1,gQ2)

                    dFY(idx1,idx1) = dFY(idx1,idx1) - eps1*(Ne - Y(idx1)*dNeY1)
                    dFY(idx1,idx2) = dFY(idx1,idx2) + eps2*(Ne**2 + Y(idx2)*2*Ne*dNeY2)

                    dFY(idx2,idx2) = dFY(idx2,idx2) - eps2*(Ne**2 + Y(idx2)*2*Ne*dNeY2)
                    dFY(idx2,idx1) = dFY(idx2,idx1) + eps1*(Ne - Y(idx1)*dNeY1)
                endif


                ! - - - - - radiational processes
                eps1 = 0.d0
                eps2 = 0.d0
                if(fRbf.gt.0) then
                    eps1 = fRbf * R_ION_lim(fRad,TeV,dEeV,ilev)
                    eps2 = fRbf * R_REC_lim(fRad,TeV,dEeV,ilev,gQ1,gQ2)

                    dFY(idx1,idx1) = dFY(idx1,idx1) - eps1
                    dFY(idx1,idx2) = dFY(idx1,idx2) + eps2*(Ne + Y(idx2)*dNeY2)

                    dFY(idx2,idx2) = dFY(idx2,idx2) - eps2*(Ne + Y(idx2)*dNeY2)
                    dFY(idx2,idx1) = dFY(idx2,idx1) + eps1
                endif
    enddo


    ! excitation/deexcitation part
    do i = 1,AtomData%NQt_bb
        idx1 = AtomData%QtbbList(i)%iLQd1
        idx2 = AtomData%QtbbList(i)%iLQd2
        iE1  = AtomData%QdList(idx1)%iE
        iI1  = AtomData%QdList(idx1)%iI
        iE2  = AtomData%QdList(idx2)%iE
        iI2  = AtomData%QdList(idx2)%iI

        dEeV = AtomData%QtbbList(i)%dE12
        gf12 = AtomData%QtbbList(i)%gf12
        gQ1  = AtomData%QtbbList(i)%gQd1
        gQ2  = AtomData%QtbbList(i)%gQd2

        dNeY1 = AtomData%elem(iE1)%ion(iI1)%ZQ
        dNeY2 = AtomData%elem(iE1)%ion(iI2)%ZQ

        ! - - - - - collisional processes
        eps1 = 0.d0
        eps2 = 0.d0
        if(fCbb.eq.1) then
            eps1 = fCbb * C_EX(TeV,dEeV,gf12)
            eps2 = fCbb * C_DEX(TeV,dEeV,gf12,gQ1,gQ2)

            dFY(idx1,idx1) = dFY(idx1,idx1) - eps1*(Ne - Y(idx1)*dNeY1)
            dFY(idx1,idx2) = dFY(idx1,idx2) + eps2*(Ne + Y(idx2)*dNeY2)

            dFY(idx2,idx2) = dFY(idx2,idx2) - eps2*(Ne + Y(idx2)*dNeY2)
            dFY(idx2,idx1) = dFY(idx2,idx1) + eps1*(Ne - Y(idx1)*dNeY1)
        endif


        ! - - - - - radiational processes
        eps1 = 0.d0
        eps2 = 0.d0
        if(fRbb.eq.1) then
            eps1 = fRbb * R_EX_lim(fRad,TeV,dEeV,gf12)
            eps2 = fRbb * R_DEX_lim(fRad,TeV,dEeV,gf12,gQ1,gQ2)

            dFY(idx1,idx1) = dFY(idx1,idx1) - eps1
            dFY(idx1,idx2) = dFY(idx1,idx2) + eps2

            dFY(idx2,idx2) = dFY(idx2,idx2) - eps2
            dFY(idx2,idx1) = dFY(idx2,idx1) + eps1
        endif


    enddo

    return
end subroutine KinJAC


! ======= collisional rates =================================
! ---------------- СКОРОСТЬ СТОЛКНОВИТЕЛЬНОГО ВОЗБУЖДЕНИЯ -----------------
REAL(8) FUNCTION C_EX(TeV,dEeV,gf)
      use module_Constants
      !      Вайнштейн, Собельман, Юков: Возбуждение атомов и уширение спектральных
      !      линий (стр. 118-119, 127)
      REAL(8) TeV,dEeV,gf

      C_EX = 0.0e00
      if(abs(dEeV/TeV).lt.80.e0) then
          C_EX = 4.36e-6*gf*dexp(-dEeV/TeV)/(dEeV*TeV**0.5)
          C_EX = C_EX*0.3e0*dlog(1.017 + 0.462*TeV/dEeV)      ! исправлен кулоновский лог.
      endif
      RETURN
END

! ---------------- СКОРОСТЬ СТОЛКНОВИТЕЛЬНОГО ДЕВОЗБУЖДЕНИЯ -----------------
REAL(8) FUNCTION C_DEX(TeV,dEeV,gf,J1,J2)
      use module_Constants
      !      Вайнштейн, Собельман, Юков: Возбуждение атомов и уширение спектральных
      !      линий (стр. 118-119, 127)
      REAL(8) TeV,dEeV,gf,J1,J2

      C_DEX = 4.36e-6*J1/J2*gf/(dEeV*TeV**0.5)
      C_DEX = C_DEX*0.3e0*dlog(1.017 + 0.462*TeV/dEeV)      ! исправлен кулоновский лог.
      RETURN
END

! -------------- СКОРОСТЬ СТОЛКНОВИТЕЛЬНОЙ ИОНИЗАЦИИ -----------------
REAL(8) FUNCTION C_ION(TeV,dEeV)
      use module_Constants
      !      Вайнштейн, Собельман, Юков: Возбуждение атомов и уширение спектральных
      !      линий (стр. 118-119, 127) формула Лотца:
      REAL(8) TeV,dEeV,eint

      C_ION = 0.0e00
      if (dEeV/TeV.lt.80.e0) then
            C_ION = 3.02e-6*eint(dEeV/TeV)/(dEeV*TeV**0.5)
      else
!            C_ION = 3.02e-6*dexp(-dEeV/TeV)*TeV**0.5/(dEeV**2)
            C_ION = 0.02e+00
      endif
      C_ION = 1.e0*C_ION
      RETURN
END

! -------------- СКОРОСТЬ СТОЛКНОВИТЕЛЬНОЙ РЕКОМБИНАЦИИ ---------------
REAL(8) FUNCTION C_REC(TeV,dEeV,J1,J2)
      use module_Constants
      !      Вайнштейн, Собельман, Юков: Возбуждение атомов и уширение спектральных
      !      линий (стр. 118-119, 127)
      REAL(8) TeV,dEeV,J1,J2,C0,eint
!      REAL(8) ExpEint
      !c1 = 0.5d0*(h_Pl**2/(2*pi*eV_erg*m_e))**1.5d0
      C0 = 1.65641518d-22
!      C_REC = 3.02e-6*C0*J1/J2*ExpEint(dEeV/TeV)/(dEeV*TeV**2)
      if (dEeV/TeV.lt.80.e0) then
            C_REC = 3.02e-6*C0*J1/J2*dexp(dEeV/TeV)*eint(dEeV/TeV)/(dEeV*TeV**2)
      else
            C_REC = 3.02e-6*C0*J1/J2/(dEeV**2*TeV)
      endif
      C_REC = 1.e0*C_REC
      RETURN
END


! ======= radiation rates limits: LTE, koronal ================
! ------------- photo emission (excitation) --------------
REAL(8) FUNCTION R_EX_lim(fRad,TeV,dEeV,gf12)
      use module_Constants
      implicit none
      REAL(8) TeV,dEeV,gf12
      REAL(8) :: A0,R0
      REAL(8) :: t,Ft
      REAL(8) :: fRad

      A0 = 6.33d+11
!      R0 = 6.71d-5*A0
!      R0 = 8*pi**2*q_e**2/(m_e*vc**3*h_Pl**2)*eV_erg**2
      R0 = 4.34515489d+07

      R_EX_lim = 0.d0
      t = dEeV/TeV
      Ft = 0.d0
      if(t.lt.10.d0) then
            Ft = 1.d0/(dexp(t)-1.d0)
      else
            Ft = dexp(-t)
      endif

      R_EX_lim = 0.d0
      if(fRad.eq.0) then
            R_EX_lim = 0.d0
      elseif(fRad.gt.0) then
            R_EX_lim = R0*gf12*dEeV**2
            R_EX_lim = fRad * R_EX_lim * Ft
      endif

      RETURN
END

! ------------- photo absorption (deexcitation) --------------
REAL(8) FUNCTION R_DEX_lim(fRad,TeV,dEeV,gf12,gQ1,gQ2)
      use module_Constants
      implicit none
      REAL(8) :: TeV,dEeV,gf12,gQ1,gQ2
      REAL(8) :: A0,R0
      REAL(8) :: t,Ft
      REAL(8) :: fRad

      A0 = 6.33d+11
!      R0 = 6.71d-5*A0
!      R0 = 8*pi**2*q_e**2/(m_e*vc**3*h_Pl**2)*eV_erg**2
      R0 = 4.34515489d+07
      R_DEX_lim = 0.d0
      t = dEeV/TeV
      Ft = 0.d0
      if(t.lt.10.d0) then
            Ft = dexp(t)/(dexp(t)-1.d0)
      else
            Ft = 1.d0
      endif

      R_DEX_lim = 0.d0
      if(fRad.eq.0) then
            R_DEX_lim = R0*gf12*dEeV**2*gQ1/gQ2
      elseif(fRad.gt.0) then
            R_DEX_lim = R0*gf12*dEeV**2*gQ1/gQ2
            R_DEX_lim = fRad * R_DEX_lim * Ft
      endif

      RETURN
END

! ------------- photo ionization --------------
REAL(8) FUNCTION R_ION_lim(fRad,TeV,dEeV,K)
      use module_Constants
      implicit none
      INTEGER :: K
      REAL(8) :: TeV,dEeV
!      REAL(8) :: A0,S0,a0_Bor
      REAL(8) :: R0
      REAL(8) :: t,Ft,rint
      REAL(8) :: fRad

!!      A0     = 2*h_Pl/vc**2
!      A0      = 1.47449944d-47
!!      a0_Bor = h_Pl/(2*pi*m_e*vc*alpha)
!      a0_Bor  = 5.2917721d-9                 ! cm
!!      S0     = 64*alpha*pi*a0_Bor**2/3**1.5e0   ! sigma0 - ioniz. crossection
!      S0      = 7.90706994d-18                   ! cm**2

!      R0 = 8*pi*S0*eV_erg**3/h_Pl**3/vc**2
      R0 = 3.12591480d+06

      R_ION_lim = 0.d0
      t = dEeV/TeV
      Ft = 0.d0
      if(t.lt.80.e0) then
            Ft = rint(t)
      else
            Ft = dexp(-t)/t
      endif

      R_ION_lim = 0.d0
      if(fRad.eq.0) then
            R_ION_lim = 0.d0
      elseif(fRad.gt.0) then
            R_ION_lim = fRad * K*R0*dEeV**3 * Ft
      endif

      RETURN
END

! ------------- photo recombination --------------
REAL(8) FUNCTION R_REC_lim(fRad,TeV,dEeV,K,gQ1,gQ2)
      use module_Constants
      implicit none
      INTEGER :: K
      REAL(8) :: TeV,dEeV,gQ1,gQ2
!      REAL(8) :: A0,S0,a0_Bor
      REAL(8) :: R0,C0,C1
      REAL(8) :: t,F0t,F1t,eint,rint
      REAL(8) :: fRad

!!      A0     = 2*h_Pl/vc**2
!      A0      = 1.47449944d-47
!!      a0_Bor = h_Pl/(2*pi*m_e*vc*alpha)
!      a0_Bor  = 5.2917721d-9                 ! cm
!!      S0     = 64*alpha*pi*a0_Bor**2/3**1.5e0   ! sigma0 - ioniz. crossection
!      S0      = 7.90706994d-18                   ! cm**2

!      C0 = 0.5d0*(h_Pl**2/(2*pi*eV_erg*m_e))**1.5d0
      C0 = 1.65641518d-22
!      R0 = 8*pi*S0*eV_erg**3/h_Pl**3/vc**2
      R0 = 3.12591480d+06
      C1 = C0*gQ1/gQ2/TeV**1.5e0

      t = dEeV/TeV
      F0t = 0.d0
      F1t = 0.d0
      if(t.lt.80.e0) then
            F0t = dexp(t)*eint(t)
            F1t = dexp(t)*rint(t)
      else
            F0t = 1.0d0/t
            F1t = 1.0d0/t
      endif

      R_REC_lim = 0.d0
      if(fRad.eq.0) then
            R_REC_lim = K*R0*C1*dEeV**3 * F0t
      elseif(fRad.gt.0) then
            R_REC_lim = fRad * K*R0*C1*dEeV**3 * F1t
      endif

      RETURN
END


! ------------- photo absorption (excitation) --------------
FUNCTION R_EX_int(cUx,Ne,TeV,iQtbb) result (R_EX)
      use module_Constants
      use module_AtomData
      use module_Spectrum
      implicit none
      ! input:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: cUx
      real(8)                                    :: Ne
      real(8)                                    :: TeV
      integer                                    :: iQtbb
      ! output:
      real(8)                                    :: R_EX
      ! local:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: prof
      integer :: i,imin,imax
      real(8) :: dX,Xi
      real(8) :: dEeV,gf12
      real(8) :: A0,R0

      R_EX = 0.e0
      A0 = 6.33d+11
!      R0 = 6.71d-5*A0
!      R0 = 8*pi**2*q_e**2/(m_e*vc**3*h_Pl**2)*eV_erg**2
      R0 = 4.34515489d+07
      dEeV = AtomData%QtbbList(iQtbb)%dE12
      gf12 = AtomData%QtbbList(iQtbb)%gf12

      call prof_Voigt_arr(Ne,TeV,iQtbb,prof)
      imin = SpectrData%QtBBOpts(iQtbb)%iGr_min
      imax = SpectrData%QtBBOpts(iQtbb)%iGr_max
      do i = imin,imax
        dX  =       (SpectrData%EnGrid%X(i+1)-SpectrData%EnGrid%X(i))
        Xi  = 0.5e0*(SpectrData%EnGrid%X(i+1)+SpectrData%EnGrid%X(i))
        R_EX = R_EX + cUx(i)*prof(i)*dX/Xi
      end do

      R_EX = 6.71e-05 * gf12 * R_EX
END

! ------------- photo emission (deexcitation) --------------
FUNCTION R_DEX_int(cUx,Ne,TeV,iQtbb) result (R_DEX)
      use module_Constants
      use module_AtomData
      use module_Spectrum
      implicit none
      ! input:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: cUx
      real(8)                                    :: Ne
      real(8)                                    :: TeV
      integer                                    :: iQtbb
      ! output:
      real(8)                                    :: R_DEX
      ! local:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: prof
      integer :: i,imin,imax
      real(8) :: dX,Xi
      real(8) :: dEeV,gf12,gQ1,gQ2
      real(8) :: A0,R0
      real(8) :: eps

      R_DEX = 0.e0
      A0 = 6.33d+11
!      R0 = 6.71d-5*A0
!      R0 = 8*pi**2*q_e**2/(m_e*vc**3*h_Pl**2)*eV_erg**2
      R0 = 4.34515489d+07
      dEeV = AtomData%QtbbList(iQtbb)%dE12
      gf12 = AtomData%QtbbList(iQtbb)%gf12
      gQ1  = AtomData%QtbbList(iQtbb)%gQd1
      gQ2  = AtomData%QtbbList(iQtbb)%gQd2

      call prof_Voigt_arr(Ne,TeV,iQtbb,prof)
      imin = SpectrData%QtBBOpts(iQtbb)%iGr_min
      imax = SpectrData%QtBBOpts(iQtbb)%iGr_max
      do i = imin,imax
        dX  =       (SpectrData%EnGrid%X(i+1)-SpectrData%EnGrid%X(i))
        Xi  = 0.5e0*(SpectrData%EnGrid%X(i+1)+SpectrData%EnGrid%X(i))
        eps = (dEeV-Xi)/TeV
        if(abs(eps).lt.80.d0) then
            R_DEX = R_DEX + (cUx(i)+6.59e+11*Xi**3)*prof(i)*exp((dEeV-Xi)/TeV)*dX/Xi
        endif
      end do

      R_DEX = 6.71e-05 * gf12*gQ1/gQ2 * R_DEX
END


! ------------- photo ionization --------------
FUNCTION R_ION_int(cUx,dEeV,ilev) result (R_ION)
      use module_Constants
      use module_AtomData
      use module_Spectrum
      implicit none
      ! input:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: cUx
      real(8)                                    :: dEeV
      integer                                    :: ilev
      ! output:
      real(8)                                    :: R_ION
      ! local:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: sigma
      integer :: i
      real(8) :: dX,Xi

      R_ION = 0.e0

      call sigma_pi_arr(ilev,dEeV,sigma)
      do i = 1,SpectrData%EnGrid%NX-1
        dX  =       (SpectrData%EnGrid%X(i+1)-SpectrData%EnGrid%X(i))
        Xi  = 0.5e0*(SpectrData%EnGrid%X(i+1)+SpectrData%EnGrid%X(i))
        R_ION = R_ION + cUx(i)*sigma(i)*dX/Xi
      end do

      R_ION = R_ION / eV_erg
END


! ------------- photo ionization --------------
FUNCTION R_REC_int(cUx,TeV,dEeV,ilev,gQ1,gQ2) result (R_REC)
      use module_Constants
      use module_AtomData
      use module_Spectrum
      implicit none
      ! input:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: cUx
      real(8)                                    :: TeV
      real(8)                                    :: dEeV
      integer                                    :: ilev
      real(8)                                    :: gQ1,gQ2
      ! output:
      real(8)                                    :: R_REC
      ! local:
      real(8), dimension(SpectrData%EnGrid%NX-1) :: sigma
      integer :: i
      real(8) :: dX,Xi
      real(8) :: eps

      R_REC = 0.e0

      call sigma_pi_arr(ilev,dEeV,sigma)
      do i = 1,SpectrData%EnGrid%NX-1
        dX  =       (SpectrData%EnGrid%X(i+1)-SpectrData%EnGrid%X(i))
        Xi  = 0.5e0*(SpectrData%EnGrid%X(i+1)+SpectrData%EnGrid%X(i))
        eps = (dEeV-Xi)/TeV
        if(abs(eps).lt.80.d0) then
            R_REC = R_REC + (cUx(i)+6.59e+11*Xi**3)*sigma(i)*exp((dEeV-Xi)/TeV)*dX/Xi
        endif
      end do

      R_REC = 1.66e-22/eV_erg/TeV**1.5e0 * gQ1/gQ2 * R_REC
END



! ------------------ integral exponent ------------------------
function eint(z)
      implicit doubleprecision (a-h,o-z)
      parameter (a1=8.5733287401,a2=18.0590169730,a3=8.6347608925)
      parameter (a4=0.2677737343,b1=9.5733223454,b2=25.6329561486)
      parameter (b3=21.0996530827,b4=3.9584969228)
      parameter (aa0=-0.57721566,aa1=0.99999193,aa2=-0.24991055)
      parameter (aa3=0.05519968,aa4=-0.00976004,aa5=0.00107857)

      eint=0.d0

      if(z .gt. 80.d0) then
            eint = exp(-z)/z
            return
      end if

      if(z .lt. 1.0d-10) then
            eint = aa0-log(z)
            return
      end if

      if(z.lt.1.0) then
            eint=aa0-log(z)+(((((((((aa5*z)+aa4)*z)+aa3)*z)+aa2)*z)+aa1)*z)
      else
            part1=((((((z+a1)*z)+a2)*z)+a3)*z)+a4
            part2=((((((z+b1)*z)+b2)*z)+b3)*z)+b4
            part3=z*exp(z)
            eint=part1/(part2*part3)
      endif
      return
end

! ------------------ integral function in rad. rates ------------------------
function rint(z)
      implicit none
      real(8) :: rint,z

      rint = dexp(-z)/z

      return
end

! ------------------ integral function in rad. rates ------------------------
function ExpEint(z)
      implicit doubleprecision (a-h,o-z)
      parameter (a1=8.5733287401,a2=18.0590169730,a3=8.6347608925)
      parameter (a4=0.2677737343,b1=9.5733223454,b2=25.6329561486)
      parameter (b3=21.0996530827,b4=3.9584969228)
      parameter (aa0=-0.57721566,aa1=0.99999193,aa2=-0.24991055)
      parameter (aa3=0.05519968,aa4=-0.00976004,aa5=0.00107857)
      ExpEint = 0.d0
      if(z.lt.1e-10) return
      if(z.lt.1.0) then
            eps=aa0-log(z)+(((((((((aa5*z)+aa4)*z)+aa3)*z)+aa2)*z)+aa1)*z)
            ExpEint = eps*dexp(z)
      else
            part1=((((((z+a1)*z)+a2)*z)+a3)*z)+a4
            part2=((((((z+b1)*z)+b2)*z)+b3)*z)+b4
            part3=z
            ExpEint = part1/(part2*part3)
      endif
      return
end



