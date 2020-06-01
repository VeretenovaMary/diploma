! ---------------------------------------------------
subroutine SahaKin_TeV_Nt(TeV,Ntot,SahaNi,Ne,IonDeg)
    use module_Settings
    use module_Constants
    use module_TxtRead
    use module_AtomData

    implicit none
    ! - - - input
    type(TConcData),dimension(AtomData%NElems) :: SahaData
    real(8), dimension(AtomData%NQds)        :: SahaNi

    real(8), dimension(AtomData%NElems)   :: Ntot
    real(8)                     :: TK,TeV
    real(8),parameter           :: rtol = 1.d-8
    real(8)                     :: Ne,IonDeg

    ! - - - local
    real(8)              :: EQ,EQd,dEQd
    real(8)              :: gQ,gQd
    real(8)              :: SWe,dEi
    integer              :: i,j,k,nx
    integer              :: NIons
    real(8)              :: N_all
    character(100)       :: fpath,strbuf

    real(8) :: eps,yi,xPk
    real(8) :: EdT
    real(8) :: x1,x2,dx
    real(8) :: xl,xr,xc
    real(8) :: fl,fr,fc,dfc
    real(8) :: f1,f2
    real(8) :: df1,df2
    real(8) :: ShSigm12,dShSigm12
    integer :: ires,istp,idx
    real(8) :: dCi,Sml,Smu
    integer,dimension(AtomData%NElems) :: MajIonIdx
!    logical :: dirExists                    ! KV

!    fpath = trim(OutputDir)//'levels/'
!    call slash_win2lin(fpath)
!    inquire( file=trim(fpath)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
!    if(.not.dirExists) call system('mkdir '//trim(fpath),ios)

    do i = 1,AtomData%NElems
        k = AtomData%elem(i)%NIons
        call SetSahaDataSize(SahaData(i),k)
    enddo

    !   calculate free electrons stat. weight
    TK = TeV*eV_erg/k_Bol
    call electron_stat_weight(TK,SWe)
    !   set Ci coefficients

    N_all = 0.d0
    do i = 1,AtomData%NElems
        N_all = N_all + Ntot(i)
    enddo


    do i = 1,AtomData%NElems
        NIons = AtomData%elem(i)%NIons
        !write(*,*) 'Saha: ',AtomData%elem(i)%title
        SahaData(i)%title = AtomData%elem(i)%title
        SahaData(i)%NIons = AtomData%elem(i)%NIons

        ! calculate average ions parameters: average energy, statsum
        do j=1,NIons
            eps = AtomData%elem(i)%ion(j)%EQ
            do k = 1,AtomData%elem(i)%ion(j)%Nqd               ! find lower energy level
                EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                eps = min(eps,EQd)                             ! -> lower energy value
            enddo
            gQ  = 0.d0
            EQ  = 0.d0
            do k = 1,AtomData%elem(i)%ion(j)%Nqd               !... find average energy, statsum
                EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                dEQd = EQd - eps
                gQd  = AtomData%elem(i)%ion(j)%Qd(k)%gQd
                EdT = dEQd/TeV
                if(EdT.lt.80.d0) then
                    EQ = EQ + EQd*gQd*dexp(-dEQd/TeV)       ! -> average energy value
                    gQ = gQ +     gQd*dexp(-dEQd/TeV)       ! -> statsum
                endif
            enddo
            EQ = EQ/gQ
            SahaData(i)%EQ(j) = AtomData%elem(i)%ion(j)%EQ    ! v.1 - average ion energy
            !SahaData(i)%EQ(j) = eps                          ! v.2 - average ion energy
            !SahaData(i)%EQ(j) = EQ                           ! v.3 - average ion energy
            SahaData(i)%gQ(j) = gQ
            !write(*,'(a,1pE12.4e3)') 'Saha: EQ = ',SahaData(i)%EQ(j)
        enddo

        ! calculate ions Ci-coefficients
        SahaData(i)%Ci(:) = 1.d0
        do j=1,NIons-1
            dEi = SahaData(i)%EQ(j+1)-SahaData(i)%EQ(j)
            !dEi = AtomData%elem(i)%ion(j)%dEion
            !dEi = AtomData%elem(i)%ion(j+1)%Qd(1)%EQd-SahaData(i)%EQ(j)
            !write(*,'(a,2(1x,1pE12.4e3))') 'Saha -> dEion: ',dEi,AtomData%elem(i)%ion(j)%dEion
            EdT = dEi/TeV
            SahaData(i)%Ci(j) = 0.d0
            if(EdT.lt.80.d0) then
                SahaData(i)%Ci(j) = SWe*SahaData(i)%gQ(j+1)/SahaData(i)%gQ(j)*dexp(-dEi/TeV)
            endif
            SahaData(i)%Ci(j) = SahaData(i)%Ci(j)/N_all
        enddo
    enddo

    do i=1,AtomData%NElems
        dCi = SahaData(i)%Ci(1)
        NIons = AtomData%elem(i)%NIons
        do j=1,NIons-1
            dCi = max(dCi,SahaData(i)%Ci(j))
        enddo
        if(dCi.eq.0.d0) then
            dCi = 1.d0
        end if
        SahaData(i)%dCi = dCi
        SahaData(i)%Ci = SahaData(i)%Ci/dCi
    enddo

    !    find interval [x1,x2] where solution apears
    NIons = AtomData%elem(1)%NIons
    do i = 1,AtomData%NElems
        NIons = max(NIons,AtomData%elem(i)%NIons)
    enddo

    x1 = 0
    f1 = 0.d0
    do i = 1,AtomData%NElems
        yi  = Ntot(i)/N_all
        f1 = f1 + yi*(AtomData%elem(i)%NIons-1)
    enddo
    f1   = f1 - x1
    x2   = 0
    do k = 2,NIons
        x2 = 1.d0*(k-1)
        f2 = 0.d0
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            yi  = Ntot(i)/N_all
            f2 = f2 + yi*ShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,x2/dCi)
        enddo
        f2 = f2 - x2
        if(f1*f2.le.0.d0) then
            exit
        endif
        x1 = x2
        f1 = f2
    enddo


    !write(*,*) 'find interval [x1,x2] where solution apears: end'
    !   - - - plot function near solution - - -
!    nx = 1001
    nx = 101
    dx = (x2-x1)/(nx-1)
    strbuf = trim(fpath)//'fcn_Saha.dat'
!    open(2,file=strbuf)
!    write(2,'(a,1x,a)') ' data:','function'
    do j=1,nx
        xc = x1 + dx*(j-1)
        !xc = 1.d-80*10.d0**(j/5.0)
        fc  = 0.d0
        dfc = 0.d0
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            yi  = Ntot(i)/N_all
            fc  = fc  + yi* ShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,xc/dCi)
            dfc = dfc + yi*dShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,xc/dCi)/dCi
        enddo
        fc  = fc  - xc
        dfc = dfc - 1.d0
!        write(2,'(3(1pE17.9e3,2x))') xc,fc,dfc
    enddo
!    close(2)

    !   - - - derivatives - - -
    f1  = 0.d0
    df1 = 0.d0
    do i=1,AtomData%NElems
        dCi = SahaData(i)%dCi
        yi  = Ntot(i)/N_all
        f1  = f1  + yi* ShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,x1/dCi)
        df1 = df1 + yi*dShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,x1/dCi)/dCi
    enddo
    f1  = f1  - x1
    df1 = df1 - 1.d0

    f2  = 0.d0
    df2 = 0.d0
    do i=1,AtomData%NElems
        dCi = SahaData(i)%dCi
        yi  = Ntot(i)/N_all
        f2  = f2  + yi* ShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,x2/dCi)
        df2 = df2 + yi*dShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,x2/dCi)/dCi
    enddo
    f2  = f2  - x2
    df2 = df2 - 1.d0

!    write(*,'(a,1pE12.4e3)')    'dCi : ',dCi
!    write(*,'(a,2(1pE12.4e3))') 'x_i : ',x1,x2
!    write(*,'(a,2(1pE12.4e3))') 'f_i : ',f1,f2
!    write(*,'(a,2(1pE12.4e3))') 'df_i: ',df1,df2

    ! start solution search
    if(abs(df1).gt.abs(df2)) then
        xc = x1              ! current x (start for search solution)
        istp = 0
    else
        xc = x2              ! current x (start for search solution)
        istp = 1
    endif
    eps = 1.e0
    fc  = 0.d0
    dfc = 0.d0
    do i=1,AtomData%NElems
        dCi = SahaData(i)%dCi
        yi  = Ntot(i)/N_all
        fc  = fc  + yi* ShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,xc/dCi)
        dfc = dfc + yi*dShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,xc/dCi)/dCi
    enddo
    fc  = fc  - xc
    dfc = dfc - 1.d0
    xl = x1
    fl = f1
    xr = x2
    fr = f2

    ires = 0
!    open(2,file=strbuf,access='append')
!    write(2,*)
!    write(2,*)
!    write(2,'(a,1x,a)') ' data:','solution_steps'
!    write(2,'(4(1pE17.9e3,2x),i3)') xc,fc,dfc,eps
    j = 0
    do while(eps.gt.rtol)
        j = j + 1
        xc = xc - fc/dfc

        if(xc.le.xl .or. xc.ge.xr) then
            ires = 2                               ! out of range [xl,xr] -> midle solution
            xc = 0.5d0*(xl+xr)
        endif

        eps = (xr-xl)/(xc-xl)*(1-istp) + (xr-xl)/(xr-xc)*istp
!        if(eps.gt.1.d2 .and. ires.eq.0) then
        if(eps.gt.1.d2) then
            ires = 3                               ! too huge derivative -> midle solution
            xc = 0.5d0*(xl+xr)
        endif

        fc  = 0.d0
        dfc = 0.d0
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            yi  = Ntot(i)/N_all
            fc  = fc  + yi* ShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,xc/dCi)
            dfc = dfc + yi*dShSigm12(AtomData%elem(i)%NIons,SahaData(i)%Ci,xc/dCi)/dCi
        enddo
        fc  = fc  - xc
        dfc = dfc - 1.d0

        if(fl*fc.le.0.d0) then
            xr = xc
            fr = fc
            istp = 1
        endif
        if(fr*fc.le.0.d0) then
            xl = xc
            fl = fc
            istp = 0
        endif
        ires = 1                                   ! tangent solution

        eps = 2.0*abs(fc)/(abs(f1)+abs(f2))
!        write(2,'(4(1pE17.9e3,2x),i3,4(1pE17.9e3,2x))') xc,fc,dfc,eps,ires,xl,xr
    enddo
!    close(2)

    ! calculate electrons number
    istp = 0
    Ne = xc*N_all                ! initial value (before correction inerations)
    do istp=1,10        ! number of iterations

        ! find major ions indexes (ions that have higher populations)
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            ires = 1
            if(xc.gt.0.d0) then
                do k=2,AtomData%elem(i)%NIons
                    if(SahaData(i)%Ci(k-1)*dCi/xc.ge.1.d0) then
                        ires = k
                    endif
                enddo
            endif
            SahaData(i)%MajIonIdx = ires
        enddo

        ! calculate major ions number
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            Sml = 0.d0
            Smu = 0.d0
            eps = 1.d0

            NIons = AtomData%elem(i)%NIons
            ires = SahaData(i)%MajIonIdx
            ! lower energy ions
            do k=1,ires-1
                xPk = 1.d0
                do j=1,k
                    if(SahaData(i)%Ci(ires-j).gt.0.d0) then
                        xPk = xPk*(xc/SahaData(i)%Ci(ires-j)/dCi)
                    else
                        xPk = 0.d0
                    endif
                enddo
                Sml = Sml + xPk
            enddo
            ! higher energy ions
            do k=1,NIons-ires
                xPk = 1.d0
                do j=0,k-1
                    if(xc.gt.0.d0) then
                        xPk = xPk*(SahaData(i)%Ci(ires+j)*dCi/xc)
                    else
                        xPk = 0.d0
                    endif
                enddo
                Smu = Smu + xPk
            enddo
            SahaData(i)%Ni(ires) = Ntot(i)/(1.d0+Sml+Smu)
        enddo

!        ! calculate other ions number
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            NIons = AtomData%elem(i)%NIons
            ires = SahaData(i)%MajIonIdx
            ! lower then MajIonIdx ions:
            do k=1,ires-1
                if(SahaData(i)%Ci(ires-k).gt.0.d0) then
                    SahaData(i)%Ni(ires-k) = SahaData(i)%Ni(ires-k+1)*xc/dCi/SahaData(i)%Ci(ires-k)
                else
                    SahaData(i)%Ni(ires-k) = 0.d0
                endif
            enddo
            ! higher then MajIonIdx ions:
            do k=1,NIons-ires
                if(xc.gt.0.d0) then
                    SahaData(i)%Ni(ires+k) = SahaData(i)%Ni(ires+k-1)*SahaData(i)%Ci(ires+k-1)*dCi/xc
                else
                    SahaData(i)%Ni(ires+k) = 0.d0
                endif
            enddo
        enddo

        ! recalculate Ne (value correction in case of Ne->0)
        Ne = 0.d0
        do i=1,AtomData%NElems
            do k=1,AtomData%elem(i)%NIons
                Ne = Ne + SahaData(i)%Ni(k)*(k-1)
            enddo
        enddo
        do i=1,AtomData%NElems
            SahaData(i)%Ne = Ne
        enddo
        xc = Ne/N_all
    enddo


    ! find other populations: redistribute values between different ion configurations
    do i = 1,AtomData%NElems
        do j = 1,AtomData%elem(i)%NIons

            ires = 1
            EQd = AtomData%elem(i)%ion(j)%Qd(1)%EQd
            gQd = AtomData%elem(i)%ion(j)%Qd(1)%gQd
            do k = 2,AtomData%elem(i)%ion(j)%Nqd
                eps = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                if(eps.lt.EQd) then
                    ires = k
                    EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                    gQd = AtomData%elem(i)%ion(j)%Qd(k)%gQd
                endif
            enddo

            eps = 0.d0
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                dEi = AtomData%elem(i)%ion(j)%Qd(k)%EQd - EQd
                EdT = dEi/TeV
                if(EdT.lt.80.d0) then
                    eps = eps + AtomData%elem(i)%ion(j)%Qd(k)%gQd/gQd*dexp(-dEi/TeV)
                endif
            enddo

            idx = AtomData%elem(i)%ion(j)%Qd(ires)%iList
            if(eps.gt.0.d0) then
                SahaNi(idx) = SahaData(i)%Ni(j)/eps
            else
                SahaNi(idx) = 0.d0
            endif

            eps = SahaNi(idx)
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                dEi = AtomData%elem(i)%ion(j)%Qd(k)%EQd - EQd
                EdT = dEi/TeV
                SahaNi(idx) = 0.d0
                if(EdT.lt.80.d0) then
                    SahaNi(idx) = eps * AtomData%elem(i)%ion(j)%Qd(k)%gQd/gQd * dexp(-dEi/TeV)
                endif
            enddo
        enddo
    enddo

    MajIonIdx = SahaData(:)%MajIonIdx


    IonDeg = 0.d0
    do i = 1,AtomData%NElems
        do j=2,AtomData%elem(i)%NIons
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                IonDeg = IonDeg + SahaNi(idx)
            enddo
        enddo
    enddo
    IonDeg = IonDeg/N_all

!    ! screen output
!    write(*,'(a)')              'Saha results: '
!    write(*,'(a,1pE12.4e3)')    '  Ne       =  ',Ne
!    write(*,'(a,1pE12.4e3)')    '  IonDeg   =  ',IonDeg
!    write(*,"(a)",advance='no') '  Maj.Ion id: '
!    do i=1,AtomData%NElems
!        write(*,'(i5,5x)',advance='no') MajIonIdx(i)
!    enddo
!    write(*,*)
!    write(*,*)
!    write(*,*)

    do i = 1,AtomData%NElems
        call UnsetSahaDataSize(SahaData(i))
    enddo

    return
end subroutine SahaKin_TeV_Nt


! ---------------------------------------------------
subroutine SahaKin_TeV_Ne(TeV,Ne,SahaNi,Ntot,IonDeg)
    use module_Settings
    use module_Constants
    use module_AtomData
    use module_TxtRead

    implicit none
    ! - - - input
    type(TConcData),dimension(AtomData%NElems) :: SahaData
    real(8), dimension(AtomData%NQds)        :: SahaNi

    real(8), dimension(AtomData%NElems)   :: Ntot
    real(8)                     :: TK,TeV
    real(8)                     :: Ne

    ! - - - local
    real(8)              :: EQ,EQd,dEQd
    real(8)              :: gQ,gQd
    real(8)              :: SWe,dEi
    integer              :: i,j,k
    integer              :: NIons
    character(100)       :: fpath

    real(8) :: eps
    integer :: ires,istp,idx
    real(8) :: dCi,N_all,IonDeg
    integer,dimension(AtomData%NElems) :: MajIonIdx
    !logical :: dirExists
    !integer :: ios

    fpath = trim(OutputDir)//'levels/'
    call slash_correct(fpath,OSType)
    !inquire( file=trim(fpath)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
    !if(.not.dirExists) call system('mkdir '//trim(fpath),ios)
!    call system('mkdir '//trim(fpath),ios)

    do i = 1,AtomData%NElems
        k = AtomData%elem(i)%NIons
        call SetSahaDataSize(SahaData(i),k)
    enddo


    !   calculate free electrons stat. weight
    TK  = TeV*eV_erg/k_Bol
    call electron_stat_weight(TK,SWe)

    !   set Ci coefficients
    do i = 1,AtomData%NElems
        NIons = AtomData%elem(i)%NIons
        !write(*,*) 'Saha: ',AtomData%elem(i)%title
        SahaData(i)%title = AtomData%elem(i)%title
        SahaData(i)%NIons = AtomData%elem(i)%NIons

        ! calculate average ions parameters: average energy, statsum
        do j=1,NIons
            eps = AtomData%elem(i)%ion(j)%EQ
            do k = 1,AtomData%elem(i)%ion(j)%Nqd               ! find lower energy level
                EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                eps = min(eps,EQd)                       ! -> lower energy value
            enddo
            gQ  = 0.d0
            EQ  = 0.d0
            do k = 1,AtomData%elem(i)%ion(j)%Nqd               !... find average energy, statsum
                EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                dEQd = EQd - eps
                gQd  = AtomData%elem(i)%ion(j)%Qd(k)%gQd
                EQ = EQ + EQd*gQd*dexp(-dEQd/TeV)       ! -> average energy value
                gQ = gQ +     gQd*dexp(-dEQd/TeV)       ! -> statsum
            enddo
            EQ = EQ/gQ
            SahaData(i)%EQ(j) = AtomData%elem(i)%ion(j)%EQ    ! v.1 - average ion energy
            !SahaData(i)%EQ(j) = eps                     ! v.2 - average ion energy
            !SahaData(i)%EQ(j) = EQ                      ! v.3 - average ion energy
            SahaData(i)%gQ(j) = gQ
            !write(*,'(a,1pE12.4e3)') 'Saha: EQ = ',SahaData(i)%EQ(j)
        enddo

        ! calculate ions Ci-coefficients
        SahaData(i)%Ci(:) = 1.d0
        do j=1,NIons-1
            dEi = SahaData(i)%EQ(j+1)-SahaData(i)%EQ(j)
            !dEi = SahaData(i)%EQ(j+1)-AtomData%elem(i)%ion(j)%Qd(1)%EQd
            SahaData(i)%Ci(j) = SWe*SahaData(i)%gQ(j+1)/SahaData(i)%gQ(j)*dexp(-dEi/TeV)
        enddo
    enddo

    do i=1,AtomData%NElems
        dCi = SahaData(i)%Ci(1)
        NIons = AtomData%elem(i)%NIons
        do j=1,NIons-1
            dCi = max(dCi,SahaData(i)%Ci(j))
        enddo
        SahaData(i)%dCi = dCi
        SahaData(i)%Ci = SahaData(i)%Ci/dCi
    enddo


    ! calculate ions populations
    do istp=1,1        ! number of iterations

        ! find major ions indexes (ions that have higher populations)
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            ires = 1
            do k=2,AtomData%elem(i)%NIons
                if(SahaData(i)%Ci(k-1)*dCi.ge.1.d0) then
                    ires = k
                endif
            enddo
            SahaData(i)%MajIonIdx = ires
        enddo

        ! first set major ions number -> Ne (then normalise)
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            ires = SahaData(i)%MajIonIdx
            SahaData(i)%Ni(ires) = Ne
        enddo

!        ! calculate other ions number
        do i=1,AtomData%NElems
            dCi = SahaData(i)%dCi
            NIons = AtomData%elem(i)%NIons
            ires = SahaData(i)%MajIonIdx
            ! lower then MajIonIdx ions:
            do k=1,ires-1
                if(SahaData(i)%Ci(ires-k).gt.0.d0) then
                    SahaData(i)%Ni(ires-k) = SahaData(i)%Ni(ires-k+1)*Ne/dCi/SahaData(i)%Ci(ires-k)
                else
                    SahaData(i)%Ni(ires-k) = 0.d0
                endif
            enddo
            ! higher then MajIonIdx ions:
            do k=1,NIons-ires
                if(Ne.gt.0.d0) then
                    SahaData(i)%Ni(ires+k) = SahaData(i)%Ni(ires+k-1)*SahaData(i)%Ci(ires+k-1)*dCi/Ne
                else
                    SahaData(i)%Ni(ires+k) = 0.d0
                endif
            enddo
        enddo

        ! recalculate Ne (eps) for Ni values renomalise
        eps = 0.d0
        do i=1,AtomData%NElems
            do k=1,AtomData%elem(i)%NIons
                eps = eps + SahaData(i)%Ni(k)*(k-1)
            enddo
        enddo

        if(eps.gt.0.d0) then
            do i=1,AtomData%NElems
                do k=1,AtomData%elem(i)%NIons
                    SahaData(i)%Ni(k) = SahaData(i)%Ni(k)*Ne/eps
                enddo
            enddo
        endif

        ! calculate Ntot
        Ntot = 0.d0
        do i=1,AtomData%NElems
            do k=1,AtomData%elem(i)%NIons
                Ntot(i) = Ntot(i) + SahaData(i)%Ni(k)
            enddo
        enddo
    enddo


    ! find other populations: redistribute values between different ion configurations
    do i = 1,AtomData%NElems
        do j = 1,AtomData%elem(i)%NIons

            ires = 1
            EQd = AtomData%elem(i)%ion(j)%Qd(1)%EQd
            gQd = AtomData%elem(i)%ion(j)%Qd(1)%gQd
            do k = 2,AtomData%elem(i)%ion(j)%Nqd
                eps = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                if(eps.lt.EQd) then
                    ires = k
                    EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                    gQd = AtomData%elem(i)%ion(j)%Qd(k)%gQd
                endif
            enddo

            eps = 0.d0
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                dEi = AtomData%elem(i)%ion(j)%Qd(k)%EQd - EQd
                eps = eps + AtomData%elem(i)%ion(j)%Qd(k)%gQd/gQd*dexp(-dEi/TeV)
            enddo

            idx = AtomData%elem(i)%ion(j)%Qd(ires)%iList
            if(eps.gt.0.d0) then
                SahaNi(idx) = SahaData(i)%Ni(j)/eps
            else
                SahaNi(idx) = 0.d0
            endif

            eps = SahaNi(idx)
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                dEi = AtomData%elem(i)%ion(j)%Qd(k)%EQd - EQd
                SahaNi(idx) = eps * AtomData%elem(i)%ion(j)%Qd(k)%gQd/gQd * dexp(-dEi/TeV)
            enddo
        enddo
    enddo

    MajIonIdx = SahaData(:)%MajIonIdx

    Ntot  = 0.d0
    N_all = 0.d0
    do i=1,AtomData%NElems
        do j=1,AtomData%elem(i)%NIons
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                Ntot(i) = Ntot(i) + SahaNi(idx)
            enddo
        enddo
        N_all = N_all + Ntot(i)
    enddo

    IonDeg = 0.d0
    do i = 1,AtomData%NElems
        do j=2,AtomData%elem(i)%NIons
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                IonDeg = IonDeg + SahaNi(idx)
            enddo
        enddo
    enddo
    IonDeg = IonDeg/N_all

!    ! screen output
!    write(*,'(a)')              'Saha results: '
!    write(*,'(a,1pE12.4e3)')    '  Ne       =  ',Ne
!    write(*,"(a)",advance='no') '  Maj.Ion id: '
!    do i=1,AtomData%NElems
!        write(*,'(i5,5x)',advance='no') MajIonIdx(i)
!    enddo
!    write(*,*)
!    write(*,*)
!    write(*,*)

    do i = 1,AtomData%NElems
        call UnsetSahaDataSize(SahaData(i))
    enddo

    return
end subroutine SahaKin_TeV_Ne




! **************************************************
subroutine SetSahaDataSize(A,n)
    use module_AtomData
    implicit none
    type(TConcData)          :: A
    integer          :: n
    allocate(A%EQ(n))
    allocate(A%gQ(n))
    allocate(A%Ni(n))
    allocate(A%Ci(n))
    return
end subroutine SetSahaDataSize

! **************************************************
subroutine UnsetSahaDataSize(A)
    use module_AtomData
    implicit none
    type(TConcData)          :: A

    deallocate(A%EQ)
    deallocate(A%gQ)
    deallocate(A%Ni)
    deallocate(A%Ci)
    return
end subroutine UnsetSahaDataSize

! ********************************************************************
subroutine SahaOutput(TeV,Ntot,SahaNi,Ne)
    use module_Settings
    use module_AtomData
    use module_TxtRead
    use module_Constants

    implicit none
    ! - - - input
    real(8), dimension(AtomData%NQds)        :: SahaNi
    real(8), dimension(AtomData%NElems)      :: Ntot
    real(8)                     :: TK,TeV
    real(8)                     :: Ne

    ! - - - local
    integer         :: i,j,k
    character(100)  :: OutputPath
    real(8)         :: EQ,EQd,dEi
    integer         :: NQsum,idx,num
    real(8), allocatable :: SahaQsum(:)
    logical         :: dirExists

    TK = TeV*eV_erg/k_Bol

    ! output all configurations
    OutputPath  = trim(OutputDir)//'/levels/'
    call slash_correct(OutputPath,OSType)
    inquire(file=trim(OutputPath)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
    if(.not.dirExists) call system('mkdir '//trim(OutputPath))

    open(2,file=trim(OutputPath)//'Nij_abs.dat')
    write(2,'(a)')               '# Saha populations calculation '
    write(2,'(a,10(2x,a10))')    '# components:',AtomData%elem(:)%title
    write(2,'(a,10(1pE12.4e3))') '# Nk_tot    = ',Ntot
    write(2,'(a,1pE12.4e3)')     '# TK        = ',TK
    write(2,'(a,1pE12.4e3)')     '# TeV       = ',TeV
    write(2,'(a,1pE12.4e3)')     '# -> Ne     = ',Ne
    write(2,'(a)')               '# - - - - - '
    write(2,*)
    write(2,*)


    NQsum = 0
    do i = 1,AtomData%NElems
        NQsum = NQsum + AtomData%elem(i)%NIons
        write(2,'(i3,a11,1x,a)') i,'_component: ',AtomData%elem(i)%title
        write(2,'(a13,1x,a12)')    '# E_ion, eV  ',' populations'
        EQ = AtomData%elem(i)%ion(1)%EQ
        do j = 1,AtomData%elem(i)%NIons
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                dEi = EQd - EQ
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                write(2,'(2(1pE12.4e3,2x),i4)') dEi,SahaNi(idx),idx
            enddo
        enddo
        write(2,*)
        write(2,*)
    enddo
    close(2)

    open(2,file=trim(OutputPath)//'Nij_rel.dat')
    write(2,'(a)')               '# Saha populations calculation '
    write(2,'(a,10(2x,a10))')    '# components:',AtomData%elem(:)%title
    write(2,'(a,10(1pE12.4e3))') '# Nk_tot    = ',Ntot
    write(2,'(a,1pE12.4e3)')     '# TK        = ',TK
    write(2,'(a,1pE12.4e3)')     '# TeV       = ',TeV
    write(2,'(a,1pE12.4e3)')     '# -> Ne     = ',Ne
    write(2,'(a)')               '# - - - - - '
    write(2,*)
    write(2,*)


    do i = 1,AtomData%NElems
        write(2,'(i3,a11,1x,a)') i,'_component: ',AtomData%elem(i)%title
        write(2,'(a13,1x,a12)')    '# E_ion, eV  ',' populations'
        EQ = AtomData%elem(i)%ion(1)%EQ
        do j = 1,AtomData%elem(i)%NIons
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                EQd = AtomData%elem(i)%ion(j)%Qd(k)%EQd
                dEi = EQd - EQ
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                if(Ntot(i).gt.0) then
                    write(2,'(2(1pE12.4e3,2x),i4)') dEi,SahaNi(idx)/Ntot(i),idx
                else
                    write(2,'(2(1pE12.4e3,2x),i4)') dEi,0.d0,idx
                endif
            enddo
        enddo
        write(2,*)
        write(2,*)
    enddo
    close(2)

    ! output sums of ions with same charge
    allocate(SahaQsum(NQsum))
    SahaQsum = 0.d0
    num = 0
    do i = 1,AtomData%NElems
        do j = 1,AtomData%elem(i)%NIons
            num = num + 1
            do k = 1,AtomData%elem(i)%ion(j)%Nqd
                idx = AtomData%elem(i)%ion(j)%Qd(k)%iList
                SahaQsum(num) = SahaQsum(num) + SahaNi(idx)
            enddo
        enddo
    enddo


    open(2,file=trim(OutputPath)//'Ni_abs.dat')
    write(2,'(a)')               '# Saha populations calculation '
    write(2,'(a,10(2x,a10))')    '# components:',AtomData%elem(:)%title
    write(2,'(a,10(1pE12.4e3))') '# Nk_tot    = ',Ntot
    write(2,'(a,1pE12.4e3)')     '# TK        = ',TK
    write(2,'(a,1pE12.4e3)')     '# TeV       = ',TeV
    write(2,'(a,1pE12.4e3)')     '# -> Ne     = ',Ne
    write(2,'(a)')               '# - - - - - '
    write(2,*)
    write(2,*)

    num = 0
    do i = 1,AtomData%NElems
        write(2,'(i3,a11,1x,a)') i,'_component: ',AtomData%elem(i)%title
        write(2,'(a13,1x,a12)')    '# E_ion, eV  ',' populations'
        EQ = AtomData%elem(i)%ion(1)%EQ
        do j = 1,AtomData%elem(i)%NIons
            num = num + 1
            dEi = AtomData%elem(i)%ion(j)%EQ - EQ
            write(2,'(2(1pE12.4e3,2x))') dEi,SahaQsum(num)
        enddo
        write(2,*)
        write(2,*)
    enddo
    close(2)

    open(2,file=trim(OutputPath)//'Ni_rel.dat')
    write(2,'(a)')               '# Saha populations calculation '
    write(2,'(a,10(2x,a10))')    '# components:',AtomData%elem(:)%title
    write(2,'(a,10(1pE12.4e3))') '# Nk_tot    = ',Ntot
    write(2,'(a,1pE12.4e3)')     '# TK        = ',TK
    write(2,'(a,1pE12.4e3)')     '# TeV       = ',TeV
    write(2,'(a,1pE12.4e3)')     '# -> Ne     = ',Ne
    write(2,'(a)')               '# - - - - - '
    write(2,*)
    write(2,*)


    num = 0
    do i = 1,AtomData%NElems
        write(2,'(i3,a11,1x,a)') i,'_component: ',AtomData%elem(i)%title
        write(2,'(a13,1x,a12)')    '# E_ion, eV  ',' populations'
        EQ = AtomData%elem(i)%ion(1)%EQ
        do j = 1,AtomData%elem(i)%NIons
            num = num + 1
            dEi = AtomData%elem(i)%ion(j)%EQ - EQ
            if(Ntot(i).gt.0) then
                write(2,'(2(1pE12.4e3,2x))') dEi,SahaQsum(num)/Ntot(i)
            else
                write(2,'(2(1pE12.4e3,2x))') dEi,0.d0
            endif
        enddo
        write(2,*)
        write(2,*)
    enddo
    close(2)


    deallocate(SahaQsum)
    return
end subroutine SahaOutput

! ********************************************************************
function ShSigm12(N,Ci,x)
    implicit none
    integer                :: N
    real(8),dimension(N)   :: Ci
    real(8)                :: x,xmin
    real(8)                :: ShSigm12
    real(8)                :: sigm1,sigm2

    xmin = 1.d-80
    if(x.ge.xmin) then
        ShSigm12 = sigm1(N,Ci,x)/sigm2(N,Ci,x)
!        write(*,'(a,i3,10(1pE12.4e3,2x))') 'fc: ',N,x,sigm1(N,Ci,x),sigm2(N,Ci,x),ShSigm12
!        pause
    else
        ShSigm12 = (N-1)
    endif

    return
end function ShSigm12

! ---------------------------------------------------------
function dShSigm12(N,Ci,x)
    implicit none
    integer                :: N
    real(8),dimension(N)   :: Ci
    real(8)                :: x,xmin
    real(8)                :: dShSigm12
    real(8)                :: sigm1,dsigm1,sigm2,dsigm2
    real(8)      :: eps1,eps2

    xmin = 1.d-80
    if(x.ge.xmin) then
        eps1 = dsigm1(N,Ci,x)/sigm2(N,Ci,x)
        eps2 = sigm1(N,Ci,x)/sigm2(N,Ci,x)
        eps2 = eps2*(dsigm2(N,Ci,x)/sigm2(N,Ci,x))
        dShSigm12 = (eps1-eps2)
    elseif(Ci(N-1).gt.0.d0) then
        dShSigm12 = -1.d0/Ci(N-1)               ! <!> Ci = 0 in some cases
    else
        dShSigm12 = 0.d0
    endif

    return
end function dShSigm12


! ---------------------------------------------------------
function sigm1(N,Ci,x)
    implicit none
    integer                :: N
    real(8),dimension(N)   :: Ci
    real(8)                :: x,xmin
    real(8)                :: sigm1
    real(8)          :: xPk
    integer          :: k,j

    xmin = 1.d-10
    sigm1 = 0.d0
    if(x.lt.xmin) then
        continue
        sigm1 = 0.d0
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)
            enddo
            sigm1 = sigm1 + (k-1)*x**(N-k)*xPk
        enddo
    else
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)/x
            enddo
            sigm1 = sigm1 + (k-1)*xPk
        enddo
    endif

    return
end function sigm1

! ---------------------------------------------------------
function dsigm1(N,Ci,x)
    implicit none
    integer                :: N
    real(8),dimension(N)   :: Ci
    real(8)                :: x,xmin
    real(8)                :: dsigm1
    real(8)          :: xPk
    integer          :: k,j

    xmin = 1.d-10
    dsigm1 = 0.d0
    if(x.lt.xmin) then
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)
            enddo
            dsigm1 = dsigm1 + (k-1)*(N-k)*x**(N-k-1)*xPk
        enddo
    else
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)/x
            enddo
            dsigm1 = dsigm1 - (k-1)**2/x*xPk
        enddo
    endif


    return
end function dsigm1

! ---------------------------------------------------------
function sigm2(N,Ci,x)
    implicit none
    integer                :: N
    real(8),dimension(N)   :: Ci
    real(8)                :: x
    real(8)                :: sigm2
    real(8)          :: xPk,xmin
    integer          :: k,j

    xmin = 1.d-10
    sigm2 = 1.d0
    if(x.lt.xmin) then
        sigm2 = x**(N-1)
        continue
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)
            enddo
            sigm2 = sigm2 + x**(N-k)*xPk
        enddo
    else
        sigm2 = 1.d0
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)/x
            enddo
            sigm2 = sigm2 + xPk
        enddo
    endif

    return
end function sigm2

! ---------------------------------------------------------
function dsigm2(N,Ci,x)
    implicit none
    integer                :: N
    real(8),dimension(N)   :: Ci
    real(8)                :: x,xmin
    real(8)                :: dsigm2
    real(8)          :: xPk
    integer          :: k,j

    xmin = 1.d-10
    dsigm2 = 0.d0
    if(x.lt.xmin) then
        dsigm2 = (N-1)*x**(N-2)
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)
            enddo
            dsigm2 = dsigm2 + (N-k)*x**(N-k-1)*xPk
        enddo
    else
        do k=2,N
            xPk = 1.d0
            do j=1,k-1
                xPk = xPk*Ci(j)/x
            enddo
            dsigm2 = dsigm2 - (k-1)/x*xPk
        enddo
    endif


    return
end function dsigm2

! ---------------------------------------------------------
subroutine electron_stat_weight(TK,SWe)
    use module_constants

    implicit none
    real(8) :: TK,SWe
    real(8) :: lam,lam3

    !lam  = dsqrt(h_Pl**2/(2.d0*pi*m_e*k_Bol*TK))
    !lam3 = dsqrt(h_Pl**2/(2.d0*pi*m_e*k_Bol))**3
    !write(*,*) lam3 * 0.5d0*dsqrt(k_Bol/eV_erg)**3        ! = c1

    lam  = 7.454e-08/TK**0.5d0   ! de Broil wave length
    lam3 = 4.141053e-16/TK**1.5d0   ! lam**3
    SWe = 2.d0/lam3

    return
end subroutine electron_stat_weight
