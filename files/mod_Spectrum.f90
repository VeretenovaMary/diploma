MODULE module_Spectrum
! <!>  - mark places for check
    use module_Settings

    type TSpectrGrid
        integer(i4p)       :: NX
        real(r8p), pointer :: X(:) => NULL()
        character(8)       :: units
    end type TSpectrGrid

    type TGridIntervalData
        integer   :: GrType
        real(r8p) :: Xmin       ! grid interval range
        real(r8p) :: Xmax       ! grid interval range
    end type TGridIntervalData

    type TGridIntervals
        integer                                  :: NX
        type(TGridIntervalData),dimension(10000) :: Xint
    end type TGridIntervals

    type TQtbbGrData                   ! spectral grid groups range for bb-transitions
        integer(i4p) :: iGr_min        ! left group index
        integer(i4p) :: iGr_max        ! right group index
        integer(i2p) :: TpProf         ! profile type (box/Dopler/Voigt)
        real(r8p)    :: mult           ! normalization coefficient
    end type TQtbbGrData

    type TSpectrData
        type(TSpectrGrid)          :: EnGrid
        type(TQtbbGrData), pointer :: QtBBOpts(:) => NULL()
    end type TSpectrData

    type(TSpectrData) :: SpectrData

    public TSpectrData
    public SpectrData

CONTAINS

! ------------------- Planck function ----------------------------------
    subroutine close_spectral_data()
        implicit none

        deallocate(SpectrData%EnGrid%X)
        deallocate(SpectrData%QtBBOpts)

        return
    end subroutine close_spectral_data
! -------------------------------------------------------------------------

END MODULE module_Spectrum


! ------------------- Planck function ----------------------------------
subroutine set_QtBBOptions(SpGridOpts,SpRange)
    use module_Settings
    use module_AtomData
    use module_Spectrum
    implicit none
    ! input
    type TSpRange
        sequence
        real(r8p)               :: Xmin,Xmax
        character(16)           :: GridType
        character(16)           :: ProfType
        integer(i4p)            :: NGrLine  ! number of groups for detailed line
        integer(i4p)            :: NGrCntr  ! number of groups for detailed line center
        real(r8p)               :: NDWLine  ! number of Dopler widths per detailed line
        real(r8p)               :: NDWCntr  ! number of Dopler widths per detailed line center
    end type TSpRange

    type TSpGridOpts
        sequence
        character(4)                  :: SpUnits
        integer(i4p)                  :: NSpRange             ! number of spectral intervals
        real(r8p)                     :: Xmin,Xmax
        real(r8p)                     :: dXmin,dXmax
        real(r8p)                     :: SpTemp               ! characteristic temperature for line width calculation
    end type TSpGridOpts

    type(TSpGridOpts)             :: SpGridOpts
    type(TSpRange), dimension(20) :: SpRange         ! spectral intervals

    ! local
    integer(i4p) :: i,j
    integer(i4p) :: iflag,NLnGr
    real(r8p)    :: X1,X2,dE12

    allocate(SpectrData%QtBBOpts(AtomData%NQt_bb))
    do i = 1,AtomData%NQt_bb
        dE12 = AtomData%QtbbList(i)%dE12
        SpectrData%QtBBOpts(i)%TpProf  = 0
        SpectrData%QtBBOpts(i)%iGr_min = 0
        SpectrData%QtBBOpts(i)%iGr_max = 0

        iflag = 0
        NLnGr = 1
        do j = 1,SpGridOpts%NSpRange
            X1 = SpRange(j)%Xmin
            X2 = SpRange(j)%Xmax
            if(dE12.gt.X1 .and. dE12.le.X2) then
                iflag = 1                                      ! line is inside detailed spectral range
                NLnGr = SpRange(j)%NGrLine
                selectcase(trim(SpRange(j)%ProfType))
                case("box")
                    SpectrData%QtBBOpts(i)%TpProf  = 0         ! 0/1/2/3 = box/Dopler/Lorentz/Voigt
                    SpectrData%QtBBOpts(i)%iGr_min = 0
                    SpectrData%QtBBOpts(i)%iGr_max = 0
                case("Dopler")
                    SpectrData%QtBBOpts(i)%TpProf  = 1         ! 0/1/2/3 = box/Dopler/Lorentz/Voigt
                    SpectrData%QtBBOpts(i)%iGr_min = -NLnGr
                    SpectrData%QtBBOpts(i)%iGr_max = +NLnGr
                case("Lorentz")
                    SpectrData%QtBBOpts(i)%TpProf  = 2         ! 0/1/2/3 = box/Dopler/Lorentz/Voigt
                    SpectrData%QtBBOpts(i)%iGr_min = -NLnGr
                    SpectrData%QtBBOpts(i)%iGr_max = +NLnGr
                case("Voigt")
                    SpectrData%QtBBOpts(i)%TpProf  = 3         ! 0/1/2/3 = box/Dopler/Lorentz/Voigt
                    SpectrData%QtBBOpts(i)%iGr_min = -NLnGr
                    SpectrData%QtBBOpts(i)%iGr_max = +NLnGr
                end select
            end if
        end do

        do j = 1,SpectrData%EnGrid%NX-1
            X1 = SpectrData%EnGrid%X(j)
            X2 = SpectrData%EnGrid%X(j+1)
            if(dE12.gt.X1 .and. dE12.le.X2) then
                ! <!>
                AtomData%QtbbList(i)%dE12 = 0.5d0*(X1+X2)        ! <!> energy shift to the center of group
                ! <!>
                selectcase(SpectrData%QtBBOpts(i)%TpProf)
                case(0)
                    SpectrData%QtBBOpts(i)%iGr_min = j
                    SpectrData%QtBBOpts(i)%iGr_max = j
                case(1)
                    SpectrData%QtBBOpts(i)%iGr_min = max(j-NLnGr,1)
                    SpectrData%QtBBOpts(i)%iGr_max = min(j+NLnGr,SpectrData%EnGrid%NX-1)
                case(2)
                    SpectrData%QtBBOpts(i)%iGr_min = max(j-NLnGr,1)
                    SpectrData%QtBBOpts(i)%iGr_max = min(j+NLnGr,SpectrData%EnGrid%NX-1)
                case(3)
                    SpectrData%QtBBOpts(i)%iGr_min = max(j-NLnGr,1)
                    SpectrData%QtBBOpts(i)%iGr_max = min(j+NLnGr,SpectrData%EnGrid%NX-1)
                endselect
            end if
        end do
    end do


    return
end subroutine set_QtBBOptions


! ------------------- Planck function ----------------------------------
subroutine Planck(TeV,Xk,F_Pl)
    use module_Settings
    use module_Constants
    implicit none
    ! input
    real(r8p) :: TeV,Xk
    ! output
    real(r8p) :: F_Pl
    ! local
    real(r8p) :: c1
    !F_Pl = 1.d0*eV_erg**4/(4*pi**3*v_c**2*h_Pl**3)*Xk**3/(dexp(Xk/TeV)-1)
    !c1 = 2.d0*eV_erg**4/(v_c**2*h_Pl**3)

    c1 = 5.0403474695329d+010

    F_Pl = 0.d0
    if(Xk/TeV.le.80.d0) then
!        F_Pl = c1*Xk**3/(dexp(Xk/TeV)-1)                   ! erg/(eV*cm2*s)
        F_Pl = c1*Xk**3*dexp(-Xk/TeV)/(1.d0-dexp(-Xk/TeV))  ! erg/(eV*cm2*s)
    endif
    return
end subroutine Planck
! -------------------------------------------------------------------------


! ------------------- ionization cross section (Castor (7.68)) ---------------------------------------
subroutine sigma_pi(ilev,Epi,Xk,s_pi)
    use module_Settings
    use module_Constants
    implicit none
    ! input
    integer   :: ilev
    real(r8p) :: Epi
    real(r8p) :: Xk
    ! output
    real(r8p) :: s_pi
    ! local
    !real(r8p) :: c1
    ! n    - electron level number
    ! Epi  - photoionization energy
    ! S_pi - cross section value, cm2
    ! c1 = 64.d0/(3.d0*3**0.5d0)*alpha*pi*a_0**2
    ! c1 = 7.90707093084d-018
    s_pi = 7.9071d-018*0.5d0*(sign(1.d0,Xk-Epi)+1)*ilev*(Epi/Xk)**3
    s_pi = s_pi/2.d0**0.5
    if(s_pi.lt.0.d0) then
        s_pi = 0.d0
    end if

    return
end subroutine sigma_pi

! ------------------- ionization cross section (Castor (7.68)) ---------------------------------------
subroutine sigma_pi_arr(ilev,Epi,s_pi)
    use module_Settings
    use module_Constants
    use module_Spectrum
    implicit none
    ! input
    integer   :: ilev
    real(r8p) :: Epi
    ! output
    real(r8p),dimension(SpectrData%EnGrid%NX-1) :: s_pi
    ! local
    integer   :: k
    real(r8p) :: Xk
    !real(r8p) :: c1
    ! n    - electron level number
    ! Epi  - photoionization energy
    ! S_pi - cross section value, cm2
    ! c1 = 64.d0/(3.d0*3**0.5d0)*alpha*pi*a_0**2
    ! c1 = 7.90707093084d-018
    do k = 1,SpectrData%EnGrid%NX-1
        Xk   = 0.5*(SpectrData%EnGrid%X(k)+SpectrData%EnGrid%X(k+1))
        s_pi(k) = (sign(1.d0,Xk-Epi)+1)*(Epi/Xk)**3
        if(s_pi(k).lt.0.d0) then            ! <!> - ?
            s_pi(k) = 0.d0
        end if
    enddo
    s_pi = s_pi / 2.d0**0.5
    s_pi = s_pi * 7.9071d-018*0.5d0*ilev

    return
end subroutine sigma_pi_arr

! ------------------- box profile (dimension) -----------------------------------
subroutine prof_box(dX,prof)
    use module_Settings
    use module_Constants
    implicit none
    ! input
    real(r8p) :: dX
    ! output
    real(r8p) :: prof
    ! local

    ! --------------
    prof = 1.d0/dX

    return
end subroutine prof_box

! ------------------- Voigt profile (value) -----------------------------------
subroutine prof_Voigt(Ne,TeV,iQtbb,Xk,prof)
    use module_Settings
    use module_Constants
    use module_AtomData
    implicit none
    ! input
    real(r8p) :: Ne,TeV
    real(r8p) :: Xk
    integer   :: iQtbb
    ! output
    real(r8p) :: prof
    ! local
    integer   :: j
    real(r8p) :: dEeV
    real(r8p) :: Dw,DeV,Lw
    real(r8p) :: y
    real(r8p), parameter :: A = 1.d0          ! g/mole
    real(r8p), dimension(4), parameter :: si = [0.133776446996, 0.62432469019, &
                                                1.342537825640, 2.26266447700]
    real(r8p), dimension(4), parameter :: wi = [0.325302999700, 0.42110710185, &
                                                0.133344250036,0.006374323490]
    real(r8p) :: w,v_T
    real(r8p) :: gf,J1,J2
    real(r8p) :: C_EX,C_DEX

    ! --------------
    dEeV = AtomData%QtbbList(iQtbb)%dE12
    J1   = AtomData%QtbbList(iQtbb)%gQd1
    J2   = AtomData%QtbbList(iQtbb)%gQd2
    gf   = AtomData%QtbbList(iQtbb)%gf12

    w   = dEeV*eV_erg/h_Pl                                    ! frequency
    v_T = sqrt(2*TeV*eV_erg*N_Av/A)                           ! termal velocity
    Dw  = w*v_T/v_c                                           ! Dopler width, [Dw] = s
    DeV = Dw*h_Pl/eV_erg
    Lw  = Ne*C_EX(TeV,dEeV,gf) + Ne*C_DEX(TeV,dEeV,gf,J1,J2)  ! Lorentz width, [Lw] = s

    y = Lw/Dw

    prof = 0.d0
    do j = 1,4
        prof = prof + 2.71828/pi*wi(j)*                                   &
           ((1.d0+y)*dcos(2*si(j))-((Xk-dEeV)/DeV-si(j))*dsin(2*si(j)))/  &
           (((Xk-dEeV)/DeV-si(j))**2+(1.d0+y)**2)
        prof = prof + 2.71828/pi*wi(j)*                                   &
           ((1.d0+y)*dcos(2*si(j))+((Xk-dEeV)/DeV+si(j))*dsin(2*si(j)))/  &
           (((Xk-dEeV)/DeV+si(j))**2+(1.d0+y)**2)
    end do

    if(abs(Xk-dEeV).gt.5.d0) then
        prof = 0.d0            !
    endif
!    prof = 1.0d0/(Dw*pi**0.5) * prof
    prof = 1.0d0/(DeV*pi**0.5) * prof

    return
end subroutine prof_Voigt


! ------------------- Voigt profile (value) -----------------------------------
subroutine prof_Voigt_arr(Ne,TeV,iQtbb,prof)
    use module_Settings
    use module_Constants
    use module_AtomData
    use module_Spectrum
    implicit none
    ! input
    real(r8p) :: Ne,TeV
    integer   :: iQtbb
    ! output
    real(r8p), dimension(SpectrData%EnGrid%NX-1) :: prof
    ! local
    real(r8p) :: Xk,dX
    integer   :: i,j,imin,imax,NXgr
    real(r8p) :: dEeV
    real(r8p) :: Dw,DeV,Lw
    real(r8p) :: y,mult
    real(r8p), parameter :: A = 1.d0          ! g/mole
    real(r8p), dimension(4), parameter :: si = [0.133776446996, 0.62432469019, &
                                                1.342537825640, 2.26266447700]
    real(r8p), dimension(4), parameter :: wi = [0.325302999700, 0.42110710185, &
                                                0.133344250036,0.006374323490]
    real(r8p) :: w,v_T
    real(r8p) :: gf,J1,J2
    real(r8p) :: C_EX,C_DEX

    ! --------------
    dEeV = AtomData%QtbbList(iQtbb)%dE12
    J1   = AtomData%QtbbList(iQtbb)%gQd1
    J2   = AtomData%QtbbList(iQtbb)%gQd2
    gf   = AtomData%QtbbList(iQtbb)%gf12

    w   = dEeV*eV_erg/h_Pl                                    ! frequency
    v_T = sqrt(2*TeV*eV_erg*N_Av/A)                           ! termal velocity
    Dw  = w*v_T/v_c                                           ! Dopler width, [Dw] = s
    DeV = Dw*h_Pl/eV_erg
    Lw  = Ne*C_EX(TeV,dEeV,gf) + Ne*C_DEX(TeV,dEeV,gf,J1,J2)  ! Lorentz width, [Lw] = s

    y = Lw/Dw

    prof(:) = 1.0d-40
    mult    = 0.0d0
    NXgr = SpectrData%EnGrid%NX-1
    imin = SpectrData%QtBBOpts(iQtbb)%iGr_min
    imax = SpectrData%QtBBOpts(iQtbb)%iGr_max
    do i = imin,imax
        dX  = (SpectrData%EnGrid%X(i+1)-SpectrData%EnGrid%X(i))
        Xk  = 0.5*(SpectrData%EnGrid%X(i)+SpectrData%EnGrid%X(i+1))
        do j = 1,4
            prof(i) = prof(i) + 2.71828/pi*wi(j)*                             &
               ((1.d0+y)*dcos(2*si(j))-((Xk-dEeV)/DeV-si(j))*dsin(2*si(j)))/  &
               (((Xk-dEeV)/DeV-si(j))**2+(1.d0+y)**2)
            prof(i) = prof(i) + 2.71828/pi*wi(j)*                             &
               ((1.d0+y)*dcos(2*si(j))+((Xk-dEeV)/DeV+si(j))*dsin(2*si(j)))/  &
               (((Xk-dEeV)/DeV+si(j))**2+(1.d0+y)**2)
        end do
        if(abs(Xk-dEeV).gt.5.d0) then          ! <!> - ?
            prof(i) = 0.d0            !
        endif
        mult = mult + prof(i)*dX      ! normlization coefficient
    enddo

!    prof = 1.0d0/(Dw*pi**0.5) * prof
    prof = 1.0d0/(DeV*pi**0.5) * prof

    ! normalzation: int(-inf:+inf){prof(X)*dX} = 1
    mult = mult / (DeV*pi**0.5)
    prof = prof / mult

    return
end subroutine prof_Voigt_arr


! --------- radiation transport coefficients ------------
subroutine kappa_etta(TeV,Ne,N,KAPPA,ETTA)
    use module_Settings
    use module_Constants
    use module_AtomData
    use module_Spectrum

    implicit none

    integer   :: ilev
    real(r8p) :: Ne,TeV
    real(r8p) :: Zion,Z2Ni
    real(r8p) :: dE12,gf12,gQd1,gQd2
    real(r8p) :: eps
    real(r8p) :: c1,c2
    real(r8p) :: X1,X2,Xk,dXk,prof,s_pi
    real(r8p), dimension(SpectrData%EnGrid%NX-1) :: KAPPA,ETTA
    real(r8p) :: dkappa,detta
    real(r8p), dimension(AtomData%NQds)          :: N
    integer :: NEnGroups
    integer :: i,j
    integer(i4p) :: iGr_min,iGr_max
    integer      :: iLQd1,iLQd2
    integer(i2p) :: TpProf
    integer      :: iE1,iE2,iI1,iI2,iQ1,iQ2

    NEnGroups = SpectrData%EnGrid%NX-1
    KAPPA(:) = 0.d0
    ETTA(:)  = 0.d0
    dkappa = 0.d0
    detta  = 0.d0

! ------------------- B-B -------------------------
    !c1 = pi*q_e**2/(m_e*v_c)
    !c1 = c1*h_Pl/eV_erg/(2*pi)        ! correction due to profile units: 1/nu -> 1/EeV
    c1 = 1.749298d-17
    !c2 = c1 * 2*eV_erg**4/(h_Pl**3*v_c**2)
    c2 = 8.817102d-07

    do i = 1,AtomData%NQt_bb
        iGr_min = SpectrData%QtBBOpts(i)%iGr_min
        iGr_max = SpectrData%QtBBOpts(i)%iGr_max
        iLQd1    = AtomData%QtbbList(i)%iLQd1
        iLQd2    = AtomData%QtbbList(i)%iLQd2
        dE12     = AtomData%QtbbList(i)%dE12
        gf12     = AtomData%QtbbList(i)%gf12
        gQd1     = AtomData%QtbbList(i)%gQd1
        gQd2     = AtomData%QtbbList(i)%gQd2
        dkappa = 0.d0
        TpProf = SpectrData%QtBBOpts(i)%TpProf
        !write(*,*) i,AtomData%NQt_bb," TpProf: ",TpProf
        if(iGr_min.gt.0) then       ! group is inside spectral range
            do j = iGr_min,iGr_max
                X1 = SpectrData%EnGrid%X(j)
                X2 = SpectrData%EnGrid%X(j+1)
                Xk  = 0.5d0*(X1+X2)
                dXk = X2-X1
                prof = 0.d0
                select case(TpProf)
                    case(0)           ! box
                        call prof_box(dXk,prof)
                    case(3)           ! Voigt
                        call prof_Voigt(Ne,TeV,i,Xk,prof)
                end select
                dkappa = c1*      gf12*prof*(N(iLQd1) - N(iLQd2)*(gQd1/gQd2))
                detta  = c2*Xk**3*gf12*prof*N(iLQd2)*(gQd1/gQd2)
                if(dkappa.lt.0.e0) then
                    dkappa = c1*      gf12*prof*N(iLQd1)*(1.0e0_r8p - exp(-dE12/TeV)) ! <!> <?> - exp(-Xk/TeV)
                    detta  = c2*Xk**3*gf12*prof*N(iLQd1)*exp(-dE12/TeV)               ! <!> <?> - exp(-Xk/TeV)
                end if
                if(detta.lt.0.e0 .or. dkappa.lt.0.e0) then
                    write(*,*) "dkappa(bb) = ",dkappa
                    write(*,*) "detta(bb)  = ",detta
                    write(*,*) "prof       = ",prof
                    write(*,*) "N(iLQd1)   = ",N(iLQd1)
                    write(*,*) "N(iLQd2)   = ",N(iLQd2)
                    write(*,*) "dE12       = ",dE12
                    write(*,*) "TeV        = ",TeV
                    iE1 = AtomData%QdList(iLQd1)%iE
                    iE2 = AtomData%QdList(iLQd2)%iE
                    iI1 = AtomData%QdList(iLQd1)%iI
                    iI2 = AtomData%QdList(iLQd2)%iI
                    iQ1 = AtomData%QdList(iLQd1)%iQ
                    iQ2 = AtomData%QdList(iLQd2)%iQ
                    write(*,*) "E1         = ",AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%EQd
                    write(*,*) "E2         = ",AtomData%elem(iE2)%ion(iI2)%Qd(iQ2)%EQd

                    detta  = 0.d0
                    write(*,*) " new -> dkappa(bb) = ",dkappa
                    write(*,*) " new -> detta(bb)  = ",detta
                    stop
                end if

                KAPPA(j) = KAPPA(j) + dkappa
                ETTA(j)  = ETTA(j)  + detta
            end do
        endif
    end do


! ------------------- B-F -------------------------
!    c1 = 0.5*(2*pi/m_e*h_Pl**2/eV_erg)**1.5
    c1 = 1.656421282064654d-22
    c1 = c1/TeV**1.5
    !c2 = 2*eV_erg**4/(h_Pl**3*v_c**2)
    c2 = 5.0403474695329d+10
    do i = 1,AtomData%NQt_bf
        iLQd1    = AtomData%QtbfList(i)%iLQd1
        iLQd2    = AtomData%QtbfList(i)%iLQd2
        dE12     = AtomData%QtbfList(i)%dE12
        gQd1     = AtomData%QtbfList(i)%gQd1
        gQd2     = AtomData%QtbfList(i)%gQd2
        ilev     = AtomData%QtbfList(i)%ilev
!        write(*,"(a,3(1pE12.5e2))") "test: ",c1*gQd1/gQd2*Ne*N(iLQd2),N(iLQd1)*dexp(-dE12/TeV)
        do j = 1,NEnGroups
            X1 = SpectrData%EnGrid%X(j)
            X2 = SpectrData%EnGrid%X(j+1)
            Xk  = 0.5d0*(X1+X2)
            dXk = X2-X1
            call sigma_pi(ilev,dE12,Xk,s_pi)
            dkappa = 0.d0
            detta  = 0.d0
            eps = (Xk-dE12)/TeV
            if(abs(eps).lt.80.d0) then
                ! c1*gQd1/gQd2*Ne*N(iLQd2)*dexp((dE12-Xk)/TeV) = N(iLQd1)*dexp(-Xk/TeV)
                dkappa = s_pi*(N(iLQd1)  - c1*gQd1/gQd2*Ne*N(iLQd2)*dexp(-eps))
                detta  = s_pi*c2*(Xk**3.0)*c1*gQd1/gQd2*Ne*N(iLQd2)*dexp(-eps)
            else
                dkappa = s_pi*N(iLQd1)
                detta  = 0.d0
            endif
            if(dkappa.lt.0.e0) then
                dkappa = s_pi*N(iLQd1)*(1.0e0_r8p - dexp(-Xk/TeV))
                detta  = s_pi*N(iLQd1)*c1*c2*Xk**3*dexp(-Xk/TeV)
            end if

            if(detta.lt.0.e0 .or. dkappa.lt.0.e0) then
                write(*,*) "dkappa(bf) = ",dkappa
                write(*,*) "detta(bf)  = ",detta
                write(*,*) "s_pi       = ",s_pi
                write(*,*) "Ne         = ",Ne
                write(*,*) "N(iLQd1)   = ",N(iLQd1)
                write(*,*) "N(iLQd2)   = ",N(iLQd2)
                write(*,*) "dE12       = ",dE12
                write(*,*) "TeV        = ",TeV
                iE1 = AtomData%QdList(iLQd1)%iE
                iE2 = AtomData%QdList(iLQd2)%iE
                iI1 = AtomData%QdList(iLQd1)%iI
                iI2 = AtomData%QdList(iLQd2)%iI
                iQ1 = AtomData%QdList(iLQd1)%iQ
                iQ2 = AtomData%QdList(iLQd2)%iQ
                write(*,*) "E1         = ",AtomData%elem(iE1)%ion(iI1)%Qd(iQ1)%EQd
                write(*,*) "E2         = ",AtomData%elem(iE2)%ion(iI2)%Qd(iQ2)%EQd
                detta  = 0.d0
                write(*,*) " new -> dkappa(bf) = ",dkappa
                write(*,*) " new -> detta(bf)  = ",detta
                stop
            end if
            !dkappa = 0.d0
            !detta  = 0.d0
            KAPPA(j) = KAPPA(j) + dkappa
            ETTA(j)  = ETTA(j)  + detta
        end do
    end do

!
! ------------ free-free (Castor p. 168) --------------------------------
    !c1 = 4.d0/3.d0*sqrt(2*pi/3)*q_e**6*h_Pl**2/(v_c*m_e**1.5d0*eV_erg**3.5d0)
    c1 = 2.43451927d-37
    !c2 = 2*eV_erg**4/(h_Pl**3*v_c**2)
    c2 = 5.04036611d+10

    Z2Ni = 0.d0
    do i = 1,AtomData%NQds
        Zion = AtomData%QdList(i)%iI - 1
        Z2Ni = Z2Ni + Zion**2*N(i)
    end do

    do j = 1,NEnGroups
        X1 = SpectrData%EnGrid%X(j)
        X2 = SpectrData%EnGrid%X(j+1)
        Xk  = 0.5d0*(X1+X2)
        dkappa = 0.d0
        detta  = 0.d0
        eps = Xk/TeV
        if(abs(eps).lt.1.d-2) then
            dkappa = c1    * (Ne*Z2Ni)/   &
                      Xk**3/TeV**0.5 * (0.d0+Xk/TeV)
            detta  = c1*c2 * (Ne*Z2Ni)/   &
                            TeV**0.5 * (1.d0-Xk/TeV)
        elseif(abs(eps).lt.80.d0) then
            dkappa = c1    * (Ne*Z2Ni)/   &
                      Xk**3/TeV**0.5 * (1.d0-dexp(-Xk/TeV))
            detta  = c1*c2 * (Ne*Z2Ni)/   &
                            TeV**0.5 * dexp(-Xk/TeV)
        else
            dkappa = c1    * (Ne*Z2Ni)/   &
                      Xk**3/TeV**0.5 * (1.d0-0.d0)
            detta  = 0.d0
        endif
!                if(dkappa.le.0.e0) then
!                    write(*,*) "dkappa(ff) = ",dkappa
!                    write(*,*) "Xk/TeV     = ",eps
!                    write(*,*) "Z2Ni       = ",Z2Ni
!                    pause
!                end if
        KAPPA(j) = KAPPA(j) + dkappa
        ETTA(j)  = ETTA(j)  + detta
    enddo


!    open(34,file=trim(OutputDir)//"/spectr/fPl_kp_et.txt")
!    do j = 1,SpectrData%EnGrid%NX-1
!        X1 = SpectrData%EnGrid%X(j)
!        X2 = SpectrData%EnGrid%X(j+1)
!        Xk  = 0.5d0*(X1+X2)
!        call Planck(TeV,Xk,eps)
!        write(34,"(4(1pE12.5e2))") Xk,eps,KAPPA(j),ETTA(j)
!    enddo
!    close(34)


    return
end subroutine kappa_etta


! ============================ subroutines =============================
subroutine init_mod_spectrum(cfgName)
    use module_Settings
    use module_Spectrum
    use module_AtomData
    use module_Constants
    use module_TxtRead

    implicit none
    ! input:
    character(*) :: cfgName

    ! local variables
    type TSpRange
        sequence
        real(r8p)               :: Xmin,Xmax
        character(16)           :: GridType
        character(16)           :: ProfType
        integer(i4p)            :: NGrLine  ! number of groups for detailed line
        integer(i4p)            :: NGrCntr  ! number of groups for detailed line center
        real(r8p)               :: NDWLine  ! number of Dopler widths per detailed line
        real(r8p)               :: NDWCntr  ! number of Dopler widths per detailed line center
    end type TSpRange

    type TSpGridOpts
        sequence
        character(4)                  :: SpUnits
        integer(i4p)                  :: NSpRange             ! number of spectral intervals
        real(r8p)                     :: Xmin,Xmax
        real(r8p)                     :: dXmin,dXmax
        real(r8p)                     :: SpTemp               ! characteristic temperature for line width calculation
    end type TSpGridOpts

    ! temporal arrays for detailed spectral ranges:
    real(r8p), dimension(20)     :: DetXmin,DetXmax
    character(16), dimension(20) :: DetGridType
    character(16), dimension(20) :: DetProfType
    integer(i4p), dimension(20)  :: DetNGrLine  ! number of groups for detailed line
    integer(i4p), dimension(20)  :: DetNGrCntr  ! number of groups for detailed line center
    real(r8p), dimension(20)     :: DetNDWLine  ! number of Dopler widths per detailed line
    real(r8p), dimension(20)     :: DetNDWCntr  ! number of Dopler widths per detailed line center

    type(TSpGridOpts)             :: SpGridOpts
    type(TSpRange), dimension(20) :: SpRange         ! spectral intervals

    real(8), dimension(AtomData%NQt_bb) :: Xbb
    real(8), dimension(10000) :: XGrid
!    real(8), dimension(AtomData%NQt_bf) :: Xbf
    integer :: NX,Nbb
!    integer :: Nbb,Nbf
    ! output parameters
    real(r8p) :: X1,X2,dX,eps
    real(r8p) :: Xmin,Xmax
    real(r8p) :: dXmin,dXmax
    integer :: i,j,k
    integer(i4p) :: GetUnit,cfgUnit

    ! ........................................................
    type(TGridIntervals) :: SpGridPart
    ! ........................................................


    ! local
    real(r8p) :: Xact,dXact
    real(r8p) :: TeV,Ne
    real(r8p) :: TK,N_all,IonDeg
    real(r8p), allocatable ::  KAPPA(:), ETTA(:)
    real(8), dimension(AtomData%NQds)   :: SahaNi
    real(8), dimension(AtomData%NElems) :: N_tot
    ! file output
    character(128) :: OutputPath
    integer        :: ios,ipos
    integer        :: irange
    integer        :: NDetRange,NSpRange
    logical        :: dirExists

    ! calculate populations for spectr/kappa-etta output
!    N_tot(:) = 3.0d16          ! for hydrogen
    N_tot(:) = 1.0d15          ! for nitrogen
    N_all    = sum(N_tot)
!    TeV = 1.0d0
    cfgUnit = GetUnit()
    open(cfgUnit,file=cfgName)
    ! read atom names
    !call txt_read_value(cfgUnit,"OutputDir",OutputDir,ios)
    call txt_read_struct_value(cfgUnit,"SpectrData","TeV",TeV,ios)
    write(*,*) "TeV = ",TeV
    TK = TeV*eV_erg/k_Bol

    ! set and sorting Xbb (energy points)
    Xbb(:) = AtomData%QtbbList(:)%dE12
    Nbb = AtomData%NQt_bb
    do i = 1,AtomData%NQt_bb-1
        do j = 1,AtomData%NQt_bb-1
            X1 = Xbb(j)
            X2 = Xbb(j+1)
            if(X1.gt.X2) then
                Xbb(j)   = X2
                Xbb(j+1) = X1
            end if
        end do
    enddo
    ! ... output results
    OutputPath = trim(OutputDir)//"/spectr/"
    call slash_correct(OutputPath,OSType)
    inquire(file=trim(OutputPath)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
    if(.not.dirExists) call system('mkdir '//trim(OutputPath))
    open(407,file=trim(OutputPath)//"Xbb.txt")
    dXact = 1.d-3
    do i = 1,AtomData%NQt_bb
        write(407,"(2(1x,1pE12.4e2))") Xbb(i)-dXact,0.d0
        write(407,"(2(1x,1pE12.4e2))") Xbb(i)-dXact,1.d-1/dXact
        write(407,"(2(1x,1pE12.4e2))") Xbb(i)+dXact,1.d-1/dXact
        write(407,"(2(1x,1pE12.4e2))") Xbb(i)+dXact,0.d0
        write(407,"(2(1x,1pE12.4e2))")
    enddo
    close(407)


    ! set multigroup intervals:
    call txt_read_struct_value(cfgUnit,"SpectrData","TeV",   SpGridOpts%SpTemp,ios)
    call txt_read_struct_value(cfgUnit,"SpectrData","Xunits",SpGridOpts%SpUnits,ios)
    call txt_read_struct_value(cfgUnit,"SpectrData","Xmin",  SpGridOpts%Xmin,ios)
    call txt_read_struct_value(cfgUnit,"SpectrData","Xmax",  SpGridOpts%Xmax,ios)
    call txt_read_struct_value(cfgUnit,"SpectrData","dXmin", SpGridOpts%dXmin,ios)
    call txt_read_struct_value(cfgUnit,"SpectrData","dXmax", SpGridOpts%dXmax,ios)
    ! set detailet intervals:
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","Xmin",      DetXmin,NDetRange,ios)
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","Xmax",      DetXmax,NDetRange,ios)
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","TpLineGrid",DetGridType,NDetRange,ios)
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","TpLineProf",DetProfType,NDetRange,ios)
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","NGrPerLine",DetNGrLine,NDetRange,ios)
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","NGrPerCntr",DetNGrCntr,NDetRange,ios)
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","NDWPerLine",DetNDWLine,NDetRange,ios)
    call txt_read_struct_array(cfgUnit,"SpectrData%DetRange","NDWPerCntr",DetNDWCntr,NDetRange,ios)
    close(cfgUnit)

    write(*,*) " NDetRange = ",NDetRange

    ! read and sort spectral grid intervals
    Xmin  = -1.0d-5
    Xmax  = SpGridOpts%Xmax
    Xact  = Xmax
    X2    = SpGridOpts%Xmin

    irange = 1
    SpRange(irange)%Xmin = SpGridOpts%Xmin         ! 1st interval left bound

    do i = 1,NDetRange
        k = 1
        do j = 1,NDetRange
            X1 = DetXmin(j)
            X2 = DetXmax(j)
            if(Xact.gt.X1 .and. X1.gt.Xmin) then
                Xact = 0.5*(X1+X2)
                k = j
            end if
        end do
        X1 = DetXmin(k)
        X2 = DetXmax(k)

        if(SpRange(irange)%Xmin.lt.X1 .and. X1.lt.X2) then        ! close draft interval
            SpRange(irange)%Xmax = X1
            SpRange(irange)%GridType = "draft"
            SpRange(irange)%ProfType = "box"
            SpRange(irange)%NDWLine = 1
            SpRange(irange)%NDWCntr = 1
            SpRange(irange)%NGrLine = 1
            SpRange(irange)%NGrCntr = 1
            irange = irange + 1
        end if

        if(X1.lt.X2) then
            SpRange(irange)%Xmin = X1
            SpRange(irange)%Xmax = X2
            SpRange(irange)%GridType = DetGridType(k)
            SpRange(irange)%ProfType = DetProfType(k)
            SpRange(irange)%NDWLine = DetNDWLine(k)
            SpRange(irange)%NDWCntr = DetNDWCntr(k)
            SpRange(irange)%NGrLine = DetNGrLine(k)
            SpRange(irange)%NGrCntr = DetNGrCntr(k)
        endif

        ! open draft interval
        if(SpGridOpts%Xmax.gt.X2 .and. X1.lt.X2) then
            irange = irange + 1
            SpRange(irange)%Xmin = X2
        endif

        Xmin = Xact
        Xact = Xmax

    end do
    if(SpGridOpts%Xmax.gt.X2) then              ! close final interval
        SpRange(irange)%Xmax = SpGridOpts%Xmax
        SpRange(irange)%GridType = "draft"
        SpRange(irange)%ProfType = "box"
        SpRange(irange)%NDWLine = 1
        SpRange(irange)%NDWCntr = 1
        SpRange(irange)%NGrLine = 1
        SpRange(irange)%NGrCntr = 1
    endif
    NSpRange = irange
    SpGridOpts%NSpRange = NSpRange

    ! output intervals
!    write(*,*) " NSpRange = ",NSpRange
!    do i = 1,NSpRange
!        X1 = SpRange(i)%Xmin
!        X2 = SpRange(i)%Xmax
!        write(*,"(i3,2(1x,1pE12.4e3),4x,a)") i,X1,X2,trim(SpRange(i)%GridType)
!    end do

    open(409,file=trim(OutputPath)//"XGrid_BB_Draft.txt")
    open(410,file=trim(OutputPath)//"XGrid_BB_Uniform.txt")
    open(411,file=trim(OutputPath)//"XGrid_BB_Progress.txt")


    ! create BB grid
    dXmin = SpGridOpts%dXmin
    dXmax = SpGridOpts%dXmax
    NX       = 1
    XGrid(:) = 0.d0
    X1 = SpGridOpts%Xmin
    XGrid(NX) = X1
    do i = 1,NSpRange
        selectcase(trim(SpRange(i)%GridType))
        case("draft")
            call set_SpGrid_BB_Draft(Xbb,SpGridOpts,SpRange(i),SpGridPart)
        case("uniform")
            call set_SpGrid_BB_Uniform(Xbb,SpGridOpts,SpRange(i),SpGridPart)
        case("progression")
            call set_SpGrid_BB_Progression(Xbb,SpGridOpts,SpRange(i),SpGridPart)
        endselect

        if(SpGridPart%NX.gt.0) then
            do j = 1,SpGridPart%NX
                X1 = XGrid(NX)
                X2 = SpGridPart%Xint(j)%Xmin
                eps = abs(X2-X1)/dXmin
                if(eps.gt.1.0d-8) then
                    NX = NX + 1
                    XGrid(NX) = SpGridPart%Xint(j)%Xmin
                end if

                X1 = XGrid(NX)
                X2 = SpGridPart%Xint(j)%Xmax
                eps = abs(X2-X1)/dXmin
                if(eps.gt.1.0d-8) then
                    NX = NX + 1
                    XGrid(NX) = SpGridPart%Xint(j)%Xmax
                end if
            end do
        end if
    end do
    SpectrData%EnGrid%NX = NX

    close(409)
    close(410)
    close(411)

    ! create background grid (split wide intervals):
    ipos = 2
    do i = 2,SpectrData%EnGrid%NX
        dX = XGrid(ipos) - XGrid(ipos-1)
        eps =       dX*(1.d0+1.d-3)/dXmax
        j   = floor(dX*(1.d0+1.d-3)/dXmax + 0.5d0)
        if (j.ge.1) then
            dXact = dX/j
            do k = 1,NX-ipos+1     ! shift elements to the right (create free space)
                XGrid(NX-k+j) = XGrid(NX-k+1)
            end do
            do k = 1,j-1           ! paste elements into free space
                XGrid(ipos+k-1) = XGrid(ipos-1) + k*dXact
            end do
        else
            j = 1
        end if
        ipos = ipos + j
        NX   = NX + (j-1)
    end do
    ! last (right bound) interval:
    Xmax  = SpGridOpts%Xmax
    dX = Xmax - XGrid(NX)
    eps =       dX*(1.d0+1.d-3)/dXmax
    j   = floor(dX*(1.d0+1.d-3)/dXmax + 0.5d0)
    if (j.ge.1) then
        dXact = dX/j
        do k = 1,j            ! paste elements into free space
            NX = NX + 1
            XGrid(NX) = XGrid(NX-1) + dXact
        end do
    end if
    XGrid(NX) = Xmax
    SpectrData%EnGrid%NX = NX


    open(412,file=trim(OutputPath)//"XGrid.txt")
    do i = 1,NX-1
        dXact = XGrid(i+1) - XGrid(i)
        write(412,"(2(1x,1pE12.4e2))") XGrid(i),0.d0
        write(412,"(2(1x,1pE12.4e2))") XGrid(i),1.d0/dXact
        write(412,"(2(1x,1pE12.4e2))") XGrid(i+1),1.d0/dXact
        write(412,"(2(1x,1pE12.4e2))") XGrid(i+1),0.d0
        write(412,"(2(1x,1pE12.4e2))")
    enddo
    close(412)

    allocate(SpectrData%EnGrid%X(NX))
    SpectrData%EnGrid%X(1:NX) = XGrid(1:NX)
    SpectrData%EnGrid%NX = NX


    ! Qtbb energy location on spectral grid
    call set_QtBBOptions(SpGridOpts,SpRange)


    call SahaKin_TeV_Nt(TeV,N_tot,SahaNi,Ne,IonDeg)
    call SahaOutput(TeV,N_tot,SahaNi,Ne)
    write(*,'(a,1pE12.4E3)') ' IonDeg = ',IonDeg


    open(407,file=trim(OutputPath)//"kppa_etta_inuse.txt")
    allocate(KAPPA(NX),ETTA(NX))
    call kappa_etta(TeV,Ne,SahaNi,KAPPA,ETTA)
    do i = 1,NX-1
        write(407,"(3(1pE12.4e3,1x))") XGrid(i)  ,KAPPA(i),ETTA(i)
        write(407,"(3(1pE12.4e3,1x))") XGrid(i+1),KAPPA(i),ETTA(i)
    enddo
    close(407)

!    open(407,file=trim(sbuf)//"kppa_etta_line.txt")
!    call kappa_etta(TeV,Ne,SahaNi,KAPPA,ETTA)
!
!    open(408,file="config.txt")
!    call txt_read_value(408,"ActLine",i,ios)
!    close(408)
!
!    j = SpectrData%QtBBOpts(i)%iGr_min
!    k = SpectrData%QtBBOpts(i)%iGr_max
!    write(*,*) " iGr_min = ",j
!    write(*,*) " iGr_max = ",k
!    do i = j,k
!        write(407,"(3(1pE12.4e3,1x))") XGrid(i)  ,KAPPA(i),ETTA(i)
!        write(407,"(3(1pE12.4e3,1x))") XGrid(i+1),KAPPA(i),ETTA(i)
!    enddo
!    close(407)

!    pause

    deallocate(KAPPA,ETTA)

    return
end subroutine init_mod_spectrum


! --------------------------------------
subroutine set_SpGrid_BB_Draft(Xbb,SpGridOpts,SpRange,SpIntervals)
    use module_Settings
    use module_AtomData
    use module_Spectrum
    implicit none
    ! input
    real(8), dimension(AtomData%NQt_bb) :: Xbb
    type TSpRange
        sequence
        real(r8p)               :: Xmin,Xmax
        character(16)           :: GridType
        character(16)           :: ProfType
        integer(i4p)            :: NGrLine  ! number of groups for detailed line
        integer(i4p)            :: NGrCntr  ! number of groups for detailed line center
        real(r8p)               :: NDWLine  ! number of Dopler widths per detailed line
        real(r8p)               :: NDWCntr  ! number of Dopler widths per detailed line center
    end type TSpRange
    type TSpGridOpts
        sequence
        character(4)            :: SpUnits
        integer(i4p)            :: NSpRange             ! number of spectral intervals
        real(r8p)               :: Xmin,Xmax
        real(r8p)               :: dXmin,dXmax
        real(r8p)               :: SpTemp               ! characteristic temperature for line width calculation
    end type TSpGridOpts

    type(TSpGridOpts)             :: SpGridOpts
    type(TSpRange)                :: SpRange         ! spectral intervals
    ! output
    ! ........................................................
    type(TGridIntervals) :: SpIntervals

    ! local
    real(r8p) :: Xmin,Xmax
    real(r8p) :: dXmin,dXmax,dXact,dX
    real(r8p) :: Xl,Xr
    real(r8p) :: eps
    integer   :: i,j,k,ipos,NX


    Xmin = SpRange%Xmin
    Xmax = SpRange%Xmax
    dXmin = SpGridOpts%dXmin
    dXmax = SpGridOpts%dXmax

    SpIntervals%Xint(:)%Xmin  = 0.d0
    SpIntervals%Xint(:)%Xmax  = 0.d0
    SpIntervals%NX            = 0


    ! ======================== 1.a create 1st ( left boundary) interval: =======================
    i = 1
    do while(Xbb(i).le.Xmin .and. i.lt.AtomData%NQt_bb)
        i = i + 1
    end do
    NX = 0
    if(Xbb(i).ge.Xmin .and. Xbb(i).lt.Xmax) then
        NX = NX + 1
        if(Xbb(i)-0.5d0*dXmin.gt.Xmin) then
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXmin
        else
            SpIntervals%Xint(NX)%Xmin = Xmin
        end if
        if(Xbb(i)+0.5d0*dXmin.lt.Xmax) then
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXmin
        else
            SpIntervals%Xint(NX)%Xmax = Xmax
        end if
    end if

    if(i.eq.AtomData%NQt_bb) then
        return
    end if
    ! 1.b set other intervals
    i = i + 1
    do while(Xbb(i).le.Xmax .and. i.lt.AtomData%NQt_bb)
        Xl = SpIntervals%Xint(NX)%Xmin
        Xr = SpIntervals%Xint(NX)%Xmax
        if(Xbb(i)-1.5d0*dXmin.gt.Xr) then
            NX = NX + 1
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXmin
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXmin
        else
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXmin
        end if
        i = i + 1
    end do
    if(i.eq.AtomData%NQt_bb) then         ! last Xbb spectral point
        Xl = SpIntervals%Xint(NX)%Xmin
        Xr = SpIntervals%Xint(NX)%Xmax
        if(Xbb(i)-1.5d0*dXmin.gt.Xr) then
            NX = NX + 1
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXmin
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXmin
        else
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXmin
        end if
    end if
    SpIntervals%NX = NX


    ! ========================== 2. correct multigroup intervals (interval length = n * dXmin): ===================
    do i = 1,SpIntervals%NX        ! skip 1st - boundary interval
        dXact = SpIntervals%Xint(i)%Xmax - SpIntervals%Xint(i)%Xmin
        j = floor(dXact*(1.d0+1.d-3)/dXmin)
        eps = 0.5d0*(dXact - j*dXmin)
        SpIntervals%Xint(i)%Xmin = SpIntervals%Xint(i)%Xmin + eps
        SpIntervals%Xint(i)%Xmax = SpIntervals%Xint(i)%Xmax - eps
    end do
    if(NX.gt.0) then
        if(SpIntervals%Xint(1)%Xmin.lt.Xmin) then
            SpIntervals%Xint(1)%Xmin = Xmin
        end if
        if(SpIntervals%Xint(NX)%Xmax.gt.Xmax) then
            SpIntervals%Xint(NX)%Xmax = Xmax
        end if
    endif
    SpIntervals%NX = NX


    ! =========================== 3. split wide intervals: ===============================
    ipos = 1
    do i = 1,SpIntervals%NX
        dX = SpIntervals%Xint(ipos)%Xmax - SpIntervals%Xint(ipos)%Xmin
        j = floor(dX*(1.d0+1.d-3)/dXmin)
        if (j.gt.1) then
            eps = dX/j
            do k = 1,NX-ipos     ! shift elements to the right (create free space)
                SpIntervals%Xint(NX-k+j)%Xmin = SpIntervals%Xint(NX-k+1)%Xmin
                SpIntervals%Xint(NX-k+j)%Xmax = SpIntervals%Xint(NX-k+1)%Xmax
            end do
            SpIntervals%Xint(ipos)%Xmax = SpIntervals%Xint(ipos)%Xmin + eps
            do k = 1,j-1         ! paste elements into free space
                SpIntervals%Xint(ipos+k)%Xmin = SpIntervals%Xint(ipos+k-1)%Xmax
                SpIntervals%Xint(ipos+k)%Xmax = SpIntervals%Xint(ipos+k)%Xmin + eps
            end do
            SpIntervals%Xint(ipos+j-1)%Xmax = SpIntervals%Xint(ipos)%Xmin + dX
        else
            j = 1
        end if
        ipos = ipos + j
        NX   = NX + (j-1)
    end do
    SpIntervals%NX = NX


    write(409,*) "# -------------------------------------------- "
    do i = 1,NX
        dXact = SpIntervals%Xint(i)%Xmax - SpIntervals%Xint(i)%Xmin
        write(409,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmin,0.d0
        write(409,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmin,1.d0/dXact
        write(409,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmax,1.d0/dXact
        write(409,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmax,0.d0
        write(409,"(2(1x,1pE12.4e2))")
    enddo


    return
end subroutine


! ---------------------------------------------------------------
subroutine set_SpGrid_BB_Uniform(Xbb,SpGridOpts,SpRange,SpIntervals)
    use module_Settings
    use module_Constants
    use module_AtomData
    use module_Spectrum
    implicit none
    ! input
    real(8), dimension(AtomData%NQt_bb) :: Xbb
    type TSpRange
        sequence
        real(r8p)               :: Xmin,Xmax
        character(16)           :: GridType
        character(16)           :: ProfType
        integer(i4p)            :: NGrLine  ! number of groups for detailed line
        integer(i4p)            :: NGrCntr  ! number of groups for detailed line center
        real(r8p)               :: NDWLine  ! number of Dopler widths per detailed line
        real(r8p)               :: NDWCntr  ! number of Dopler widths per detailed line center
    end type TSpRange
    type TSpGridOpts
        sequence
        character(4)            :: SpUnits
        integer(i4p)            :: NSpRange             ! number of spectral intervals
        real(r8p)               :: Xmin,Xmax
        real(r8p)               :: dXmin,dXmax
        real(r8p)               :: SpTemp               ! characteristic temperature for line width calculation
    end type TSpGridOpts

    type(TSpGridOpts)             :: SpGridOpts
    type(TSpRange)                :: SpRange         ! spectral intervals
    ! output
    ! ........................................................
    type(TGridIntervals) :: SpIntervals

    ! local
    real(r8p) :: Xact,Xmin,Xmax
    real(r8p) :: dXmin,dXmax,dXact,dX
    real(r8p) :: Xl,Xr,dXl,dXr
    real(r8p) :: eps
    integer   :: i,j,k,ipos,NX
    real(r8p) :: NDWWing
    integer   :: NGrWing
    real(r8p) :: DWCoef,TeV,Mw           ! Dopler width coefficient, Dopler width in eV, molar weight
    real(8),dimension(100) :: XWingGrid
    integer                :: NWx


    Xmin = SpRange%Xmin
    Xmax = SpRange%Xmax
    dXmin = SpGridOpts%dXmin
    dXmax = SpGridOpts%dXmax
    Mw     = 1.0d0
    TeV    = SpGridOpts%SpTemp
    DWCoef = sqrt(2*TeV*eV_erg*N_Av/Mw)/v_c


    SpIntervals%Xint(:)%Xmin  = 0.d0
    SpIntervals%Xint(:)%Xmax  = 0.d0
    SpIntervals%NX            = 0


    ! ======================== 1.a create 1st ( left boundary) interval: =======================
    i = 1
    do while(Xbb(i).le.Xmin .and. i.lt.AtomData%NQt_bb)
        i = i + 1
    end do
    NX = 0
    if(Xbb(i).ge.Xmin .and. Xbb(i).lt.Xmax) then
        NX = NX + 1
        dX = Xbb(i)*DWCoef*SpRange%NDWCntr
        dXact = min(dX,dXmin)
        if(Xbb(i)-0.5d0*dXact.gt.Xmin) then
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmin = Xmin
        end if
        if(Xbb(i)+0.5d0*dXact.lt.Xmax) then
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmax = Xmax
        end if
    end if
    if(i.eq.AtomData%NQt_bb) then
        return
    end if

    ! 1.b set other intervals
    i = i + 1
    do while(Xbb(i).le.Xmax .and. i.lt.AtomData%NQt_bb)
        Xl = SpIntervals%Xint(NX)%Xmin
        Xr = SpIntervals%Xint(NX)%Xmax
        dX = Xbb(i)*DWCoef*SpRange%NDWCntr
        dXact = min(dX,dXmin)
        if(Xbb(i)-1.5d0*dXact.gt.Xr) then
            NX = NX + 1
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXact
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        end if
        i = i + 1
    end do
    if(i.eq.AtomData%NQt_bb) then         ! last Xbb spectral point
        Xl = SpIntervals%Xint(NX)%Xmin
        Xr = SpIntervals%Xint(NX)%Xmax
        dX = Xbb(i)*DWCoef*SpRange%NDWCntr
        dXact = min(dX,dXmin)
        if(Xbb(i)-1.5d0*dXact.gt.Xr) then
            NX = NX + 1
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXact
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        end if
    end if
    SpIntervals%NX = NX


    ! ========================== 2. correct multigroup intervals (interval length = n * dXmin): ===================
    do i = 1,SpIntervals%NX        ! skip 1st - boundary interval
        Xact  = 0.5d0*(SpIntervals%Xint(i)%Xmax+SpIntervals%Xint(i)%Xmin)
        dX    = Xact*DWCoef*SpRange%NDWCntr/SpRange%NGrCntr
        dXact = min(dX,dXmin)
        dX = SpIntervals%Xint(i)%Xmax - SpIntervals%Xint(i)%Xmin
        j = floor(dX*(1.d0+1.d-3)/dXact)
        eps = 0.5d0*(dX - j*dXact)
        SpIntervals%Xint(i)%Xmin = SpIntervals%Xint(i)%Xmin + eps
        SpIntervals%Xint(i)%Xmax = SpIntervals%Xint(i)%Xmax - eps
    end do
    if(NX.gt.0) then
        if(SpIntervals%Xint(1)%Xmin.lt.Xmin) then
            SpIntervals%Xint(1)%Xmin = Xmin
        end if
        if(SpIntervals%Xint(NX)%Xmax.gt.Xmax) then
            SpIntervals%Xint(NX)%Xmax = Xmax
        end if
    endif
    SpIntervals%NX = NX


    ! =========================== 3. split wide central intervals: ===============================
    ipos = 1
    do i = 1,SpIntervals%NX
        Xact  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)
        dX    = Xact*DWCoef*SpRange%NDWCntr/SpRange%NGrCntr
        dXact = min(dX,dXmin)
        dX = SpIntervals%Xint(ipos)%Xmax - SpIntervals%Xint(ipos)%Xmin
        j = floor(dX*(1.d0+1.d-3)/dXact)
        if (j.gt.1) then
            eps = dX/j
            do k = 1,NX-ipos     ! shift elements to the right (create free space)
                SpIntervals%Xint(NX-k+j)%Xmin = SpIntervals%Xint(NX-k+1)%Xmin
                SpIntervals%Xint(NX-k+j)%Xmax = SpIntervals%Xint(NX-k+1)%Xmax
                SpIntervals%Xint(NX-k+j)%GrType = 1    ! center
            end do
            SpIntervals%Xint(ipos)%Xmax = SpIntervals%Xint(ipos)%Xmin + eps
            SpIntervals%Xint(ipos)%GrType = 1    ! center
            do k = 1,j-1         ! paste elements into free space
                SpIntervals%Xint(ipos+k)%Xmin = SpIntervals%Xint(ipos+k-1)%Xmax
                SpIntervals%Xint(ipos+k)%Xmax = SpIntervals%Xint(ipos+k)%Xmin + eps
                SpIntervals%Xint(ipos+k)%GrType = 1    ! center
            end do
            SpIntervals%Xint(ipos+j-1)%Xmax = SpIntervals%Xint(ipos)%Xmin + dX
            SpIntervals%Xint(ipos+j-1)%GrType = 1    ! center
        else
            j = 1
        end if
        ipos = ipos + j
        NX   = NX + (j-1)
    end do
    SpIntervals%NX = NX

    ! =========================== 4. create wings: ===============================
    NDWWing = (SpRange%NDWLine-SpRange%NDWCntr)/2.d0
    NGrWing = (SpRange%NGrLine-SpRange%NGrCntr)/2
    ipos = 1
    Xl   = SpRange%Xmin
    dXl  = 0.d0
    do i = 1,SpIntervals%NX
        ! .............................. left wing ...............................
        Xr   = SpIntervals%Xint(ipos)%Xmin
        dXr  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)*DWCoef*NDWWing
        ! correction of Xl position
        if(Xr-dXr.lt.Xl+dXl) then
            eps  = 0.5d0*((Xr-dXr) + (Xl+dXl))
            Xact = min(eps,Xl+dXl)
            if(Xact.ge.Xl) then
                Xl = Xact
            end if
        else
            Xl = Xr - dXr
        end if
        ! set dX after Xl correction
        dX = Xr - Xl

        Xact  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)
        dXact = Xact*DWCoef*NDWWing/NGrWing
        call WingMeshLeft(1.d0,Xl,Xr,dXact,dXmax,XWingGrid,NWx)
        j = NWx-1          ! "number of intervals" = "number of points" - 1
        if (j.gt.0) then
            do k = 1,NX-ipos+1     ! shift elements to the right (create free space)
                SpIntervals%Xint(NX-k+j+1)%Xmin = SpIntervals%Xint(NX-k+1)%Xmin
                SpIntervals%Xint(NX-k+j+1)%Xmax = SpIntervals%Xint(NX-k+1)%Xmax
            end do

            do k = 1,j              ! paste elements into free space
                SpIntervals%Xint(ipos+k-1)%Xmin = XWingGrid(k)
                SpIntervals%Xint(ipos+k-1)%Xmax = XWingGrid(k+1)
            end do
        else
            j = 0
        end if
        ipos = ipos + j              ! stay on current group for right wing building
        NX   = NX + j

        ! .............................. right wing .......................................
        Xl   = SpIntervals%Xint(ipos)%Xmax
        dXl  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)*DWCoef*NDWWing
        if(ipos+1.le.NX) then
            Xr  = SpIntervals%Xint(ipos+1)%Xmin
            dXr = 0.5d0*(SpIntervals%Xint(ipos+1)%Xmax+SpIntervals%Xint(ipos+1)%Xmin)*DWCoef*NDWWing
        else
            Xr  = Xmax
            dXr = 0.d0
        end if

        ! correction of Xr position
        if(Xl+dXl.gt.Xr-dXr) then
            eps = 0.5d0*((Xr-dXr) + (Xl+dXl))
            Xact = max(eps,Xr-dXr)
            if(Xact.le.Xr) then
                Xr = Xact
            end if
        else
            Xr = Xl + dXl
        end if
        ! set dX after Xr correction
        dX = Xr - Xl

        Xact  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)
        dXact = Xact*DWCoef*NDWWing/NGrWing
        call WingMeshRight(1.d0,Xl,Xr,dXact,dXmax,XWingGrid,NWx)
        j = NWx-1          ! "number of intervals" = "number of points" - 1
        if (j.gt.0) then
            do k = 1,NX-ipos       ! shift elements to the right (create free space)
                SpIntervals%Xint(NX-k+j+1)%Xmin = SpIntervals%Xint(NX-k+1)%Xmin
                SpIntervals%Xint(NX-k+j+1)%Xmax = SpIntervals%Xint(NX-k+1)%Xmax
            end do

            do k = 1,j              ! paste elements into free space
                SpIntervals%Xint(ipos+k)%Xmin = XWingGrid(k)
                SpIntervals%Xint(ipos+k)%Xmax = XWingGrid(k+1)
            end do
        else
            j = 0
        end if
        ipos = ipos + j
        NX   = NX + j

        if(ipos-1.gt.0) then
            Xl   = SpIntervals%Xint(ipos-1)%Xmax
            dXl  = 0.5d0*(SpIntervals%Xint(ipos-1)%Xmax+SpIntervals%Xint(ipos-1)%Xmin)*DWCoef*NDWWing
        endif

        ipos = ipos + 1       ! "+ 1" to go the next group
    enddo
    SpIntervals%NX = NX


    write(409,*) "# -------------------------------------------- "
    do i = 1,NX
        dXact = SpIntervals%Xint(i)%Xmax - SpIntervals%Xint(i)%Xmin
        write(410,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmin,0.d0
        write(410,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmin,1.d0/dXact
        write(410,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmax,1.d0/dXact
        write(410,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmax,0.d0
        write(410,"(2(1x,1pE12.4e2))")
    enddo


    return
end subroutine


! ---------------------------------------------------------------
subroutine set_SpGrid_BB_Progression(Xbb,SpGridOpts,SpRange,SpIntervals)
    use module_Settings
    use module_Constants
    use module_AtomData
    use module_Spectrum
    implicit none
    ! input
    real(8), dimension(AtomData%NQt_bb) :: Xbb
    type TSpRange
        sequence
        real(r8p)               :: Xmin,Xmax
        character(16)           :: GridType
        character(16)           :: ProfType
        integer(i4p)            :: NGrLine  ! number of groups for detailed line
        integer(i4p)            :: NGrCntr  ! number of groups for detailed line center
        real(r8p)               :: NDWLine  ! number of Dopler widths per detailed line
        real(r8p)               :: NDWCntr  ! number of Dopler widths per detailed line center
    end type TSpRange
    type TSpGridOpts
        sequence
        character(4)            :: SpUnits
        integer(i4p)            :: NSpRange             ! number of spectral intervals
        real(r8p)               :: Xmin,Xmax
        real(r8p)               :: dXmin,dXmax
        real(r8p)               :: SpTemp               ! characteristic temperature for line width calculation
    end type TSpGridOpts

    type(TSpGridOpts)             :: SpGridOpts
    type(TSpRange)                :: SpRange         ! spectral intervals
    ! output
    ! ........................................................
    type(TGridIntervals) :: SpIntervals

    ! local
    real(r8p) :: Xact,Xmin,Xmax
    real(r8p) :: dXmin,dXmax,dXact,dX
    real(r8p) :: Xl,Xr,dXl,dXr
    real(r8p) :: LineWidth,CntrWidth
    real(r8p) :: eps
    integer   :: i,j,k,ipos,NX
    real(r8p) :: NDWWing
    integer   :: NGrLine,NGrCntr,NGrWing
    real(r8p) :: StpInc,DWCoef,TeV,Mw           ! Dopler width coefficient, Dopler width in eV, molar weight
    real(8),dimension(100) :: XWingGrid
    integer                :: NWx


    Xmin = SpRange%Xmin
    Xmax = SpRange%Xmax
    dXmin = SpGridOpts%dXmin
    dXmax = SpGridOpts%dXmax
    Mw     = 1.0d0
    TeV    = SpGridOpts%SpTemp
    DWCoef = sqrt(2*TeV*eV_erg*N_Av/Mw)/v_c


    SpIntervals%Xint(:)%Xmin  = 0.d0
    SpIntervals%Xint(:)%Xmax  = 0.d0
    SpIntervals%NX            = 0


    ! ======================== 1.a create 1st ( left boundary) interval: =======================
    i = 1
    do while(Xbb(i).le.Xmin .and. i.lt.AtomData%NQt_bb)
        i = i + 1
    end do
    NX = 0
    if(Xbb(i).ge.Xmin .and. Xbb(i).lt.Xmax) then
        NX = NX + 1
        dX = Xbb(i)*DWCoef*SpRange%NDWCntr
        dXact = min(dX,dXmin)
        if(Xbb(i)-0.5d0*dXact.gt.Xmin) then
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmin = Xmin
        end if
        if(Xbb(i)+0.5d0*dXact.lt.Xmax) then
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmax = Xmax
        end if
    end if
    if(i.eq.AtomData%NQt_bb) then
        return
    end if

    ! 1.b set other intervals
    i = i + 1
    do while(Xbb(i).le.Xmax .and. i.lt.AtomData%NQt_bb)
        Xl = SpIntervals%Xint(NX)%Xmin
        Xr = SpIntervals%Xint(NX)%Xmax
        dX = Xbb(i)*DWCoef*SpRange%NDWCntr
        dXact = min(dX,dXmin)
        if(Xbb(i)-1.5d0*dXact.gt.Xr) then
            NX = NX + 1
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXact
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        end if
        i = i + 1
    end do
    if(i.eq.AtomData%NQt_bb) then         ! last Xbb spectral point
        Xl = SpIntervals%Xint(NX)%Xmin
        Xr = SpIntervals%Xint(NX)%Xmax
        dX = Xbb(i)*DWCoef*SpRange%NDWCntr
        dXact = min(dX,dXmin)
        if(Xbb(i)-1.5d0*dXact.gt.Xr) then
            NX = NX + 1
            SpIntervals%Xint(NX)%Xmin = Xbb(i)-0.5d0*dXact
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        else
            SpIntervals%Xint(NX)%Xmax = Xbb(i)+0.5d0*dXact
        end if
    end if
    SpIntervals%NX = NX


    ! ========================== 2. correct multigroup intervals (interval length = n * dXmin): ===================
    do i = 1,SpIntervals%NX        ! skip 1st - boundary interval
        Xact  = 0.5d0*(SpIntervals%Xint(i)%Xmax+SpIntervals%Xint(i)%Xmin)
        dX    = Xact*DWCoef*SpRange%NDWCntr/SpRange%NGrCntr
        dXact = min(dX,dXmin)
        dX = SpIntervals%Xint(i)%Xmax - SpIntervals%Xint(i)%Xmin
        j = floor(dX*(1.d0+1.d-3)/dXact)
        eps = 0.5d0*(dX - j*dXact)
        SpIntervals%Xint(i)%Xmin = SpIntervals%Xint(i)%Xmin + eps
        SpIntervals%Xint(i)%Xmax = SpIntervals%Xint(i)%Xmax - eps
    end do
    if(NX.gt.0) then
        if(SpIntervals%Xint(1)%Xmin.lt.Xmin) then
            SpIntervals%Xint(1)%Xmin = Xmin
        end if
        if(SpIntervals%Xint(NX)%Xmax.gt.Xmax) then
            SpIntervals%Xint(NX)%Xmax = Xmax
        end if
    endif
    SpIntervals%NX = NX


    ! =========================== 3. split wide central intervals: ===============================
    ipos = 1
    do i = 1,SpIntervals%NX
        Xact  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)
        dX    = Xact*DWCoef*SpRange%NDWCntr/SpRange%NGrCntr
        dXact = min(dX,dXmin)
        dX = SpIntervals%Xint(ipos)%Xmax - SpIntervals%Xint(ipos)%Xmin
        j = floor(dX*(1.d0+1.d-3)/dXact)
        if (j.gt.1) then
            eps = dX/j
            do k = 1,NX-ipos     ! shift elements to the right (create free space)
                SpIntervals%Xint(NX-k+j)%Xmin = SpIntervals%Xint(NX-k+1)%Xmin
                SpIntervals%Xint(NX-k+j)%Xmax = SpIntervals%Xint(NX-k+1)%Xmax
                SpIntervals%Xint(NX-k+j)%GrType = 1    ! center
            end do
            SpIntervals%Xint(ipos)%Xmax = SpIntervals%Xint(ipos)%Xmin + eps
            SpIntervals%Xint(ipos)%GrType = 1    ! center
            do k = 1,j-1         ! paste elements into free space
                SpIntervals%Xint(ipos+k)%Xmin = SpIntervals%Xint(ipos+k-1)%Xmax
                SpIntervals%Xint(ipos+k)%Xmax = SpIntervals%Xint(ipos+k)%Xmin + eps
                SpIntervals%Xint(ipos+k)%GrType = 1    ! center
            end do
            SpIntervals%Xint(ipos+j-1)%Xmax = SpIntervals%Xint(ipos)%Xmin + dX
            SpIntervals%Xint(ipos+j-1)%GrType = 1    ! center
        else
            j = 1
        end if
        ipos = ipos + j
        NX   = NX + (j-1)
    end do
    SpIntervals%NX = NX

    ! =========================== 4. create wings: ===============================
    NGrLine = SpRange%NGrLine
    NGrCntr = SpRange%NGrCntr
    NDWWing = (SpRange%NDWLine-SpRange%NDWCntr)/2.d0
    NGrWing = (SpRange%NGrLine-SpRange%NGrCntr)/2
    ipos = 1
    Xl   = SpRange%Xmin
    dXl  = 0.d0
    do i = 1,SpIntervals%NX
        ! .............................. left wing ...............................
        Xr   = SpIntervals%Xint(ipos)%Xmin
        dXr  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)*DWCoef*NDWWing
        ! correction of Xl position
        if(Xr-dXr.lt.Xl+dXl) then
            eps  = 0.5d0*((Xr-dXr) + (Xl+dXl))
            Xact = min(eps,Xl+dXl)
            if(Xact.ge.Xl) then
                Xl = Xact
            end if
        else
            Xl = Xr - dXr
        end if
        ! set dX after Xl correction
        dX = Xr - Xl

        Xact  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)
        LineWidth = Xact*DWCoef*SpRange%NDWLine
        CntrWidth = Xact*DWCoef*SpRange%NDWCntr
        call GetStepIncrement(LineWidth,CntrWidth,NGrLine,NGrCntr,StpInc,k)
        dXact = (SpIntervals%Xint(ipos)%Xmax-SpIntervals%Xint(ipos)%Xmin)
        call WingMeshLeft(StpInc,Xl,Xr,dXact,dXmax,XWingGrid,NWx)

        j = NWx-1          ! "number of intervals" = "number of points" - 1
        if (j.gt.0) then
            do k = 1,NX-ipos+1     ! shift elements to the right (create free space)
                SpIntervals%Xint(NX-k+j+1)%Xmin = SpIntervals%Xint(NX-k+1)%Xmin
                SpIntervals%Xint(NX-k+j+1)%Xmax = SpIntervals%Xint(NX-k+1)%Xmax
            end do

            do k = 1,j              ! paste elements into free space
                SpIntervals%Xint(ipos+k-1)%Xmin = XWingGrid(k)
                SpIntervals%Xint(ipos+k-1)%Xmax = XWingGrid(k+1)
            end do
        else
            j = 0
        end if
        ipos = ipos + j              ! stay on current group for right wing building
        NX   = NX + j

        ! .............................. right wing .......................................
        Xl   = SpIntervals%Xint(ipos)%Xmax
        dXl  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)*DWCoef*NDWWing
        if(ipos+1.le.NX) then
            Xr  = SpIntervals%Xint(ipos+1)%Xmin
            dXr = 0.5d0*(SpIntervals%Xint(ipos+1)%Xmax+SpIntervals%Xint(ipos+1)%Xmin)*DWCoef*NDWWing
        else
            Xr  = Xmax
            dXr = 0.d0
        end if

        ! correction of Xr position
        if(Xl+dXl.gt.Xr-dXr) then
            eps = 0.5d0*((Xr-dXr) + (Xl+dXl))
            Xact = max(eps,Xr-dXr)
            if(Xact.le.Xr) then
                Xr = Xact
            end if
        else
            Xr = Xl + dXl
        end if
        ! set dX after Xr correction
        dX = Xr - Xl

        Xact  = 0.5d0*(SpIntervals%Xint(ipos)%Xmax+SpIntervals%Xint(ipos)%Xmin)
        LineWidth = Xact*DWCoef*SpRange%NDWLine
        CntrWidth = Xact*DWCoef*SpRange%NDWCntr
        call GetStepIncrement(LineWidth,CntrWidth,NGrLine,NGrCntr,StpInc,k)
        dXact = (SpIntervals%Xint(ipos)%Xmax-SpIntervals%Xint(ipos)%Xmin)
        call WingMeshRight(StpInc,Xl,Xr,dXact,dXmax,XWingGrid,NWx)
        j = NWx-1          ! "number of intervals" = "number of points" - 1
        if (j.gt.0) then
            do k = 1,NX-ipos       ! shift elements to the right (create free space)
                SpIntervals%Xint(NX-k+j+1)%Xmin = SpIntervals%Xint(NX-k+1)%Xmin
                SpIntervals%Xint(NX-k+j+1)%Xmax = SpIntervals%Xint(NX-k+1)%Xmax
            end do

            do k = 1,j              ! paste elements into free space
                SpIntervals%Xint(ipos+k)%Xmin = XWingGrid(k)
                SpIntervals%Xint(ipos+k)%Xmax = XWingGrid(k+1)
            end do
        else
            j = 0
        end if
        ipos = ipos + j
        NX   = NX + j

        if(ipos-1.gt.0) then
            Xl   = SpIntervals%Xint(ipos-1)%Xmax
            dXl  = 0.5d0*(SpIntervals%Xint(ipos-1)%Xmax+SpIntervals%Xint(ipos-1)%Xmin)*DWCoef*NDWWing
        endif

        ipos = ipos + 1       ! "+ 1" to go the next group
    enddo
    SpIntervals%NX = NX


    write(409,*) "# -------------------------------------------- "
    do i = 1,NX
        dXact = SpIntervals%Xint(i)%Xmax - SpIntervals%Xint(i)%Xmin
        write(411,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmin,0.d0
        write(411,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmin,1.d0/dXact
        write(411,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmax,1.d0/dXact
        write(411,"(2(1x,1pE12.4e2))") SpIntervals%Xint(i)%Xmax,0.d0
        write(411,"(2(1x,1pE12.4e2))")
    enddo


    return
end subroutine



! ---------------------------------------------------------------
subroutine WingMeshLeft(StpInc,Xl,Xr,dXmin,dXmax,XWingGrid,NWx)
    implicit none
    ! input:
    real(8) :: StpInc
    real(8) :: Xl,Xr
    real(8) :: dXmin,dXmax
!    integer :: NWGr
    ! output:
    real(8),dimension(100) :: XWingGrid
    integer                :: NWx
    ! local:
    real(8) :: dX,Xc,Xstop,dXstop
    real(8) :: rcor,eps           ! correctors
    integer :: i

    XWingGrid = 0.d0
    NWx    = 0

    ! 1. find number if groups in the interval [Xl,Xr]
    i = 1
    dX = dXmin*StpInc**(i-1)*(1.d0-1.d-7)        ! *(1.d0-1.d-5) - corrector, for not
    if(dX.gt.dXmax) then
        dX = dXmax
    end if
    Xc = Xr - dX
    Xstop  = Xc
    dXstop = dX
    NWx = 0
    do while (Xc-0.5d0*dX.gt.Xl)
        i = i + 1
        Xstop  = Xc
        dXstop = dX
        dX = dXmin*StpInc**(i-1)*(1.d0-1.d-7)
        if(dX.gt.dXmax) then
            dX = dXmax
        end if
        Xc = Xc - dX
    end do
!    eps = abs(Xl-Xc)/dX
    eps = (Xl-Xc)/dX               ! "<!>" if abs(Xc-Xr)/dX then wouldn't go to bound

    if(eps.lt.1.d-2) then
        Xstop  = Xc
        dXstop = dX
        rcor = (Xr-Xl)/(Xr-Xstop)                 ! distance to left bound -> 0
!        write(*,*) "case: 1"
    else
        rcor = (Xr-Xl)/(Xr-(Xstop-0.5d0*dXstop))  ! distance to left bound = 0.5*dXstop
!        write(*,*) "case: 2"
    end if

    if(i-1.le.0) then
        XWingGrid(:) = 0.d0
        NWx = 0
        return
    end if



    ! grid after correction:
    i = 0
    Xc = Xr
    XWingGrid(:) = 0.d0
    do while (Xc.ge.Xl-1.d-3*dX)
        i = i + 1
        XWingGrid(i) = Xc
        dX = dXmin*StpInc**(i-1)*(1.d0-1.d-7) * rcor
        if(dX.gt.dXmax) then
            dX = dXmax * rcor
        end if
        Xc = Xc - dX
    end do
    NWx = i                    ! NNodes
    if(XWingGrid(NWx).lt.Xl) then
        XWingGrid(NWx) = Xl
    end if

    ! reverse XWingGrid
    do i = 1,NWx/2
        Xc = XWingGrid(i)
        XWingGrid(i) = XWingGrid(NWx+1-i)
        XWingGrid(NWx+1-i) = Xc
    end do


    return
end subroutine

! ---------------------------------------------------------------
subroutine WingMeshRight(StpInc,Xl,Xr,dXmin,dXmax,XWingGrid,NWx)
    implicit none
    ! input:
    real(8) :: StpInc
    real(8) :: Xl,Xr
    real(8) :: dXmin,dXmax
    ! output:
    real(8),dimension(100) :: XWingGrid
    integer                :: NWx
    ! local:
    real(8) :: dX,Xc,Xstop,dXstop
    real(8) :: rcor,eps             ! correctors
    integer :: i

    XWingGrid = 0.d0
    NWx    = 0

    ! 1. find number if groups in the interval [Xl,Xr]
    i = 1
    dX = dXmin*StpInc**(i-1)*(1.d0-1.d-7)        ! *(1.d0-1.d-5) - corrector, for not
    if(dX.gt.dXmax) then
        dX = dXmax
    end if
    Xc = Xl + dX
    Xstop  = Xc
    dXstop = dX
    NWx = 0
    do while (Xc+0.5d0*dX.lt.Xr)
        i = i + 1
        Xstop  = Xc
        dXstop = dX
        dX = dXmin*StpInc**(i-1)*(1.d0-1.d-7)
        if(dX.gt.dXmax) then
            dX = dXmax
        end if
        Xc = Xc + dX
    end do
!    eps = abs(Xc-Xr)/dX
    eps = (Xc-Xr)/dX         ! "<!>" if abs(Xc-Xr)/dX then wouldn't go to bound

    if(eps.lt.1.d-1) then
        Xstop  = Xc
        dXstop = dX
        rcor = (Xr-Xl)/(Xstop-Xl)                 ! distance to left bound -> 0
    else
        rcor = (Xr-Xl)/((Xstop+0.5d0*dXstop)-Xl)  ! distance to left bound = 0.5*dXstop
    end if

    if(i-1.le.0) then
        XWingGrid(:) = 0.d0
        NWx = 0
        return
    end if

    ! grid after correction:
    i = 0
    Xc = Xl
    XWingGrid(:) = 0.d0
    do while (Xc.lt.Xr+1.d-2*dX)
        i = i + 1
        XWingGrid(i) = Xc
        dX = dXmin*StpInc**(i-1)*(1.d0-1.d-7) * rcor
        if(dX.gt.dXmax) then
            dX = dXmax * rcor
        end if
        Xc = Xc + dX
    end do
    NWx = i                    ! NNodes
    if(XWingGrid(NWx).gt.Xr) then
        XWingGrid(NWx) = Xr
    end if

    return
end subroutine


! ---------------------------------------------------------------
subroutine GetStepIncrement(LineWidth,CntrWidth,NLGr,NCGr,StpInc,niter)
    implicit none
    ! input:
    real(8) :: LineWidth
    real(8) :: CntrWidth
    integer :: NLGr,NCGr
    ! output:
    real(8) :: StpInc
    ! local:
    integer :: NWGr
    real(8) :: WingWidth
    real(8) :: dX1
    real(8) :: StpInc0
    real(8) :: eps,rtol
    integer :: niter

    WingWidth = 0.5d0*(LineWidth-CntrWidth)
    NWGr      = (NLGr-NCGr)/2
    StpInc0 = 2.0
    dX1 = CntrWidth
    eps = 1.d0/(NWGr+1)
    rtol = 1.0d0
    niter = 0
    do while (rtol.gt.1.0d-2)
        StpInc  = (1.d0 + WingWidth/dX1*(StpInc0-1.d0))**eps
        rtol = abs(StpInc-StpInc0)/StpInc
        StpInc0 = StpInc
        niter = niter + 1
        if(niter.gt.1000) then
            write(*,*) "error: GetStepIncrement - can't calculate step increment, iterations > ",niter
            stop
        end if
    end do

    return
end subroutine

