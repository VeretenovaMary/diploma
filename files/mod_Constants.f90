MODULE module_Constants
! <!>  - mark places for check


    ! phys constants:
    real(8), parameter :: pi       = 3.14159265d+00   ! pi
    real(8), parameter :: v_c      = 2.99792458d+10   ! cm/s, speed of light
    real(8), parameter :: h_Pl     = 6.62607004d-27   ! erg*s, Planck constant (without line)
    real(8), parameter :: k_Bol    = 1.38064852d-16   ! erg/K, Boltzman constant
    real(8), parameter :: N_Av     = 6.02214086d+23   ! 1/mole, Avogadro constant
    real(8), parameter :: eV_erg   = 1.60217662d-12   ! erg/eV, convertation: 1eV = eV_erg*1erg
    real(8), parameter :: R_ugc    = k_Bol*N_Av       ! erg/K/mole, universal gas constant
    real(8), parameter :: m_e      = 9.10938356d-28   ! g, electron mass
    real(8), parameter :: alpha    = 1/137.036d0      ! g, electron mass
    real(8), parameter :: q_e      = 4.8065d-10       ! CGSE, electron charge
    real(8), parameter :: a_0      = 0.52917725d-8    ! cm, Bohr radius


END MODULE module_Constants


