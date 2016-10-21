!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : This module converts python definitions of model's variables to FORTRAN95 derived data type.
module statevars_mod
    implicit none
    type StateVariables_type
        real*8, pointer :: y_CamKII(:)
        real*8, pointer :: PP1
        real*8, pointer :: V
        real*8, pointer :: Ca_cyt
        real*8, pointer :: o_CB1R
        real*8, pointer :: o_NMDA
        real*8, pointer :: Ca_ER
        real*8, pointer :: m_caL13
        real*8, pointer :: I1P
        real*8, pointer :: h_caL13
        real*8, pointer :: AEA
        real*8, pointer :: IP3
        real*8, pointer :: h_CICR
        real*8, pointer :: DAG
        real*8, pointer :: d_CB1R
        real*8, pointer :: twoAG
        real*8, pointer :: d_AMPA
        real*8, pointer :: fpre
        real*8, pointer :: DAGLP
        real*8, pointer :: o_AMPA
    end type StateVariables_type
    contains
    subroutine SVs_get(SVs, y)
        type(StateVariables_type) :: SVs
        real*8, target :: y(32)
        SVs%y_CamKII => y(1:13)
        SVs%PP1 => y(14)
        SVs%V => y(15)
        SVs%Ca_cyt => y(16)
        SVs%o_CB1R => y(17)
        SVs%o_NMDA => y(18)
        SVs%Ca_ER => y(19)
        SVs%m_caL13 => y(20)
        SVs%I1P => y(21)
        SVs%h_caL13 => y(22)
        SVs%AEA => y(23)
        SVs%IP3 => y(24)
        SVs%h_CICR => y(25)
        SVs%DAG => y(26)
        SVs%d_CB1R => y(27)
        SVs%twoAG => y(28)
        SVs%d_AMPA => y(29)
        SVs%fpre => y(30)
        SVs%DAGLP => y(31)
        SVs%o_AMPA => y(32)
    end subroutine SVs_get
end module statevars_mod
