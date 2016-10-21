!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : This module converts python definitions of parameters to FORTRAN95 derived data type.
module pars_mod
    implicit none
    type stimulation_type
        real*8, pointer :: regular_on
        real*8, pointer :: tsdt
        real*8, pointer :: tpost
        real*8, pointer :: post_on
        real*8, pointer :: dt_stim
        real*8, pointer :: num_stim
        real*8, pointer :: pre_on
        real*8, pointer :: Freq
        real*8, pointer :: tables_step
    end type stimulation_type
    type Glu_release_type
        real*8, pointer :: steadyrise_on
        real*8, pointer :: Glumax
        real*8, pointer :: BaseLevel
        real*8, pointer :: tauGlu
    end type Glu_release_type
    type action_type
        real*8, pointer :: APdur
        real*8, pointer :: action_as_VDCC
        real*8, pointer :: DPmax
        real*8, pointer :: APmax
        real*8, pointer :: tausbAP
    end type action_type
    type integration_type
        real*8, pointer :: t_step
        real*8, pointer :: t_start
        real*8, pointer :: t_end
        real*8, pointer :: ATOL
        real*8, pointer :: RTOL
        real*8, pointer :: MXSTEP
        real*8, pointer :: HMAX
    end type integration_type
    type common_type
        real*8, pointer :: R
        real*8, pointer :: RT
        real*8, pointer :: T
        real*8, pointer :: F
        real*8, pointer :: zS
        real*8, pointer :: Ca_out
    end type common_type
    type CB1R_type
        real*8, pointer :: on
        real*8, pointer :: Alpha
        real*8, pointer :: Beta
        real*8, pointer :: Gamma
        real*8, pointer :: Epsilon
    end type CB1R_type
    type mem_type
        real*8, pointer :: EL
        real*8, pointer :: gL
        real*8, pointer :: Cm
    end type mem_type
    type NMDA_type
        real*8, pointer :: on
        real*8, pointer :: gNMDA
        real*8, pointer :: p_ca
        real*8, pointer :: Mg
        real*8, pointer :: Alpha
        real*8, pointer :: Beta
    end type NMDA_type
    type AMPA_type
        real*8, pointer :: on
        real*8, pointer :: gAMPA
        real*8, pointer :: Epsilon
        real*8, pointer :: Beta
        real*8, pointer :: Alpha
        real*8, pointer :: Gamma
    end type AMPA_type
    type DA_type
        real*8, pointer :: gamma1DA
        real*8, pointer :: gamma2DA
        real*8, pointer :: DA
    end type DA_type
    type post_CaMKII_plast_type
        real*8, pointer :: kpkaI1
        real*8, pointer :: PP10
        real*8, pointer :: KM
        real*8, pointer :: K5
        real*8, pointer :: k12
        real*8, pointer :: k11
        real*8, pointer :: kcan0I1
        real*8, pointer :: km11
        real*8, pointer :: CaMKT
        real*8, pointer :: k7
        real*8, pointer :: k6
        real*8, pointer :: ncanI1
        real*8, pointer :: I10
        real*8, pointer :: CaMT
        real*8, pointer :: KdpkaI1
        real*8, pointer :: kcanI1
        real*8, pointer :: KdcanI1
        real*8, pointer :: kpka0I1
        real*8, pointer :: Ka3
        real*8, pointer :: Ka2
        real*8, pointer :: Ka1
        real*8, pointer :: Ka4
        real*8, pointer :: npkaI1
    end type post_CaMKII_plast_type
    type ECb_smooth_type
        real*8, pointer :: on
        real*8, pointer :: K
        real*8, pointer :: n
        real*8, pointer :: kw
        real*8, pointer :: tau
        real*8, pointer :: kadd
        real*8, pointer :: kn
    end type ECb_smooth_type
    type ECb_type
        real*8, pointer :: on
        real*8, pointer :: CB1R_on
        real*8, pointer :: kCB1R
        real*8, pointer :: alphaAEACB1
        real*8, pointer :: P1
        real*8, pointer :: P2
        real*8, pointer :: P3
        real*8, pointer :: P4
        real*8, pointer :: vATAEA
        real*8, pointer :: LTDstart
        real*8, pointer :: LTDstop
        real*8, pointer :: LTDMax
        real*8, pointer :: LTPstart
        real*8, pointer :: LTPMax
        real*8, pointer :: KFAAH
        real*8, pointer :: vFAAH
    end type ECb_type
    type KandP_on_DAGLP_type
        real*8, pointer :: nK
        real*8, pointer :: rP
        real*8, pointer :: rK
    end type KandP_on_DAGLP_type
    type DGLandDAG_type
        real*8, pointer :: KDGL
        real*8, pointer :: kMAGL
        real*8, pointer :: rDGL
        real*8, pointer :: kDAGK
    end type DGLandDAG_type
    type IP3_type
        real*8, pointer :: kappad
        real*8, pointer :: kdelta
        real*8, pointer :: r5p
        real*8, pointer :: v3k
        real*8, pointer :: kd
        real*8, pointer :: vdelta
        real*8, pointer :: k3
        real*8, pointer :: kr
        real*8, pointer :: kp
        real*8, pointer :: vbeta
        real*8, pointer :: kpi
        real*8, pointer :: n3
    end type IP3_type
    type CICR_type
        real*8, pointer :: a2
        real*8, pointer :: ver
        real*8, pointer :: d3
        real*8, pointer :: rc
        real*8, pointer :: rl
        real*8, pointer :: rhoER
        real*8, pointer :: d5
        real*8, pointer :: d2
        real*8, pointer :: ker
        real*8, pointer :: d1
    end type CICR_type
    type CaBuff_type
        real*8, pointer :: BT
        real*8, pointer :: Cab
        real*8, pointer :: tauCab
        real*8, pointer :: KdB
    end type CaBuff_type
    type I_to_Ca_flux_type
        real*8, pointer :: NMDA
        real*8, pointer :: VDCC
        real*8, pointer :: TRPV1
    end type I_to_Ca_flux_type
    type TRPV1_type
        real*8, pointer :: on
        real*8, pointer :: gTRPV1
        real*8, pointer :: p_ca
        real*8, pointer :: C
        real*8, pointer :: D
        real*8, pointer :: DH
        real*8, pointer :: KD
        real*8, pointer :: J0
        real*8, pointer :: L
        real*8, pointer :: P
        real*8, pointer :: z
        real*8, pointer :: DS
        real*8, pointer :: K
    end type TRPV1_type
    type caL13_type
        real*8, pointer :: on
        real*8, pointer :: pcaLbar
        real*8, pointer :: mslope
        real*8, pointer :: hshift
        real*8, pointer :: vm
        real*8, pointer :: mshift
        real*8, pointer :: hslope
        real*8, pointer :: kpr
        real*8, pointer :: c
        real*8, pointer :: k
        real*8, pointer :: hvhalf
        real*8, pointer :: mvhalf
        real*8, pointer :: cpr
        real*8, pointer :: htau
        real*8, pointer :: hqfact
        real*8, pointer :: qfact
    end type caL13_type
    type pars_type
        type(caL13_type) :: caL13
        type(TRPV1_type) :: TRPV1
        type(I_to_Ca_flux_type) :: I_to_Ca_flux
        type(CaBuff_type) :: CaBuff
        type(CICR_type) :: CICR
        type(IP3_type) :: IP3
        type(DGLandDAG_type) :: DGLandDAG
        type(KandP_on_DAGLP_type) :: KandP_on_DAGLP
        type(ECb_type) :: ECb
        type(ECb_smooth_type) :: ECb_smooth
        type(post_CaMKII_plast_type) :: post_CaMKII_plast
        type(DA_type) :: DA
        type(AMPA_type) :: AMPA
        type(NMDA_type) :: NMDA
        type(mem_type) :: mem
        type(CB1R_type) :: CB1R
        type(common_type) :: common
        type(integration_type) :: integration
        type(action_type) :: action
        type(Glu_release_type) :: Glu_release
        type(stimulation_type) :: stimulation
    end type pars_type
    contains
    subroutine parameters_get(pars, lv)
        type(pars_type) :: pars
        real*8, target :: lv(165)
        pars%caL13%on => lv(1)
        pars%caL13%pcaLbar => lv(2)
        pars%caL13%mslope => lv(3)
        pars%caL13%hshift => lv(4)
        pars%caL13%vm => lv(5)
        pars%caL13%mshift => lv(6)
        pars%caL13%hslope => lv(7)
        pars%caL13%kpr => lv(8)
        pars%caL13%c => lv(9)
        pars%caL13%k => lv(10)
        pars%caL13%hvhalf => lv(11)
        pars%caL13%mvhalf => lv(12)
        pars%caL13%cpr => lv(13)
        pars%caL13%htau => lv(14)
        pars%caL13%hqfact => lv(15)
        pars%caL13%qfact => lv(16)
        pars%TRPV1%on => lv(17)
        pars%TRPV1%gTRPV1 => lv(18)
        pars%TRPV1%p_ca => lv(19)
        pars%TRPV1%C => lv(20)
        pars%TRPV1%D => lv(21)
        pars%TRPV1%DH => lv(22)
        pars%TRPV1%KD => lv(23)
        pars%TRPV1%J0 => lv(24)
        pars%TRPV1%L => lv(25)
        pars%TRPV1%P => lv(26)
        pars%TRPV1%z => lv(27)
        pars%TRPV1%DS => lv(28)
        pars%TRPV1%K => lv(29)
        pars%I_to_Ca_flux%NMDA => lv(30)
        pars%I_to_Ca_flux%VDCC => lv(31)
        pars%I_to_Ca_flux%TRPV1 => lv(32)
        pars%CaBuff%BT => lv(33)
        pars%CaBuff%Cab => lv(34)
        pars%CaBuff%tauCab => lv(35)
        pars%CaBuff%KdB => lv(36)
        pars%CICR%a2 => lv(37)
        pars%CICR%ver => lv(38)
        pars%CICR%d3 => lv(39)
        pars%CICR%rc => lv(40)
        pars%CICR%rl => lv(41)
        pars%CICR%rhoER => lv(42)
        pars%CICR%d5 => lv(43)
        pars%CICR%d2 => lv(44)
        pars%CICR%ker => lv(45)
        pars%CICR%d1 => lv(46)
        pars%IP3%kappad => lv(47)
        pars%IP3%kdelta => lv(48)
        pars%IP3%r5p => lv(49)
        pars%IP3%v3k => lv(50)
        pars%IP3%kd => lv(51)
        pars%IP3%vdelta => lv(52)
        pars%IP3%k3 => lv(53)
        pars%IP3%kr => lv(54)
        pars%IP3%kp => lv(55)
        pars%IP3%vbeta => lv(56)
        pars%IP3%kpi => lv(57)
        pars%IP3%n3 => lv(58)
        pars%DGLandDAG%KDGL => lv(59)
        pars%DGLandDAG%kMAGL => lv(60)
        pars%DGLandDAG%rDGL => lv(61)
        pars%DGLandDAG%kDAGK => lv(62)
        pars%KandP_on_DAGLP%nK => lv(63)
        pars%KandP_on_DAGLP%rP => lv(64)
        pars%KandP_on_DAGLP%rK => lv(65)
        pars%ECb%on => lv(66)
        pars%ECb%CB1R_on => lv(67)
        pars%ECb%kCB1R => lv(68)
        pars%ECb%alphaAEACB1 => lv(69)
        pars%ECb%P1 => lv(70)
        pars%ECb%P2 => lv(71)
        pars%ECb%P3 => lv(72)
        pars%ECb%P4 => lv(73)
        pars%ECb%vATAEA => lv(74)
        pars%ECb%LTDstart => lv(75)
        pars%ECb%LTDstop => lv(76)
        pars%ECb%LTDMax => lv(77)
        pars%ECb%LTPstart => lv(78)
        pars%ECb%LTPMax => lv(79)
        pars%ECb%KFAAH => lv(80)
        pars%ECb%vFAAH => lv(81)
        pars%ECb_smooth%on => lv(82)
        pars%ECb_smooth%K => lv(83)
        pars%ECb_smooth%n => lv(84)
        pars%ECb_smooth%kw => lv(85)
        pars%ECb_smooth%tau => lv(86)
        pars%ECb_smooth%kadd => lv(87)
        pars%ECb_smooth%kn => lv(88)
        pars%post_CaMKII_plast%kpkaI1 => lv(89)
        pars%post_CaMKII_plast%PP10 => lv(90)
        pars%post_CaMKII_plast%KM => lv(91)
        pars%post_CaMKII_plast%K5 => lv(92)
        pars%post_CaMKII_plast%k12 => lv(93)
        pars%post_CaMKII_plast%k11 => lv(94)
        pars%post_CaMKII_plast%kcan0I1 => lv(95)
        pars%post_CaMKII_plast%km11 => lv(96)
        pars%post_CaMKII_plast%CaMKT => lv(97)
        pars%post_CaMKII_plast%k7 => lv(98)
        pars%post_CaMKII_plast%k6 => lv(99)
        pars%post_CaMKII_plast%ncanI1 => lv(100)
        pars%post_CaMKII_plast%I10 => lv(101)
        pars%post_CaMKII_plast%CaMT => lv(102)
        pars%post_CaMKII_plast%KdpkaI1 => lv(103)
        pars%post_CaMKII_plast%kcanI1 => lv(104)
        pars%post_CaMKII_plast%KdcanI1 => lv(105)
        pars%post_CaMKII_plast%kpka0I1 => lv(106)
        pars%post_CaMKII_plast%Ka3 => lv(107)
        pars%post_CaMKII_plast%Ka2 => lv(108)
        pars%post_CaMKII_plast%Ka1 => lv(109)
        pars%post_CaMKII_plast%Ka4 => lv(110)
        pars%post_CaMKII_plast%npkaI1 => lv(111)
        pars%DA%gamma1DA => lv(112)
        pars%DA%gamma2DA => lv(113)
        pars%DA%DA => lv(114)
        pars%AMPA%on => lv(115)
        pars%AMPA%gAMPA => lv(116)
        pars%AMPA%Epsilon => lv(117)
        pars%AMPA%Beta => lv(118)
        pars%AMPA%Alpha => lv(119)
        pars%AMPA%Gamma => lv(120)
        pars%NMDA%on => lv(121)
        pars%NMDA%gNMDA => lv(122)
        pars%NMDA%p_ca => lv(123)
        pars%NMDA%Mg => lv(124)
        pars%NMDA%Alpha => lv(125)
        pars%NMDA%Beta => lv(126)
        pars%mem%EL => lv(127)
        pars%mem%gL => lv(128)
        pars%mem%Cm => lv(129)
        pars%CB1R%on => lv(130)
        pars%CB1R%Alpha => lv(131)
        pars%CB1R%Beta => lv(132)
        pars%CB1R%Gamma => lv(133)
        pars%CB1R%Epsilon => lv(134)
        pars%common%R => lv(135)
        pars%common%RT => lv(136)
        pars%common%T => lv(137)
        pars%common%F => lv(138)
        pars%common%zS => lv(139)
        pars%common%Ca_out => lv(140)
        pars%integration%t_step => lv(141)
        pars%integration%t_start => lv(142)
        pars%integration%t_end => lv(143)
        pars%integration%ATOL => lv(144)
        pars%integration%RTOL => lv(145)
        pars%integration%MXSTEP => lv(146)
        pars%integration%HMAX => lv(147)
        pars%action%APdur => lv(148)
        pars%action%action_as_VDCC => lv(149)
        pars%action%DPmax => lv(150)
        pars%action%APmax => lv(151)
        pars%action%tausbAP => lv(152)
        pars%Glu_release%steadyrise_on => lv(153)
        pars%Glu_release%Glumax => lv(154)
        pars%Glu_release%BaseLevel => lv(155)
        pars%Glu_release%tauGlu => lv(156)
        pars%stimulation%regular_on => lv(157)
        pars%stimulation%tsdt => lv(158)
        pars%stimulation%tpost => lv(159)
        pars%stimulation%post_on => lv(160)
        pars%stimulation%dt_stim => lv(161)
        pars%stimulation%num_stim => lv(162)
        pars%stimulation%pre_on => lv(163)
        pars%stimulation%Freq => lv(164)
        pars%stimulation%tables_step => lv(165)
    end subroutine parameters_get
end module pars_mod
