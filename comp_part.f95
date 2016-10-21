!       -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : RHS
module comp_part

    use pars_mod
    use statevars_mod
    use general_math

    use caL13
    use TRPV1

    use subcellular
    use CaMKII_plast

    use AMPA
    use NMDA

    use stims

    use CB1R

    implicit none
    type(pars_type), save :: pars

    type currents_type
        real*8 :: caL13
        real*8 :: TRPV1
        real*8 :: action
        real*8 :: AMPA, NMDA
    end type currents_type

    type conductance_type
        real*8 :: TRPV1, NMDA
    end type conductance_type

    type calcum_fluxes
        real*8 :: tot, CaER, Ca_ch, IP3R, serca, leak
    end type calcum_fluxes

    contains

    real*8 function dV_func(V, I, pars)
        implicit none
        type(currents_type) :: I
        type(pars_type) :: pars
        real*8 :: Itotal, Ileak, V
        Ileak = pars%mem%gL * (V - pars%mem%EL)
        Itotal = -Ileak -I%caL13 -I%TRPV1 -I%AMPA-I%NMDA -I%action
        dV_func = Itotal/pars%mem%Cm
    end function dV_func

    subroutine tables_make
        implicit none
        real*8, parameter :: xst = -100, xfin = 100
        integer, parameter :: n_x = 401
        call stims_tables_make(pars, stims_tabs)
        call caL13_tables_make(xst,xfin,n_x,pars, caL13_tabs)
        call NMDA_tables_make(xst,xfin,n_x,pars, NMDA_tabs)
    end subroutine tables_make

    subroutine tables_clean
        implicit none
        call stims_tables_clean(stims_tabs)
        call caL13_tables_clean(caL13_tabs)
        call NMDA_tables_clean(NMDA_tabs)
        subcellular_compute_once = .true.
    end subroutine tables_clean

    subroutine RHS(NEQ, t, y, dy)
        implicit none
        integer, intent(in) :: NEQ
        real*8, intent(in) :: y(NEQ), t
        real*8, intent(out) :: dy(NEQ)
        type(StateVariables_type) :: v, d
        type(currents_type) :: I
        type(calcum_fluxes) :: J
        type(conductance_type) :: G

        real*8 :: CaM, CaMKIIact

        real*8 :: Glu
        real*8 :: vglu, vplcg, vip3prod, v3k
        real*8 :: ctrl1, ctrl2

        dy = 0

        call SVs_get(v, y)
        call SVs_get(d, dy)

        ! To avoid negative Ca flux
        if (v%Ca_cyt < 0) then
            v%Ca_cyt = 0
        end if

        if (pars%caL13%on == 1) then
            I%caL13 = ical_caL13_func(v%V, v%h_caL13, v%m_caL13, pars%common%Ca_out, v%Ca_cyt, pars)
            call dh_dm_caL13(v%V, v%h_caL13, v%m_caL13, pars, d%h_caL13, d%m_caL13)
        else
            I%caL13 = 0
        end if

        ! stims
        Glu = pars%Glu_release%BaseLevel

        if (pars%stimulation%pre_on == 1) then
            Glu = Glu + Glu_func(t, stims_tabs, pars)
        end if

        if (pars%stimulation%post_on == 1) then
            I%action = Iact_func(t, stims_tabs, pars)
        else
            I%action = 0
        end if

        ! synapse

        ! AMPA
        if (pars%AMPA%on == 1) then
            I%AMPA = i_AMPA_func(V%V, v%o_AMPA, pars%AMPA%gAMPA)
            call do_dd_AMPA(Glu, v%o_AMPA, v%d_AMPA,pars, d%o_AMPA, d%d_AMPA)
        else
            I%AMPA=0
        endif

        ! NMDA
        if (pars%NMDA%on == 1) then
            G%NMDA = g_NMDA_func(v%V, v%o_NMDA)
            I%NMDA = pars%NMDA%gNMDA * v%V * G%NMDA
            call do_NMDA(Glu, v%o_NMDA, pars, d%o_NMDA)
        else
            G%NMDA=0
            I%NMDA=0
        endif

        ! TRPV1
        if (pars%TRPV1%on == 1) then
            G%TRPV1 = g_TRPV1_func(v%AEA, v%V, pars)
            I%TRPV1 = pars%TRPV1%gTRPV1 * v%V * G%TRPV1
        else
            G%TRPV1=0
            I%TRPV1=0
        endif


        ! CaM and CaMKII plasticity
        CaM = CaM_conc(v%Ca_cyt, pars)
        call dy_CaMKII(v%y_CaMKII, v%PP1, CaM, pars,  d%y_CaMKII, CaMKIIact)
        call d_PP1_I1P(v%PP1, v%I1P, CaM, pars, d%PP1, d%I1P)


        ! subcellular calcium, IP3, DAG and 2-AG

        ! Ca, for NMDA's and TRPV1's need to compute conductances G in advance

        d%h_CICR = dh_CICR(v%Ca_cyt, v%IP3, v%h_CICR, pars)

        J%IP3R = JIP3R_CICR_func(v%IP3, v%Ca_cyt, v%Ca_ER, v%h_CICR, pars)
        J%serca = Jserca_CICR_func(v%Ca_cyt, pars)
        J%leak = Jleak_CICR_func(v%Ca_cyt, v%Ca_ER, pars)
        J%CaER = J%IP3R-J%serca+J%leak ! Ca from ER

        J%Ca_ch = -pars%I_to_Ca_flux%VDCC * I%caL13

        J%Ca_ch = J%Ca_ch - pars%I_to_Ca_flux%NMDA*I%NMDA
        J%Ca_ch = J%Ca_ch - pars%I_to_Ca_flux%TRPV1*I%TRPV1

        J%tot = J%CaER + J%Ca_ch

        d%Ca_ER = dCa_ER_func(J%CaER, v%Ca_ER, pars)
        d%Ca_cyt = dCa_cyt_func(J%tot, v%Ca_cyt, pars)

        ! IP3, DAG, ECb
        vglu = vglu_func(Glu,v%Ca_cyt,pars)
        vplcg = vplcg_func(v%IP3,v%Ca_cyt,pars)
        vip3prod = vglu + vplcg

        v3k=v3k_func(v%IP3,CaMKIIact,pars)

        d%IP3 = dIP3_func(v%IP3, vip3prod, pars, v3k)
        d%DAG = dDAG_func(v%DAG, v%DAGLP, vip3prod, pars)

        d%DAGLP = dDAGLP_simple_func(v%DAGLP, v%Ca_cyt, pars)

        if (pars%ECb%on == 1) then
            call dtwoAG_dAEA_ECb(v%twoAG, v%AEA, v%DAG, v%DAGLP, v%Ca_cyt, pars,  d%twoAG, d%AEA)
            if (pars%ECb%CB1R_on == 1) then
                call ctrl1_ctrl2_ECb(pars%ECb%kCB1R*v%o_CB1R, pars, ctrl1, ctrl2)
            else
                call ctrl1_ctrl2_ECb(v%twoAG+pars%ECb%alphaAEACB1*v%AEA, pars,  ctrl1, ctrl2)
            end if
            if (pars%ECb_smooth%on == 0) then
                call dfpre_ECb(Sharp_Om_ECb, ctrl1, ctrl2, v%fpre, pars, d%fpre)
            else
                call dfpre_ECb(Smooth_Om_ECb, ctrl1, ctrl2, v%fpre, pars, d%fpre)
            end if
        end if

        ! membrane potential

        d%V = dV_func(v%V, I, pars)

        ! presynaptic CB1R

        if (pars%CB1R%on == 1) then
            call do_dd_CB1R(v%twoAG+pars%ECb%alphaAEACB1*v%AEA, v%o_CB1R, v%d_CB1R,pars,  d%o_CB1R,d%d_CB1R)
        end if

    end subroutine RHS
end module comp_part
