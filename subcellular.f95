!        -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : This module defines models of mGluR-dependent signaling, CICR, calcium buffering, production of endocannabinoids,
! phenomenological model of endocannabinoid-dependent plasticity
module subcellular

    use pars_mod
    use general_math

    implicit none

    logical, save :: subcellular_compute_once = .true.

    contains

    real*8 function tauCa(C, pars)
        implicit none
        real*8 :: C
        type(pars_type) :: pars
        tauCa = (1+pars%CaBuff%BT/pars%CaBuff%KdB/(1+C/pars%CaBuff%KdB)**2)
    end function tauCa

    real*8 function JIP3R_CICR_func(IP3, Ca_cyt, Ca_ER, h, pars)
        implicit none
        real*8 :: IP3, Ca_cyt, Ca_ER, h
        type(pars_type) :: pars
        real*8 :: minf, ninf
        minf = IP3/(IP3+pars%CICR%d1)
        ninf = Ca_cyt/(Ca_cyt+pars%CICR%d5)
        JIP3R_CICR_func = pars%CICR%rc * (minf*ninf*h)**3 * (Ca_ER-Ca_cyt)
    end function JIP3R_CICR_func

    real*8 function Jserca_CICR_func(Ca_cyt, pars)
        implicit none
        real*8 :: Ca_cyt
        type(pars_type) :: pars
        Jserca_CICR_func = pars%CICR%ver * Hill(Ca_cyt,pars%CICR%ker,2)
    end function Jserca_CICR_func

    real*8 function Jleak_CICR_func(Ca_cyt, Ca_ER, pars)
        implicit none
        real*8 :: Ca_cyt, Ca_ER
        type(pars_type) :: pars
        Jleak_CICR_func = pars%CICR%rl * (Ca_ER-Ca_cyt)
    end function Jleak_CICR_func

    real*8 function dh_CICR(Ca_cyt, IP3, h, pars)
        implicit none
        real*8 :: Ca_cyt, IP3, h
        type(pars_type) :: pars
        type(CICR_type) :: p
        p = pars%CICR
        dh_CICR = (p%a2 * p%d2 * (IP3+p%d1)/(IP3+p%d3))*(1-h) - p%a2*Ca_cyt*h
    end function dh_CICR

    real*8 function dCa_cyt_func(JCa, Ca_cyt, pars)
        implicit none
        real*8 :: JCa, Ca_cyt
        type(pars_type) :: pars
        dCa_cyt_func = ( JCa - ( Ca_cyt - pars%CaBuff%Cab )/pars%CaBuff%tauCab )/tauCa(Ca_cyt, pars)
    end function dCa_cyt_func

    real*8 function dCa_ER_func(JCaER, Ca_ER, pars)
        implicit none
        real*8 :: JCaER, Ca_ER
        ! JCaER = JIP3R-Jserca+Jleak
        type(pars_type) :: pars
        dCa_ER_func = -JCaER*pars%CICR%rhoER/tauCa(Ca_ER, pars)
    end function dCa_ER_func

    real*8 function vglu_func(Glu,Ca_cyt,pars)
        implicit none
        real*8 :: Glu,Ca_cyt
        type(pars_type) :: pars
        real*8 :: Hill1
        Hill1=Ca_cyt/(Ca_cyt+pars%IP3%kpi)
        vglu_func = pars%IP3%vbeta * Glu / (Glu+( pars%IP3%kr * (1+pars%IP3%kp/pars%IP3%kr*Hill1) ))
    end function vglu_func

    real*8 function vplcg_func(IP3,Ca_cyt,pars)
        implicit none
        real*8 :: IP3,Ca_cyt
        type(pars_type) :: pars
        vplcg_func = pars%IP3%vdelta/(1+IP3/pars%IP3%kappad)*Hill(Ca_cyt,pars%IP3%kdelta,2)
    end function vplcg_func

    real*8 function v3k_func(IP3,CaMKIIact,pars)
        implicit none
        real*8 :: IP3,CaMKIIact
        type(pars_type) :: pars
        v3k_func = pars%IP3%v3k*CaMKIIact * Hill(IP3,pars%IP3%k3,int(pars%IP3%n3))
    end function v3k_func

    real*8 function dIP3_func(IP3, vip3prod, pars,    v3k)
        implicit none
        real*8 :: IP3, vip3prod,    v3k
        type(pars_type) :: pars
        dIP3_func = vip3prod - v3k - pars%IP3%r5p*IP3
    end function dIP3_func

    real*8 function dDAG_func(DAG, DAGLP, vip3prod, pars)
        implicit none
        real*8 :: DAG, DAGLP, vip3prod
        type(pars_type) :: pars
        type(DGLandDAG_type) :: p
        p = pars%DGLandDAG
        dDAG_func = vip3prod-p%rDGL*DAGLP*DAG/(DAG+p%KDGL)-p%kDAGK*DAG
    end function dDAG_func

    real*8 function dDAGLP_simple_func(DAGLP, Ca_cyt, pars)
        implicit none
        real*8 :: DAGLP, Ca_cyt
        type(pars_type) :: pars
        type(KandP_on_DAGLP_type) :: p
        p = pars%KandP_on_DAGLP
        dDAGLP_simple_func = p%rK*Ca_cyt**int(p%nK)*(1-DAGLP) - p%rP*DAGLP
    end function dDAGLP_simple_func

    subroutine dtwoAG_dAEA_ECb(twoAG, AEA, DAG, DAGLP, Ca_cyt, pars,    dtwoAG, dAEA)
        implicit none
        real*8, intent(in) :: twoAG, AEA, DAG, DAGLP, Ca_cyt
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: dtwoAG, dAEA
        type(ECb_type) :: pECb
        type(DGLandDAG_type) :: parDnD
        pECb = pars%ECb
        parDnD = pars%DGLandDAG
        !d2-AG/dt
        dtwoAG = parDnD%rDGL * DAGLP * DAG / (DAG + parDnD%KDGL) - parDnD%kMAGL * twoAG
        !dAEA/dt
        dAEA = pECb%vATAEA * Ca_cyt - pECb%vFAAH * AEA / (pECb%KFAAH + AEA)
    end subroutine dtwoAG_dAEA_ECb

    subroutine ctrl1_ctrl2_ECb(x, pars,    ctrl1, ctrl2)
        implicit none
        real*8, intent(in) :: x
        real*8, intent(out) :: ctrl1, ctrl2
        type(pars_type), intent(in) :: pars
        ctrl1 = x + pars%DA%gamma1DA*pars%DA%DA
        ctrl2 = x + pars%DA%gamma2DA*pars%DA%DA
    end subroutine ctrl1_ctrl2_ECb

    subroutine Smooth_Om_ECb(x, fpre, pars, Omega)
        implicit none
        real*8, intent(in) :: x, fpre
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: Omega
        ! {to tabulate LTPwin}
        integer, parameter :: n_x=100
        real*8, save :: xst, xstep, LTPwin_tab(n_x)
        real*8 :: v(n_x), ka
        integer :: i
        ! end {to tabulate LTPwin}
        real*8 :: x0, LTDw, w, LTDwin, LTPwin
        x0 = 0.5*(pars%ECb%LTDstart + pars%ECb%LTDstop)
        LTDw = x0-pars%ECb%LTDstart
        w = pars%ECb_smooth%kw*LTDw
        ka = pars%ECb_smooth%K + pars%ECb_smooth%kadd * (1-abs(x-x0)/LTDw)**int(pars%ECb_smooth%kn) * Rect(x, x0-LTDw, x0+LTDw)
        LTDwin = (Hill(x-x0, ka, int(pars%ECb_smooth%n))-1)*Rect(x,x0-w,x0+w)

        ! tabulate for speed
        if (subcellular_compute_once) then
            subcellular_compute_once = .false.
            xst = pars%ECb%LTPstart-6*pars%ECb_smooth%tau
            xstep = (2*6*pars%ECb_smooth%tau)/n_x
            forall (i=1:n_x)
                v(i) = xst+xstep*(i-1)
                LTPwin_tab(i) = 1/(1+exp(-(v(i)-pars%ECb%LTPstart)/pars%ECb_smooth%tau))
            end forall
        end if
        LTPwin = lin_interpTable_brds(LTPwin_tab, n_x, (x-xst)/xstep)

        Omega = 1.0 + pars%ECb%LTDMax*LTDwin + pars%ECb%LTPMax*LTPwin
    end subroutine Smooth_Om_ECb

    subroutine Sharp_Om_ECb(x, fpre, pars, Omega)
        implicit none
        real*8, intent(in) :: x, fpre
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: Omega
        Omega = 1 - pars%ECb%LTDMax*(heaviside(x-pars%ECb%LTDstart)-heaviside(x-pars%ECb%LTDstop)) + heaviside(x-pars%ECb%LTPstart)*pars%ECb%LTPMax
    end subroutine Sharp_Om_ECb

    subroutine dfpre_ECb(Om_func, ctrl1, ctrl2, fpre, pars,    dfpre)
        implicit none
        real*8, intent(in) :: ctrl1, ctrl2, fpre
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: dfpre
        type(ECb_type) :: pECb
        real*8 :: Omega, taufpre
        interface
            subroutine Om_func(x, fpre, pars, Omega)
                use pars_mod
                implicit none
                real*8, intent(in) :: x, fpre
                type(pars_type), intent(in) :: pars
                real*8, intent(out) :: Omega
            end subroutine Om_func
        end interface
        pECb = pars%ECb
        call Om_func(ctrl1, fpre, pars, Omega)
        !make it static
        taufpre = pECb%P1 / ( pECb%P2**pECb%P3 + ctrl2**pECb%P3 ) + pECb%P4
        !dfpre/dt
        dfpre = (Omega-fpre)/taufpre
    end subroutine dfpre_ECb

end module subcellular
