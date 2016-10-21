!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : solver of ODE system with lsoda and interface with python

subroutine lsoda_calc(y, iStatus, y0, lv, t, n)
    use comp_part
    implicit none
    external :: FEX
    integer, parameter :: NEQ=32
    integer, parameter :: JT = 2
    integer, parameter :: LRW = 22 + NEQ * max(16, NEQ + 9)
    integer, parameter :: LIW = 20 + NEQ
    real*8 :: RWORK(LRW)
    integer :: IWORK(LIW), JDUM
    integer :: IOPT=0, ISTATE, ITASK, iCRIT
    integer, parameter :: ITOL = 1
    integer, intent(in) :: n
    real*8, intent(in) :: lv(165)
    real*8, intent(in) :: y0(NEQ), t(n)
    real*8, intent(out) :: y(NEQ,n)
    integer, intent(out) :: iStatus
    integer, parameter :: ncrit = 7
    real*8 :: ynew(NEQ), tcrit(ncrit), tcritlast
    real*8 :: t_in
    integer :: i

    ITASK = 1
    ISTATE = 1
    iCRIT=1

    call parameters_get(pars, lv)
    call tables_make
    call stims_first_tcrit_set(pars, stims_tabs, tcrit)
    tcritlast = tcrit(ncrit)+pars%stimulation%num_stim/pars%stimulation%Freq

    ! setting optional integration parameters
    IOPT=1
    RWORK(5:10)=0.0
    IWORK(5:10)=0
    RWORK(6)=pars%integration%HMAX
    IWORK(6)=int(pars%integration%MXSTEP) ! number of interally defined steps for each integration step

    ynew=y0
    y(:,1)=y0
    t_in=t(1)
    do i=2,n
        call tcrit_per_stim_set_helper(t(i), tcritlast, 1d0/pars%stimulation%Freq, ITASK, RWORK, LRW, tcrit, ncrit, iCRIT)
        call DLSODA(FEX,NEQ,ynew,t_in,t(i),ITOL, pars%integration%RTOL, pars%integration%ATOL, ITASK,ISTATE, IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
        if (ISTATE<0) exit
        y(:,i)=ynew
    end do
    iStatus = ISTATE
    call tables_clean
end subroutine lsoda_calc

subroutine FEX(NEQ, T, Y, YDOT)
    use comp_part
    implicit none
    integer, intent(in) :: NEQ
    real*8, intent(in) :: T, Y(NEQ)
    real*8, intent(out) :: YDOT(NEQ)
    call RHS(NEQ, T, Y, YDOT)
end subroutine FEX


! --- helpers to be called from python ---


! python interface to functions from comp_part
module extras
    use comp_part
    implicit none
    contains

    subroutine CaM_conc_func(Ca_cyt, n, lv, CaM)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: Ca_cyt(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: CaM(n)
        integer :: i
        !f2py depend(n) Ca_cyt :: n=len(Ca_cyt)
        !f2py depend(n) CaM
        call parameters_get(pars, lv)
        do i=1,n
            CaM(i) = CaM_conc(Ca_cyt(i), pars)
        end do
    end subroutine CaM_conc_func

    subroutine CaMKIIpho(y, n, CaMKIIpho_answer)
        integer, intent(in) :: n
        real*8, intent(in) :: y(n,13)
        real*8, intent(out) :: CaMKIIpho_answer(n)
        integer :: i
        !f2py depend(n) y :: n=len(y)
        !f2py depend(n) CaMKIIpho_answer
        forall(i=1:n)
            CaMKIIpho_answer(i)=CaMKIIpho_func(y(i,:)) ! CaMKIIpho_func(y(i,:)) - Probably wrong! Not wrong for this definition
        end forall
    end subroutine CaMKIIpho

    subroutine CaMKIIpho_4F(y, n, CaMKIIpho_answer)
        integer, intent(in) :: n
        real*8, intent(in) :: y(13, n)
        real*8, intent(out) :: CaMKIIpho_answer(n)
        integer :: i
        !f2py depend(n) y :: n=len(y)
        !f2py depend(n) CaMKIIpho_answer
        forall(i=1:n)
            CaMKIIpho_answer(i)=CaMKIIpho_func(y(:,i))
        end forall
    end subroutine CaMKIIpho_4F

    subroutine ical_caL13(v, h, m, calo, cali, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: v(n), h(n), m(n), calo, cali(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        integer :: i
        call parameters_get(pars, lv)
        call caL13_tables_make(-100d0,100d0,401,pars, caL13_tabs)
        do i=1,n
            answer(i) = ical_caL13_func(v(i), h(i), m(i), calo, cali(i), pars)
        end do
        call caL13_tables_clean(caL13_tabs)
    end subroutine ical_caL13

    subroutine i_TRPV1(AEA, V, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: AEA(n), V(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        integer :: i
        call parameters_get(pars, lv)
        do i=1,n
            answer(i) = pars%TRPV1%gTRPV1 * V(i) * g_TRPV1_func(AEA(i), V(i), pars)
        end do
    end subroutine i_TRPV1

    subroutine j_ca_TRPV1(AEA, Ca_cyt, V, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: AEA(n), Ca_cyt(n), V(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        real*8 :: G_TRPV1
        integer :: i
        call parameters_get(pars, lv)
        do i=1,n
            G_TRPV1 = g_TRPV1_func(AEA(i), V(i), pars)
            answer(i) = -pars%I_to_Ca_flux%TRPV1 * pars%TRPV1%gTRPV1 * V(i) * G_TRPV1
        end do
    end subroutine j_ca_TRPV1

    subroutine i_AMPA(V, o_AMPA, gAMPAmax, n, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: V(n), o_AMPA(n), gAMPAmax
        real*8, intent(out) :: answer(n)
        integer :: i
        do i=1,n
            answer(i) = i_AMPA_func(V(i), o_AMPA(i), gAMPAmax)
        end do
    end subroutine i_AMPA

    subroutine i_NMDA(V, o_NMDA, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: V(n), o_NMDA(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        integer :: i
        call parameters_get(pars, lv)
        call NMDA_tables_make(-100d0,100d0,401,pars, NMDA_tabs)
        do i=1,n
            answer(i) = pars%NMDA%gNMDA * V(i) * g_NMDA_func(V(i), o_NMDA(i))
        end do
        call NMDA_tables_clean(NMDA_tabs)
    end subroutine i_NMDA

    subroutine j_ca_NMDA(V, Ca_cyt, o_NMDA, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: V(n), Ca_cyt(n), o_NMDA(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        real*8 :: G_NMDA
        integer :: i
        call parameters_get(pars, lv)
        call NMDA_tables_make(-100d0,100d0,401,pars, NMDA_tabs)
        do i=1,n
            G_NMDA=g_NMDA_func(V(i), o_NMDA(i))
            answer(i) = -pars%I_to_Ca_flux%NMDA * pars%NMDA%gNMDA * V(i) * G_NMDA
        end do
        call NMDA_tables_clean(NMDA_tabs)
    end subroutine j_ca_NMDA

    subroutine JIP3R_CICR(IP3, Ca_cyt, Ca_ER, h, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: IP3(n), Ca_cyt(n), Ca_ER(n), h(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        integer :: i
        call parameters_get(pars, lv)
        do i=1,n
            answer(i) = JIP3R_CICR_func(IP3(i), Ca_cyt(i), Ca_ER(i), h(i), pars)
        end do
    end subroutine JIP3R_CICR

    subroutine Jserca_CICR(Ca_cyt, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: Ca_cyt(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        integer :: i
        call parameters_get(pars, lv)
        do i=1,n
            answer(i) = Jserca_CICR_func(Ca_cyt(i), pars)
        end do
    end subroutine Jserca_CICR

    subroutine Jleak_CICR(Ca_cyt, Ca_ER, n, lv, answer)
        integer, intent(in) :: n
        real*8, intent(in) :: Ca_cyt(n), Ca_ER(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: answer(n)
        integer :: i
        call parameters_get(pars, lv)
        do i=1,n
            answer(i) = Jleak_CICR_func(Ca_cyt(i), Ca_ER(i), pars)
        end do
    end subroutine Jleak_CICR

end module extras

!!! steady state solutions

module steady_states
    use comp_part
    implicit none
    contains
    ! Ca channels
    subroutine h0_m0_caL13(v, n, lv, h, m)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: v(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: h(n), m(n)
        integer :: i
        call parameters_get(pars, lv)
        call caL13_tables_make(-100d0,100d0,401,pars, caL13_tabs)
        do i=1,n
            h(i) = caL13_hinf(v(i))
            m(i) = caL13_minf(v(i))
        end do
        call caL13_tables_clean(caL13_tabs)
    end subroutine h0_m0_caL13
    ! NMDA
    subroutine o0_NMDA(Glu, n, lv, o_NMDA)
        ! o_NMDA - probability of the open state
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: Glu(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: o_NMDA(n)
        call parameters_get(pars, lv)
        call NMDA_tables_make(-100d0,100d0,401,pars, NMDA_tabs)
        o_NMDA = pars%NMDA%Alpha*Glu/( pars%NMDA%Alpha*Glu + pars%NMDA%Beta )
        call NMDA_tables_clean(NMDA_tabs)
    end subroutine o0_NMDA
    ! CICR
    subroutine h0_CICR(Ca_cyt, IP3, n, lv, h_CICR)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: Ca_cyt(n), IP3(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: h_CICR(n)
        real*8 :: H1(n)
        call parameters_get(pars, lv)
        H1 = (pars%CICR%a2 * pars%CICR%d2 * (IP3+pars%CICR%d1)/(IP3+pars%CICR%d3))
        h_CICR = H1/(H1+pars%CICR%a2*Ca_cyt)
    end subroutine h0_CICR
    subroutine o0_d0_AMPA(Glu, n, lv, o_AMPA0,d_AMPA0)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: Glu(n)
        real*8, intent(in) :: lv(165)
        real*8, intent(out) :: o_AMPA0(n), d_AMPA0(n)
        integer :: i
        call parameters_get(pars, lv)
        do i=1,n
            call o_d_AMPA_IC_setup(Glu(i), pars,  o_AMPA0(i), d_AMPA0(i))
        end do
    end subroutine o0_d0_AMPA
end module steady_states
