!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : file defines functions describing periodical stimulation,
! tabulating functions needed for stimulation and helpers defining
! singularity time points for integration
module stims

    use pars_mod
    use general_math

    implicit none

    type stims_tabs_type
        real*8, allocatable :: Glu_tab(:), Iact_tab(:)
        integer :: Glu_tab_li, Iact_tab_li
        real*8 :: Glu_tab_step, Iact_tab_step
    end type stims_tabs_type

    type(stims_tabs_type), save :: stims_tabs

    contains

    subroutine stims_first_tcrit_set(pars, stims_tabs, stims_first_tcrit)
        use qsort_c_module
        implicit none
        type(pars_type), intent(in) :: pars
        type(stims_tabs_type), intent(in) :: stims_tabs
        real*8, intent(out) :: stims_first_tcrit(5)
        ! Glu
        stims_first_tcrit(1)=pars%integration%t_start + pars%stimulation%tpost-pars%stimulation%tsdt-pars%stimulation%dt_stim
        stims_first_tcrit(2)=stims_first_tcrit(1)+stims_tabs%Glu_tab_li*pars%stimulation%tables_step
        ! Iact
        stims_first_tcrit(3)=pars%integration%t_start + pars%stimulation%tpost-2.0*pars%stimulation%tsdt
        stims_first_tcrit(4) = stims_first_tcrit(3) + pars%stimulation%tsdt
        stims_first_tcrit(5)=stims_first_tcrit(3)+stims_tabs%Iact_tab_li*pars%stimulation%tables_step
        call QsortC(stims_first_tcrit)
    end subroutine stims_first_tcrit_set

    subroutine tcrit_per_stim_set_helper(t, tcritlast, period, ITASK, RWORK, LRW, tcrit, ncrit, icrit)
        ! to correctly set TCRIT for the solvers from ODEPACK
        ! pediod is a period of periodic stimulation = 1d0/pars%stimulation%Freq
        ! icrit should be declared in main program and set icrit=1 if all values in tcrit should be considered
        ! call in a loop before a solver. Ex. with LSODA:
        !    ! initializations...
        !    do i=2,n
        !    call tcrit_set_helper(t(i), tcrit, ncrit, ITASK, RWORK, LRW, iCRIT)
        !    call DLSODA(FEX,NEQ,ynew,t_in,t(i),ITOL,RTOL,ATOL,ITASK,ISTATE, IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
        !    y(:,i)=ynew
        !    end do
        implicit none
        integer, intent(in) :: LRW, ncrit
        real*8, intent(in) :: t, tcritlast, period
        real*8, intent(inout) :: RWORK(LRW)
        integer, intent(inout) :: ITASK
        integer, intent(inout) :: icrit
        real*8, intent(inout) :: tcrit(ncrit)
        do
            if (t>tcritlast) then
                ITASK=1
                exit
            end if
            if (icrit<=ncrit) then
                if (tcrit(icrit)>=t) then
                    ITASK=4
                    RWORK(1)=tcrit(icrit)
                    exit
                else
                    icrit=icrit+1
                end if
            else
                icrit = 1
                tcrit = tcrit + period
            end if
        end do
    end subroutine tcrit_per_stim_set_helper

    subroutine stims_tables_make(pars, stims_tabs)
        implicit none
        type(pars_type), intent(in) :: pars
        type(stims_tabs_type), intent(out) :: stims_tabs

        stims_tabs%Glu_tab_step = pars%stimulation%tables_step
        call Glu_pulse(stims_tabs%Glu_tab_step, pars, stims_tabs%Glu_tab)
        stims_tabs%Glu_tab_li = size(stims_tabs%Glu_tab)

        stims_tabs%Iact_tab_step = pars%stimulation%tables_step
        call Iact_pulse(stims_tabs%Iact_tab_step, pars, stims_tabs%Iact_tab)
        stims_tabs%Iact_tab_li = size(stims_tabs%Iact_tab)
    end subroutine stims_tables_make

    subroutine stims_tables_clean(stims_tabs)
        type(stims_tabs_type) :: stims_tabs
        deallocate(stims_tabs%Glu_tab)
        deallocate(stims_tabs%Iact_tab)
    end subroutine stims_tables_clean

    subroutine Glu_pulse(table_step, pars, Tabl)
        implicit none
        real*8, intent(in) :: table_step
        type(pars_type), intent(in) :: pars
        integer :: i
        integer :: li
        real*8, allocatable, intent(out) :: Tabl(:)
        real*8 :: t
        li = pars%Glu_release%tauGlu*30/table_step+1
        allocate(Tabl(li))
        do i = 1,li
            t=(i-1)*table_step
            Tabl(i) = pars%Glu_release%Glumax*exp(-t/pars%Glu_release%tauGlu)
        end do
    end subroutine Glu_pulse

    subroutine Iact_pulse(table_step, pars, Tabl)
        implicit none
        real*8, intent(in) :: table_step
        type(pars_type), intent(in) :: pars
        integer :: i
        integer :: li
        real*8, allocatable, intent(out) :: Tabl(:)
        real*8 :: t
        li = int(pars%action%APdur/table_step)+1
        allocate(Tabl(li))
        do i = 1,li
            t=(i-1)*table_step - pars%stimulation%tsdt
            Tabl(i) = -pars%action%DPmax -pars%action%APmax * exp(-t/pars%action%tausbAP)*heaviside(t)
        end do
    end subroutine Iact_pulse

    ! functions for periodical stims

    real*8 function Glu_func(t, stims_tabs, pars)
        implicit none
        type(pars_type) :: pars
        type(stims_tabs_type) :: stims_tabs
        real*8 :: t, t0, tp
        t0=pars%integration%t_start + pars%stimulation%tpost-pars%stimulation%tsdt-pars%stimulation%dt_stim
        tp=t-t0
        if (tp>0) then
            Glu_func = kernel_repeat(tp, pars%stimulation%Freq, int(pars%stimulation%num_stim), stims_tabs%Glu_tab, stims_tabs%Glu_tab_li, stims_tabs%Glu_tab_step)
        else
            Glu_func = 0.0
        end if
    end function Glu_func

    real*8 function Iact_func(t, stims_tabs, pars)
        implicit none
        type(pars_type) :: pars
        type(stims_tabs_type) :: stims_tabs
        real*8 :: t,t0, tp
        t0=pars%integration%t_start + pars%stimulation%tpost-2.0*pars%stimulation%tsdt
        tp=t-t0
        if (tp>0) then
            Iact_func = kernel_repeat(tp, pars%stimulation%Freq, int(pars%stimulation%num_stim), stims_tabs%Iact_tab, stims_tabs%Iact_tab_li, stims_tabs%Iact_tab_step)
        else
            Iact_func = 0.0
        end if
    end function Iact_func

end module stims
