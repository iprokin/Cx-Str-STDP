!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : NMDAR model
module NMDA

    use pars_mod
    use general_math
    use ghk_flux
    implicit none

    type NMDA_tabs_type
    real*8, allocatable :: B(:)
    integer :: n
    real*8 :: xst, xfin, xstep
    end type NMDA_tabs_type

    type(NMDA_tabs_type), save :: NMDA_tabs

    contains

    subroutine NMDA_tables_make(xst,xfin,n_x,pars, NMDA_tabs)
        integer, intent(in) :: n_x
        real*8, intent(in) :: xst,xfin
        type(pars_type), intent(in) :: pars
        type(NMDA_tabs_type), intent(out) :: NMDA_tabs
        real*8 :: v(n_x), xstep
        integer :: i
        xstep = (xfin-xst)/n_x
        NMDA_tabs%n=n_x
        NMDA_tabs%xst=xst
        NMDA_tabs%xfin=xfin
        NMDA_tabs%xstep=xstep
        allocate(NMDA_tabs%B(n_x))
        forall (i=1:n_x)
            v(i) = xst+xstep*(i-1)
            NMDA_tabs%B(i) = 1.0/(1+pars%NMDA%Mg/3.57*exp(-0.062*v(i)))
        end forall
    end subroutine NMDA_tables_make

    subroutine NMDA_tables_clean(NMDA_tabs)
        type(NMDA_tabs_type) :: NMDA_tabs
        deallocate(NMDA_tabs%B)
    end subroutine NMDA_tables_clean

    real*8 function NMDA_B(v)
        real*8, intent(in) :: v
        NMDA_B=lin_interpTable_brds(NMDA_tabs%B, NMDA_tabs%n, (v-NMDA_tabs%xst)/NMDA_tabs%xstep)
    end function NMDA_B

    real*8 function g_NMDA_func(V, o_NMDA)
        implicit none
        real*8 :: o_NMDA,V
        g_NMDA_func =  o_NMDA * NMDA_B(V)
    end function g_NMDA_func

    subroutine do_NMDA(Glu,o_NMDA,pars, d_o_NMDA)
        ! o_NMDA - probability of the open state
        implicit none
        real*8, intent(in) :: Glu,o_NMDA
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: d_o_NMDA
        d_o_NMDA = pars%NMDA%Alpha*Glu*(1-o_NMDA)-pars%NMDA%Beta*o_NMDA
    end subroutine do_NMDA

end module NMDA
