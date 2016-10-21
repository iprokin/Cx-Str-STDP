!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : Voltage-Dependent Calcium Channel L-type (1.3)
! converted from https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=112834&file=%2Fnacb_msp%2FcaL13.mod (Wolf et al 2005)

module caL13

    use pars_mod
    use general_math
    use ghk_flux
    implicit none

    type caL13_tabs_type
        real*8, allocatable :: mtau(:), hinf(:), minf(:)
        integer :: n
        real*8 :: xst, xfin, xstep
    end type caL13_tabs_type

    type(caL13_tabs_type), save :: caL13_tabs

    contains

    subroutine caL13_tables_make(xst,xfin,n_x,pars, caL13_tabs)
        integer, intent(in) :: n_x
        real*8, intent(in) :: xst,xfin
        type(pars_type), intent(in) :: pars
        type(caL13_tabs_type), intent(out) :: caL13_tabs
        real*8 :: v, xstep
        real*8 :: malpha, mbeta
        integer :: i
        xstep = (xfin-xst)/n_x
        caL13_tabs%n=n_x
        caL13_tabs%xst=xst
        caL13_tabs%xfin=xfin
        caL13_tabs%xstep=xstep
        allocate(caL13_tabs%mtau(n_x))
        allocate(caL13_tabs%hinf(n_x))
        allocate(caL13_tabs%minf(n_x))
        do i=1,n_x
            v = xst+xstep*(i-1)
            caL13_tabs%minf(i) = 1  /  ( 1 + exp( (v-pars%caL13%mvhalf-pars%caL13%mshift) / pars%caL13%mslope) )
            caL13_tabs%hinf(i) = 1  /  ( 1 + exp( (v-pars%caL13%hvhalf-pars%caL13%hshift) / pars%caL13%hslope) )
            ! to match caL13_tabs%hinf(i) data from Bell 2001, fig 2

            malpha = pars%caL13%c * (v-pars%caL13%vm) / ( exp((v-pars%caL13%vm)/pars%caL13%k) - 1 )
            mbeta = pars%caL13%cpr * exp(v/pars%caL13%kpr)! Kasai 1992, fig 15
            caL13_tabs%mtau(i) = 1 / (malpha + mbeta)
        end do
    end subroutine caL13_tables_make

    subroutine caL13_tables_clean(caL13_tabs)
        type(caL13_tabs_type) :: caL13_tabs
        deallocate(caL13_tabs%mtau)
        deallocate(caL13_tabs%hinf)
        deallocate(caL13_tabs%minf)
    end subroutine caL13_tables_clean

    real*8 function caL13_mtau(v)
        real*8, intent(in) :: v
        caL13_mtau=lin_interpTable_brds(caL13_tabs%mtau, caL13_tabs%n, (v-caL13_tabs%xst)/caL13_tabs%xstep)
    end function caL13_mtau

    real*8 function caL13_hinf(v)
        real*8, intent(in) :: v
        caL13_hinf=lin_interpTable_brds(caL13_tabs%hinf, caL13_tabs%n, (v-caL13_tabs%xst)/caL13_tabs%xstep)
    end function caL13_hinf

    real*8 function caL13_minf(v)
        real*8, intent(in) :: v
        caL13_minf=lin_interpTable_brds(caL13_tabs%minf, caL13_tabs%n, (v-caL13_tabs%xst)/caL13_tabs%xstep)
    end function caL13_minf

    real*8 function ical_caL13_func(v, h, m, calo, cali, pars)
        implicit none
        real*8 :: h, calo, m, cali, v
        type(pars_type) :: pars
        ical_caL13_func = ghk(v,cali,calo, pars) * pars%caL13%pcaLbar * m * m * h  ! Kasai 92, Brown 93
    end function ical_caL13_func

    subroutine dh_dm_caL13(v, h, m, pars, dh, dm)
        implicit none
        real*8, intent(in) :: v, h, m
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: dh, dm
        dh = (caL13_hinf(v) - h) / (pars%caL13%htau/pars%caL13%hqfact)
        dm = (caL13_minf(v) - m) / (caL13_mtau(v)/pars%caL13%qfact)
    end subroutine dh_dm_caL13

end module caL13
