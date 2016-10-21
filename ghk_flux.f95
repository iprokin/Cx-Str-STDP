!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : GHK flux
module ghk_flux
    implicit none
    contains

    real*8 function efun(x)
        implicit none
        real*8 :: x
        if (abs(x)<1e-4) then
            efun = 1-x/2
        elseif (x<=85) then
            efun = x/(exp(x)-1)
        else
            efun = 0
        endif
    end function efun

    real*8 function ghk(v, ci, co, pars)
        use pars_mod, only : pars_type
        implicit none
        real*8 :: v, ci, co
        type(pars_type) :: pars
        real*8 :: x, eco, eci
        ! (v(mV), ci(mM), co(mM)) (.001 coul/cm3)
        x = pars%common%zS*pars%common%F*v*1e-3/pars%common%RT
        ! 1e-3 mV to V
        eco = co*efun(x)
        eci = ci*efun(-x)
        ! high cao charge moves inward
        ! negative potential charge moves inward
        ghk = pars%common%zS*pars%common%F*(eci - eco)
    end function ghk

end module ghk_flux
