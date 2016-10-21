!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : A simple three state kinetic model of activation of cannabinoid receptor type 1 (CB1R)
! by 2-arachidonoylglycerol (2-AG) and anandamide (AEA)
! The modeling approach is similar the one of Destexhe et al 1994.

module CB1R

    use pars_mod
    !use pars_CB1R_mod

    implicit none

    contains

    subroutine o_d_CB1R_IC_setup(eCB_bl, pars,  o_CB1R0, d_CB1R0)
        implicit none
        real*8, intent(in) :: eCB_bl
        type(pars_type), intent(in) :: pars
        real*8 :: gameps, Ab
        real*8, intent(out) :: o_CB1R0, d_CB1R0
        gameps=pars%CB1R%Gamma/pars%CB1R%Epsilon
        Ab=(pars%CB1R%Beta+pars%CB1R%Gamma)/pars%CB1R%Alpha
        o_CB1R0 = eCB_bl/(Ab+(1+gameps)*eCB_bl)
        d_CB1R0 = o_CB1R0*gameps
    end subroutine o_d_CB1R_IC_setup

    subroutine do_dd_CB1R(eCB,o_CB1R,d_CB1R,pars,  do_CB1R,dd_CB1R)
        implicit none
        real*8, intent(in) :: eCB,o_CB1R,d_CB1R
        type(pars_type), intent(in) :: pars
        real*8 :: c_CB1R
        real*8, intent(out) :: do_CB1R, dd_CB1R
        ! o_CB1R - open state probatility
        ! d_CB1R - the fraction of desensetisized or moved to extracellular compart state probatility
        ! c_CB1R - closed state probatility
        ! C <-alpha,beta-> O -gamma-> D -Epsilon-> C
        c_CB1R = 1-o_CB1R-d_CB1R
        do_CB1R = pars%CB1R%Alpha*eCB*c_CB1R - (pars%CB1R%Beta+pars%CB1R%Gamma)*o_CB1R
        dd_CB1R = -pars%CB1R%Epsilon*d_CB1R + pars%CB1R%Gamma*o_CB1R
    end subroutine do_dd_CB1R

end module CB1R
