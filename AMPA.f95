!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : An implementation of simple AMPAR model of Destexhe, A, Z F Mainen, and T J Sejnowski. “Synthesis of Models for
! Excitable Membranes, Synaptic Transmission and Neuromodulation Using a Common Kinetic Formalism.” Journal of Computational
! Neuroscience 1, no. 3 (August 1994): 195–230.

module AMPA

    use pars_mod
    implicit none

    contains

    real*8 function i_AMPA_func(V, o_AMPA, gAMPAmax)
        implicit none
        real*8 :: gAMPAmax,o_AMPA,V
        i_AMPA_func = gAMPAmax*o_AMPA*V
    end function i_AMPA_func

    subroutine o_d_AMPA_IC_setup(glu_bl, pars,  o_AMPA0, d_AMPA0)
        implicit none
        real*8, intent(in) :: glu_bl
        type(pars_type), intent(in) :: pars
        real*8 :: gameps, Ab
        real*8, intent(out) :: o_AMPA0, d_AMPA0
        gameps = pars%AMPA%Gamma/pars%AMPA%Epsilon
        Ab = (pars%AMPA%Beta+pars%AMPA%Gamma)/pars%AMPA%Alpha
        o_AMPA0 = glu_bl/(Ab+(1+gameps)*glu_bl)
        d_AMPA0 = o_AMPA0*gameps
    end subroutine o_d_AMPA_IC_setup

    subroutine do_dd_AMPA(Glu,o_AMPA,d_AMPA,pars,  do_AMPA,dd_AMPA)
        implicit none
        real*8, intent(in) :: Glu,o_AMPA,d_AMPA
        type(pars_type), intent(in) :: pars
        real*8 :: c_AMPA
        real*8, intent(out) :: do_AMPA, dd_AMPA
        ! o_AMPA - open state probatility
        ! d_AMPA - desensetisized state probatility
        ! c_AMPA - closed state probatility
        ! C <-alpha,beta-> O -gamma-> D -Epsilon-> C
        c_AMPA=1-o_AMPA-d_AMPA
        do_AMPA = pars%AMPA%Alpha*Glu*c_AMPA - (pars%AMPA%Beta+pars%AMPA%Gamma)*o_AMPA
        dd_AMPA = -pars%AMPA%Epsilon*d_AMPA + pars%AMPA%Gamma*o_AMPA
    end subroutine do_dd_AMPA

end module AMPA
