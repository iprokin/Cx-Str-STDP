!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : TRPV1 channel
! An implementation of the allosteric TRPV1 model by Matta & Ahern, 2007; Calcium flux modeled by GHK.

module TRPV1

    use pars_mod
    use ghk_flux

    implicit none

    contains

    ! Matta & Ahern, 2007 (their eq 3)
    real*8 function g_TRPV1_func(AEA,V,pars)
        implicit none
        real*8 :: AEA,V
        type(pars_type) :: pars
        real*8 :: L,D,C,P,K,Q,x,J
        !con0stants
        L=pars%TRPV1%L
        D=pars%TRPV1%D
        C=pars%TRPV1%C
        P=pars%TRPV1%P
        K=pars%TRPV1%K
        Q=AEA/pars%TRPV1%KD
        x=pars%TRPV1%z*pars%common%F*V/pars%common%RT
        if (x <= 85.0) then
            J=pars%TRPV1%J0*exp(x)
            g_TRPV1_func = 1/(1+(1+J+K+Q+J*K+J*Q+K*Q+J*K*Q)/(L*(1+J*D+K*C+Q*P+J*K*C*D+J*Q*D*P+K*Q*C*P+J*K*Q*D*C*P)))
        else
            g_TRPV1_func = 1/(1+(1+K+Q+K*Q)/(L*(D+K*C*D+Q*D*P+K*Q*D*C*P)))
        endif
    end function g_TRPV1_func

    real*8 function j_ca_TRPV1_func(g, V, calo, cali, pars)
        implicit none
        ! g = i/V
        real*8, intent(in) :: g, V, calo, cali
        type(pars_type), intent(in) :: pars
        j_ca_TRPV1_func = -pars%TRPV1%p_ca * g * ghk(V,cali,calo, pars)
    end function j_ca_TRPV1_func

end module TRPV1
