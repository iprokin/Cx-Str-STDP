!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : An implementation of bistable CaMKII model of Graupner, Michael, and Nicolas Brunel. “STDP in a Bistable Synapse Model Based on CaMKII and
! Associated Signaling Pathways.” Edited by Karl J Friston. PLoS Computational Biology 3, no. 11 (November 2007): e221–e221.
! doi:10.1371/journal.pcbi.0030221.

module CaMKII_plast

    use pars_mod

    implicit none

    contains

    real*8 function CaM_conc(Ca_cyt, pars) 
        implicit none
        real*8 :: Ca_cyt
        type(pars_type) :: pars
        type(post_CaMKII_plast_type) :: p
        p = pars%post_CaMKII_plast
        CaM_conc = p%CaMT/(1 + p%Ka4/Ca_cyt + p%Ka3*p%Ka4/(Ca_cyt**2) + p%Ka2*p%Ka3*p%Ka4/(Ca_cyt**3) + p%Ka1*p%Ka2*p%Ka3*p%Ka4/(Ca_cyt**4))
    end function CaM_conc

    subroutine dy_CaMKII(y, PP1, CaM, pars,   dy, phossum)
        implicit none
        real*8, intent(in) :: y(13), PP1, CaM
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: phossum
        real*8, intent(out) :: dy(13)
        type(post_CaMKII_plast_type) :: p
        real*8 :: B0, sum_y23, sum_y24, sum_y57, sum_y58, sum_y911, rr
        real*8 :: k10, gamma, gamma2, k6gamma2, k7gamma
        dy = 0.0
        p = pars%post_CaMKII_plast
        rr=sum(y)
        ! B0 is whats left from total
        B0=2*p%CaMKT-rr
        ! kinetic equations
        sum_y23=sum(y(2:3))
        sum_y24=sum_y23+y(4)
        sum_y57=sum(y(5:7))
        sum_y58=sum_y57+y(8)
        sum_y911=sum(y(9:11))
        phossum=y(1) + 2*sum_y24 + 3*sum_y58 + 4*sum_y911 + 5*y(12) + 6*y(13)
        k10=p%k12*PP1/(p%KM + phossum)
        gamma=CaM/(p%K5+CaM)

        !dBi/dt
        gamma2=gamma*gamma
        k6gamma2=p%k6*gamma2
        k7gamma=p%k7*gamma
        dy(1) = 6*k6gamma2*B0 - (4*k6gamma2 + k7gamma + k10)*y(1) + 2*k10*sum_y24
        dy(2) = (k7gamma + k6gamma2)*y(1) - (3*k6gamma2 + k7gamma + 2*k10)*y(2) + k10*(y(5) + sum_y57)
        dy(3) = 2*k6gamma2*y(1) - 2*(k7gamma + k6gamma2 + k10)*y(3) + k10*(sum_y57 + 3*y(8)) 
        dy(4) = k6gamma2*y(1) - 2*(k7gamma + k6gamma2 + k10)*y(4) + k10*(y(6) + y(7))
        dy(5) = k7gamma*(sum_y23 - y(5)) + k6gamma2*(y(2) - 2*y(5)) + k10*(2*y(9) + y(10) - 3*y(5))
        dy(6) = k6gamma2*(sum_y23 - y(6))  + k7gamma*(2*y(4)  - 2*y(6)) +k10*(-3*y(6) + sum_y911 + y(11))
        dy(7) = k6gamma2*(y(2) + 2*y(4) - y(7)) + k7gamma*(y(3) - 2*y(7)) +k10*(-3*y(7) + y(9) + y(10) + 2*y(11))
        dy(8) = k6gamma2*y(3) - 3*k7gamma*y(8) + k10*(y(10)- 3*y(8))
        dy(9) = k7gamma*(sum_y57 - y(9)) + k6gamma2*(y(5) - y(9)) +k10*(-4*y(9) + 2*y(12))
        dy(10)=  k6gamma2*y(5) + k6gamma2*y(6) + k7gamma*(y(7) + 3*y(8) - 2*y(10))  + k10*(2*y(12)- 4*y(10))
        dy(11)= k7gamma*(y(6)- 2*y(11)) +  k6gamma2*y(7)  + k10*(y(12)- 4*y(11))
        dy(12)= k6gamma2*y(9) +k7gamma*(2*sum_y911-y(9) - y(12))  + k10*(6*y(13)- 5*y(12))
        dy(13)= k7gamma*y(12) - 6*k10*y(13)
    end subroutine dy_CaMKII

    pure real*8 function CaMKIIpho_func(y)
        implicit none
        real*8, intent(in) :: y(13)
        real*8 :: sum_y24, sum_y58, sum_y911
        sum_y24=sum(y(2:4))
        sum_y58=sum(y(5:8))
        sum_y911=sum(y(9:11))
        CaMKIIpho_func = y(1) + 2*sum_y24 + 3*sum_y58 + 4*sum_y911 + 5*y(12) + 6*y(13)
    end function CaMKIIpho_func

    subroutine d_PP1_I1P(PP1, I1P, CaM, pars,   dPP1, dI1P)
        implicit none
        real*8, intent(in) :: PP1, I1P, CaM
        type(pars_type), intent(in) :: pars
        real*8, intent(out) :: dPP1, dI1P
        type(post_CaMKII_plast_type) :: p
        real*8 :: vPKA, vCaN, k11, km11
        p = pars%post_CaMKII_plast
        vPKA = p%kpka0I1 + p%kpkaI1/(1 + (p%KdpkaI1/CaM)**p%npkaI1)
        vCaN = p%kcan0I1 + p%kcanI1/(1 + (p%KdcanI1/CaM)**p%ncanI1)
        k11=p%k11
        km11=p%km11
        dPP1= -k11*I1P*PP1 + km11*(p%PP10 - PP1) !RHS PP1
        dI1P= dPP1 + vPKA*p%I10 - vCaN*I1P !RHS I1P
    end subroutine d_PP1_I1P

end module CaMKII_plast
