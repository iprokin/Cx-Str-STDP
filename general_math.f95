!    -*- f95 -*-
! (c) 2016 - Ilya Prokin - isprokin@gmail.com - https://sites.google.com/site/ilyaprokin
! INRIA Rhone-Alpes
! STDP model : math functions
module general_math
    implicit none

    contains

    real*8 function heaviside(x)
        implicit none
        real*8 :: x
        if (x == 0) then
            heaviside = 0.5
        elseif (x < 0) then
            heaviside = 0
        else
            heaviside = 1.0
        endif
    end function heaviside

    real*8 function Rect(x, xlo, xhi)
        implicit none
        real*8 :: x, xlo, xhi
        if (x < xlo) then
            Rect = 0
        elseif (x > xhi) then
            Rect = 0
        else
            Rect = 1.0
        endif
    end function Rect

    real*8 function Hill(x,k,n)
        implicit none
        real*8 :: x,k
        integer :: n
        real*8 :: xn
        xn=x**n
        Hill = xn/(xn+k**n)
    end function Hill

    real*8 function lin_interpTable(Table,li,i_point, default)
        implicit none
        real*8 :: i_point
        integer :: ip, li
        real*8 :: Table(li), default
        ip=floor(i_point)
        if (ip<1) then
            lin_interpTable = default
        elseif (ip<li) then
            lin_interpTable = Table(ip)+(Table(ip+1)-Table(ip))*(i_point-ip)
        else
            lin_interpTable = default
        endif
    end function lin_interpTable

    real*8 function lin_interpTable_brds(Table,li,i_point)
        implicit none
        real*8 :: i_point
        integer :: ip, li
        real*8 :: Table(li)
        ip=floor(i_point)
        if (ip<1) then
            lin_interpTable_brds = Table(1)
        elseif (ip<li) then
            lin_interpTable_brds = Table(ip)+(Table(ip+1)-Table(ip))*(i_point-ip)
        else
            lin_interpTable_brds = Table(li)
        endif
    end function lin_interpTable_brds

    subroutine search_xa_inds_around_x(xa,n,x, klo, khi)
        implicit none
        integer, intent(in) :: n
        real*8, intent(in) :: xa(n), x
        integer, intent(out) :: khi,klo
        integer :: k
        klo=1
        khi=n
        do
            if (khi-klo <= 1) exit
            k=(khi+klo)/2
            if(xa(k) > x) then
                khi=k
            else
                klo=k
            endif
        end do
    end subroutine search_xa_inds_around_x

    real*8 function lin_interpTable_xy(x,y,n, x_new)
        implicit none
        integer :: n, klo, khi
        real*8 :: x(n), y(n), x_new
        if (x_new<=x(1)) then
            lin_interpTable_xy = y(1)
        elseif (x_new>=x(n)) then
            lin_interpTable_xy = y(n)
        else
            call search_xa_inds_around_x(x,n,x_new, klo, khi)
            if (y(khi)/=y(klo)) then
                lin_interpTable_xy = y(klo)+(y(khi)-y(klo))*(x_new-x(klo))/(x(khi)-x(klo))
            else
                lin_interpTable_xy = y(klo)
            end if
        endif
    end function lin_interpTable_xy

    real*8 function kernel_repeat(t, Freq, n, kernel, l_ker, step_ker)
        !t -grid starts from 0.0 (t=t'-t0)
        implicit none
        real*8 :: t, Freq
        integer :: n, l_ker
        real*8 :: kernel(l_ker), step_ker
        integer :: ic, idepth, k
        ic = floor(t*Freq)
        idepth = floor(l_ker*step_ker*Freq)
        kernel_repeat = 0
        if (t>l_ker*step_ker+(n-1)/Freq) then
            return
        end if
        do k=0,idepth
            kernel_repeat = kernel_repeat + lin_interpTable(kernel, l_ker, (t-(k+ic)/Freq)/step_ker, 0d0)
        end do
    end function kernel_repeat

end module general_math
