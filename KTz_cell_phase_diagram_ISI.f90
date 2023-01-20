    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
 
        !KTz model
        subroutine ktz_model(Xr, T, X_o, Y_o, Z_o, X_f, Y_f, Z_f)
            implicit none
        
            ! Declaration of variables:
            real(kind=8), intent(in) :: Xr, T, X_o, Y_o, Z_o
            real(kind=8), intent(out) :: X_f, Y_f, Z_f
            real(kind=8) :: arg, k, del, lamb, Ie
    
            ! KTz model parameters:
            k = 0.6d0
            del = 0.001d0
            lamb = 0.001d0
            Ie = 0.0d0
    
            ! KTz equations:
            arg = (X_o - k * Y_o + Z_o + Ie) / T
            X_f = arg / (1.0d0 + dabs(arg)) !Logistic KTz
        !   X_f = tanh(arg)                   !Hyperbolic KTz
            Y_f = X_o
            Z_f = (1.0d0 - del) * Z_o - lamb * (X_o - Xr)
            
        end subroutine ktz_model

        ! Iterate time and determine maximum and minimum ISI for a point (Xr,T)
        subroutine measure_ISI(Xr, T, max_ISI, min_ISI, n_ISIs)
            implicit none
            
            ! Declaration of variables:
            real(kind=8), intent(in) :: Xr, T
            integer, intent(out) :: max_ISI, min_ISI, n_ISIs
            integer :: t1, t_max, t_trans, n_ISIs_proxy, t_o, ISI
            integer, allocatable, dimension(:) :: V_ISI
            real(kind=8) :: Zo, Yo, Xo, Z, Y, X

            ! Additional parameters:
            t_max = 100000 !Total time steps
            t_trans = 50000 !Transient  

            ! Allocate memory:
            allocate(V_ISI(t_max-t_trans))
        
            ! Vector to store the series of ISIs found in the time series:
            allocate(V_ISI(t_max-t_trans)) !Number of ISIs always <= to the time steps considered
            V_ISI = -1
            
            n_ISIs = 0 !Count the number of ISIs in the time series
            t_o = 0 !Store the start of the ISI interval
            
            ! Iterate time:
            do t1 = 1, t_max
                ! Initial conditions:
                if (t1 == 1) then
                    Xo = 0.0d0
                    Yo = 0.0d0
                    Zo = 0.0d0
                ! Update variables:
                else if (t1 > 1) then
                    Xo = X
                    Yo = Y
                    Zo = Z
                end if
                ! Calculate the KTz equations:
                call ktz_model(Xr, T, Xo, Yo, Zo, X, Y, Z)
                ! Calculate the ISI:
                if (t1 >= t_trans) then
                    if ((X*Xo < 0.0d0) .and. (X < Xo)) then
                        if (t_o == 0) then
                            t_o = t1
                        else if (t_o > 0) then
                            n_ISIs = n_ISIs + 1
                            ISI = t1 - t_o
                            t_o = t1
                            V_ISI(n_ISIs) = ISI
                        end if
                    end if
                end if
            end do
            
            ! Get maximum and minimum ISI in the series of ISIs:
            if (n_ISIs /= 0) then !Only when ISIs were found
                max_ISI = maxval(V_ISI)
                do n_ISIs_proxy = 1, (t_max-t_trans)
                    if (V_ISI(n_ISIs_proxy) == -1) then
                        V_ISI(n_ISIs_proxy) = max_ISI
                    end if
                end do
                min_ISI = minval(V_ISI)
            else
                max_ISI = 0
                min_ISI = 0
            end if
            
        end subroutine measure_ISI
        
    end module funcs_subrs
    
    ! Main program:
    program KTz_cell_phase_diagram_ISI_m
        use funcs_subrs
        implicit none
        
        !Calculates the phase diagram of the KTz model. The interspike interval (ISI)
        !is measured in the time series for each point (Xr,T). The maximum and minimum
        !values are obtained. The maximum ISI for each (Xr,T) is plotted associated
        !with a color scale
    
        ! Declaration of variables:
        integer :: n_T, n_T_max, n_Xr, n_Xr_max, number_ISIs, ISI_max, ISI_min
        real(kind=8) :: T_p, T_0, T_f, Xr_p, Xr_0, Xr_f, d_T, dXr
        
        ! Start and end of the parameter T range:
        T_0 = 1.0e-8 
        T_f = 0.6d0
        
        ! Start and end of the parameter Xr range:
        Xr_0 = -0.6d0
        Xr_f = 0.0d0
        
        ! Interval between points:
        d_T = 0.01d0
        dXr = 0.01d0
          
        ! Number of points in the T axis of the diagram:
        n_T_max = nint((T_f-T_0)/d_T)

        ! Number of points in the Xr axis of the diagram:
        n_Xr_max = nint((Xr_f-Xr_0)/dXr)

        ! Open output files:
        open(unit=11, file="diag_ISI_max_m.dat", status="replace")
    !   open(unit=12, file="diag_ISI_min.dat", status="replace")

        ! Iterate paramater Xr:
        Xr_p = Xr_0 - dXr
        do n_Xr = 1, n_Xr_max
            Xr_p = Xr_p + dXr
            ! Iterate parameter T:
            T_p = T_0 - d_T 
            do n_T = 1, n_T_max
                T_p = T_p + d_T
                ! Determine maximum and minimum ISI in the time series:
                call measure_ISI(Xr_p, T_p, ISI_max, ISI_min, number_ISIs)
                ! Write to files (only when ISIs were found in the time series):
                if (number_ISIs /= 0) then
                    write(11,*) T_p, Xr_p, ISI_max
    !               write(12,*) T_p, Xr_p, ISI_min
                end if
            end do
        end do

        ! Close output files:
        close(unit=11)
        close(unit=12)
        
        ! Gnuplot scripts
        ! Plot the phase diagram with the maximum ISI in a png file:
        open(unit=13, file="diag_ISI_max.plt", status="replace")
            write(13, '(a)') "set terminal png enhanced size 1920,1080"
            write(13, '(a)') "set size square"
            write(13, '(a)') "unset key"
            write(13, '(a)') "set tics font ',20'"
            write(13, '(a)') "set cbrange[0:5000]"
            write(13, '(a)') "set yrange[-0.6:0]"
            write(13, '(a)') "set xrange[0:0.6]"
            write(13, '(a)') "set xlabel 'T' font ',30' offset '3,0,0'"
            write(13, '(a)') "set ylabel 'x_{r}' font ',30'"
            write(13, '(a)') "set cblabel 'ISI_{max}' font ',30'"
            write(13, '(a)') "set output 'diag_ISI_max_m.png'"
            write(13, '(a)') "plot 'diag_ISI_max.dat' using 1:2:3 with points pointtype 7 pointsize 4 palette"
        close(unit=13)

        call system("gnuplot -p diag_ISI_max_m.plt")

    end program KTz_cell_phase_diagram_ISI