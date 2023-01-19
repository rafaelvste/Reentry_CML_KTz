    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
        
        !KTz model
        subroutine ktz_model(X_o, Y_o, Z_o, I_t, X_f, Y_f, Z_f)
            implicit none
        
			! Declaration of variables:
            real(kind=8), intent(in) :: X_o, Y_o, Z_o, I_t
            real(kind=8), intent(out) :: X_f, Y_f, Z_f
			real(kind=8) :: arg, k, T, del, lamb, Xr, Ie
    
			! KTz model parameters:
			k = 0.6d0
			T = 0.154d0
			del = 0.001d0
			lamb = 0.001d0
			Xr = -0.48d0
			Ie = 0.0d0
	
			! KTz equations:
            arg = (X_o - k * Y_o + Z_o + Ie + I_t) / T
            X_f = arg / (1.0d0 + dabs(arg)) !Logistic KTz
        !  	X_f = tanh(arg) 				  !Hyperbolic KTz
            Y_f = X_o
            Z_f = (1.0d0 - del) * Z_o - lamb * (X_o - Xr)
            
        end subroutine ktz_model
		        
    end module funcs_subrs
    
    ! Main program:
    program KTz_cell_bifurc_APD
        use funcs_subrs
        implicit none
        
		!Calculates the action potential duration (APD) bifurcation diagram for the single
		!KTz cell. APD is measured in the time series for different pacing periods P.
		!All the durations found are plotted in the corresponding P.
    
        ! Declaration of variables:
        integer :: dt, t1, St, t_max, stim_dur, t_trans, t_apd, APD
        real(kind=8) :: Zo, Yo, Xo, Z, Y, X, I, I_stim

        ! Additional parameters:
        t_max = 1000000 !Total time steps
		t_trans = 950000 !Transient		
        I_stim = 0.1d0 !Stimulus current value
		stim_dur = 10 !Stimulus duration
		
		! Open output file:
		open(unit=11, file="bifurc_cell.dat", status="replace")

		! Iterate the pacing period:
		do dt = 11, 350 !Pacing period (P in the paper)
			t_apd = 0 !Store the start of the APD interval
			! Iterate time:
			do t1 = 1, t_max
				! Initial conditions:
				if (t1 == 1) then
					Xo = 0.0d0
					Yo = 0.0d0
					Zo = 0.0d0
					I = I_stim 
					St = 1 !Count the time steps for stimulation
				! Update variables:
				else if (t1 > 1) then
					Xo = X
					Yo = Y
					Zo = Z
					I = 0.0d0
					St = St + 1
					! Periodic stimulus:
					if (St <= stim_dur) then
						I = I_stim
					end if
					if (St == dt) then
						St = 0
					end if
				end if
				! Calculate the KTz equations:
				call ktz_model(Xo, Yo, Zo, I, X, Y, Z)
				! Calculate the APD:
				if (t1 >= t_trans) then
					if ((X*Xo < 0.0d0) .and. (X > Xo)) then
						t_apd = t1
					else if ((X*Xo < 0.0d0) .and. (X < Xo)) then
						if (t_apd > 0) then
							APD = t1 - t_apd
							t_apd = 0
							write(11,*) dt, APD
						end if
					end if
				end if
			end do
		end do

        ! Close output files:
        close(unit=11)
		
		! Gnuplot script
		! Plot the APD bifucation diagram against pacing P:
		open(unit=12, file="bifurc_cell.plt", status="replace")
			write(12, '(a)') "set terminal wxt"
			write(12, '(a)') "set grid"
			write(12, '(a)') "unset key"		   
			write(12, '(a)') "set tics font ', 15'"
			write(12, '(a)') "set pointsize 0.4"
			write(12, '(a)') "set xrange[11:350]"
			write(12, '(a)') "set yrange[0:65]"		   
			write(12, '(a)') "set xlabel 'P' font ', 20'"
			write(12, '(a)') "set ylabel 'APD' font ', 20'"
			write(12, '(a)') "plot 'bifurc_cell.dat' pointtype 21 lt 8"
		close(unit=12)

		call system("gnuplot -p bifurc_cell.plt")

    end program KTz_cell_bifurc_APD