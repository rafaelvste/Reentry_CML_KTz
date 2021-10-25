      program KTz_cell_bifurc_APD
      implicit none
	  
	  !Calculates the action potential duration (APD) bifurcation diagram for the single
	  !KTz cell. APD is measured in the time series for different pacing periods P.
	  !All the durations found are plotted in the corresponding P.
		  
		  ! Declaration of variables:
		  integer :: dt, t1, St, t_max, t_trans, t_apd
		  real*8 :: k, T, del, lamb, Xr, Ie, I_stim, APD
		  real*8 :: Zo, Yo, Xo, Z, Y, X, I, arg
			
		  ! KTz model parameters:
		  k = 0.6d0
		  T = 0.154d0
		  Ie = 0.0d0
		  del = 0.001d0
		  lamb = 0.001d0
		  Xr = -0.48d0
		  I_stim = 0.1d0 !Stimulus current value		  
		  t_max = 1000000 !Total time steps
		  t_trans = 950000 !Transient
		  
		  ! Open output file:
		  open(unit=11, file="bifurc_cell.dat", status="replace")		  

		  ! Iterate the pacing period:
		  do dt = 11, 350 !Variable P in the paper
			  t_apd = 0
			  ! Iterate time:
			  do t1 = 1, t_max
				 ! Initial conditions:
				 if (t1 == 1) then
					 Zo = 0.0d0
					 Yo = 0.0d0
					 Xo = 0.0d0
					 I = I_stim
					 St = 1
				 else if (t1 > 1) then
					 ! Update variables:
					 Zo = Z
					 Yo = Y
					 Xo = X
					 I = 0.0d0
					 ! Periodic stimulus:
					 St = St + 1
					 if (St <= 10) then !Adjust the duration of the stimulus
						I = I_stim
					 end if
					 if (St == dt) then
						St = 0
					 end if
				 end if
				 ! Calculate the KTz equations:
				 arg = (Xo - k * Yo + Zo + Ie + I) / T
				 X = arg / (1.0d0 + dabs(arg)) !Logistic KTz
!				 X = tanh(arg)				   !Hyperbolic KTz				 
				 Y = Xo
				 Z = (1.0d0 - del) * Zo - lamb * (Xo - Xr)
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

		 ! Close output file:
	     close(unit=11)
	  
		 ! Gnuplot script
		 ! Plot the APD bifucation diagram against pacing P:
		 open(unit=11, file="bifurc_cell.plt", status="replace")
		   write(11, '(a)') "set terminal wxt"
		   write(11, '(a)') "set grid"
		   write(11, '(a)') "unset key"		   
		   write(11, '(a)') "set tics font ', 15'"
		   write(11, '(a)') "set pointsize 0.4"
		   write(11, '(a)') "set xrange[11:350]"
		   write(11, '(a)') "set yrange[0:65]"		   
		   write(11, '(a)') "set xlabel 'P' font ', 20'"
		   write(11, '(a)') "set ylabel 'APD' font ', 20'"
		   write(11, '(a)') "plot 'bifurc_cell.dat' pointtype 21 lt 8"
		 close(unit=11)
		  
		 call system("gnuplot -p bifurc_cell.plt")
		  
      end program