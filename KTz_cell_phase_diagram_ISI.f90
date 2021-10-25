      program KTz_phase_diagram_ISI
      implicit none
			
	  !Calculates the phase diagram of the KTz model. The interspike interval (ISI)
	  !is measured in the time series for each point (Xr,T). The maximum and minimum
	  !values are obtained. The maximum ISI for each (Xr,T) is plotted associated
	  !with a color scale.
			
		  ! Declaration of variables:		  
		  integer :: t1, St, t_trans, t_max, t_o, n_T, n_T_max, n_Xr, n_Xr_max, ISI, k1, k2, min_ISI, max_ISI
		  integer, allocatable, dimension(:) :: V_ISI
		  real*8 :: k, T, T_0, T_f, del, lamb, Xr, Xr_0, Xr_f, Ie, d_T, dXr
		  real*8 :: Zo, Yo, Xo, Z, Y, X, I, arg

		  ! KTz model parameters:
		  k = 0.6d0
		  T_0 = 1.0e-8 
		  T_f = 0.6d0	
		  del = 0.001d0
		  lamb = 0.001d0
		  Xr_0 = -0.6d0
		  Xr_f = 0.0d0
		  Ie = 0.0d0
		  d_T = 0.01d0
		  dXr = 0.01d0
		  t_max = 100000 !Total time steps
		  t_trans = 50000 !Transient	  
		  
		  ! Number of points in the T axis of the diagram:
		  n_T_max = nint((T_f-T_0)/d_T)
		  
		  ! Number of points in the Xr axis of the diagram:
		  n_Xr_max = nint((Xr_f-Xr_0)/dXr)

		  ! Allocate memory:
		  allocate(V_ISI(t_max-t_trans)) !Store the ISIs found in a time series

		  ! Open output files:
		  open(unit=11, file="diag_ISI_max.dat", status="replace")
	!	  open(unit=12, file="diag_ISI_min.dat", status="replace")
		  
		  ! Iterate paramater Xr:
		  Xr = Xr_0 - dXr
		  do n_Xr = 1, n_Xr_max
			Xr = Xr + dXr
			! Iterate parameter T:
			T = T_0 - d_T	
			do n_T = 1, n_T_max
			  T = T + d_T
			  do k1 = 1, (t_max-t_trans)
				V_ISI(k1) = -1
			  end do
			  k2 = 0
			  t_o = 0				  
			  do t1 = 1, t_max
				 ! Initial conditions:
				 if (t1 == 1) then
					 Zo = 0.0d0
					 Yo = 0.0d0
					 Xo = 0.0d0
					 I = 0.0d0
				 ! Update variables:	
				 else if (t1 > 1) then
					 Zo = Z
					 Yo = Xo
					 Xo = X
				 end if
				 ! Calculate the KTz equations:
				 arg = (Xo - k * Yo + Zo + Ie + I) / T
				 X = arg / (1.0d0 + dabs(arg)) !Logistic KTz
!				 X = tanh(arg)				   !Hyperbolic KTz
				 Y = Xo
				 Z = (1.0d0 - del) * Zo - lamb * (Xo - Xr)	
				 ! Measure the ISI in the time series after transient:
				 if (t1 >= t_trans) then
					if ((X*Xo < 0.0d0) .and. (X < Xo)) then
						if (t_o == 0) then
							t_o = t1
						else if (t_o > 0) then
							k2 = k2 + 1
							ISI = t1 - t_o
							t_o = t1
							V_ISI(k2) = ISI
						end if
					end if
				 end if
			  end do
			  ! Obtain the maximum and minimum ISI for the point (Xr,T):
			  if (k2 /= 0) then
				  max_ISI = maxval(V_ISI)
				  do k1 = 1, (t_max-t_trans)
					  if (V_ISI(k1) == -1) then
						  V_ISI(k1) = max_ISI
					  end if
				  end do
				  min_ISI = minval(V_ISI)
				  ! Write to files:
				  write(11,*) T, Xr, max_ISI
	!			  write(12,*) T, Xr, min_ISI
			  end if
			end do
		  end do	
		
		  ! Close output files:
		  close(unit=11)
	!	  close(unit=12)
		  
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
			write(13, '(a)') "set output 'diag_ISI_max.png'"
			write(13, '(a)') "plot 'diag_ISI_max.dat' using 1:2:3 with points pointtype 7 pointsize 4 palette"
		  close(unit=13)
		  
		  call system("gnuplot -p diag_ISI_max.plt")
		  
      end program
	  	  
