      program KTz_cell
      implicit none
	  
	  !Simulates the single KTz cell. An external periodic current pulse
	  !of adjustable duration is also implemented. Outputs the membrane
	  !potential x(t), external current I(t) and the first return map.
	  
		  ! Declaration of variables:
		  integer :: dt, t1, St, t_max
		  real*8 :: k, T, del, lamb, Xr, Ie, I_stim
		  real*8 :: Zo, Yo, Xo, Z, Y, X, I, arg

		  ! KTz model parameters:
		  k = 0.6d0
		  T = 0.154d0
		  del = 0.001d0
		  lamb = 0.001d0
		  Xr = -0.48d0
		  Ie = 0.0d0
		  I_stim = 0.1d0 !Stimulus current value
		  t_max = 10000 !Total time steps
		  dt = 232 !Pacing period (P in the paper)

		  ! Open output files:
		  open(unit=11, file="potential.dat", status="replace")
		  open(unit=12, file="current.dat", status="replace")
		  open(unit=13, file="1st_return.dat", status="replace")

		  ! Iterate time:
		  do t1 = 1, t_max
			 ! Initial conditions:
			 if (t1 == 1) then
				 Zo = 0.0d0
				 Yo = -0.5d0
				 Xo = -0.5d0
				 I = I_stim 
				 St = 1
			 ! Update variables:
			 else if (t1 > 1) then
				 Zo = Z
				 Yo = Xo
				 Xo = X
				 I = 0.0d0
				 St = St + 1
				 ! Periodic stimulus:
				 if (St <= 10) then !Adjust the duration of the stimulus
					I = I_stim
				 end if
				 if (St == dt) then
					St = 0
				 end if
			 end if
			 ! Output membrane potential and current:
			 write(11,*) t1, Xo
			 write(12,*) t1, I
			 ! Calculate the KTz equations:
			 arg = (Xo - k * Yo + Zo + Ie + I) / T
			 X = arg / (1.0d0 + dabs(arg)) !Logistic KTz
		!	 X = tanh(arg) 				   !Hyperbolic KTz
			 Y = Xo
			 Z = (1.0d0 - del) * Zo - lamb * (Xo - Xr)			 
			 ! Output 1st return map (x(t),x(t+1)):
			 if (t1 >= 2000) then !Ignore transient
				write(13,*) Xo, X
			 end if
		  end do		 

	      ! Close output files:
		  close(unit=11)
		  close(unit=12)
		  close(unit=13)
		  
		  ! Gnuplot scripts
		  ! Plot the membrane potential x(t) and stimulus current I(t):
		  open(unit=14, file="potential.plt", status="replace")
			  write(14, '(a)') "set terminal wxt"
			  write(14, '(a)') "set multiplot layout 2,1"
			  write(14, '(a)') "set grid"
			  write(14, '(a)') "unset key"			  
			  write(14, '(a)') "set ytics font ', 12'"
			  write(14, '(a)') "set xtics font ', 12'"			  
			  write(14, '(a)') "unset xlabel"
			  write(14, '(a)') "set ylabel 'x(t)' font ', 16'"
			  write(14, '(a)') "set xrange[8100:8700]"
			  write(14, '(a)') "set size 1,0.7"
			  write(14, '(a)') "set origin 0,0.3"
			  write(14, '(a)') "plot 'potential.dat' w linespoints pointtype 21 lt 6"
			  write(14, '(a)') "set ytics font ', 12' 0.1"
			  write(14, '(a)') "set xlabel 't' font ', 16'"
			  write(14, '(a)') "set ylabel 'I(t)' font ', 16'"
			  write(14, '(a)') "set size 1,0.35"
			  write(14, '(a)') "plot 'current.dat' w linespoints pointtype 21 lt 7"
		  close(unit=14)
		  
		  call system("gnuplot -p potential.plt")

		  ! Plot the 1st return map (x(t),x(t+1)):
		  open(unit=15, file="1st_return.plt", status="replace")
			  write(15, '(a)') "set term wxt"
			  write(15, '(a)') "set grid"		  
			  write(15, '(a)') "set tics font ', 12'"			  
			  write(15, '(a)') "set xlabel 'x(t)' font ', 16'"
			  write(15, '(a)') "set ylabel 'x(t+1)' font ', 16'"
			  write(15, '(a)') "plot '1st_return.dat' pt 21 lc 8"
		  close(unit=15)	
		  
		  call system("gnuplot -p 1st_return.plt")

      end program
	  	  