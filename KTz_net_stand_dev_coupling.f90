    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
	
        !KTz model
        subroutine ktz_model(N, X_o, Y_o, Z_o, I_t, X_f, Y_f, Z_f)
            implicit none
        
			! Declaration of variables:
			integer, intent(in) :: N
            real(kind=8), dimension(N), intent(in) :: X_o, Y_o, Z_o, I_t
            real(kind=8), dimension(N), intent(out) :: X_f, Y_f, Z_f
			real(kind=8) :: k, T, del, lamb, Xr, Ie
			real(kind=8), dimension(N) :: arg
			
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
		
		!Diffusive coupling (electrical synapse)
		subroutine diff_coupl(L, N, J, X_s, I_s)
			implicit none
			
			! Declaration of variables:
			integer, intent(in) :: L, N
			real(kind=8), intent(in) :: J
			real(kind=8), dimension(N), intent(in) :: X_s
			real(kind=8), dimension(N), intent(out) :: I_s
			integer :: a
			real(kind=8), dimension(N,4) :: O, G, C
		
			! Matrix filled with the differences between potentials of first nearest neighbors
			O = 0.0d0
			do a = 1, N !Numbering from left to right and from bottom to top
				! Connection with the right neighbor:
				if (a < N) then
					O(a,1) = X_s(a+1) - X_s(a)
				end if
				! Connection with the left neighbor:
				if (a > 1) then
					O(a,2) = X_s(a-1) - X_s(a)
				end if
				! Connection with the neighbor above:
				if (a <= (N-L)) then
					O(a,3) = X_s(a+L) - X_s(a)
				end if
				! Connection with the neighbor below:
				if (a > L) then
					O(a,4) = X_s(a-L) - X_s(a)
				end if
			end do
			
			! Fix the remaining boundary conditions of the left and right columns:
			do a = 1, (L-1)
				O(a*L,1) = 0.0d0
				O(a*L+1,2) = 0.0d0
			end do
			
			! Diffusive coupling equation:
			G = J
			C = G * O
			I_s = sum(C, DIM=2)
			
		end subroutine diff_coupl
		
		!Iterate time and calculate the time average of the network standard deviation
		subroutine stand_dev_x(L, N, J, I_stim, rand_y, rand_x, avg_sigma)
			implicit none

			! Declaration of variables:
			integer, intent(in) :: L, N
			real(kind=8), intent(in) :: J
			real(kind=8), dimension(N), intent(in) :: I_stim, rand_y, rand_x
			real(kind=8), intent(out) :: avg_sigma
			integer :: a, a_o, t_max, t_trans, stim_dur, dt, t1, St
			real(kind=8) :: sum_sigma
			real(kind=8), dimension(N) :: Zo, Yo, Xo, Z, Y, X, I, I_snps, sum_x, sum_x2, sigma
			
			! Additional parameters:
			t_max = 100000 !Total time steps
			t_trans = 50000 !Transient
			dt = 92 !Pacing period (P in the paper)
			stim_dur = 10 !Stimulus duration
			
			sum_sigma = 0.0d0 !Store standard deviation values for time average
			
			! Iterate time:
			do t1 = 1, t_max
				! Initial conditions:
				if (t1 == 1) then
					Xo(:) = -0.5d0 + rand_x(:)
					Yo(:) = -0.5d0 + rand_y(:)
					Zo = 0.0d0
					! Stimulus:
					St = 1 !Count the time steps for stimulation
					I = I_stim
				! Update variables:
				else if (t1 > 1) then
					Xo = X
					Yo = Y
					Zo = Z
					! Periodic stimulus:
					I = 0.0d0
					St = St + 1
					if (St <= stim_dur) then
						I = I_stim
					end if
					if (St == dt) then
						St = 0
					end if
				end if
				! Diffusive coupling:
				call diff_coupl(L, N, J, Xo, I_snps)
				I = I + I_snps
				! Calculate the KTz equations:
				call ktz_model(N, Xo, Yo, Zo, I, X, Y, Z)
				! Calculate the column-based network standard deviation at time t:
				if (t1 > t_trans) then
					sum_x = 0.0d0
					sum_x2 = 0.0d0
					! Sum X and X**2 in each column:
					do a = 1, L !Column index
						do a_o = 0, (L-1)
							sum_x(a) = sum_x(a) + X(a+a_o*L)
							sum_x2(a) = sum_x2(a) + (X(a+a_o*L))**2
						end do
					end do
					! Obtain the standard deviation of each column:
					sigma = sqrt(dabs(sum_x2/(L*1.0d0) - (sum_x/(L*1.0d0))**2))
					! Average per column gives the network SD at time t, sum for time average:
					sum_sigma = sum_sigma + sum(sigma)/(L*1.0d0)
				end if
			end do
			
			! Time average:
			avg_sigma = sum_sigma/((t_max-t_trans)*1.0d0)
			
		end subroutine stand_dev_x
		
    end module funcs_subrs
	
	! Main program
	program KTz_net_stand_dev_coupling
		use funcs_subrs
	!	use ifport !use to compile with ifort
		implicit none
		
		!Simulates the network and calculates the column-based standard deviation.
		!Calculates the standad deviation of x in each column and the average of all columns.
		!Then averages all time steps in a realization and averages the realizations.
		!This is calculated for different values of coupling J.
		
		! Declaration of variables:
		integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran
		integer :: L_side, N_cells, a
		integer :: n_J_max, n_J, realiz_max, realiz
		real(kind=8) :: J_p, J_0, J_f, dJ, sum_realiz, sigma_avg, sigma_x
		real(kind=8), allocatable, dimension(:) :: I_stimulus
		real(kind=8), allocatable, dimension(:,:) :: random_y, random_x
		
		! Random seed:
	!	call srand(666) !use to compile with ifort	  
		call srand(seed) !use to compile with gfortran
		
		! Start and end of the parameter J (coupling) range:		  
		J_0 = 0.01d0
		J_f = 0.04d0
		dJ = 0.0001d0 !Interval between points

		! Discretize J parameter:
		n_J_max = nint((J_f - J_0) / dJ)
		J_p = J_0 - dJ
		
		! Network parameters:
		L_side = 10
		N_cells = L_side**2		
		
		!Number of realizations for each J:
		realiz_max = 10
	
		! Allocate memory:
		allocate(random_y(N_cells,realiz_max))
		allocate(random_x(N_cells,realiz_max))	
		allocate(I_stimulus(N_cells))
	
		! Generate random numbers for the IC for various realizations:
		random_y = 0.0d0
		random_x = 0.0d0
		do realiz = 1, realiz_max
			do a = 1, L_side
		!		random_y(a*L_side,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
		!		random_x(a*L_side,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
				random_y(a*L_side,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
				random_x(a*L_side,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7
			end do
		end do
		
		! Generate the stimulus vector (cells of the rightmost column):
		I_stimulus = 0.0d0
		do a = 1, L_side
			I_stimulus(a*L_side) = 0.1d0 !Stimulus current value
		end do
		
		! Open output files:
		open(unit=10, file="sigma_x_J.dat", status="replace")

		! Iterate parameter J:
		do n_J = 1, n_J_max
			J_p = J_p + dJ
			sum_realiz = 0.0d0
			! Iterate realizations:
			do realiz = 1, realiz_max
				! Iterate time and calculate the average standard deviation of the network:
				call stand_dev_x(L_side, N_cells, J_p, I_stimulus, random_y(:,realiz), random_x(:,realiz), sigma_avg)
				sum_realiz = sum_realiz + sigma_avg
			end do
			! Average standard deviation by realizations:
			sigma_x = sum_realiz / (realiz_max*1.0d0)
			write(10,*) J_p, sigma_x
		end do

		! Close output files:
		close(unit=10)
		
		! Gnuplot script
		! Plot the standard deviation against coupling J:
		open(unit=11, file="sigma_x_J.plt", status="replace")
			write(11, '(a)') "set terminal wxt"
			write(11, '(a)') "set grid"
			write(11, '(a)') "unset key"
			write(11, '(a)') "set ytics font ', 15'"
			write(11, '(a)') "set xtics 0.01 font ', 15'"
			write(11, '(a)') "set pointsize 0.6"
			write(11, '(a)') "set xrange[0.01:0.04]"
			write(11, '(a)') "set xlabel 'J' font ', 20'"
			write(11, '(a)') "set ylabel '<~{/Symbol s}{0.4-}_{x}>' font ', 20'"
			write(11, '(a)') "plot 'sigma_x_J.dat' w linespoints pt 21 lt 8 lw 1.5"
		close(unit=11)

		! Call Gnuplot script for plot:		  
		call system("gnuplot -p sigma_x_J.plt")
		
	end program KTz_net_stand_dev_coupling