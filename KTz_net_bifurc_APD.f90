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
		
		!Iterate time and measure APD for different coupling J
		subroutine measure_APD(L, N, J, I_stim, rand_y, rand_x)
			implicit none

			! Declaration of variables:
			integer, intent(in) :: L, N
			real(kind=8), intent(in) :: J
			real(kind=8), dimension(N), intent(in) :: I_stim, rand_y, rand_x
			integer :: cell_sample, a, t_max, t_trans, stim_dur, dt, t1, St, t_apd, APD
			real(kind=8), dimension(N) :: Zo, Yo, Xo, Z, Y, X, I, I_snps

			! Cell for measures:
			cell_sample = (N-L) / 2
			
			! Additional parameters:
			t_max = 100000 !Total time steps
			t_trans = 90000 !Transient
			dt = 92 !Pacing period (P in the paper)
			stim_dur = 10 !Stimulus duration
				
			t_apd = 0 !Store the start of the APD interval
			
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
				! Calculate the APD:
				if (t1 >= t_trans) then
					if ((X(cell_sample)*Xo(cell_sample) < 0.0d0) .and. (X(cell_sample) > Xo(cell_sample))) then
						t_apd = t1
					else if ((X(cell_sample)*Xo(cell_sample) < 0.0d0) .and. (X(cell_sample) < Xo(cell_sample))) then
						if (t_apd > 0) then
							APD = t1 - t_apd
							t_apd = 0
							write(10,*) J, APD
						end if
					end if
				end if
			end do
			
		end subroutine measure_APD
		
    end module funcs_subrs
	
	! Main program
	program KTz_net_bifurc_APD
		use funcs_subrs
	!	use ifport !use to compile with ifort
		implicit none
		
		!Simulates the lattice of diffusively coupled KTz cells and calculates the action potential
		!duration (APD) bifurcation diagram for one of the cells of the network. APD is 
		!measured in the time series for different couplings J. All the durations found are 
		!plotted in the corresponding J.
		
		! Declaration of variables:
		integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran		  		  		
		integer :: L_side, N_cells, a
		integer :: n_J_max, n_J, realiz_max, realiz
		real(kind=8) :: J_p, J_0, J_f, dJ
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
		L_side = 2
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
		open(unit=10, file="bifurc_net.dat", status="replace")

		! Iterate parameter J:
		do n_J = 1, n_J_max
			J_p = J_p + dJ
			! Iterate realizations:
			do realiz = 1, realiz_max
				! Iterate time and calculate the APD:
				call measure_APD(L_side, N_cells, J_p, I_stimulus, random_y(:,realiz), random_x(:,realiz))
			end do
		end do
		
		! Close output files:
		close(unit=10)
		
		! Gnuplot script
		! Plot the APD bifucation diagram against coupling J:
		open(unit=10, file="bifurc_net.plt", status="replace")
			write(10, '(a)') "set terminal wxt"
			write(10, '(a)') "set grid"
			write(10, '(a)') "unset key"		   
			write(10, '(a)') "set tics font ', 15'"
			write(10, '(a)') "set pointsize 0.4"
			write(10, '(a)') "set xrange[0.01:0.04]"
			write(10, '(a)') "set yrange[0:100]"		   
			write(10, '(a)') "set xlabel 'J' font ', 20'"
			write(10, '(a)') "set ylabel 'APD' font ', 20'"
			write(10, '(a)') "plot 'bifurc_net.dat' pointtype 21 lt 8"
		close(unit=10)

		call system("gnuplot -p bifurc_net.plt")
		
	end program KTz_net_bifurc_APD