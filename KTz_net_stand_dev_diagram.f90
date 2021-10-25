	  module ktz_LxL
	  implicit none
	  
	  contains
	  
		  ! Simulates de network and returns the column average of the standard deviation of each column (averaged by n realizations):
		  function sigma_x_avg(L,J) result(sigma_x)
		!  use ifport !use to compile with ifort	  
		  implicit none

			  ! Declaration of variables:
			  integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran
			  integer, intent(in) :: L
			  integer :: N, b, a, t_max, t_trans, t1, q, qo, dt, St, a_o, realiz, realiz_max
			  real*8, intent(in) :: J
			  real*8 :: k, T, del, lamb, Xr, Ie, sum_realiz, sum_sigma, sigma_x
			  real*8, allocatable, dimension(:) :: Zo, Yo, Xo, Z, Y, X, I, I_stim, arg
			  real*8, allocatable, dimension(:) :: sum_x, sum_x2, sigma
			  real*8, allocatable, dimension(:,:) :: rand_y, rand_x, O, G, C

			  ! Random seed:		  		  
		!	  call srand(666) !use to compile with ifort	  
			  call srand(seed) !use to compile with gfortran

			  ! Parameters for the KTz model:			  			  
			  k = 0.6d0
			  T = 0.154d0
			  del = 0.001d0
			  lamb = 0.001d0
			  Xr = -0.48d0
			  Ie = 0.0d0
			  
			  ! Network parameters:	  
			  N = L**2
			  dt = 92 !Pacing period (P in the paper)
			  
			  !Total time steps and transient:
			  t_max = 100000
			  t_trans = 50000
			  
			  !Number of realizations:
			  realiz_max = 10
			  
			  ! Allocate memory:
			  allocate(rand_y(N,realiz_max))
			  allocate(rand_x(N,realiz_max))
			  allocate(Xo(N))
			  allocate(Yo(N))
			  allocate(Zo(N))
			  allocate(X(N))
			  allocate(Y(N))
			  allocate(Z(N))
			  allocate(I(N))
			  allocate(I_stim(N))
			  allocate(arg(N))
			  allocate(O(N,4))
			  allocate(G(N,4))
			  allocate(C(N,4))
			  ! SD variables:
			  allocate(sum_x(L))
			  allocate(sum_x2(L))
			  allocate(sigma(L))
			  rand_y = 0.0d0
			  rand_x = 0.0d0
			  Xo = 0.0d0
			  Yo = 0.0d0
			  Zo = 0.0d0
			  X = 0.0d0
			  Y = 0.0d0
			  Z = 0.0d0
			  I = 0.0d0
			  I_stim = 0.0d0
			  arg = 0.0d0
			  O = 0.0d0
			  G = 0.0d0
			  C = 0.0d0
			  sigma = 0.0d0
			  
			  ! Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
			  do realiz = 1, realiz_max
				  do a = 1, L
			!		rand_y(a*L,realiz_max)) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
			!		rand_x(a*L,realiz_max)) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
					rand_y(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
					rand_x(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7			
					if (realiz == 1) then
						I_stim(a*L) = 0.1d0
					end if
					! Use to stimulate more columns (count from right to left):
			!		 do a_o = 1, 4
			!			rand_y(a*L-a_o,realiz)) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
			!			rand_x(a*L-a_o,realiz)) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
			!			rand_y(a*L,realiz)) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran			
			!			rand_x(a*L,realiz)) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7
			!			if (realiz == 1) then
			!				I_stim(a*L-a_o) = 0.1d0
			!			end if	
			!		 end do		
				  end do
			  end do
			  
			  sum_realiz = 0.0d0
			  do realiz = 1, realiz_max
				  sum_sigma = 0.0d0
				  ! Iterate time:			 		
				  do t1 = 1, t_max
					 ! Initial conditions for the KTz model:			 
					 if (t1 == 1) then
						 St = 1
						 Zo = 0.0d0
						 Yo = -0.5d0
						 Xo = -0.5d0
						 ! Random IC and stimulus:
						 Yo = Yo + rand_y(:,realiz)
						 Xo = Xo + rand_x(:,realiz)
						 I = I_stim 				 
					 ! Update variables of the KTz model:
					 else if (t1 > 1) then
						 St = St + 1
						 Zo = Z
						 Yo = Y
						 Xo = X
						 ! Periodic stimulus:
						 if (St <= 10) then !Adjust the duration of the stimulus
							 I = I + I_stim
						 end if
						 if (St == dt) then
							St = 0
						 end if
					 end if
					 ! Calculate the KTz equations:			 
					 arg = (Xo - k * Yo + Zo + Ie + I) / T
					 X = arg / (1.0d0 + dabs(arg)) !Logistic KTz
			!		 X = tanh(arg)                 !Hyperbolic KTz
					 Y = Xo
					 Z = (1.0d0 - del) * Zo - lamb * (Xo - Xr)
					 ! Calculate the column average SD at time t:
					 if (t1 >= t_trans) then
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
						 ! Average by column gives the network SD at time t, sum for the time average:
						 sum_sigma = sum_sigma + sum(sigma)/(L*1.0d0)
					 end if				 
					 ! Matrix filled with the differences between potentials of first nearest neighbors:
					 O(1,1) = X(2) - X(1)
					 O(1,2) = 0.0d0
					 O(1,3) = X(1+L) - X(1)
					 O(1,4) = 0.0d0
					 do a = 2, L
						O(a,1) = X(a+1) - X(a)
						O(a,2) = X(a-1) - X(a)
						O(a,3) = X(a+L) - X(a)
						O(a,4) = 0.0d0
					 end do
					 do a = (L+1), (N-L)
						O(a,1) = X(a+1) - X(a)
						O(a,2) = X(a-1) - X(a)
						O(a,3) = X(a+L) - X(a)
						O(a,4) = X(a-L) - X(a)
					 end do
					 do a = ((N-L)+1), (N-1)
						O(a,1) = X(a+1) - X(a)
						O(a,2) = X(a-1) - X(a)
						O(a,3) = 0.0d0
						O(a,4) = X(a-L) - X(a)
					 end do
					 O(N,1) = 0.0d0
					 O(N,2) = X(N-1) - X(N)
					 O(N,3) = 0.0d0
					 O(N,4) = X(N-L) - X(N)				 
					 do a = 1, (L-1)
						O(a*L,1) = 0.0d0
						O(a*L+1,2) = 0.0d0
					 end do		 
					 ! Diffusive coupling:
					 G = J
					 C = G * O
					 I = sum(C, DIM=2)
				  end do
				  ! Time average of the network SD, sum for the realizations average:
				  sum_realiz = sum_realiz + sum_sigma/((t_max-t_trans)*1.0d0)
			  end do
			  
			  ! Average by realizations
			  sigma_x = sum_realiz/(realiz_max*1.0d0)

		  end function
		  
	  end module
	  
	  program KTz_net_stand_dev_diagram
	  use ktz_LxL
	  implicit none
	  
	  !Iterates the coupling and network size and calls the function to calculate
	  !the standard deviation.	  
	  
		integer :: L, L_max, n_J_max, n_J
		real*8 :: sigma_x, J, J_0, J_f, dJ
				  
		! Diffusive coupling constant:		  
		J_0 = 0.01d0
		J_f = 0.04d0
		dJ = 0.001d0

		! Discretize J variable:
		n_J_max = nint((J_f - J_0) / dJ)
		
		! Maximum network size:
		L_max = 20

		! Open output file:								
		open(unit=10, file="diag_sigma.dat", status="replace")
		
		! Iterates the network size and coupling and calls the previous function:
		do L = 2, L_max
		  J = J_0 - dJ 			  
		  do n_J = 1, n_J_max
			  J = J + dJ	
			  sigma_x = sigma_x_avg(L,J)
			  write(10,*) L, J, sigma_x
		  end do
		end do
		
		! Close output file:		  		
		close(unit=10)
		
		! Gnuplot script
		! Plot the standard deviation against coupling J:
		open(unit=11, file="diag_sigma.plt", status="replace")
		  write(11, '(a)') "set terminal wxt"
		  write(11, '(a)') "set cbrange[0:0.3]"
		  write(11, '(a)') "set yrange[0.01:0.04]"
		  write(11, '(a)') "set xrange[2:20]"
		  write(11, '(a)') "set cblabel '<~{/Symbol s}{0.4-}_{x}>' font ', 30' offset -3,0,0"
		  write(11, '(a)') "set xlabel 'L' font ', 30' offset 3,1.4,0"
		  write(11, '(a)') "set ylabel 'J' font ', 30' offset 4,0,0"
		  write(11, '(a)') "set cbtics 0.1 font ', 20'"
		  write(11, '(a)') "set ytics 0.01 font ', 20'"
		  write(11, '(a)') "set xtics 20 font ', 20'"
		  write(11, '(a)') "set palette"
		  write(11, '(a)') "plot 'diag_sigma.dat' u 1:2:3 with image"
		close(unit=11)
		  
		! Call Gnuplot script for plot:		  
		call system("gnuplot -p diag_sigma.plt")
		
	  end program
