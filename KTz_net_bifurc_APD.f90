      program KTz_net_bifurc_APD
    !  use ifport !use to compile with ifort	  
      implicit none
	  
	  !Simulates the lattice of diffusively coupled KTz cells and calculates the action potential
	  !duration (APD) bifurcation diagram for one of the cells of the network. APD is 
	  !measured in the time series for different couplings J. All the durations found are 
	  !plotted in the corresponding J.
	  
          ! Declaration of variables:
		  integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran		  		  
		  integer :: L, N, b, a, t_max, t_trans, t1, q, qo, dt, St, a_o
		  integer :: n_J_max, n_J, t_apd, APD, realiz_max, realiz
		  real*8 :: k, T, del, lamb, Xr, Ie, J, J_0, J_f, dJ
		  real*8, allocatable, dimension(:) :: Zo, Yo, Xo, Z, Y, X, I, I_stim, arg
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

		  ! Diffusive coupling constant:		  
		  J_0 = 0.01d0
		  J_f = 0.04d0
		  dJ = 0.0001d0
		  
		  ! Discretize J variable:
		  n_J_max = nint((J_f - J_0) / dJ)
		  J = J_0 - dJ 		  
		  
		  ! Network parameters:	  
		  L = 2
		  N = L**2
		  dt = 92 !Pacing period (P in the paper)
		  
		  ! Cell for measures:
		  b = (N-L)/2
		  
		  !Total time steps and transient:
		  t_max = 100000
		  t_trans = 90000
		  
		  !Number of realizations for each J:
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

		  ! Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
		  do realiz = 1, realiz_max
			  do a = 1, L
	!			rand_y(a*L,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
	!			rand_x(a*L,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
				rand_y(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
				rand_x(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7
				if (realiz == 1) then
					I_stim(a*L) = 0.1d0
				end if
				! Use to stimulate more columns (count from right to left):
	!			 do a_o = 1, 4
	!				rand_y(a*L-a_o,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
	!				rand_x(a*L-a_o,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
	!				rand_y(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran			
	!				rand_x(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7
	!				if (realiz == 1) then
	!					I_stim(a*L-a_o) = 0.1d0
	!				end if					
	!			 end do					
			  end do
		  end do
		  	  
		  ! Open output files:
		  open(unit=10, file="bifurc_net.dat", status="replace")
		  
		  !Iterate parameter J:
		  do n_J = 1, n_J_max
			  J = J + dJ
			  do realiz = 1, realiz_max	  
				  t_apd = 0
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
						 ! Specify the duration of the periodic (period dt) stimulus:
						 if (St <= 10) then
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
					 ! Calculate the APD (analogous condition used to ISI in the KTz Log paper Plos One)
					 if (t1 >= t_trans) then
						if ((X(b)*Xo(b) < 0.0) .and. (X(b) > Xo(b))) then
							t_apd = t1
						else if ((X(b)*Xo(b) < 0.0) .and. (X(b) < Xo(b))) then
							if (t_apd > 0) then
								APD = t1 - t_apd
								t_apd = 0								
								write(10,*) J, APD
							end if
						end if
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
					 ! Fix the boundary conditions of the left and right columns:
					 O(1,2) = 0.0d0
					 O(N,1) = 0.0d0
					 do a = 1, (L-1)
						O(a*L,1) = 0.0d0
						O(a*L+1,2) = 0.0d0
					 end do					 
					 ! Diffusive coupling:
					 G = J
					 C = G * O
					 I = sum(C, DIM=2)
				  end do
			  end do
		  end do
		  
		  ! Close output files:
		  close(unit=10)

		  ! Gnuplot script
		  ! Plot the APD bifucation diagram against coupling J:
		  open(unit=11, file="bifurc_net.plt", status="replace")
		    write(11, '(a)') "set terminal wxt"
		    write(11, '(a)') "set grid"
		    write(11, '(a)') "unset key"		   
		    write(11, '(a)') "set tics font ', 15'"
		    write(11, '(a)') "set pointsize 0.4"
		    write(11, '(a)') "set xrange[0.01:0.04]"
		    write(11, '(a)') "set yrange[0:100]"		   
		    write(11, '(a)') "set xlabel 'J' font ', 20'"
		    write(11, '(a)') "set ylabel 'APD' font ', 20'"
		    write(11, '(a)') "plot 'bifurc_net.dat' pointtype 21 lt 8"
		  close(unit=11)
		  
		 call system("gnuplot -p bifurc_net.plt")
		  
      end program 
