	  program KTz_log_net_Lyap_coupling
       !   use ifport !use to compile with ifort
	  implicit none
	  
	  !Simulates the lattice of diffusively coupled KTz cells and calculates the Lyapunov spectrum.
	  !The exponents are approximated using the Eckmann-Ruelle method, where the Jacobian
	  !matrix is triangularized by LU decompostion at each iteration using Lapack/MKL routines.
	  !Outputs the main exponent against coupling J.

		  ! Declaration of variables:
		  integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran
		  integer :: L, N, b, a, t_max, t_trans, t1, q, qo, dt, St, a_o, u, v, u_o, v_o, u_f, info, n_J, n_J_max
		  real*8 :: k, T, del, lamb, Xr, Ie, J, J_0, J_f, dJ
		  integer, allocatable, dimension(:)  :: ipiv  
		  real*8, allocatable, dimension(:) :: rand_y, rand_x, Zo, Yo, Xo, Z, Y, X, I, I_stim, arg
		  real*8, allocatable, dimension(:) :: n_s, a0, a11, a12, a13, b11, diag_O2, lambda0, row, exp_lyap
		  real*8, allocatable, dimension(:,:) :: O, G, C, M, M_dgemm, M_0, Per, Low, O1, O2
		  real*8, allocatable, dimension(:,:,:) :: M_ii, M_ij
		  external :: DGEMM, DGETRF2 !External routines for LU decomposition
		  
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
		  
		  !Total time steps and transient:
		  t_max = 100000
		  t_trans = 1
		  
		  ! Allocate memory:
		  allocate(rand_y(N))
		  allocate(rand_x(N))		  		  
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
		  ! Lyapunov variables:
		  allocate(a0(N))
		  allocate(a11(N))
		  allocate(a12(N))
		  allocate(a13(N))
		  allocate(b11(N))
		  allocate(M_ii(3,3,N))
		  allocate(M_ij(3,3,N))
		  allocate(M_0(3,3))
		  allocate(M(3*N,3*N))
		  allocate(M_dgemm(3*N,3*N))		  
		  allocate(Per(3*N,3*N))
		  allocate(Low(3*N,3*N))
		  allocate(ipiv(3*N))
		  allocate(row(3*N))
		  allocate(O1(3*N,3*N))
		  allocate(O2(3*N,3*N))
		  allocate(diag_O2(3*N))
		  allocate(lambda0(3*N))
		  ! Main output:
		  allocate(exp_lyap(3*N))		  
		  
		  ! Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
		  do a = 1, L
	!		rand_y(a*L) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
	!		rand_x(a*L) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
			rand_y(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
			rand_x(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7			
			I_stim(a*L) = 0.1d0
			! Use to stimulate more columns (count from right to left):
	!		 do a_o = 1, 4
	!			I_stim(a*L-a_o) = 0.1d0
	!			rand_y(a*L-a_o) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
	!			rand_x(a*L-a_o) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
	!			rand_y(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran			
	!			rand_x(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7			
	!		 end do		
		  end do
		  
		  ! Number of synapses for each cell:
		  allocate(n_s(N))
		  n_s = 0.0d0
		  if (L == 2) then
			  n_s(:) = 2.0d0
		  end if
		  if (L > 2) then
			  n_s(1) = 2.0d0
			  do a = 2, (L-1)
				  n_s(a) = 3.0d0
			  end do
			  n_s(L) = 2.0d0
			  n_s(L+1) = 3.0d0
			  do a = (L+2), (N-L-1)
				  n_s(a) = 4.0d0
			  end do
			  n_s(N-L) = 3.0d0
			  n_s(N-L+1) = 2.0d0    
			  do a = ((N-L)+2), (N-1)
				  n_s(a) = 3.0d0
			  end do
			  n_s(N) = 2.0d0
		  end if
		  ! Fix the number of synapses of the left and right columns:
		  if (L > 3) then
			  do a = 2, (L-2)
				  n_s(a*L) = 3.0d0
				  n_s(a*L+1) = 3.0d0
			  end do
		  end if
		  
		  ! Open output file:
		  open(unit=10, file="exp_lyap1.dat", status="replace")

		  ! Iterate parameter J:
		  do n_J = 1, n_J_max
			  J = J + dJ
			  lambda0 = 0.0d0
			  ! Iterate time:			 		
			  do t1 = 1, t_max 
				 ! Initial conditions for the KTz model:			 
				 if (t1 == 1) then
					 St = 1
					 Zo = 0.0d0
					 Yo = -0.5d0
					 Xo = -0.5d0
					 I = 0.0d0
					 ! Random IC and stimulus:
					 Yo = Yo + rand_y
					 Xo = Xo + rand_x
					 I = I_stim 				 
					 ! Update variables of the KTz model:
					 else if (t1 > 1) then
						 St = St + 1
						 Zo = Z
						 Yo = Y
						 Xo = X
						 ! Periodic stimulus:
						 if (St <= 10) then
							 I = I + I_stim !Adjust the duration of the stimulus
						 end if
						 if (St == dt) then
							St = 0
						 end if
				 end if
				 ! Calculate the KTz equations:			 
				 arg = (Xo - k * Yo + Zo + Ie + I) / T
				 X = arg / (1.0d0 + dabs(arg))
				 Y = Xo
				 Z = (1.0d0 - del) * Zo - lamb * (Xo - Xr)
				 ! Lyapunov exponent (vectorized):
				 if (t1 >= t_trans) then
					a0 = 1.0d0 / (T*(1.0d0 + dabs(arg))**2)
					! Values in the single element Jacobian matrices M_ii (vectorized):
					a11 = (1.0d0 - n_s*J) * a0
					a12 = -k * a0
					a13 = a0
					! Single value in the cross elements Jacobian matrices M_ij (vectorized):
					b11 = J * a0
					! Place the values in the M_ii and M_ij tensors (the 3rd dimension refer to network cell):
					M_ii(1,1,:) = a11
					M_ii(1,2,:) = a12
					M_ii(1,3,:) = a13
					M_ii(2,1,:) = 1.0d0
					M_ii(3,1,:) = -lamb
					M_ii(3,3,:) = 1.0d0 - del
					M_ij(1,1,:) = b11
					! Fill the Jacobian matrix of the network with M_ii and M_ij:
					u_o = -2
					do u = 1, N
						u_o = u_o + 3
						v_o = -2
						do v = 1, N
							v_o = v_o + 3
							if (v == u) then
								M(u_o:u_o+2, v_o:v_o+2) = M_ii(:,:,u)
							else if ((v == (u-1)) .or. (v == (u+1))) then
								M(u_o:u_o+2, v_o:v_o+2) = M_ij(:,:,u)
							else if ((v == (u-L)) .or. (v == (u+L))) then
								M(u_o:u_o+2, v_o:v_o+2) = M_ij(:,:,u)
							else
								M(u_o:u_o+2, v_o:v_o+2) = M_0(:,:)
							end if
						end do
					end do
					! Fix the matrices of the left and right columns:
					u_f = 1
					do u = 1, (L-1)
						u_o = u_f + 3*(L-1)
						u_f = u_o + 3
						M(u_o:u_f-1, u_f:u_f+2) = M_0(:,:)
						M(u_f:u_f+2, u_o:u_f-1) = M_0(:,:)        
					end do
					! Perform the LU decomposition:
					M_dgemm = M
					if (t1 > t_trans) then
						call DGEMM('N','N',3*N,3*N,3*N,1.0d0,M_dgemm,3*N,O1,3*N,0.0d0,M,3*N) !DGEMM for matrix product (faster than matmul)
					end if
					!Call LAPACK routine and obtain the matrices P(er), L(ow), U (O2) and P'L (O1) from the answer:
					Low = M
					call DGETRF2(3*N,3*N,Low,3*N,ipiv,info)
					if (info /= 0) stop "lufact: info /= 0"
					O2 = 0.0d0
					Per = 0.0d0
					do u = 1, 3*N
						O2(u,u:3*N) = Low(u,u:3*N)
						Low(u,u:3*N) = 0.0d0
						Low(u,u) = 1.0d0
						Per(u,u) = 1.0d0
					end do
					do u = 1, 3*N
						row = Per(u,:)
						Per(u,:) = Per(ipiv(u),:)
						Per(ipiv(u),:) = row
					end do
					call DGEMM('T','N',3*N,3*N,3*N,1.0d0,Per,3*N,Low,3*N,0.0d0,O1,3*N)
					do u = 1, 3*N
						diag_O2(u) = O2(u,u)
					end do
					lambda0 = lambda0 + log(dabs(diag_O2))
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
			  exp_lyap = lambda0 / (t_max - t_trans)*1.0d0			  
			  write(10,*) J, exp_lyap(1)
		  end do  
		  
		  ! Close output file:		  
		  close(unit=10)
		  
		  ! Gnuplot script
		  ! Plot the main Lyapunov exponent against coupling J:
		  open(unit=11, file="exp_lyap1.plt", status="replace")
		    write(11, '(a)') "set terminal wxt"
		    write(11, '(a)') "set grid"
		    write(11, '(a)') "unset key"		   
		    write(11, '(a)') "set ytics font ', 15'"
		    write(11, '(a)') "set xtics 0.01 font ', 15'"			
		    write(11, '(a)') "set pointsize 0.6"
		    write(11, '(a)') "set xrange[0.01:0.04]"
		    write(11, '(a)') "set yrange[-0.0045:0.0025]"		   
		    write(11, '(a)') "set xlabel 'J' font ', 20'"
		    write(11, '(a)') "set ylabel '{/Symbol l}_{L}' font ', 20' offset 0,2,0"
		    write(11, '(a)') "plot 'exp_lyap1.dat' w linespoints pt 21 lt 8 lw 1.5"
		  close(unit=11)
		  
		  ! Call Gnuplot script for plot:
		  call system("gnuplot -p exp_lyap1.plt")		  
		  
	  end program	  
