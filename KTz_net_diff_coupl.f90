      program KTz_net_diff_coupl
   !   use ifport !use to compile with ifort	  
      implicit none
	  
	  !Simulates the lattice of diffusively coupled KTz cells and returns the network actvity
	  !(number of cells with x > 0), membrane potential and current of one of the cells.

          ! Declaration of variables:
		  integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran		  		  
		  integer :: L, N, b, a, t_max, t1, dt, St, a_o
		  integer, allocatable, dimension(:) :: s
		  real*8 :: k, T, del, lamb, Xr, Ie, J
		  real*8, allocatable, dimension(:) :: rand_y, rand_x, Zo, Yo, Xo, Z, Y, X, I, I_stim, arg
		  real*8, allocatable, dimension(:,:) :: O, G, C
		  
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
		  J = 0.0257d0
		  
		  ! Network parameters:	  
		  L = 80
		  N = L**2
		  dt = 92 !Pacing period
		  
		  ! Cell for measures:
		  b = (N-L)/2
		  
		  !Total time steps:
		  t_max = 50000
		  
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
		  allocate(s(t_max))
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
		  s = 0
		  
		  ! Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
		  do a = 1, L
	!		rand_y(a*L) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
	!		rand_x(a*L) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
			rand_y(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
			rand_x(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7			
			I_stim(a*L) = 0.1d0
			! Use to stimulate more columns (count from right to left):
	!		 do a_o = 1, 4
	!			rand_y(a*L-a_o) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
	!			rand_x(a*L-a_o) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
	!			rand_y(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran			
	!			rand_x(a*L) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7
	!			I_stim(a*L-a_o) = 0.1d0
	!		 end do		
		  end do
		  
		  ! Open output files:
		  open(unit=10, file="activity.dat", status="replace")
		  open(unit=11, file="potential.dat", status="replace")
		  open(unit=12, file="current.dat", status="replace")	  	  

		  ! Iterate time:			 		
		  do t1 = 1, t_max
			 ! Initial conditions for the KTz model:			 
			 if (t1 == 1) then
				 St = 1
				 Zo = 0.0d0
				 Yo = -0.5d0
				 Xo = -0.5d0
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
			 ! Output the membrane potential x(t) and current I(t) of cell b:
			 write(11,*) t1, Xo(b)			 
			 write(12,*) t1, I(b)
			 ! Calculate the KTz equations:			 
			 arg = (Xo - k * Yo + Zo + Ie + I) / T
			 X = arg / (1.0d0 + dabs(arg)) !Logistic KTz
	!		 X = tanh(arg)                 !Hyperbolic KTz
			 Y = Xo
			 Z = (1.0d0 - del) * Zo - lamb * (Xo - Xr)
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
			 !Calculate activity:
			 do a = 1, N
				if (Xo(a) > 0.0) then
					s(t1) = s(t1) + 1
				end if
			 end do
			 ! Output the network activity:
			 write(10,*) t1, s(t1)			 
			 ! Diffusive coupling:
			 G = J
			 C = G * O
			 I = sum(C, DIM=2)
		  end do
		  
		  ! Close output files:
		  close(unit=10)
		  close(unit=11)
		  close(unit=12)	
		  
		  ! Gnuplot scripts
		  ! Plot the network activity s(t):
		  open(unit=13, file="activity.plt", status="replace")
		    write(13, '(a)') "set terminal wxt"
		    write(13, '(a)') "set grid"
		    write(13, '(a)') "unset key"		   
		    write(13, '(a)') "set ytics font ', 10'"
		    write(13, '(a)') "set xtics font ', 10'"			
		    write(13, '(a)') "set pointsize 0.3"
		    write(13, '(a)') "set xlabel 't' font ', 16'"
		    write(13, '(a)') "set ylabel 's(t)' font ', 16'"
		    write(13, '(a)') "plot 'activity.dat' w linespoints pt 21 lt 8 lw 1.5"
		  close(unit=13)

		  call system("gnuplot -p activity.plt")
		  
		  ! Plot the membrane potential x(t) and current I(t) of one of the central cells:
		  open(unit=14, file="potential.plt", status="replace")
			  write(14, '(a)') "set terminal wxt"
			  write(14, '(a)') "set grid"
			  write(14, '(a)') "set ytics font ', 10'"
			  write(14, '(a)') "set xtics font ', 10'"			  
			  write(14, '(a)') "unset ylabel"
			  write(14, '(a)') "set xlabel 't' font ', 16'"
			  write(14, '(a)') "plot 'potential.dat' title 'Potential x(t)' w linespoints pointtype 21 lt 6"
			  write(14, '(a)') "replot 'current.dat' title 'Current I(t)' w linespoints pointtype 21 lt 7"
		  close(unit=14)
		  
		  call system("gnuplot -p potential.plt")
		  
      end program 
