      program KTz_net_diff_coupl_video
   !   use ifport !use to compile with ifort	  
      implicit none

	  !Simulates the lattice of diffusively coupled KTz cells and creates a video of the simulation.
	  !The program writes and calls a Gnuplot script that plots the spatiotemporal portrait of each 
	  !time step in png format inside a folder. Each png constitutes the frames of the video. At the end, 
	  !FFmpeg is called to glue all the pngs into a MP4 video.

          ! Declaration of variables:
		  integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran		  		  
		  integer :: L, N, b, a, t_max, t_trans, t1, q, qo, dt, St, a_o
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
		  dt = 92 !Pacing period (P in the paper)
		  
		  ! Cell for measures:
		  b = (N-L)/2
		  
		  !Total time steps and transient:
		  t_max = 50000
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
			 !Create and write the specs of the gnuplot script to generate the png frames:
			 if (t1 >= t_trans) then
				 open(unit=16, file="frame_script.plt", status="replace")			  
				 write(16,'(a)') 'set term png enhanced'
				 write(16,'(a)') 'set size square'
		!		 write(16,'(a)') 'set xtics 1'
		!		 write(16,'(a)') 'set ytics 1'
				 write(16,'(a)') 'set tic'
				 write(16,'(a)') 'unset key'
				 write(16,'(a)') 'set format x ""'
				 write(16,'(a)') 'set format y ""'
				 write(16,'(a,I0,a)') 'set xrange[-1:',L,']'
				 write(16,'(a,I0,a)') 'set yrange[-1:',L,']'		 
				 write(16,'(a)') 'set cblabel "Membrane Potential" font ",20" offset "2,0,0"'
				 write(16,'(a)') "set palette defined (1.'dark-blue',2.'blue',3.'cyan',5.'green',6.'yellow',7.'orange',8.'red',9.'dark-red')"
				 write(16,'(a)') 'set cbrange[-1.0:1.0]' 
			 end if
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
	!		 X = tanh(arg)               !Hyperbolic KTz
			 Y = Xo
			 Z = (1.0d0 - del) * Zo - lamb * (Xo - Xr)
			 ! Write the membrane potential of each cell after t1 iterations and call the
			 !script to create the png frames in the designated folder:
			 if (t1 >= t_trans) then			 
				 write(16,'(a,I0,a,I0,a)') 'set title "L = ',L,',  t = ',t1,'" offset 0,-1.2 font ",20"'			 
	!			 write(16,'(a,I0,a)') "set output '[INSERT_FOLDER_PATH]\net_frame",t1,".png'" !use for Windows
				 write(16,'(a,I0,a)') "set output '[INSERT_FOLDER_PATH]/net_frame",t1,".png'" !use for Linux	
				 write(16,'(a)') '$data << EOD'
				 do q = 1, L
					qo = q - 1
					write(16,*) Xo((qo*L+1):(q*L))
				 end do
				 write(16,'(a)') 'EOD'
				 write(16,'(a)') "plot '$data' matrix with image"
				 write(16,'(a)') 'replot'
				 close(unit=16)
				 call system("gnuplot frame_script.plt")
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
		  
		  ! Call ffmpeg to 'glue' all the pngs in mp4 format:
	!	  call system("ffmpeg -f image2 -r 100.0 -start_number 1 -i [INSERT_FOLDER_PATH]\net_frame%d.png -qscale 1 net.mp4") !use for Windows
		  call system("ffmpeg -f image2 -r 100.0 -start_number 1 -i [INSERT_FOLDER_PATH]/net_frame%d.png -qscale 1 net.mp4") !use for Linux		  
	!	  call system("./net.mp4")
		  
      end program 
