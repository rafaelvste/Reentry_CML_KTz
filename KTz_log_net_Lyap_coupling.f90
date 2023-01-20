    ! Functions and subroutines:
    module funcs_subrs
        implicit none
        
    contains
    
        !KTz model and network Jacobian
        subroutine ktz_log_net_jacob(L, N, n_s, time, trans, J, X_o, Y_o, Z_o, I_t, X_f, Y_f, Z_f, Jacob)
            implicit none
        
            ! Declaration of variables:
            integer, intent(in) :: L, N, time, trans
            real(kind=8), intent(in) :: J
            real(kind=8), dimension(N), intent(in) :: n_s, X_o, Y_o, Z_o, I_t
            real(kind=8), dimension(N), intent(out) :: X_f, Y_f, Z_f
            real(kind=8), dimension(3*N,3*N), intent(out) :: Jacob
            integer :: u, u_o, u_f, v, v_o
            real(kind=8) :: k, T, del, lamb, Xr, Ie
            real(kind=8), dimension(N) :: arg, a0, a11, a12, a13, b11
            real(kind=8), dimension(3,3) :: Jacob_0
            real(kind=8), dimension(3,3,N) :: Jacob_ii, Jacob_ij

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
            Y_f = X_o
            Z_f = (1.0d0 - del) * Z_o - lamb * (X_o - Xr)
 
            ! Calculate the Jacobian:
            if (time >= trans) then !Only after transient
                a0 = 1.0d0 / (T*(1.0d0 + dabs(arg))**2)
                ! Values in submatrix Jacob_ii:
                a11 = (1.0d0 - n_s*J) * a0
                a12 = -k * a0
                a13 = a0
                ! Value in submatrix Jacob_ij:
                b11 = J * a0
                ! Place the values in the Jacob_ii and Jacob_ij tensors (the 3rd dimension refer to network cell):
                Jacob_ii(1,1,:) = a11
                Jacob_ii(1,2,:) = a12
                Jacob_ii(1,3,:) = a13
                Jacob_ii(2,1,:) = 1.0d0
                Jacob_ii(3,1,:) = -lamb
                Jacob_ii(3,3,:) = 1.0d0 - del
                Jacob_ij(1,1,:) = b11
                ! Fill the Jacobian matrix of the network with Jacob_ii and Jacob_ij:
                Jacob_0 = 0.0d0
                u_o = -2
                do u = 1, N
                    u_o = u_o + 3
                    v_o = -2
                    do v = 1, N
                        v_o = v_o + 3
                        if (v == u) then
                            Jacob(u_o:u_o+2, v_o:v_o+2) = Jacob_ii(:,:,u)
                        else if ((v == (u-1)) .or. (v == (u+1))) then
                            Jacob(u_o:u_o+2, v_o:v_o+2) = Jacob_ij(:,:,u)
                        else if ((v == (u-L)) .or. (v == (u+L))) then
                            Jacob(u_o:u_o+2, v_o:v_o+2) = Jacob_ij(:,:,u)
                        else
                            Jacob(u_o:u_o+2, v_o:v_o+2) = Jacob_0(:,:)
                        end if
                    end do
                end do
                ! Fix the matrices of the left and right columns:
                u_f = 1
                do u = 1, (L-1)
                    u_o = u_f + 3*(L-1)
                    u_f = u_o + 3
                    Jacob(u_o:u_f-1, u_f:u_f+2) = Jacob_0(:,:)
                    Jacob(u_f:u_f+2, u_o:u_f-1) = Jacob_0(:,:)        
                end do
            end if
 
        end subroutine ktz_log_net_jacob
        
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

        !LU decomposition (similar to SciPy and Matlab lu() function)
        subroutine lu(N, M_input, O_1, O_2)
            implicit none
            
            integer, intent(in) :: N
            real(kind=8), dimension(3*N,3*N), intent(in) :: M_input
            real(kind=8), dimension(3*N,3*N), intent(out) :: O_1, O_2
            integer :: u, info
            real(kind=8), dimension(3*N) :: row
            real(kind=8), dimension(3*N,3*N) :: Low, Per
            integer, dimension(3*N) :: ipiv
            external :: DGEMM, DGETRF2 !External subroutines for matrix product and LU decomposition
            
            ! LAPACK subroutine that receives a matrix in Low, decomposes it into lower
            ! and upper triangular (O_2) and returns the results in Low, ignoring the 1s
            ! on the diagonal of the lower. ipiv is a vector with the indices of the
            ! permuted rows:
            Low = M_input
            call DGETRF2(3*N,3*N,Low,3*N,ipiv,info)
            
            if (info /= 0) stop "lufact: info /= 0"
            
            ! Removes the upper triangular (O_2) matrix from Low and restores the 1s
            ! to the diagonal in Low. It also creates an identity matrix in Per
            ! for the permutations:
            O_2 = 0.0d0
            Per = 0.0d0
            do u = 1, 3*N
                O_2(u,u:3*N) = Low(u,u:3*N)
                Low(u,u:3*N) = 0.0d0
                Low(u,u) = 1.0d0
                Per(u,u) = 1.0d0
            end do
            
            ! Use the ipiv vector to rewrite the identity as the permutation matrix (Per):
            do u = 1, 3*N
                row = Per(u,:)
                Per(u,:) = Per(ipiv(u),:)
                Per(ipiv(u),:) = row
            end do
            
            ! Calculate the remaining (permuted lower) matrix O_1 from the LU decomposition as
            ! O_1 = Per_transpose*Lower_triangular
            call DGEMM('T','N',3*N,3*N,3*N,1.0d0,Per,3*N,Low,3*N,0.0d0,O_1,3*N)
            
        end subroutine lu
        
        ! Iterate time and determine the Lyapunov spectrum for different coupling J:
        subroutine lyapunov(L, N, n_s, J, I_stim, rand_y, rand_x, exp_lyap)
            implicit none

            ! Declaration of variables:
            integer, intent(in) :: L, N
            real(kind=8), intent(in) :: J
            real(kind=8), dimension(N), intent(in) :: n_s, I_stim, rand_y, rand_x
            real(kind=8), dimension(3*N), intent(out) :: exp_lyap           
            integer :: t_max, t_trans, stim_dur, dt, t1, St, u
            real(kind=8), dimension(N) :: Zo, Yo, Xo, Z, Y, X, I, I_snps
            real(kind=8), dimension(3*N) :: diag_O2, sum_log_diag
            real(kind=8), dimension(3*N,3*N) :: M, M_dgemm, O1, O2          
            external :: DGEMM !External subroutine for matrix poduct
            
            ! Additional parameters:
            t_max = 100000 !Total time steps
            t_trans = 1 !Transient
            dt = 92 !Pacing period (P in the paper)
            stim_dur = 10 !Stimulus duration
            
            sum_log_diag = 0.0d0
            
            ! Iterate time:
            do t1 = 1, t_max
                ! Initial conditions:
                if (t1 == 1) then
                    Xo = -0.5d0 + rand_x
                    Yo = -0.5d0 + rand_y
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
                ! Calculate the KTz equations and Jacobian:
                call ktz_log_net_jacob(L, N, n_s, t_max, t_trans, J, Xo, Yo, Zo, I, X, Y, Z, M)
                ! Calculate the Lyapunov exponents using the Eckmann-Ruelle method:
                if (t1 >= t_trans) then
                    ! Multiply the Jacobian by the remaining matrix from the LU
                    ! decomposition:
                    M_dgemm = M
                    if (t1 > t_trans) then
                        ! DGEMM for matrix product (faster than matmul):
                        ! M = M*O1
                        call DGEMM('N','N',3*N,3*N,3*N,1.0d0,M_dgemm,3*N,O1,3*N,0.0d0,M,3*N)
                    end if
                    ! LU decomposition (retuns the upper triangular matrix O2 and
                    ! the remaining matrix O1):
                    call lu(N, M, O1, O2)                   
                    ! Get the diagonal elements of the upper triangular matrix:
                    do u = 1, 3*N
                        diag_O2(u) = O2(u,u)
                    end do
                    ! Store in a sum the log of each diagonal element:
                    sum_log_diag = sum_log_diag + log(dabs(diag_O2))
                end if
            end do

            ! Divide the sum of the diagonal elements for the the total time steps:
            exp_lyap = sum_log_diag / (t_max - t_trans + 1.0d0)*1.0d0

        end subroutine lyapunov
        
    end module funcs_subrs
    
    ! Main program
    program KTz_log_net_Lyap_coupling_m
        use funcs_subrs
    !   use ifport !use to compile with ifort
        implicit none
        
        !Simulates the lattice of diffusively coupled KTz cells and calculates the Lyapunov spectrum.
        !The exponents are approximated using the Eckmann-Ruelle method, where the Jacobian
        !matrix is triangularized by LU decompostion at each iteration using Lapack/MKL routines.
        !Outputs the main exponent against coupling J.
        
        ! Declaration of variables:
        integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran                        
        integer :: n_J_max, n_J, L_side, N_cells, a
        real(kind=8) :: J_p, J_0, J_f, dJ
        real(kind=8), allocatable, dimension(:) :: I_stimulus, random_y, random_x, n_snps, lyap_exp
        
        ! Random seed:
    !   call srand(666) !use to compile with ifort    
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
    
        ! Allocate memory:
        allocate(I_stimulus(N_cells))
        allocate(random_y(N_cells))
        allocate(random_x(N_cells))
        allocate(n_snps(N_cells))
        allocate(lyap_exp(3*N_cells))
        
        ! Generate random numbers for the IC for various realizations and stimulus
        ! vector (cells of the rightmost column):
        I_stimulus = 0.0d0
        random_y = 0.0d0
        random_x = 0.0d0
        do a = 1, L_side
            I_stimulus(a*L_side) = 0.1d0 !Stimulus current value
    !       random_y(a*L_side) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
    !       random_x(a*L_side) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
            random_y(a*L_side) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
            random_x(a*L_side) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7     
        end do

        ! Number of 'synapses'/couplings for each cell (coordination number):
        n_snps = 0.0d0
        if (L_side == 2) then
            n_snps = 2.0d0
        end if
        if (L_side > 2) then
            n_snps(1) = 2.0d0
            do a = 2, (L_side-1)
                n_snps(a) = 3.0d0
            end do
            n_snps(L_side) = 2.0d0
            n_snps(L_side+1) = 3.0d0
            do a = (L_side+2), (N_cells-L_side-1)
                n_snps(a) = 4.0d0
            end do
            n_snps(N_cells-L_side) = 3.0d0
            n_snps(N_cells-L_side+1) = 2.0d0    
            do a = ((N_cells-L_side)+2), (N_cells-1)
                n_snps(a) = 3.0d0
            end do
            n_snps(N_cells) = 2.0d0
        end if
        ! Fix the number of synapses of the left and right columns:
        if (L_side > 3) then
            do a = 2, (L_side-2)
                n_snps(a*L_side) = 3.0d0
                n_snps(a*L_side+1) = 3.0d0
            end do
        end if
        
        ! Open output files:
        open(unit=10, file="exp_lyap1.dat", status="replace")

        ! Iterate parameter J:
        do n_J = 1, n_J_max
            J_p = J_p + dJ
            ! Iterate time and obtain the Lyapunov spectrum for J:
            call lyapunov(L_side, N_cells, n_snps, J_p, I_stimulus, random_y, random_x, lyap_exp)
            ! Write to file the main Lyapunov exponent:
            write(10,*) J_p, lyap_exp(1)
        end do
        
        ! Close output files:
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
        
    end program KTz_log_net_Lyap_coupling