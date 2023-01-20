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
        !   X_f = tanh(arg)                   !Hyperbolic KTz
            Y_f = X_o
            Z_f = (1.0d0 - del) * Z_o - lamb * (X_o - Xr)
            
        end subroutine ktz_model
        
        !Diffusive coupling (electrical synapse)
        subroutine diff_coupl(L, N, X_s, I_s)
            implicit none
            
            ! Declaration of variables:
            integer, intent(in) :: L, N
            real(kind=8), dimension(N), intent(in) :: X_s
            real(kind=8), dimension(N), intent(out) :: I_s
            integer :: a
            real(kind=8) :: J
            real(kind=8), dimension(N,4) :: O, G, C

            ! Diffusive coupling constant:        
            J = 0.0264d0
        
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
        
        !Iterate time and calculate the network standard deviation averaged for several realizations
        subroutine stand_dev_x_avg(L, N, sigma_x_avg)
    !       use ifport !use to compile with ifort

            implicit none

            ! Declaration of variables:
            integer, intent(in) :: L, N
            real(kind=8), intent(out) :: sigma_x_avg
            integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran            
            integer :: realiz_max, realiz, a, a_o, t_max, t_trans, stim_dur, dt, t1, St
            real(kind=8) :: sum_sigma, sum_realiz
            real(kind=8), dimension(N) :: Zo, Yo, Xo, Z, Y, X, I, I_snps, I_stim, sum_x, sum_x2, sigma
            real(kind=8), allocatable, dimension(:,:) :: rand_y, rand_x
            
            ! Random seed:
    !       call srand(666) !use to compile with ifort
            call srand(seed) !use to compile with gfortran
            
            ! Additional parameters:
            t_max = 100000 !Total time steps
            t_trans = 50000 !Transient
            dt = 92 !Pacing period (P in the paper)
            stim_dur = 10 !Stimulus duration

            !Number of realizations for each L:
            realiz_max = 10
        
            ! Allocate memory:
            allocate(rand_y(N,realiz_max))
            allocate(rand_x(N,realiz_max))  
        
            ! Generate random numbers for the IC for various realizations:
            rand_y = 0.0d0
            rand_x = 0.0d0
            do realiz = 1, realiz_max
                do a = 1, L
            !       random_y(a*L,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
            !       random_x(a*L,realiz) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
                    rand_y(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
                    rand_x(a*L,realiz) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7
                end do
            end do
            
            ! Generate the stimulus vector (cells of the rightmost column):
            I_stim = 0.0d0
            do a = 1, L
                I_stim(a*L) = 0.1d0 !Stimulus current value 
            end do
            
            sum_realiz = 0.0d0
            ! Iterate realizations:
            do realiz = 1, realiz_max
                sum_sigma = 0.0d0 !Store standard deviation values for time average
                ! Iterate time:
                do t1 = 1, t_max
                    ! Initial conditions:
                    if (t1 == 1) then
                        Xo(:) = -0.5d0 + rand_x(:,realiz)
                        Yo(:) = -0.5d0 + rand_y(:,realiz)
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
                    call diff_coupl(L, N, Xo, I_snps)
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
                ! Time average and sum for realizations average:
                sum_realiz = sum_realiz + sum_sigma/((t_max-t_trans)*1.0d0)
            end do
            
            ! Average standard deviation by realizations:
            sigma_x_avg = sum_realiz / (realiz_max*1.0d0)
            
        end subroutine stand_dev_x_avg
        
    end module funcs_subrs
    
    ! Main program
    program KTz_net_stand_dev_size_m
        use funcs_subrs
        implicit none
        
        !Simulates the network and calculates the column-based standard deviation.
        !Calculates the standad deviation of x in each column and the average of all columns.
        !Then averages all time steps in a realization and averages the realizations.
        !This is calculated for different values of the network size L.
        
        ! Declaration of variables:
        integer :: L_max, L_side, N_cells
        real(kind=8) :: sigma_x, tic, toc

        call cpu_time(tic)
        
        ! Network parameters:
        L_max = 40
        
        ! Open output files:
        open(unit=10, file="sigma_x_L.dat", status="replace")

        ! Iterate parameter L:
        do L_side = 2, L_max
            N_cells = L_side**2
            ! Iterate time and calculate the average standard deviation of the network:
            call stand_dev_x_avg(L_side, N_cells, sigma_x)
            write(10,*) L_side, sigma_x
        end do

        ! Close output files:
        close(unit=10)

        call cpu_time(toc)

        print*, toc-tic
        
        ! Gnuplot script
        ! Plot the standard deviation against size L:
        open(unit=10, file="sigma_x_L.plt", status="replace")
            write(10, '(a)') "set terminal wxt"
            write(10, '(a)') "set grid"
            write(10, '(a)') "unset key"
            write(10, '(a)') "set ytics font ', 15'"
            write(10, '(a)') "set xtics font ', 15'"
            write(10, '(a)') "set pointsize 0.6"
            write(10, '(a)') "set xrange[2:40]"
            write(10, '(a)') "set xlabel 'L' font ', 20'"
            write(10, '(a)') "set ylabel '<~{/Symbol s}{0.4-}_{x}>' font ', 20'"
            write(10, '(a)') "plot 'sigma_x_L.dat' w linespoints pt 21 lt 8 lw 1.5"
        close(unit=10)

        ! Call Gnuplot script for plot:       
        call system("gnuplot -p sigma_x_L.plt")
        
    end program KTz_net_stand_dev_size_m