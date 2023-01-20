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
            J = 0.0257d0
        
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
        
    end module funcs_subrs
    
    ! Main program
    program KTz_net_diff_coupl
        use funcs_subrs
    !   use ifport !use to compile with ifort
        implicit none
        
        !Simulates the lattice of diffusively coupled KTz cells and returns the network actvity
        !(number of cells with x > 0), membrane potential and current of one of the cells.
        
        ! Declaration of variables:
        integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran
        integer :: L_side, N_cells, cell_sample, a, t_max, t1, dt, St, stim_dur
        integer, allocatable, dimension(:) :: s
        real(kind=8), allocatable, dimension(:) :: rand_y, rand_x, Zo, Yo, Xo, Z, Y, X, I, I_stim, I_snps
        
        ! Random seed:
    !   call srand(666) !use to compile with ifort    
        call srand(seed) !use to compile with gfortran
        
        ! Network parameters:
        L_side = 80
        N_cells = L_side**2
        
        ! Cell for measures:
        cell_sample = (N_cells-L_side) / 2
        
        ! Additional parameters:
        t_max = 50000 !Total time steps
        dt = 92 !Pacing period (P in the paper)
        stim_dur = 10 !Stimulus duration
        
        ! Allocate memory:
        allocate(rand_y(N_cells))
        allocate(rand_x(N_cells))
        allocate(Xo(N_cells))
        allocate(Yo(N_cells))
        allocate(Zo(N_cells))
        allocate(X(N_cells))
        allocate(Y(N_cells))
        allocate(Z(N_cells))
        allocate(I(N_cells))
        allocate(I_stim(N_cells))               
        allocate(I_snps(N_cells))
        allocate(s(t_max))
        s = 0
        
        ! Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
        rand_y = 0.0d0
        rand_x = 0.0d0
        I_stim = 0.0d0
        do a = 1, L_side
    !       rand_y(a*L_side) = (rand(0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with ifort
    !       rand_x(a*L_side) = (rand(0) - 1.0d0/2.0d0)*1.0e-7
            rand_y(a*L_side) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7 !use to compile with gfortran
            rand_x(a*L_side) = (rand(seed0) - 1.0d0/2.0d0)*1.0e-7
            I_stim(a*L_side) = 0.1d0 !Stimulus current value
        end do
        
        ! Open output files:
        open(unit=10, file="activity.dat", status="replace")
        open(unit=11, file="potential.dat", status="replace")
        open(unit=12, file="current.dat", status="replace")
                
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
            call diff_coupl(L_side, N_cells, Xo, I_snps)
            ! If (t1 == 1) then
                ! I_snps = 0.0d0
            ! end if
            I = I + I_snps
            ! Output the membrane potential x(t) and current I(t) from sampling cell:
            write(11,*) t1, Xo(cell_sample)
            write(12,*) t1, I(cell_sample)
            !Calculate and output the network activity:
            do a = 1, N_cells
                if (Xo(a) > 0.0d0) then
                    s(t1) = s(t1) + 1
                end if
            end do
            write(10,*) t1, s(t1)
            ! Calculate the KTz equations:
            call ktz_model(N_cells, Xo, Yo, Zo, I, X, Y, Z)
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
        
    end program KTz_net_diff_coupl
