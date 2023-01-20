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
            
            ! Diffusive coupling equation:
            G = J
            C = G * O
            I_s = sum(C, DIM=2)
            
        end subroutine diff_coupl
        
        !Create png frames in Gnuplot
        subroutine gnuplot_frame(L, N, time, X_o)
            implicit none
            
            ! Declaration of variables:
            integer, intent(in) :: L, N, time
            real(kind=8), dimension(N), intent(in) :: X_o
            integer :: q, qo
            
            ! Gnuplot script to create the frame:
            open(unit=16, file="frame_script.plt", status="replace")
            
                ! Specs for the plot:
                write(16,'(a)') 'set term png enhanced'
                write(16,'(a)') 'set size square'
        !       write(16,'(a)') 'set xtics 1'
        !       write(16,'(a)') 'set ytics 1'
                write(16,'(a)') 'set tic'
                write(16,'(a)') 'unset key'
                write(16,'(a)') 'set format x ""'
                write(16,'(a)') 'set format y ""'
                write(16,'(a,I0,a)') 'set xrange[-1:',L,']'
                write(16,'(a,I0,a)') 'set yrange[-1:',L,']'      
                write(16,'(a)') 'set cblabel "Membrane Potential" font ",20" offset "2,0,0"'
                write(16,'(a)') "set palette defined (1.'dark-blue',2.'blue',3.'cyan',5.'green',6.'yellow',7.'orange',8.'red',9.'dark-red')"
                write(16,'(a)') 'set cbrange[-1.0:1.0]'
                write(16,'(a,I0,a,I0,a)') 'set title "L = ',L,',  t = ',time,'" offset 0,-1.2 font ",20"'

                ! Name the png file using the respective frame time:
        !       write(16,'(a,I0,a)') "set output '[INSERT_FOLDER_PATH]\net_frame",time,".png'" !use for Windows
                write(16,'(a,I0,a)') "set output '[INSERT_FOLDER_PATH]/net_frame",time,".png'" !use for Linux

                ! Pass data to the script and plot:
                write(16,'(a)') '$data << EOD'
                do q = 1, L
                    qo = q - 1
                    write(16,*) X_o((qo*L+1):(q*L))
                end do
                write(16,'(a)') 'EOD'
                write(16,'(a)') "plot '$data' matrix with image"
                write(16,'(a)') 'replot'
                
            close(unit=16)
            
            ! Call the written script:
            call system("gnuplot frame_script.plt")

        end subroutine gnuplot_frame
        
    end module funcs_subrs
    
    ! Main program
    program KTz_net_diff_coupl_video
        use funcs_subrs
    !   use ifport !use to compile with ifort
        implicit none
        
        !Simulates the lattice of diffusively coupled KTz cells and creates a video of the simulation.
        !The program writes and calls a Gnuplot script that plots the spatiotemporal portrait of each 
        !time step in png format inside a folder. Each png constitutes the frames of the video. At the end, 
        !FFmpeg is called to paste all the pngs into a MP4 video.
        
        ! Declaration of variables:
        integer (kind = 4), parameter :: seed = 666, seed0 = 0 !use to compile with gfortran
        integer :: L_side, N_cells, cell_sample, a, t_max, t_trans, t1, dt, St, stim_dur
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
        t_trans = 1 !Transient
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
    !       If (t1 == 1) then
    !           I_snps = 0.0d0
    !       end if
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
            ! Write and call a Gnuplot script to generate the png frame:
            if (t1 >= t_trans) then
                call gnuplot_frame(L_side, N_cells, t1, Xo)
            end if
            ! Calculate the KTz equations:
            call ktz_model(N_cells, Xo, Yo, Zo, I, X, Y, Z)
        end do
        
        ! Close output files:
        close(unit=10)
        close(unit=11)
        close(unit=12)
        
        ! Call ffmpeg to paste all pngs in mp4 format:
    !   call system("ffmpeg -f image2 -r 100.0 -start_number 1 -i [INSERT_FOLDER_PATH]\net_frame%d.png -qscale 1 net.mp4") !use for Windows
        call system("ffmpeg -f image2 -r 100.0 -start_number 1 -i [INSERT_FOLDER_PATH]/net_frame%d.png -qscale 1 net.mp4") !use for Linux         
    !   call system("./net.mp4")        
        
    end program KTz_net_diff_coupl_video
