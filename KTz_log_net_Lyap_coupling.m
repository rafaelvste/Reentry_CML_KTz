%Simulates the lattice of diffusively coupled KTz cells and calculates the Lyapunov spectrum.
%The exponents are approximated using the Eckmann-Ruelle method, where the Jacobian
%matrix is triangularized by LU decompostion at each iteration using the LU function.
%Outputs the main exponent against coupling J.

clear
format long

    tic

    % Parameters for the KTz model:			  			  
    k = 0.6;
    T = 0.154;
    del = 0.001;
    lamb = 0.001;
    Xr = -0.48;
    Ie = 0.0;

    % Diffusive coupling constant:		  
    J_0 = 0.01;
    J_f = 0.04;
    dJ = 1.0e-4;

    % Network parameters	  
    L = 2;	
    N = L^2;
    dt = 92; %Pacing period (P in the paper)

    %Total time steps and transient:
    t_max = 100000;	
    t_trans = 1;
    
    % Allocate memory:
    rand_y = zeros(N,1);
    rand_x = zeros(N,1);      
    Xo = zeros(N,1);
    Yo = zeros(N,1);
    Zo = zeros(N,1);      
    X = zeros(N,1);
    Y = zeros(N,1);
    Z = zeros(N,1);
    I = zeros(N,1);
    I_stim = zeros(N,1);            
    arg = zeros(N,1);
    O = zeros(N,4);
    eps = zeros(N,4);
    G = zeros(N,4);
    C = zeros(N,4);
    % Lyapunov variables:
    a0 = zeros(N,1);
    a11 = zeros(N,1);
    a12 = zeros(N,1);
    a13 = zeros(N,1);
    b11 = zeros(N,1);
    M_ii = zeros(3,3,N);
    M_ij = zeros(3,3,N);
    M_0 = zeros(3);
    M_cell = cell(N);
    M = zeros(3*N);
    % Discretize J variable:
    n_J_max = int16((J_f - J_0) / dJ);
    J = J_0 - dJ;     
    % Outputs:
    exp_lyap = zeros(n_J_max,3*N);
    values_J = zeros(n_J_max,1);  
    
    % Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
    for a = 1:L
        rand_y(a*L,1) = (rand - 1.0/2.0)*1.0e-7;
        rand_x(a*L,1) = (rand - 1.0/2.0)*1.0e-7;
        I_stim(a*L,1) = 0.1;
        % Use to stimulate more columns (count from right to left):
    %    for a_o = 1:4
    %       I(a*L-a_o) = 0.1;
    %       rand_y(a*L-a_o,1) = (rand - 1.0/2.0)*1.0e-7;
    %       rand_x(a*L-a_o,1)= (rand - 1.0/2.0)*1.0e-7;         
    end       

    % Number of synapse for each cell:
    n_s = zeros(N,1);
    if L == 2
        n_s(:,1) = 2;
    end
    if L > 2
        n_s(1,1) = 2.0;
        for a = 2:(L-1)
            n_s(a,1) = 3.0;
        end
        n_s(L,1) = 2.0;
        n_s(L+1,1) = 3.0;
        for a = (L+2):(N-L-1)
            n_s(a,1) = 4.0;
        end
        n_s(N-L,1) = 3.0;
        n_s(N-L+1,1) = 2.0;    
        for a = ((N-L)+2):(N-1)
            n_s(a,1) = 3.0;
        end
        n_s(N,1) = 2.0;
    end
    % Fix the number of synapses of the left and right columns:
    if L > 3
        for a = 2:(L-2)
            n_s(a*L,1) = 3.0;
            n_s(a*L+1,1) = 3.0;
        end
    end
    
    % Iterate parameter J:
    for n_J = 1:n_J_max
        J = J + dJ;
        lambda0 = zeros(3*N,1);       
        % Iterate time:			 		
        for t1 = 1:t_max
         % Initial conditions for the KTz model:			 
         if t1 == 1
             St = 1;
             Zo(:,1) = 0.0;
             Yo(:,1) = -0.5;
             Xo(:,1) = -0.5;
             I(:,1) = 0.0;
             % Random IC and stimulus:
             Yo = Yo + rand_y;
             Xo = Xo + rand_x;
             I = I_stim; 
         % Update variables of the KTz model:
         elseif t1 > 1
             St = St + 1;
             Zo = Z;
             Yo = Y;
             Xo = X;
             % Periodic stimulus:
             if St <= 10 %Adjust the duration of the stimulus
                I = I + I_stim;
             end
             if St == dt
                St = 0;
             end
         end   
         % Calculate the KTz equations:			 
         arg = (Xo - k * Yo + Zo + Ie + I) / T;
         X = arg ./ (1.0 + abs(arg)); %Logistic KTz
    %    X = tanh(arg)                %Hyperbolic KTz
         Y = Xo;
         Z = (1.0 - del) * Zo - lamb * (Xo - Xr);
         % Lyapunov exponent (vectorized):
         if t1 >= t_trans
            a0 = 1.0 ./ (T*(1.0 + abs(arg)).^2);
            % Values in the single element Jacobian matrices M_ii (vectorized):
            a11 = (1.0 - n_s*J) .* a0;
            a12 = -k * a0;
            a13 = a0;
            % Single value in the cross elements Jacobian matrices M_ij (vectorized):
            b11 = J * a0;
            % Place the values in the M_ii and M_ij tensors (the 3rd dimension refer to network cell):
            M_ii(1,1,:) = a11;
            M_ii(1,2,:) = a12;
            M_ii(1,3,:) = a13;
            M_ii(2,1,:) = 1.0;
            M_ii(3,1,:) = -lamb;
            M_ii(3,3,:) = 1.0 - del;
            M_ij(1,1,:) = b11;
            % Fill the Jacobian matrix of the network with M_ii and M_ij:        
            for u = 1:N
                for v = 1:N
                    if v == u
                        M_cell(u,v) = {M_ii(:,:,u)};
                    elseif (v == (u-1)) || (v == (u+1))
                        M_cell(u,v) = {M_ij(:,:,u)};
                    elseif (v == (u-L)) || (v == (u+L))
                        M_cell(u,v) = {M_ij(:,:,u)};
                    else
                        M_cell(u,v) = {M_0};
                    end
                end
            end
            % Fix the matrices of the left and right columns:
            for u = 1:(L-1)
                M_cell(u*L,u*L+1) = {M_0};
                M_cell(u*L+1,u*L) = {M_0};
            end
            % Perform the LU decomposition:
            M = cell2mat(M_cell);
            if t1 > t_trans
                M = M*O1;
            end
            [O1,O2] = lu(M);        
            lambda0 = lambda0 + log(abs(diag(O2)));
         end
         % Matrix filled with the differences between potentials of first nearest neighbors:
         O(1,1) = X(2,1) - X(1,1);
         O(1,2) = 0.0;
         O(1,3) = X(1+L,1) - X(1,1);
         O(1,4) = 0.0;
         for a = 2:L
            O(a,1) = X(a+1,1) - X(a,1);
            O(a,2) = X(a-1,1) - X(a,1);
            O(a,3) = X(a+L,1) - X(a,1);
            O(a,4) = 0.0;
         end
         for a = (L+1):(N-L)
            O(a,1) = X(a+1,1) - X(a,1);
            O(a,2) = X(a-1,1) - X(a,1);
            O(a,3) = X(a+L,1) - X(a,1);
            O(a,4) = X(a-L,1) - X(a,1);
         end
         for a = ((N-L)+1):(N-1)
            O(a,1) = X(a+1,1) - X(a,1);
            O(a,2) = X(a-1,1) - X(a,1);
            O(a,3) = 0.0;
            O(a,4) = X(a-L,1) - X(a,1);
         end
         O(N,1) = 0.0;
         O(N,2) = X(N-1,1) - X(N,1);
         O(N,3) = 0.0;
         O(N,4) = X(N-L,1) - X(N,1);
         % Fix the boundary conditions of the left and right columns:
         for a = 1:(L-1)
            O(a*L,1) = 0.0;
            O(a*L+1,2) = 0.0;
         end
         % Diffusive coupling:
         eps(:,:) = 0.0;
         G(:,:) = J; %* (1.0 + eps*1.0e-5);
         C = G .* O;
         I = sum(C,2);
        end
        exp_lyap(n_J,1:3*N) = lambda0 ./ (t_max - t_trans);
        values_J(n_J,1) = J;
    end
    
    toc
    
    % Save main exponent to files:
    matriz = [values_J exp_lyap(:,1)];
    save('exp_lyap1.dat', 'matriz', '-ascii')
    
    % Save spectrum to files:
    matriz = [values_J exp_lyap];
    save('exp_lyap.dat', 'matriz', '-ascii')    
      
    % Plot the exponents:
    plot(values_J, exp_lyap, '-o')
    xlabel('J')
    ylabel('Lyapunov Exponent')
    grid on
    axis auto   
