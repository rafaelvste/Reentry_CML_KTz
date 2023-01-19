clear
format long

% Random seed:
rng(666)

% Network parameters:
L_side = 80;	
N_cells = L_side^2;
dt = 92; %Pacing period

% Cell for measures:
cell_sample = (N_cells-L_side)/2;

%Additional parameters:
t_max = 50000; %Total time steps
dt = 92; %Pacing period (P in the paper)
stim_dur = 10; %Stimulus duration

% Create arrays:
rand_y = zeros(N_cells,1);
rand_x = zeros(N_cells,1);
Xo = zeros(N_cells,1);
Yo = zeros(N_cells,1);
Zo = zeros(N_cells,1);
X = zeros(N_cells,1);
Y = zeros(N_cells,1);
Z = zeros(N_cells,1);
I = zeros(N_cells,1);
I_stim = zeros(N_cells,1);
I_snps = zeros(N_cells,1);

% Outputs:
time = zeros(t_max,1);
s = zeros(t_max,1);
potential = zeros(t_max,1);
current = zeros(t_max,1);

% Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
for a = 1:L_side
    rand_y(a*L_side,1) = (rand - 1.0/2.0)*1.0e-7;
    rand_x(a*L_side,1) = (rand - 1.0/2.0)*1.0e-7;
    I_stim(a*L_side,1) = 0.1;
end

% Start the time loop:
for t1 = 1:t_max
    % Initial conditions:
    if t1 == 1
        Xo(:,1) = -0.5 + rand_y(:,1);
        Yo(:,1) = -0.5 + rand_x(:,1);
        Zo(:,1) = 0.0;
        % Stimulus:
        St = 1; %Count the time steps for stimulation
        I = I_stim;
    % Update variables:
    elseif t1 > 1
        Xo = X;
        Yo = Y;
        Zo = Z;
        % Periodic stimulus:
        I = 0.0;
        St = St + 1;
        if St <= stim_dur
            I = I_stim;
        end
        if St == dt
            St = 0;
        end
    end
    % Diffusive coupling:
    I_snps = diff_coupl(L_side, N_cells, Xo);
%    if t1 == 1
%        I_snps = 0.0;
%    end
    I = I + I_snps;   
    % Output the membrane potential x(t) and current I(t) from sampling cell:
    time(t1,1) = t1;
    potential(t1,1) = Xo(cell_sample);
    current(t1,1) = I(cell_sample);
    % Calculate and output the network activity:
    for a = 1:N_cells
        if Xo(a) > 0.0
            s(t1,1) = s(t1,1) + 1;
        end
    end    
    % Calculate the KTz equations:
    [X, Y, Z] = ktz_model_net(N_cells, Xo, Yo, Zo, I);
end

% Save network activity to file:
matriz = [time s];
save('activity.dat', 'matriz', '-ascii', '-double')

% Plot the activity:
plot(time, s, '-o')
xlabel('t')
ylabel('s(t)')
legend({'Network activity'},'Location','northeast')
grid on
axis auto

% Plot the membrane potential and current of cell b:
plot(time, potential, '-o', time, current, '-o')
xlabel('t')
ylabel('X(t)')
legend({'Potential','Current'},'Location','northeast')
grid on
axis auto

%For MATLAB R2016b or newer; otherwise save the functions in a separate file.
%KTz model
function [X_f, Y_f, Z_f] = ktz_model_net(N, X_o, Y_o, Z_o, I_t)

   % Create array:
   arg = zeros(N_cells,1);
     
   % Model parameters:
   k = 0.6;
   T = 0.154;
   del = 0.001;
   lamb = 0.001;
   Xr = -0.48;
   Ie = 0.0;

   % Model equations:
   arg = (X_o - k * Y_o + Z_o + Ie + I_t) / T;
   X_f = arg ./(1.0 + abs(arg)); % Logistic KTz
   Y_f = X_o;
   Z_f = (1.0 - del)*Z_o - lamb*(X_o - Xr);

end

%Diffusive coupling (electrical synapse)
function I_s = diff_coupl(L, N, X_s)

   % Create arrays:
   O = zeros(N,4);
   G = zeros(N,4);
   C = zeros(N,4);

   % Diffusive coupling constant:
   J = 0.0257;
    
   % Matrix filled with the differences between potentials of first nearest neighbors
    
   for a = 1:N %Numbering from left to right and from bottom to top
       % Connection with the right neighbor:
       if a < N
           O(a,1) = X_s(a+1,1) - X_s(a,1);
       end
       % Connection with the left neighbor:
       if a > 1
           O(a,2) = X_s(a-1,1) - X_s(a,1);
       end
       % Connection with the neighbor above:
       if a <= (N-L)
           O(a,3) = X_s(a+L,1) - X_s(a,1);
       end
       % Connection with the neighbor below:
       if a > L
           O(a,4) = X_s(a-L,1) - X_s(a,1);
       end
   end

   % Fix the remaining boundary conditions of the left and right columns:
   for a = 1:(L-1)
       O(a*L,1) = 0.0;
       O(a*L+1,2) = 0.0;
   end
    
    % Diffusive coupling equation:
   G(:,:) = J;
   C = G .* O;
   I_s = sum(C,2);
    
end