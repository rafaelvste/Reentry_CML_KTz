clear
format long

% Additional parameters:
t_max = 10000; %Total time steps
dt = 232; %Pacing period
I_stim = 0.1; %Stimulus current value
stim_dur = 10; %Stimulus duration

% Outputs:
potential = zeros(t_max,1);
current = zeros(t_max,1);
time = zeros(t_max,1);
time(:,1) = 1:1:t_max;

% Iterate time:
for t1 = 1:t_max
    % Initial conditions:
    if t1 == 1
        Xo = -0.5;
        Yo = -0.5;
        Zo = 0.0;
        I = I_stim;
        St = 1; %Count the timesteps for stimulation
    % Update variables:
    else
        Xo = X;
        Yo = Y;
        Zo = Z;
        I = 0.0;
        St = St + 1;
        % Periodic stimulus:
        if St <= stim_dur % adjust the duration of the stimulus
            I = I_stim;
        end
        if St == dt
            St = 0;
        end
    end
    potential(t1,:) = Xo;
    current(t1,:) = I;
    % Calculate the KTz equations:
    [X, Y, Z] = ktz_model(Xo, Yo, Zo, I);
end

% Save to files:
matriz = [time potential(1:end,1)];
save('fig3_ins1_matlab_m.dat', 'matriz', '-ascii', '-double')

%matriz = [time current(1:end,1)];
%save('current.dat', 'matriz', '-ascii')

% Plot the membrane potential of cell 1:
plot(time, potential(:,1), 'o')
xlabel('t')
ylabel('X(t)')
grid on
axis auto

%For MATLAB R2016b or newer; otherwise save the function in a separate file.
%KTz model: 
function [X_f, Y_f, Z_f] = ktz_model(X_o, Y_o, Z_o, I_t)

    % Model parameters:
   k = 0.6;
   T = 0.154;
   del = 0.001;
   lamb = 0.001;
   Xr = -0.48;
   Ie = 0.0;

    % Model equations:
   arg = (X_o - k * Y_o + Z_o + Ie + I_t) / T;
   X_f = arg /(1.0 + abs(arg)); % Logistic KTz
   Y_f = X_o;
   Z_f = (1.0 - del)*Z_o - lamb*(X_o - Xr);

end