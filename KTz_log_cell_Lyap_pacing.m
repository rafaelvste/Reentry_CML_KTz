clear
format long

%Calculates the Lyapunov spectrum for the periodically paced logistic KTz cell.
%The exponents are approximated using the Eckmann-Ruelle method, where the Jacobian
%matrix is triangularized by LU decompostion at each iteration using the MATLAB/Octave
%function 'lu'. Outputs the three exponents plotted for different pacings P.

% Additional parameters:    
t_max = 1e5; %Total time steps
t_trans = 1; %Transient
I_stim = 0.1; %Stimulus current value
stim_dur = 10; %Stimulus duration
dt_max = 350; %Maximum pacing period

% Matrices for the LU decompostion:
M = zeros(3,3);
O1 = zeros(3,3);
O2 = zeros(3,3);

% Outputs:
exp_lyap = zeros(dt_max,3);
period(1:dt_max,1) = 1:dt_max;

% Iterate the pacing period:
for dt = (stim_dur+1):dt_max %Variable P in the paper
    sum_log_diag = zeros(3,1);
    % Iterate time:
    for t1 = 1:t_max
        % Initial conditions:
        if t1 == 1
            St = 1; %Count the time steps for stimulation
            Xo = 0.0;
            Yo = 0.0;
            Zo = 0.0;
            I = I_stim; 
            % Update variables:
            else
                Xo = X;
                Yo = Y;
                Zo = Z;
                I = 0.0;
                St = St + 1;
                % Periodic stimulus:
                if St <= stim_dur
                    I = I_stim;
                end
                if St == dt
                    St = 0;
                end
        end
        % Calculate the KTz equations and the Jacobian:
        [X, Y, Z, M] = ktz_log_jacob(Xo, Yo, Zo, I, t1, t_trans);
        % Calculate the Lyapunov exponents using the Eckmann-Ruelle method:
        if t1 >= t_trans
            % Multiply the Jacobian by the remaing matrix from the LU
            % decomposition:
            if t1 > t_trans
                M = M*O1;
            end
            % LU decomposition:
            [O1,O2] = lu(M);
            % Store the log of the diagonal elements of the triangular matrix:
            sum_log_diag = sum_log_diag + log(abs(diag(O2)));
        end
    end
    % Average the sum to obtain the Lyapunov spectrum:
    exp_lyap(dt,1:3) = sum_log_diag ./ (t_max - t_trans + 1);
end

% Save to files:
matriz = [period exp_lyap(1:end,1)];
save('exp_lyap1_m.dat', 'matriz', '-ascii', '-double')

matriz = [period exp_lyap(1:end,2)];
save('exp_lyap2.dat', 'matriz', '-ascii', '-double')

matriz = [period exp_lyap(1:end,3)];
save('exp_lyap3.dat', 'matriz', '-ascii', '-double')       

% Plot the exponents:
plot(period, exp_lyap, '-o')
xlabel('P')
ylabel('Lyapunov Exponent')
legend({'Lambda 1','Lambda 2','Lambda 3'},'Location','northeast')
grid on
axis auto

%For MATLAB R2016b or newer; otherwise save the function in a separate file.
%KTz model and Jacobian:
function [X_f, Y_f, Z_f, Jacob] = ktz_log_jacob(X_o, Y_o, Z_o, I_t, time, trans)

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
   
   Jacob = zeros(3,3);
   
    % Calculate the Jacobian:
   if time >= trans %Only after transient
       % Elements of the Jacobian matrix:
       a11 = 1.0 / (T*(1.0 + abs(arg))^2);
       a12 = -k*a11;
       a31 = -lamb;
       a33 = 1.0 - del;
       % Jacobian matrix:
       Jacob = [a11 a12 a11;
               1.0 0.0 0.0;
               a31 0.0 a33];
   end
   
end