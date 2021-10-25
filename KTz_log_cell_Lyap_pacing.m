clear
format long

%Calculates the Lyapunov spectrum for the periodically paced logistic KTz cell.
%The exponents are approximated using the Eckmann-Ruelle method, where the Jacobian
%matrix is triangularized by LU decompostion at each iteration using the MATLAB/Octave
%function 'lu'. Outputs the three exponents plotted for different pacings P.

    % Parameters of the KTz model:
    k = 0.6;
    T = 0.154;
    del = 0.001;
    lamb = 0.001;
    Xr = -0.48;
    Ie = 0.0;
    I_stim = 0.1;
	
    % Total time steps and transient:    
    t_max = 1e5;
    t_trans = 1;
    
    % Maximum pacing period:
    dt_max = 350;
    
    % Outputs:
    exp_lyap = zeros(dt_max,3);
    period(1:dt_max,1) = 1:dt_max;
    
    % Iterate the pacing period:
    for dt = 11:dt_max %Variable P in the paper
        lambda0 = zeros(3,1);
        % Iterate time:
        for t1 = 1:t_max
            % Initial conditions:
            if t1 == 1
                St = 1;
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
                    % Periodic stimulus:
                    St = St + 1;
                    if St <= 10 % adjust the duration of the stimulus
                        I = I_stim;
                    end                   
                    if St == dt
                        St = 0;
                    end
            end
            % Calculate the KTz equations:
            arg = (Xo - k * Yo + Zo + Ie + I) / T;       
            X = arg /(1.0 + abs(arg)); % Logistic KTz
            Y = Xo;
            Z = (1.0 - del)*Zo - lamb*(Xo - Xr);
            % Calculate the Lyapunov exponents:
            if t1 >= t_trans
                % Elements of the Jacobian matrix:
                a11 = 1.0 / (T*(1.0 + abs(arg))^2);
                a12 = -k*a11;
                a31 = -lamb;
                a33 = 1.0 - del;
                % Jacobian matrix:
                M = [a11 a12 a11;
                    1.0 0.0 0.0;
                    a31 0.0 a33];
                if t1 > t_trans
                    M = M*O1;
                end
                % LU decomposition/ Eckmann-Ruelle method:
                [O1,O2] = lu(M);        
                lambda0 = lambda0 + log(abs(diag(O2)));
            end            
        end
        exp_lyap(dt,1:3) = lambda0 ./ (t_max - t_trans);
    end
    
    % Save to files:
    matriz = [period exp_lyap(1:end,1)];
    save('exp_lyap1.dat', 'matriz', '-ascii')
    
    matriz = [period exp_lyap(1:end,2)];
    save('exp_lyap2.dat', 'matriz', '-ascii')
    
    matriz = [period exp_lyap(1:end,3)];
    save('exp_lyap3.dat', 'matriz', '-ascii')       
    
    % Plot the exponents:
    plot(period, exp_lyap, '-o')
    xlabel('P')
    ylabel('Lyapunov Exponent')
    legend({'Lambda 1','Lambda 2','Lambda 3'},'Location','northeast')
    grid on
    axis auto 