% CFD assignment 2
clear; close all; clc;

%% Spatial discretization
nx   = 1000;
x_min = -1; x_max = 1;
dx   = (x_max - x_min)/nx;
x_center = linspace(x_min + dx/2, x_max - dx/2, nx);

%% Parameters
g = 9.81;
flux_type = 'hll';     % options: 'rusanov','roe','hll','hllroe'

% Initial conditions
HL = 50; HR = 1;

% Time snapshots 
t_list = [0.00 0.0067 0.0133 0.0200];

% Courant number values
Cmax_list = [0.9 1.2 2.0];

%% Loop over time integration schemes and Courant numbers
for time_scheme = ["fe","rk4"]
    for Cmax = Cmax_list

        H_snapshot = zeros(numel(t_list), nx);
        u_snapshot = zeros(numel(t_list), nx);

        % Apply initial conditions
        H = HL*ones(1, nx);
        H(x_center > 0) = HR;
        Q = zeros(1, nx);
        U = [H; Q];

        t = 0;
        step = 0;

        for k = 1:numel(t_list)
            t_target = t_list(k);

            while t < t_target

                dt = compute_dt(U, dx, Cmax, g);

                if t + dt > t_target
                    dt = t_target - t;
                end

                switch lower(time_scheme)
                    case 'fe'
                        U = FE(U, dx, dt, flux_type, g);
                    case 'rk4'
                        U = RK4(U, dx, dt, flux_type, g);
                end
                t = t + dt;
                % If the water depth becomes non-physical, flag the failure and record the current state.
                Hnow = U(1,:);
                if any(~isfinite(Hnow)) || any(~isreal(Hnow)) || any(Hnow <= 0)
                    warning('Non-physical at t=%.4e (scheme=%s, Cmax=%.1f)', t, time_scheme, Cmax);
                    failed = true;
                    fail_k = k;
                    fail_tau = t;
                    break;
                end
            end

            % Store snapshot
            Hk = real(U(1,:));
            Qk = real(U(2,:));
            uk = Qk ./ max(Hk, 1e-12);

            H_snapshot(k,:) = Hk;
            u_snapshot(k,:) = uk;
        end

        % Plots
        figure('Color','w');
        sgtitle(sprintf('%s scheme | flux=%s | Cmax=%.1f', upper(time_scheme), upper(flux_type), Cmax));

        for k = 1:numel(t_list)
            % H-x
            subplot(numel(t_list), 2, 2*k-1);
            plot(x_center, H_snapshot(k,:), 'LineWidth', 1.2);
            grid on; ylabel('H');
            title(sprintf('t = %.4f s', t_list(k)));
            if k == numel(t_list), xlabel('x'); end

            % u-x
            subplot(numel(t_list), 2, 2*k);
            plot(x_center, u_snapshot(k,:), 'LineWidth', 1.2);
            grid on; ylabel('u');
            title(sprintf('t = %.4f s', t_list(k)));
            if k == numel(t_list), xlabel('x'); end
        end
    end
end

%% Helper functions
% Compute the time step from the CFL condition
function dt = compute_dt(U, dx, Cmax, g)
    H = U(1,:);
    Q = U(2,:);
    u = Q ./ H;
    speed = abs(u) + sqrt(g * H); % Local maximum characteristic speed
    speed_max = max(speed);
    dt = Cmax * dx / speed_max;
end

% Forward Euler
function U_euler = FE(U, dx, dt, flux_type, g)
    R = horiki_paragogos(U, dx, flux_type, g);
    U_euler = U + dt * R;
end

% Fourth-order Runge-Kutta
function U_RK = RK4(U, dx, dt, flux_type, g)
    k1 = horiki_paragogos(U,            dx, flux_type, g);
    k2 = horiki_paragogos(U + (dt/2)*k1,  dx, flux_type, g);
    k3 = horiki_paragogos(U + (dt/2)*k2,  dx, flux_type, g);
    k4 = horiki_paragogos(U + dt*k3, dx, flux_type, g);
    U_RK = U + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% Compute the spatial residual R for time integration
function R = horiki_paragogos(U, dx, flux_type, g)

    [~, Nx] = size(U);

    % Fluxes at cell interfaces
    F = zeros(2, Nx+1);

    % Left boundary, x = -1
    H1 = U(1,1);

    F(:,1) = [ ...
        0;                  % Q = 0, zero discharge
        0.5 * g * H1^2      
    ];

    % Interior interfaces
    for i = 1:Nx-1
        UL = U(:,i);
        UR = U(:,i+1);
        F(:,i+1) = riemann_flux(UL, UR, flux_type, g);
    end

    % Right boundary, x = +1
    HN = U(1,Nx);

    F(:,Nx+1) = [ ...
        0;
        0.5 * g * HN^2
    ];

    % Residual
    R = -(F(:,2:end) - F(:,1:end-1)) / dx;
end

% Compute the numerical flux using the selected Riemann solver
function Flux = riemann_flux(UL, UR, flux_type, g)
    HL = UL(1); QL = UL(2); uL = QL / HL;
    HR = UR(1); QR = UR(2); uR = QR / HR;
    % Physical fluxes as functions of the conserved variables 
    Flux_L = [QL;
          QL^2/HL + 0.5 * g * HL^2];
    Flux_R = [QR;
          QR^2/HR + 0.5 * g * HR^2];

    switch lower(flux_type)
        case 'rusanov'
            speed_L   = abs(uL) + sqrt(g*HL);
            speed_R   = abs(uR) + sqrt(g*HR);
            speed_max = max(speed_L, speed_R);
            Flux = 0.5 * (Flux_L + Flux_R) - 0.5 * speed_max * (UR - UL);

        case 'roe'
            H_pavla = 0.5 * (HL + HR);   % H~
            u_pavla = (uL*sqrt(HL) + uR*sqrt(HR)) / (sqrt(HL) + sqrt(HR));  % u~ 

            c  = sqrt(g * H_pavla);      
            Fr = u_pavla / c;            % Froude number 

            a = 0.5 * (abs(Fr + 1) + abs(Fr - 1));   
            b = 0.5 * (abs(Fr + 1) - abs(Fr - 1));   

            % Roe matrix
            A_visc = [ c * (a - b*Fr),          b;
                       b * c^2 * (1 - Fr^2),  c * (a + b*Fr) ];

            % Roe flux
            Flux = 0.5 * (Flux_L + Flux_R - A_visc * (UR - UL));

        case 'hll'
            cL = sqrt(g * HL);
            cR = sqrt(g * HR);

            SL = min([0, uL - cL, uR - cR]);   
            SR = max([0, uL + cL, uR + cR]);  

            if SL >= 0
                % All waves propagate right; use the left (upwind) flux
                Flux = Flux_L;
            elseif SR <= 0
                % All waves propagate left; use the right (upwind) flux
                Flux = Flux_R;
            else
                % HLL flux 
                Flux = (SR * Flux_L - SL * Flux_R + SR * SL * (UR - UL)) / (SR - SL);
            end

        case 'hllroe'
            % Roe-averaged values for wave-speed estimation 
            H_pavla = 0.5 * (HL + HR);
            u_pavla = (uL*sqrt(HL) + uR*sqrt(HR)) / (sqrt(HL) + sqrt(HR));

            c_pavla = sqrt(g * H_pavla);
            SL = min(0, u_pavla - c_pavla);
            SR = max(0, u_pavla + c_pavla);

            if SL >= 0
                Flux = Flux_L;
            elseif SR <= 0
                Flux = Flux_R;
            else
                % Same numerical flux formula as standard HLL 
                Flux = (SR * Flux_L - SL * Flux_R + SR * SL * (UR - UL)) / (SR - SL);
            end
    end
end


