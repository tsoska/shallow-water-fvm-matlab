% Question 2
clear; close all; clc;

%% Spatial discretization
nx   = 1000;
x_min = -1; x_max = 1;
dx   = (x_max - x_min)/nx;
x_center = linspace(x_min + dx/2, x_max - dx/2, nx);

%% Parameters
g = 9.81;
HL = 50; HR = 1;
t_final = 0.02;
Cmax = 0.9;

solvers = ["rusanov","roe","hll","hllroe"];
time_schemes = ["fe","rk4"];

%% Exact solution
[H_ananlytical, u_analytical] = swe_analytical(x_center, t_final, HL, HR, g);
H_ananlytical = H_ananlytical(:).'; 
u_analytical = u_analytical(:).';

%% Loop over time integration schemes and Riemann solvers
for time_scheme = time_schemes

    % Plot the exact solution
    figure_H = figure('Color','w'); hold on; grid on;
    title(sprintf('H(x) at t = %.3f s | %s | Cmax = %.1f', ...
        t_final, upper(time_scheme), Cmax));
    xlabel('x'); ylabel('H');
    plot(x_center, H_ananlytical, 'k--', 'LineWidth', 2, 'DisplayName', 'Analytic');

    figure_U = figure('Color','w'); hold on; grid on;
    title(sprintf('u(x) at t = %.3f s | %s | Cmax = %.1f', ...
        t_final, upper(time_scheme), Cmax));
    xlabel('x'); ylabel('u');
    plot(x_center, u_analytical, 'k--', 'LineWidth', 2, 'DisplayName', 'Analytic');

    fprintf('\nErrors at t = %.3f | scheme = %s\n', t_final, upper(time_scheme));

    for s = 1:numel(solvers)
        flux_type = solvers(s);

        % Initial conditions
        H = HL*ones(1,nx);
        H(x_center > 0) = HR;
        Q = zeros(1,nx);
        U = [H; Q];
        t = 0;

        % Time integration
        while t < t_final
            dt = compute_dt(U, dx, Cmax, g);
            if t + dt > t_final
                dt = t_final - t;
            end

            switch time_scheme
                case "fe"
                    U = FE(U, dx, dt, flux_type, g);
                case "rk4"
                    U = RK4(U, dx, dt, flux_type, g);
            end
            t = t + dt;
        end

        % Numerical solution
        H_numerical = real(U(1,:));
        u_numerical = real(U(2,:)) ./ max(H_numerical,1e-12);

        % Pointwise errors
        eH = H_numerical - H_ananlytical;
        eU = u_numerical - u_analytical;

        % Print RMSE values 
        fprintf('%8s | H_L2 = %.3e | u_L2 = %.3e\n', ...
            upper(flux_type), sqrt(mean(eH.^2)), sqrt(mean(eU.^2)));

        % Plot the numerical solution
        figure(figure_H);
        plot(x_center, H_numerical, 'LineWidth', 1.2, ...
            'DisplayName', upper(flux_type));

        figure(figure_U);
        plot(x_center, u_numerical, 'LineWidth', 1.2, ...
            'DisplayName', upper(flux_type));
    end

    figure(figure_H); legend('show','Location','best');
    figure(figure_U); legend('show','Location','best');
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

% Calculate the analytical solution 
function [H,u,cm,xA,xB,xC] = swe_analytical(x, t, HL, HR, g)

    x = x(:);  % Convert to a column vector

    % At t = 0, the initial water depth is discontinuous at x = 0
    if t == 0
        H = HL*ones(size(x));
        H(x>0) = HR;
        u = zeros(size(x));
        cm = NaN; xA = NaN; xB = NaN; xC = NaN;
        return
    end

    % Wave speeds 
    cL = sqrt(g*HL);
    cR = sqrt(g*HR);

    % Numerical calculation of cm
    f = @(cm) -8*(cR^2).*(cm.^2).*(cL - cm).^2 + (cm.^2 - cR^2).^2 .* (cm.^2 + cR^2);

    % Initial guess 
    cm0 = 0.5*(cL + cR);

    cm = NaN;
    try
        cm = fzero(f, cm0);
    catch
        % Minimize f^2 over the interval [cR, cL]
        obj = @(cm) f(cm).^2;
        cm = fminbnd(obj, cR, cL);
    end

    % xA(t), xB(t), xC(t): boundaries between solution regions
    xA = -t*cL;
    xB =  t*(2*cL - 3*cm);
    xC =  t * ( 2*(cm^2)*(cL - cm) ) / ( (cm^2) - (cR^2) );

    H = zeros(size(x));
    u = zeros(size(x));

    % x < xA: left undisturbed region, constant depth and zero velocity
    I1 = (x < xA);
    H(I1) = HL;
    u(I1) = 0;

    % xA <= x < xB: rarefaction region with parabolic depth profile
    I2 = (x >= xA) & (x < xB);
    H(I2) = (4/(9*g)) * (cL - x(I2)/(2*t)).^2;
    u(I2) = (2/3) * (x(I2)/t + cL);

    % xB <= x < xC: intermediate constant state
    I3 = (x >= xB) & (x < xC);
    H(I3) = (cm^2)/g;
    u(I3) = 2*(cL - cm);

    % x >= xC: right undisturbed region, constant depth and zero velocity
    I4 = (x >= xC);
    H(I4) = HR;
    u(I4) = 0;

end
