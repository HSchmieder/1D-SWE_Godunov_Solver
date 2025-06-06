%% 1D SWE Godunov solver
% TO DO
% Roe approximation

clc; clear;  close all;

%% Numerical Setting
nx = 20e2+1;            % number of cells
iter_max = 1e5;         % max number of time steps
CFL = 0.25;              
% available spatial schemes: 
% Lax-Friedrichs (LAXFR), Lax-Wendroff (LAXWE), FORCE
% GODUNOV, MUSCL, ENO3, WENO3, WENO5
set_spatialscheme = "MUSCL";
% available Riemann solvers: Rusanov (RUS), HLL, ROE (TODO)
set_riemann = "HLL";
% available time integrators: explicit Euler (EULER), RK2, RK4
set_timeintegrator = "RK4";

% set up
dfdx = get_spatialscheme(set_spatialscheme);
fct_rs = get_riemann(set_riemann);
fct_timeint = get_timeintegrator(set_timeintegrator);

%% Load Case
% define or chose case in get_case.m file
set_case = 1; 
[g, t, tend, dx, xL, xR, xcenter, q, b, fct_source, BC_type, fct_dirichletL, fct_dirichletR] = get_case(nx, set_case);
[fct_BCghost, fct_BCriemann] = get_BC(BC_type,set_spatialscheme,fct_dirichletL,fct_dirichletR);

%% Solve
for n=0:iter_max
    % set dt 
    smax = max(fct_eigenvalues(q,g),[],'all');
    dt = dx*CFL/smax;
    % stop exactly at tend
    if t + dt > tend
        dt = tend - t;
    end
    if t >= tend
        break
    end
    % calculate conserved variables at next time step
    qnew = fct_timeint(q, t, dx, dt, g, dfdx, fct_rs, fct_BCghost, fct_BCriemann, fct_source);

    % update variables
    q = qnew;
    t = t + dt;

    % plot every 100th time step
    if mod(n,100) == 0
        subplot(2,1,1); plot(xcenter,q(1,:)+b, xcenter, b)
        xlim([xL,xR]); ylabel('waterlevel h')
        %ylim([0,1.1])
        subplot(2,1,2); plot(xcenter,q(2,:))
        xlim([xL,xR]); ylabel('momentum hu')
        %ylim([0,1])
        sgtitle(sprintf('time = %f',t))
        drawnow
    end
end

% final plot
subplot(2,1,1); plot(xcenter,q(1,:)+b, xcenter, b.*ones(1,nx))
xlim([xL,xR]); ylabel('waterlevel h')
subplot(2,1,2); plot(xcenter,q(2,:))
xlim([xL,xR]); ylabel('momentum hu')
sgtitle(sprintf('time = %f',t))
