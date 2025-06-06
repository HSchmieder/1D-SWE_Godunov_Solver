function [g, t, tend, dx, xL, xR, xcenter, q, b, fct_source, BC_type, fct_dirichletL, fct_dirichletR] = get_case(nx,set_case)
%% set simulation case for 1D SWE
% Mandatory:
% gravitational acceleration g
% domain boundaries time: [tstart, tend], space: [xL, xR]
% initial state: q[2, nx], where h=q(1,:), hu=q(2,:) 
% B.C. type:
% 1: von Neumann: \partial_x = 0
% 2: reflective: \partial_x h = 0, hu = 0
% 3: periodic: h(xL,t)=h(xR,t), hu(xL,t)=hu(xR,t)
% 4: \partial_x h(xL,t)=0 , h(xR,t)=fct_dirichletR(t)
%    hu(xL,t)=fct_dirichletL(t) , \partial_x hu(xR,t)=0
% 5: h(xL,t) = fct_dirichletL(t) , \partial_x h(xR,t)=0 
%    \partial_x hu(xL,t)=0 , \partial_x hu(xR,t)=0
% Optional:
% bed elevation: b(1,nx) (needs a spatial derivative function of b for the source term)
% source term function: fct_source

    switch(set_case)
        case(1)
            %% Dam Break
            % domain
            g = 1;
            t = 0; tend = 1;
            xL = -1; xR = 1;
            dx = (xR-xL)/nx;
            xcenter = xL+dx/2:dx:xR-dx/2;
            x_middle = xcenter(1)+0.5*(xcenter(end)-xcenter(1));
            % I.C.
            q =zeros(2,length(xcenter));
            q(1,xcenter<=x_middle) = 0.8;
            q(1,xcenter>x_middle) = 0.4;
            % B.C.
            BC_type = 2;
        case(2)
            %% Double Dam Break
            % domain
            g = 1;
            t = 0; tend = 1;
            xL = -1; xR = 1;
            dx = (xR-xL)/nx;
            xcenter = xL+dx/2:dx:xR-dx/2;
            % I.C.
            q =zeros(2,length(xcenter));
            q(1,xcenter<=-0.33) = 0.3;
            q(1,xcenter>-0.33) = 0.1;
            q(1,xcenter>0.33) = 0.2;
            % B.C.
            BC_type = 1;
        case(3) 
            %% Initial Gaussian Water Level
            g = 1;
            t = 0; tend = 1;
            xL = -1; xR = 1;
            dx = (xR-xL)/nx;
            xcenter = xL+dx/2:dx:xR-dx/2;
            % I.C.
            q = zeros(2,length(xcenter));
            q(1,:) = fct_gaussian(xcenter,0.2,0.1,0.,0.5);
            % B.C.
            BC_type = 1;
        case(4)
            %% lake at rest
            % domain
            g = 1;
            t = 0; tend = 4;
            xL = -1; xR = 1;
            dx = (xR-xL)/nx;
            xcenter = xL+dx/2:dx:xR-dx/2;
            % source term
            b = 0.1+0.1*cos(2*pi/(xR-xL)*2.*xcenter);
            dbdx = -0.1*2*pi/(xR-xL)*2*sin(2*pi/(xR-xL)*2.*xcenter);
            fct_source = @(q) fct_bed_slope_term(q,g,dbdx);
            % I.C.
            waterlevel = 1;
            q =zeros(2,length(xcenter));
            q(1,:) = waterlevel-b;
            % B.C.
            BC_type = 2;
        case(5)
            %% gaussian flood
            % rough inclined channel with gaussian momentum inflow
            % friction term results from depth averaging the rough log-law
            % domain
            g = 9.81;
            t = 0; tend = 1000;
            xL = 0; xR = 5e3;
            dx = (xR-xL)/nx;
            xcenter = xL+dx/2:dx:xR-dx/2;
            % source term
            slope = 0.001;
            b = linspace(slope*(xR-xL),0,nx);
            dbdx = (b(2)-b(1))/dx*ones(size(xcenter));
            ks = 0.05;
            kappa = 0.41;
            fct_source = @(q) fct_bed_slope_roughness_term1(q,g,dbdx,ks,kappa);
            % I.C.
            h0 = 1;
            % initial velocity results from stationary (d/dt=0) and developed (d/dx=0) flow for given h0 
            u0 = sqrt(g*h0*abs(dbdx(1)))*((log(h0/ks)-1)/kappa+8.5);
            q = [h0; h0*u0].*ones(2,nx);
            % B.C.
            BC_type = 4;
            sigma = 100;
            m = h0*u0;
            a = 5*m;
            % shift temporal center of gaussian such that c=h(t=0,xL)/h0 - 1
            c = 0.01;
            t_shift = sqrt(-2*log(c))*sigma;
            fct_dirichletL = @(t) fct_gaussian(t,a,sigma,t_shift,m);
            fct_dirichletR = @(t) h0;
        case(6)
            %% oscillating left waterlevel
            % sinusoidal left waterlevel, gaussian bed, quadratic friction law
            % domain
            g = 9.81;
            t = 0; tend = 600;
            xL = -1e3; xR = 1e3;
            dx = (xR-xL)/nx;
            xcenter = xL+dx/2:dx:xR-dx/2;
            % source term
            % bed and bed slope
            a=2; sigma=200; shift=0; m=0;
            b = fct_gaussian(xcenter,a,sigma,shift,m);
            dbdx = -(xcenter-shift)/sigma^2.*fct_gaussian(xcenter,a,sigma,shift,m);
            % quadratic friction
            cd = 0.01;
            fct_source = @(q) fct_bed_slope_roughness_term2(q,g,dbdx,cd);
            % I.C.
            h0 = 4;
            u0 = 0;
            q = [h0-b; h0*u0*ones(1,nx)];
            % B.C.
            BC_type = 5;
            a = 1;
            f = 0.005;
            fct_dirichletL = @(t) h0+a*sin(f*2*pi*t);
    end
        
    if ~exist('b','var')
        b = zeros(size(xcenter));
    end
    if ~exist('fct_source','var')
        fct_source = @(q) [0;0];
    end
    if ~exist('fct_dirichletL','var')
        fct_dirichletL = [];
    end
    if ~exist('fct_dirichletR','var')
        fct_dirichletR = [];
    end
end

%% function Gauss (case3, case5, case6)
function y = fct_gaussian(x,a,sigma,shift,m)
    if ~(length(a)==length(sigma) && length(a)==length(shift))
        print('dimensions of Gauss features do not match')
        return
    end
    y = m*ones(size(x));
    for i = 1:length(a)
        y = y + a(i)*exp(-(x-shift(i)).^2./(2*sigma(i)^2));
    end
end

%% source term considering bed slope (case4)
function source_value = fct_bed_slope_term(q,g,dbdx)
    source_value = zeros(size(q));
    source_value(2,:) = -g*q(1,:).*dbdx;
end

%% source term considering bed slope and roughness (case5)
% friction term results from depth averaging the rough log-law
function source_value = fct_bed_slope_roughness_term1(q,g,dbdx,ks,kappa)
    source_value = zeros(size(q));
    u = fct_cons2prim(q);
    source_value(2,:) = -g*q(1,:).*dbdx - u(2,:).*abs(u(2,:)).*((log(u(1,:)./ks)-1)/kappa+8.5).^-2;
end

%% source term considering bed slope and roughness (case6)
% simple quadratic friction loss
function source_value = fct_bed_slope_roughness_term2(q,g,dbdx,cd)
    source_value = zeros(size(q));
    u = fct_cons2prim(q);
    source_value(2,:) = -g*q(1,:).*dbdx - cd*u(2,:).*abs(u(2,:));
end