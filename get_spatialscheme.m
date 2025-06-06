function dfdx = get_spatialscheme(set_spatialscheme)
    switch(set_spatialscheme)
        case("LAXFR")
            dfdx = @fct_laxfr;
        case("LAXWE")
            dfdx = @fct_laxwe;
        case("FORCE")
            dfdx = @fct_force;
        case("GODUNOV")
            dfdx = @fct_godunov;
        case("MUSCL")
            dfdx = @fct_muscl;
        case("WENO3")
            dfdx = @fct_weno3;
        case("ENO3")
            dfdx = @fct_eno3;
        case("WENO5")
            dfdx = @fct_weno5;
    end
end

%% Lax-Friedrich scheme
function dfdx = fct_laxfr(q,t,dx,dt,g, ~, fct_BCghost, ~)
    % expand by ghost cells for boundary treatment
    q_aug = fct_BCghost(q,t);
    f_aug = fct_flux(q_aug,g);
    % fluxes at left (m) and right (p) interface of each cell 
    f = 0.5*(f_aug(:,1:end-1)+f_aug(:,2:end)) - 0.5*dx/dt*(q_aug(:,2:end)-q_aug(:,1:end-1));
    fm = f(:,1:end-1);
    fp = f(:,2:end);
    % calculate spatial flux gradient
    dfdx = -1/dx*(fp-fm);
end

%% Lax-Wendroff scheme
function dfdx = fct_laxwe(q,t,dx,dt,g, ~, fct_BCghost, ~)
    % expand by ghost cells for boundary treatment
    q_aug = fct_BCghost(q,t);
    f_aug = fct_flux(q_aug,g);
    % fluxes at left (m) and right (p) interface of each cell 
    q_LW = 0.5*(q_aug(:,1:end-1)+q_aug(:,2:end))-0.5*dt/dx*(f_aug(:,2:end)-f_aug(:,1:end-1));
    f = fct_flux(q_LW,g);
    fm = f(:,1:end-1);
    fp = f(:,2:end);
    % calculate spatial flux gradient
    dfdx = -1/dx*(fp-fm);
end

%% FORCE scheme
function dfdx = fct_force(q,t,dx,dt,g, ~, fct_BCghost, ~)
    % get Lax-Wendroff and Lax-Friedrich spatial flux gradients
    dqdt_LW = fct_laxwe(q,t,dx,dt,g, [], fct_BCghost, []);
    dqdt_FR = fct_laxfr(q,t,dx,dt,g, [], fct_BCghost, []);
    % calculate spatial flux gradient
    dfdx = 0.5*(dqdt_FR+dqdt_LW);
end

%% 1st order Godunov scheme
function dfdx = fct_godunov(q,t,dx,~,g, fct_rs, fct_BCghost, ~)
    % expand by ghost cells for boundary treatment
    q_aug = fct_BCghost(q,t);
    
    % approximmate values at each side of nCell+1 interfaces
    wL = q_aug(:,1:end-1);
    wR = q_aug(:,2:end);
    
    % solve Riemann problem at interfaces
    f = fct_rs(wL,wR,g);
     % fluxes at left (m) and right (p) interface of each cell 
    fm = f(:,1:end-1);
    fp = f(:,2:end);
    % calculate spatial flux gradient
    dfdx = -1/dx*(fp-fm);
end

%% 2nd order MUSCL scheme
function dfdx = fct_muscl(q,t,dx,~,g, fct_rs, fct_BCghost, fct_BCriemann)
    % expand by ghost cells for boundary treatment
    q_aug = fct_BCghost(q,t);
    
    % calculate polynomial constants p(z)=a0*z + a1 where z=x-x_i
    a1 = q; % a1 is q_i for cell i for both stencils
    % stencil 1: i-1, i
    a01 = (q_aug(:,2:end-1)-q_aug(:,1:end-2))/dx;
    % stencil 2: i, i+1
    a02= (q_aug(:,3:end)-q_aug(:,2:end-1))/dx;
    % chose final polynomial (chose one of the 2 slopes)
    a0_f = fct_minmod(a01,a02);
    
    % approximmate values at each side of nCell+1 interfaces
    wL = a0_f*dx/2 + a1;
    wR = -a0_f*dx/2 + a1;
    % apply B.C. for left/right state of most left/right interface
    wBC = fct_BCriemann(q_aug,wL(:,end),wR(:,1));
    wL = [wBC(:,1), wL];
    wR = [wR, wBC(:,2)];

    % solve Riemann problem at interfaces
    f = fct_rs(wL,wR,g);
    % fluxes at left (m) and right (p) interface of each cell 
    fm = f(:,1:end-1);
    fp = f(:,2:end);
    % calculate spatial flux gradient
    dfdx = -1/dx*(fp-fm);
end

% minmod slope limiter
function slope = fct_minmod(a,b)
    slope = 0.5 * (sign(a) + sign(b)) .* min(abs(a), abs(b));
end

%% 3rd order WENO scheme
function dfdx = fct_weno3(q,t,dx,~,g, fct_rs, fct_BCghost, fct_BCriemann)
    % expand by ghost cells for boundary treatment
    q_aug = fct_BCghost(q,t);
    
    % calculate polynomial constants p(z)=a0*z + a1 where z=x-x_i
    a1 = q; % a1 is q_i for cell i for both stencils
    % stencil 1: i-1, i
    a01 = (q_aug(:,2:end-1)-q_aug(:,1:end-2))/dx;
    % stencil 2: i, i+1
    a02= (q_aug(:,3:end)-q_aug(:,2:end-1))/dx;
    
    % calculate oscillation indicators after Jiang and Shu 1996
    % stencil 1: i-1, i
    IS1 = (q_aug(:,2:end-1)-q_aug(:,1:end-2)).^2;
    % stencil 2: i, i+1
    IS2 = (q_aug(:,3:end)-q_aug(:,2:end-1)).^2;
    
    % approximmate values at each side of nCell+1 interfaces
    % by weighted combination of both polynomials
    epsilon = 1e-7;
    p = 2;
    % left side
    c1 = 1/3;
    alpha1 = c1./(epsilon+IS1).^p;
    c2 = 2/3;
    alpha2 = c2./(epsilon+IS2).^p;
    wL = (alpha1.*(a01*dx/2+a1) + alpha2.*(a02*dx/2+a1)) ...
        ./(alpha1 + alpha2);
    % right side
    c1 = 2/3;
    alpha1 = c1./(epsilon+IS1).^p;
    c2 = 1/3;
    alpha2 = c2./(epsilon+IS2).^p;
    wR = (alpha1.*(-a01*dx/2+a1) + alpha2.*(-a02*dx/2+a1)) ...
        ./(alpha1 + alpha2);
    % apply B.C. for left/right state of most left/right interface
    wBC = fct_BCriemann(q_aug,wL(:,end),wR(:,1));
    wL = [wBC(:,1), wL];
    wR = [wR, wBC(:,2)];
    
    % solve Riemann problem at interfaces
    f = fct_rs(wL,wR,g);
    % fluxes at left (m) and right (p) interface of each cell 
    fm = f(:,1:end-1);
    fp = f(:,2:end);
    % calculate spatial flux gradient
    dfdx = -1/dx*(fp-fm);
end

%% 3rd order ENO scheme
function dfdx = fct_eno3(q,t,dx,~,g, fct_rs, fct_BCghost, fct_BCriemann)
    [nVar, nCell] = size(q);
    % ghost cells for boundary treatment
    q_aug = fct_BCghost(q,t);
    
    % calculate polynomial constants p(z)=a0*z^2+a1*z+a0 where z = x- x_{center}
    a0 = zeros([3, nVar, nCell]);
    a1 = a0; a2 = a0; IS = a0;
    % stencil 1: i-2,i-1,i, (x_{center} = x_{i-1})
    [a0(1,:,:), a1(1,:,:), a2(1,:,:)] = fct_polyfactors(dx,q_aug(:,1:end-4),q_aug(:,2:end-3),q_aug(:,3:end-2));
    % stencil 2: i-1,i,i+1, (x_{center} = x_i)
    [a0(2,:,:), a1(2,:,:), a2(2,:,:)] = fct_polyfactors(dx,q_aug(:,2:end-3),q_aug(:,3:end-2),q_aug(:,4:end-1));
    % stencil 3: i,i+1,i+2, (x_{center} = x_{i+1})
    [a0(3,:,:), a1(3,:,:), a2(3,:,:)] = fct_polyfactors(dx,q_aug(:,3:end-2),q_aug(:,4:end-1),q_aug(:,5:end));
    
    % calculate oscillation indicators for all polynomials
    IS(1,:,:) = fct_IS(dx,0.5*dx,1.5*dx,a0(1,:,:),a1(1,:,:));
    IS(2,:,:) = fct_IS(dx,-0.5*dx,0.5*dx,a0(2,:,:),a1(2,:,:));
    IS(3,:,:) = fct_IS(dx,-1.5*dx,-0.5*dx,a0(3,:,:),a1(3,:,:));
    % find index of smoothest polynomial (least oscillating)
    [~, idx_h] = min(squeeze(IS(:,1,:)),[],1);
    [~, idx_hu] = min(squeeze(IS(:,2,:)),[],1);
    idx_linh = sub2ind([3, nVar, nCell],idx_h, ones(1,nCell),1:nCell);
    idx_linhu = sub2ind([3, nVar, nCell],idx_hu, 2*ones(1,nCell),1:nCell);
    
    % approximmate values at each side of nCell+1 interfaces
    % left side
    dxh = dx*(2.5-idx_h); % location of interface (right side from cell i point of view) w.r.t. x_{center} of chosen polynomial
    dxhu = dx*(2.5-idx_hu);
    wL = [a0(idx_linh).*dxh.^2+a1(idx_linh).*dxh+a2(idx_linh); ...
        a0(idx_linhu).*dxhu.^2+a1(idx_linhu).*dxhu+a2(idx_linhu)];
    % right side
    dxh = dx*(1.5-idx_h); % location of interface (left side from cell i point of view) w.r.t. x_{center} of chosen polynomial
    dxhu = dx*(1.5-idx_hu);
    wR = [a0(idx_linh).*dxh.^2+a1(idx_linh).*dxh+a2(idx_linh); ...
        a0(idx_linhu).*dxhu.^2+a1(idx_linhu).*dxhu+a2(idx_linhu)];
    % apply B.C. for left/right state of most left/right interface
    wBC = fct_BCriemann(q_aug,wL(:,end),wR(:,1));
    wL = [wBC(:,1), wL];
    wR = [wR, wBC(:,2)];

    % solve Riemann problem at interfaces
    f = fct_rs(wL,wR,g);
     % fluxes at left (m) and right (p) interface of each cell 
    fm = f(:,1:end-1);
    fp = f(:,2:end);
    % calculate spatial flux gradient
    dfdx = -1/dx*(fp-fm);
end

%% 5th order WENO scheme
function dfdx = fct_weno5(q,t,dx,~,g, fct_rs, fct_BCghost, fct_BCriemann)
    [nVar, nCell] = size(q);
    % expand by ghost cells for boundary treatment
    q_aug = fct_BCghost(q,t);
    
    % calculate polynomial constants p(z)=a0*z^2+a1*z+a0 (z = x- x_{center})
    a0 = zeros([3, nVar, nCell]);
    a1 = a0; a2 = a0; IS = a0;
    % stencil 1: i-2,i-1,i, (x_{center} = x_{i-1})
    [a0(1,:,:), a1(1,:,:), a2(1,:,:)] = fct_polyfactors(dx,q_aug(:,1:end-4),q_aug(:,2:end-3),q_aug(:,3:end-2));
    % stencil 2: i-1,i,i+1, (x_{center} = x_i)
    [a0(2,:,:), a1(2,:,:), a2(2,:,:)] = fct_polyfactors(dx,q_aug(:,2:end-3),q_aug(:,3:end-2),q_aug(:,4:end-1));
    % stencil 3: i,i+1,i+2, (x_{center} = x_{i+1})
    [a0(3,:,:), a1(3,:,:), a2(3,:,:)] = fct_polyfactors(dx,q_aug(:,3:end-2),q_aug(:,4:end-1),q_aug(:,5:end));
    
    % calculate oscillation indicators for all polynomials
    IS(1,:,:) = fct_IS(dx,0.5*dx,1.5*dx,a0(1,:,:),a1(1,:,:));
    IS(2,:,:) = fct_IS(dx,-0.5*dx,0.5*dx,a0(2,:,:),a1(2,:,:));
    IS(3,:,:) = fct_IS(dx,-1.5*dx,-0.5*dx,a0(3,:,:),a1(3,:,:));
    
    % approximmate values at each side of nCell+1 interfaces
    % by weighted combination of all 3 polynomials
    p = 2; 
    epsilon = 1e-7;
    % left side
    c1 = 3/10; c2 = 6/10; c3 = 2/10;
    alpha1 = c1./(epsilon + IS(1,:,:)).^p;
    alpha2 = c2./(epsilon + IS(2,:,:)).^p;
    alpha3 = c3./(epsilon + IS(3,:,:)).^p;
    wL = squeeze( ...
        (alpha1.*(a0(1,:,:)*(1.5*dx)^2 + a1(1,:,:)*(1.5*dx) + a2(1,:,:)) ...
        + alpha2.*(a0(2,:,:)*(0.5*dx)^2 + a1(2,:,:)*(0.5*dx) + a2(2,:,:)) ...
        + alpha3.*(a0(3,:,:)*(-0.5*dx)^2 + a1(3,:,:)*(-0.5*dx) + a2(3,:,:))) ...
        ./ (alpha1 + alpha2 + alpha3));
    % right side
    c1 = 1/10; c2 = 6/10; c3 = 3/10;
    alpha1 = c1./(epsilon + IS(1,:,:)).^p;
    alpha2 = c2./(epsilon + IS(2,:,:)).^p;
    alpha3 = c3./(epsilon + IS(3,:,:)).^p;
    wR = squeeze( ...
        (alpha1.*(a0(1,:,:)*(1.5*dx)^2 + a1(1,:,:)*(1.5*dx) + a2(1,:,:)) ...
        + alpha2.*(a0(2,:,:)*(-0.5*dx)^2 + a1(2,:,:)*(-0.5*dx) + a2(2,:,:)) ...
        + alpha3.*(a0(3,:,:)*(-1.5*dx)^2 + a1(3,:,:)*(-1.5*dx) + a2(3,:,:))) ...
        ./ (alpha1 + alpha2 + alpha3));
    % apply B.C. for left/right state of most left/right interface
    wBC = fct_BCriemann(q_aug,wL(:,end),wR(:,1));
    wL = [wBC(:,1), wL];
    wR = [wR, wBC(:,2)];

    % solve Riemann problem at interfaces
    f = fct_rs(wL,wR,g);
     % fluxes at left (m) and right (p) interface of each cell 
    fm = f(:,1:end-1);
    fp = f(:,2:end);
    % calculate spatial flux gradient
    dfdx = -1/dx*(fp-fm);
end

%%
% polynomial factors of 2nd order for WENO5 and ENO3 scheme
function [a0, a1, a2] = fct_polyfactors(dx,qm,q,qp)
    a0 = (qm - 2*q + qp)/(2*dx^2);
    a1 = (-qm + qp)/(2*dx);
    a2 = (-qm + 26*q -qp)/24;
end

% Oscillation Indicator after Jiang and Shu 1996 for WENO5 and ENO3 scheme
function SI = fct_IS(dx,L,R,a0,a1)
% for 2nd order reconstruction polynomial p(z)=a0*z^2+a1*z+a2
% z = x - x_center; x_center: x position of middle cell of each stencil
% L = x_{i-0.5} - x_center; R = x_{i+0.5} - x_center
SI = a0.^2*(4/3*dx*(R^3-L^3)+4*dx^3*(R-L)) + 2*a0.*a1*dx*(R^2-L^2);
end

