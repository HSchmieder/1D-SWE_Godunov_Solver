function [fct_BCghostcells, fct_BCriemann] = get_BC(BC_type, spatialscheme, fct_dirichletL, fct_dirichletR)
switch(BC_type)
    case(1) % Neumann B.C.
        if spatialscheme == "ENO3" || spatialscheme == "WENO5"
            fct_BCghostcells = @(q,~) [q(:,1), q(:,1), q, q(:,end), q(:,end)];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,1), q_aug(:,end)];
        elseif spatialscheme == "MUSCL" || spatialscheme == "WENO3"
            fct_BCghostcells = @(q,~) [q(:,1), q, q(:,end)];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,1), q_aug(:,end)];
        else
            fct_BCghostcells = @(q,~) [q(:,1), q, q(:,end)];
            fct_BCriemann = [];
        end
    case(2) % reflective/ mirrored B.C. (Neumann B.C. for h and Dirichlet B.C. hu=0)
        if spatialscheme == "ENO3" || spatialscheme == "WENO5"
            fct_BCghostcells = @(q,~) [[q(1,2), q(1,1); -q(2,2), -q(2,1)], q, [q(1,end), q(1,end-1); -q(2,end), -q(2,end-1)]];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,2), q_aug(:,end-1)];
        elseif spatialscheme == "MUSCL" || spatialscheme == "WENO3"
            fct_BCghostcells = @(q,~) [[q(1,2); -q(2,2)], q, [q(1,end-1); -q(2,end-1)]];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,1), q_aug(:,end)];
        else
            fct_BCghostcells = @(q,~) [[q(1,1); -q(2,1)], q, [q(1,end); -q(2,end)]];
            fct_BCriemann = [];
        end 
    case(3) % periodic B.C.
        if spatialscheme == "ENO3" || spatialscheme == "WENO5"
            fct_BCghostcells = @(q,~) [q(:,end-1:end), q, q(:,1:2)];
            fct_BCriemann = @(~,wL1,wRend) [wL1, wRend];
        elseif spatialscheme == "MUSCL" || spatialscheme == "WENO3"
            fct_BCghostcells = @(q,~) [q(:,end), q, q(:,1)];
            fct_BCriemann = @(~,wL1,wRend) [wL1, wRend];
        else
            fct_BCghostcells = @(q,~) [q(:,end), q, q(:,1)];
            fct_BCriemann = [];
        end
    case(4) % classic river flow B.C.
            % Inflow (L): \partial_x h(xL,t)=0 , hu(xL,t)=fct_dirichletL(t)
            % Outflow (R): h(xR,t)=fct_dirichletR(t) , \partial_x hu(xR,t)=0
        if spatialscheme == "ENO3" || spatialscheme == "WENO5"
            fct_BCghostcells = @(q,t) [[q(1,1), q(1,1); fct_dirichletL(t), fct_dirichletL(t)], q, [fct_dirichletR(t), fct_dirichletR(t); q(2,end), q(2,end)]];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,1), q_aug(:,end)];
        elseif spatialscheme == "MUSCL" || spatialscheme == "WENO3"
            fct_BCghostcells = @(q,t) [[q(1,1); fct_dirichletL(t)], q, [fct_dirichletR(t); q(2,end)]];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,1), q_aug(:,end)];
        else
            fct_BCghostcells = @(q,t) [[q(1,1); fct_dirichletL(t)], q, [fct_dirichletR(t); q(2,end)]];
            fct_BCriemann = [];
        end
    case(5) % Inflow (L): h(xL,t) = fct_dirichletL(t) , \partial_x hu(xL,t)=0
            % Outflow (R): \partial_x h(xR,t)=0 , \partial_x hu(xR,t)=0
        if spatialscheme == "ENO3" || spatialscheme == "WENO5"
            fct_BCghostcells = @(q,t) [[fct_dirichletL(t), fct_dirichletL(t); q(2,1), q(2,1)], q, q(:,end), q(:,end)];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,1), q_aug(:,end)];
        elseif spatialscheme == "MUSCL" || spatialscheme == "WENO3"
            fct_BCghostcells = @(q,t) [[fct_dirichletL(t); q(2,1)], q, q(:,end)];
            fct_BCriemann = @(q_aug,~,~) [q_aug(:,1), q_aug(:,end)];
        else
            fct_BCghostcells = @(q,t) [[fct_dirichletL(t); q(2,1)], q, q(:,end)];
            fct_BCriemann = [];
        end
end