function u = fct_cons2prim(q)
% return primitive variables for 1D-SWE given conserved variables
if any(q(1,:)<=0)
    disp('negative or zero waterlevel is found')
    return
end
u = [q(1,:); q(2,:)./q(1,:)];
