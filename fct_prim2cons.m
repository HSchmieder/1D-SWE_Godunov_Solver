function q = fct_prim2cons(u)
% return conserved variables given primitive variables for 1D-SWE
q = [u(1,:); u(1,:).*u(2,:)];