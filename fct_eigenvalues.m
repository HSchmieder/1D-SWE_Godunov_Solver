function lambda = fct_eigenvalues(q,g)
% calculates flux eigenvalues for 1D-SWE 
u = fct_cons2prim(q);
c = sqrt(g*u(1,:));
lambda =[u(2,:)-c; u(2,:)+c];