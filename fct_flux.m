function flux = fct_flux(q,g)
% calculates the flux for 1D-SWE
u = fct_cons2prim(q); % mainly to check for water levels <= 0
flux = [q(2,:); u(2,:).*q(2,:)+0.5*g*q(1,:).^2];
