function fct_rs = get_riemann(set_rs)
    switch set_rs
        case('RUS')
            fct_rs = @fct_rus;
        case('ROE')
            fct_rs = @fct_roe;
        case('HLL')
            fct_rs = @fct_hll;
    end
end

%% approximate Rusanov solver
function flux = fct_rus(wL,wR,g)
    smax = max([fct_eigenvalues(wL,g); fct_eigenvalues(wR,g)],[],1);
    flux = 0.5*(fct_flux(wL,g)+fct_flux(wR,g)) - 0.5*smax.*(wR-wL);
end

%% approximate HLL-solver (Harten-Lax-van Leer)
function flux = fct_hll(wL,wR,g)
    lambdaL = fct_eigenvalues(wL,g);
    lambdaR = fct_eigenvalues(wR,g);
    sL = min(0, min([lambdaL; lambdaR],[],1));
    sR = max(0, max([lambdaL; lambdaR],[],1));
    flux = (sR.*fct_flux(wL,g)-sL.*fct_flux(wR,g))./(sR-sL) ...
        + sR.*sL./(sR-sL).*(wR-wL);
end

%% approximate Roe-solver
function flux = fct_roe(wL,wR,g)
    % TODO
end