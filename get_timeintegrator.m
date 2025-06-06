function fct_timeint = get_timeintegrator(set_timeintegrator)
    switch set_timeintegrator
        case("EULER")
            fct_timeint = @fct_euler;
        case("RK2")
            fct_timeint = @fct_rk2;
        case("RK4")
            fct_timeint = @fct_rk4;
    end
end

% explicit euler
function qnew = fct_euler(q,t,dx,dt,g,dfdx, fct_rs, fct_BCghost, fct_BCriemann, fct_source)
    qnew = q + dt * (dfdx(q,t,dx,dt,g, fct_rs, fct_BCghost, fct_BCriemann) + fct_source(q));
end

% 2nd order Runge Kutta
function qnew = fct_rk2(q,t,dx,dt,g,dfdx, fct_rs, fct_BCghost, fct_BCriemann, fct_source)
    k1 = dfdx(q, t, dx, dt, g, fct_rs, fct_BCghost, fct_BCriemann) + fct_source(q);
    k2 = dfdx(q+dt/2*k1, t, dx, dt, g, fct_rs, fct_BCghost, fct_BCriemann) + fct_source(q+dt/2*k1);
    qnew = q + dt * k2;
end

% 4th order Runge Kutta
function qnew = fct_rk4(q,t,dx,dt,g,dfdx, fct_rs, fct_BCghost, fct_BCriemann, fct_source)
    k1 = dfdx(q,t,dx,dt,g, fct_rs, fct_BCghost, fct_BCriemann) + fct_source(q);
    k2 = dfdx(q+dt/2*k1,t,dx,dt,g, fct_rs, fct_BCghost, fct_BCriemann) + fct_source(q+dt/2*k1);
    k3 = dfdx(q+dt/2*k2,t,dx,dt,g, fct_rs, fct_BCghost, fct_BCriemann) + fct_source(q+dt/2*k2);
    k4 = dfdx(q+dt*k3,t,dx,dt,g, fct_rs, fct_BCghost, fct_BCriemann) + fct_source(q+dt*k3);
    qnew = q + dt/6 * (k1+2*k2+2*k3+k4);
end