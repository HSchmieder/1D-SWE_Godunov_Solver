[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=HSchmieder/1D-SWE_Godunov_Solver&file=https://github.com/HSchmieder/1D-SWE_Godunov_Solver/SWE1D_main.m)

# 1D-SWE_Godunov_Solver

This code solves the 1D shallow water equations using a semi-discrete finite volume scheme of Godunov-type.  
[mass](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20h%7D%7B%5Cpartial%20t%7D%20&plus;%20%5Cfrac%7B%5Cpartial%20(hu)%7D%7B%5Cpartial%20x%7D%20%3D%200)

![momentum](https://latex.codecogs.com/svg.image?%5Cfrac%7B%5Cpartial%20(hu)%7D%7B%5Cpartial%20t%7D%20&plus;%20%5Cfrac%7B%5Cpartial%7Bhu%5E2%7D%7D%7B%5Cpartial%20x%7D%20&plus;%20%5Cfrac%7B1%7D%7B2%7Dgh%5E2%20%3D%20%5Cleft(-gh%20%5Cfrac%7B%5Cpartial%20b%7D%7B%5Cpartial%20x%7D%20&plus;%20S_f%5Cright))

$$
\frac{\partial}{\partial t} h + \frac{\partial}{\partial x} hu = 0
$$
$ \frac{\partial}{\partial t} hu + \frac{\partial}{\partial x} hu^2 + \frac{1}{2} gh^2 = \left( -gh \frac{\partial}{\partial x} b + S_f \right) $

In the SWE1D_main.m file one can set the numerical settings (number of cells: nx, CFL-number: CFL), load a case from get_case.m (set_case) and chose between different Godunov-fluxes (set_spatialscheme), approximate Riemann solvers (set_riemann) and time-integratior (set_timeintegrator).

__Godunov Fluxes:__  
Lax-Friedrichs: LAXFR  
Lax-Wendroff: LAXWE  
FORCE: FORCE  
2nd-order monotonic upstream- centered scheme for conservation laws with minmod slope-limiter: MUSCL  
3rd order essentially non-oscillatory with oscillation indicator (Jiang and Shu 1996): ENO3  
3rd order weighted essentially non-oscillatory with oscillation indicator (Jiang and Shu 1996):WENO3  
5th order weigthed essentially non-oscillatory with oscillation indicator (Jiang and Shu 1996): WENO5  
  
__Approximate Riemann solver:__  
Rusanov: RUS  
Harten-Lax-van Leer: HLL  
Roe: ROE (to do)  

__Time- Integrators:__  
explicite Euler: EULER  
2nd-order Runge-Kutta: RK2  
4th-order Runge-Kutta: RK4
