SPRINT_Load_All_COde

% mass vector for stage 1
m1_v = mo1 - SPRINT_stage_1_m_dot * t_plot1;
% mass vector for stage 2
m2_v = m2_tot - mdot_2 * t_plot2;
% overall mass vector
m_v = [m1_v m2_v(2:1000)]';

% combined time vector
time = [t_plot1 t_plot2(2:1000)]';

% free stream temperature
T_inf = [T1;T2(2:1000)];

% free stream Mach number
M_inf = [M1;M2(2:1000)];

% free stream pressure
p_inf = [p1;p2(2:1000)];

% free stream stagnation temperature
T0_inf = T_inf .* (1 + 0.2 .* (M_inf .^ 2));

% free stream stagnation pressure
p0_inf = p_inf .* (1 + 0.2 .* (M_inf .^ 2)) .^ (1.4 ./ 0.4);

% local stagnation pressure and temperature vectors
p0_L = p0_inf;
T0_L = T0_inf;

% dynamic pressure
q = 0.5 .* 1.4 .* p_inf .* (M_inf .^ 2);

% coefficient of pressure
Cp = abs((p0_L - p_inf) ./ q);
Cp(1) = 0;

% local pressure
%p_L = q .* Cp .* (cosd(45) .^ 2) + p_inf;
p_L = p_inf;

% local Mach number
M_L = real(sqrt((2 ./ 0.4) .* ((p0_L ./ p_L) .^ (0.4 ./ 1.4) - 1)));

% local temperature
T_L = T0_L ./ (1 + 0.2 .* (M_L .^ 2));

% dynamic viscosity
mu = 1.458e-6 .* (T_L .^ 1.5) ./ (T_L + 110.4);

% specific heat of air
cp_air = 1000;

% thermal conductivity
%k = 5.75e-5 .* (1 + 0.00317 .* (T_L - 273.15) - 0.0000021 .* ((T_L - 273.15) .^ 2)) * 100 * 4.184;
k = 0.025;

% Prandtl number
Pr = mu .* cp_air ./ k;

% recovery parameter
r = Pr .^ (1 ./ 3);

% recovery temperature
T_R = T0_inf .* (r + ((1 - r) ./ (1 + 0.2 .* M_L .^ 2)));

% structure cp
cp_rocket = 1465.3793658;

% adiabatic wall temperature
T_AW = T_L .* (1 + 0.2 .* r .* (M_L .^ 2));

% Boltzmann constant
sigma = 1.380649e-23;

% emmisivity
epsilon = 0.79;

% density
rho_L = p_L ./ (287 .* T_L);

% Reynolds number
Re = rho_L .* M_L .* (1 + 0.5 .* 0.4 .* (M_L .^ 2)) .^ (2.4 ./ 0.8) .* sqrt(1.4 .* 287 .* T0_L) .* SPRINT_length ./ mu;

% Nusselt number
Nu = 0.0292 .* (Re .^ 0.8) .* (Pr .^ (1 ./ 3));

% convective heat transfer coefficient
h = k .* Nu ./ SPRINT_length;

% solve the energy balance
IC = T_L(1);
[time_sol,T_S] = ode45(@(t,T_S) nrg_bal(t, T_S, m_v, cp_rocket, h, T_AW, sigma, epsilon, T_R, time), time, IC);

% plot skin temperature
plot(time_sol,T_S,'linewidth',2.5);
title('Surface Tempearture');
xlabel("Time (s)");
ylabel('Surface Temperature (K)');
grid on;

function T_sdot = nrg_bal(t, T_S, m, cp_rocket, h, T_AW, sigma, epsilon, T_R, time)
    hi = interp1(time,h,t);
    T_AWi = interp1(time,T_AW,t);
    m_vi = interp1(time,m,t);
    T_Ri = interp1(time,T_R,t);
    
    T_sdot = (hi * (T_AWi - T_S) - sigma * epsilon * (T_S ^ 4 - T_Ri ^ 4)) / (m_vi * cp_rocket);
end