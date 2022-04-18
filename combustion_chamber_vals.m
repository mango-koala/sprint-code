
%% SOLVING STAGNATION CONDITIONS IN COMBUSTION CHAMBER
close all

syms t0 pc pe
c_star_eq = @(gamma, R) sqrt(R * t0 / gamma * ((gamma+1)/2)^((gamma+1)/(gamma-1)));
t0c = double(solve(FAE7_c_star_si == c_star_eq(FAE7_gamma, FAE7_R_si)));
% stag temp in combustion chamber

% since we now know t0, t, p we can find p0:
p0c = (t0c / t_c) ^ (FAE7_gamma / (FAE7_gamma-1)) * pc;

m_flux = mass_flux(FAE7_gamma, FAE7_R_si, t0c, p0c);

A_star = SPRINT_stage_1_m_dot / m_flux;

A_star_range = double(subs(A_star, pc, p_c));

pa_sea = 101325;
syms e

pe_space = linspace(pa_sea/2, pa_sea*2, 100);
e_space = linspace(1,20,100);
[e_gg, pe_gg] = meshgrid(e_space, pe_space);

cf_case = double(cf(FAE7_gamma, pe_gg, pa_sea, subs(p0c, mean(p_c)), e_gg));

% figure
% s = surf(e_gg, pe_gg/pa_sea, cf_case);
% s.EdgeAlpha = 0;
% texit("", "Expansion ratio $$\varepsilon$$", "Nozzle exit pressure fraction $$\frac{p_e}{p_a}$$ [pa]")
% zlabel("Thrust coefficient $$C_F$$")

e_val = double(solve(FAE7_isp_vac*g0 == isp_func(FAE7_gamma, FAE7_R_si, t0c, ...
    subs(p0c, mean(p_c)), 2*pa_sea, pa_sea, FAE7_c_star_si, e)));

A_e_range = 12*A_star_range;
diam_exit_range = 2*sqrt(A_e_range/pi);
SPRINT_base_diameter;

p0c_range = double(subs(p0c, p_c));
mass_flux_range = double(subs(m_flux, p_c));
