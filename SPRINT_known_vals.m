%% Load known global constants
SPRINT_global_vals

%% Martin Marietta report for contract AMMRC CTR 7447

%Combustion chamber exhaust gas (pdf page 22)
p_c = [2000, 3000] * 6894.76; %psia -> pa, pressure range
t_c = 3644.261; %6100 F -> K, combustion chamber gas temp

%FAE-7 propellant chemical composition by mass (Table II (C) pdf page 24)
FAE7_nitro_cellulose = 0.18;
FAE7_nitro_glycerine = 0.30;
FAE7_ammonium_perchlorate = 0.36;
FAE7_aluminum = 0.072; %in the form of metal staples
FAE7_triacetin = 0.067;
FAE7_resorcinol = 0.011;
FAE7_2_nitro_diphenylamine = 0.010;

%FAE-7 notable exhaust products by mass (Table II (C) pdf page 24)
FAE7_exhaust_Al2O3 = 0.134; %liquid aluminum oxide

%FAE-7 combustion chamber properties (Table II (C) pdf page 24)
FAE7_t_c = 3609.261; %6037 F -> K, combustion gas temp
FAE7_c_star_imperial = 5191; %ft/sec, characteristic exhaust velocity
FAE7_c_star_si = 1582.217; %ft/sec -> m/sec, characteristic exhaust velocity
FAE7_isp_vac = 305.3; %seconds, vacuum specific impulse
FAE7_gamma = 1.197; %ratio of specific heats
FAE7_effective_gas_constant = 52.63; %R, no units listed : - (
FAE7_effective_mol_weight = 29.356; %M_bar, no units listed : - (
FAE7_prop_density = 1710.61811; %0.0618 lb/in^3 -> kg/m^3, propellant density

% Derived quantities from this report
FAE7_R_si = Ru / FAE7_effective_mol_weight;
FAE7_c_f = FAE7_isp_vac * g0 / (FAE7_c_star_si); %coefficient of thrust in vacuum

%%  Bell labs ABM project history report (Oct 1975)

%Timeline info (pdf page 310)
TIMELINE_first_prototype_test_flight = datetime(1965, 3, 1); 
%day ambiguous, first launch of Propulsion Test Vehicle (pdf page 326)
TIMELINE_heat_shield_test_flight = datetime(1968, 6, 1); %date ambiguous (pdf page 326)
TIMELINE_first_flight_test = datetime(1965, 11, 17);
TIMELINE_first_deployment = datetime(1974, 6, 1); %day ambiguous

%Basic physical dimensions (pdf page 313 onwards)
SPRINT_length = 8.2296; %27 ft -> m, longitudinal length
SPRINT_base_diameter = 1.34112; %4.4 ft -> m, axial diameter at base
SPRINT_launch_mass = 3447.302; %7600 lb -> kg, mass at launch (approximate)

%derived system dimensions
SPRINT_cone_half_angle = rad2deg(atan(SPRINT_base_diameter / (2 * SPRINT_length))); 
%deg, half-angle of approximate cone body frame. Reported as 4 deg for
%interior motor casing
SPRINT_frontal_area = pi * (SPRINT_base_diameter/2)^2; %area of SPRINT head-on [m^2]
SPRINT_nozzle_area_limit = SPRINT_frontal_area;

SPRINT_total_volume = pi / 3 * SPRINT_length * (SPRINT_base_diameter/2)^2;
SPRINT_stage1_combustor_radius = (215 + 150) / (2*355) * SPRINT_base_diameter / 2; %mean radius. found with pixel counting
SPRINT_stage1_combustor_diameter = SPRINT_stage1_combustor_radius*2;
SPRINT_stage1_combustor_area = pi * (SPRINT_stage1_combustor_radius)^2;
SPRINT_stage1_combustor_height = 140/608*SPRINT_length;
SPRINT_stage1_combustor_volume = SPRINT_stage1_combustor_height * SPRINT_stage1_combustor_area; 
%approximating cylinder volume by pixel counting the schematic in the bell
%labs report

SPRINT_stage1_combustor_L_over_D = 140/608*SPRINT_length / (2*SPRINT_stage1_combustor_radius);

%% HIBEX Requirements from [https://www.alternatewars.com/WW3/WW3_Documents/DARPA/DARPA_II_HIBEX.htm]
% System known to use SPRINT first stage

HIBEX_burn_time = 1.124; %seconds
HIBEX_burnout_vel = 2562.758; %m/s after 1.2 sec burn
HIBEX_payload_mass = 136.078; %kg
HIBEX_STAGE1_prop_mass = 765.2103; %kg, we have reason to believe it carried a lot less propellant
HIBEX_max_acc = 377; %gs

%% HEDI Report (system known to use SPRINT stage 1 and 2 with different payload)
SPRINT_azimuth_uncertainty = 6; %deg, +/-

%% DARPA Report (includes some SPRINT info)
SPRINT_intercept_altitude = 13716; %m, intended intercept point (~45000 ft -> m)  [section 3.1]

%% Designation Systems page
STAGE1_thrust = 2900e3; %N, from Hercules X-265
STAGE1_time = 1.2; %seconds burn time
STAGE2_time = 3.8; %seconds burn time
SPRINT_flight_ceiling = 30e3; %m
SPRINT_flight_range = 40e3; %m
SPRINT_total_flight_time = 15; %s, upper limit
SPRINT_2nd_stage_burnout_mach = 10; %approximate, but quoted everywhere

%% Nuclear ABMs page [https://www.nuclearabms.info/Sprint.html]
SPRINT_heat_shield_temp = 3699.817; %k, due to atmospheric friction
SPRINT_heat_dissipation_rate = 850; %BTU/ft/sec, sorry about the shitty units lol
SPRINT_min_burnout_alt = 1500; %m, mentioned as lowest possible intercept altitude so we assume 
% burnout occurs near here

%derived atmospheric characteristics at launch and burnout, from std atmosphere in
%anderson aero text
SPRINT_ta_launch = ta; %k
SPRINT_pa_launch = pa; %pa
SPRINT_rho_launch = rhoa; %kg/m^3

SPRINT_ta_burnout = 278.41; %k
SPRINT_pa_burnout = 8.4560e4; %pa
SPRINT_rhoa_burnout = 1.0581; %kg/m^3

%derived flight conditions at burnout
SPRINT_a_burnout = 296; %[m/s] burnout speed of sound at ~10 km, we have to estimate burnout
% altitude as it's very variable
SPRINT_vel_burnout = SPRINT_2nd_stage_burnout_mach * SPRINT_a_burnout;

%% More derived values (from Liam's work)
SPRINT_stage_1_m_dot = STAGE1_thrust / (FAE7_c_f * FAE7_c_star_si);

%% Some random US patent that mentions HIBEX values [in references folder onedrive]
FAE7_burn_rate = 0.1; %[m/s] propellant burn rate
