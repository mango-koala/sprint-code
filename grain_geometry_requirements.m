%% Grain geometry params based on sizing
SPRINT_volumetric_loading = mp_stage_1 / (FAE7_prop_density * SPRINT_stage1_combustor_volume);
SPRINT_web_fraction = 2 * FAE7_burn_rate * STAGE1_time / (2*SPRINT_stage1_combustor_radius);
SPRINT_port_area = pi/4 * (2*SPRINT_stage1_combustor_radius)^2 * (1-SPRINT_volumetric_loading);
SPRINT_Ap_over_At_range = SPRINT_port_area ./ A_star_range;
