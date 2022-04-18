%% numerical integration for flight path
close all
launch_el = 70; %assumed launch angle from horizontal
tpts = 1e3; %plot point count

%integrating stage 1
ics1 = [0; 0; mo1]; tspan1 = [0,STAGE1_time]; 
v1_sol = ode45(@(t, y) rocket_eq_eom(g0, Isp1, -SPRINT_stage_1_m_dot, y, launch_el), tspan1, ics1);

t_plot1 = linspace(0,STAGE1_time,tpts);
y_plot1 = deval(v1_sol, t_plot1);
v_stage1 = y_plot1(2,:)';
acc_of_t1 = diff(y_plot1(2,:)) ./ diff(t_plot1);
h_of_t1 = sind(launch_el)*y_plot1(1,:);

[T1, rho1] = temp_density_at_alt(50, 50, h_of_t1);
p1 = rho1 .* R_air .* T1;
a1 = sqrt(gamma_air * R_air * T1);
M1 = v_stage1 ./ a1;

%integrating stage 2
ics2 = [y_plot1(1,end); v_stage1(end); m2_tot]; tspan2 = [STAGE1_time, STAGE1_time + STAGE2_time]; 
v2_sol = ode45(@(t, y) rocket_eq_eom(g0, Isp2, -mdot_2, y, launch_el), tspan2, ics2);

t_plot2 = linspace(tspan2(1),tspan2(2),tpts);
y_plot2 = deval(v2_sol, t_plot2);
v_stage2 = y_plot2(2,:)';
acc_of_t2 = diff(y_plot2(2,:)) ./ diff(t_plot2);
h_of_t2 = sind(launch_el)*y_plot2(1,:);

[T2, rho2] = temp_density_at_alt(50, 50, h_of_t2);
p2 = rho2 .* R_air .* T2;
a2 = sqrt(gamma_air * R_air * T2);
M2 = v_stage2 ./ a2;

% t=initTiles(2,2);
% legendarr = ["Stage 1", "Stage 2"];
% % f=figure;
% nexttile
% plot(t_plot1, M1, 'linewidth', 2.5)
% hold on
% plot(t_plot2, M2, 'linewidth', 2.5)
% texit("Mach Number", "Time in [$$s$$]", "Mach number", legendarr, 'northwest')
% 
% % f=figure;
% nexttile
% plot(t_plot1(1:end-1), acc_of_t1/g0, 'linewidth', 2.5)
% hold on
% plot(t_plot2(1:end-1), acc_of_t2/g0, 'linewidth', 2.5)
% texit("Acceleration", "Time in [$$s$$]", "Acceleration in [$$g_0$$]", legendarr, 'northwest')
% 
% % f=figure;
% nexttile
% plot(t_plot1, v_stage1, 'linewidth', 2.5)
% hold on
% plot(t_plot2, v_stage2, 'linewidth', 2.5)
% texit("Velocity", "Time in [$$s$$]", "Velocity in [$$\frac{m}{s}$$]", legendarr, 'northwest')
% 
% % f=figure;
% nexttile
% plot(t_plot1, h_of_t1, 'linewidth', 2.5)
% hold on
% plot(t_plot2, h_of_t2, 'linewidth', 2.5)
% texit("Altitude", "Time in [$$s$$]", "Altitude in [$$m$$]", legendarr, 'northwest')
% 
% t=initTiles(2,2);
% nexttile
% plot(t_plot1, T1, 'linewidth', 2.5)
% hold on
% plot(t_plot2, T2, 'linewidth', 2.5)
% texit("Ambient Temperature", "Time in [$$s$$]", "Temperature in [$$K$$]", legendarr, 'southwest')
% 
% nexttile
% plot(t_plot1, rho1, 'linewidth', 2.5)
% hold on
% plot(t_plot2, rho2, 'linewidth', 2.5)
% texit("Ambient Density", "Time in [$$s$$]", "Density in [$$\frac{kg}{m^3}$$]")
% 
% nexttile
% plot(t_plot1, p1, 'linewidth', 2.5)
% hold on
% plot(t_plot2, p2, 'linewidth', 2.5)
% texit("Ambient Pressure", "Time in [$$s$$]", "Pressure in [$$Pa$$]")
% 
% nexttile
% plot(t_plot1, a1, 'linewidth', 2.5)
% hold on
% plot(t_plot2, a2, 'linewidth', 2.5)
% texit("Ambient Speed of Sound", "Time in [$$s$$]", "Speed of Sound in [$$\frac{m}{s}$$]")
