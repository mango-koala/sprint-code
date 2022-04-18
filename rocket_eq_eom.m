function [dy] = rocket_eq_eom(g0, Isp, mdot, y, launch_el)
    r = y(1);
    v = y(2);
    m = y(3);
    dv = -g0*Isp*mdot / m - sind(launch_el)*g0*0; %acc due to gravity turned off for now
    dm = mdot;
    dy = [v; dv; dm];
end
