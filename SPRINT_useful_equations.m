%% Characteristic exhaust velocity for isentropic flow
c_star = @(gamma, R, t0) sqrt(R * t0 / gamma * ((gamma+1)/2)^((gamma+1)/(gamma-1)));

%% Mass flux for isentropic flow
mass_flux = @(gamma, R, t0, p0) p0 * sqrt(gamma) / sqrt(R * t0) * (2/(gamma+1))^((gamma+1)/(2*gamma-2));

%% Coefficient of thrust for isentropic flow
cf = @(gamma, pe, pa, p0, epsilon) sqrt(2*gamma^2 ./ (gamma-1) .* (2/(gamma+1)) ^ ((gamma+1)/(gamma-1)) ...
    .*(1-(pe./p0).^((gamma-1)/gamma))) + (pe-pa)./p0 .* epsilon;

isp_func = @(g, R, t0, p0, pe, pa, c_star, e) ...
    sqrt(2*g/(g-1)*R*t0.*(1-(pe./p0)^((g-1)/g))) + (pe-pa)./p0.*c_star.*e;
