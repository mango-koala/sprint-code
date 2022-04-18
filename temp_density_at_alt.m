function [T, rho] = temp_density_at_alt(lat, lon, h_of_t)
    [T,rho] = atmosnrlmsise00(h_of_t,lat,lon,2022,4,0);
    T = T(:,2); %temp at alt, not exoatmospheric
    rho = rho(:,6); %selects air density in kg/m^3
end
