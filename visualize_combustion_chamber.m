close all; clear; clc; format;

%Loading known values
SPRINT_known_vals

%Loading useful equations
SPRINT_useful_equations

%Loading known global values
SPRINT_global_vals

w = 0.55; %outer buffer
r1 = 1 - w;
r2 = 0.1; %inner buffer

N = 5; %symmetry number
zeta = pi/N * 0.8;
secant_fillet_angle = pi/N - zeta;

pix_grit = 300;
pix_range = linspace(-1, 1, pix_grit);
[pxx, ~] = meshgrid(pix_range, pix_range);
pyy = fliplr(pxx)';
pix_real_width = SPRINT_stage1_combustor_diameter/pix_grit; %width of each pixel in m
pixels_burned_per_sec = FAE7_burn_rate / pix_real_width; %pixels burned (with given burn rate) per second

pix_r = sqrt(pxx.^2 + pyy.^2);
pix_theta = atan2(pyy, pxx);
pzz = ones(size(pxx));

theta_min = 0;
theta_max = pi / N;
section_cond = abs(pix_theta) < theta_max & pix_theta > 0;
section_theta = pix_theta(section_cond);
section_r = pix_r(section_cond);

high_theta_r = r1;
low_theta_2 = r2;

pzz(pix_r > r1+w) = -1;
pzz(pix_r < r2) = 0;


%% Star-ish geometry
for i = 1:length(pzz(:))
    pix_val = pzz(i);
        
    pr = pix_r(i);
    pt = pix_theta(i);
    thinkness_factor = 8;
    symmetry_factor = 2;
    pt_mod = mod(pt*symmetry_factor, 2*thinkness_factor*pi/N);

    min_r = (r1-r2)*sin(N*pt_mod);

    if pr < min_r
        pzz(i) = 0;
    end
end


f=figure;
writerObj = initVid("burnback.mp4", 20);

tsteps = 100;
tspan = linspace(0,STAGE1_time, tsteps);
pzz_hist = zeros(pix_grit, pix_grit, tsteps);
tstep_real_time = mean(diff(tspan));
straight_pixels_burned_per_tstep = tstep_real_time * pixels_burned_per_sec;
diag_pixels_burned_per_tstep = straight_pixels_burned_per_tstep * sqrt(2)/2;

for tinds = 1:length(tspan)
    pzz_next = pzz;

    if sum(pzz == 1, 'all') ~= 0
        for i = 2:size(pzz,1)-1
            for j = 2:size(pzz,2)-1
                ns = pix_neighbors(pzz, i, j);
                if ~isnan(ns) %occurs if we're on a border
                    if ns(2,2) == 1 %if current pixel is unburned
                        ns_diag = [ns(3,3), ns(1,3), ns(3,1), ns(1,1)];
                        ns_straight = [ns(1,2), ns(2,1), ns(3,2), ns(2,3)];
                        n_all = [ns_diag, ns_straight];
    
                        if any(~ns_diag, 'all') && rand < diag_pixels_burned_per_tstep
                            %there is a burned pixel next to this
                            pzz_next(i,j) = 0; %burn this pixel
                        end
                        if any(~ns_straight, 'all') && rand < straight_pixels_burned_per_tstep
                            pzz_next(i,j) = 0; %burn this pixel
                        end
    
                        if sum(n_all) == 0 && sum(n_all < 0) == 0
                            pzz_next(i,j) = 0; %burn this pixel
                        end
                    end
                end
            end
        end
    end

    burned_area(tinds) = 1/0.7766*sum(pzz == 0, 'all') / numel(pzz) * SPRINT_stage1_combustor_area;
    burned_volume(tinds) = burned_area(tinds) * SPRINT_stage1_combustor_height;
    burned_mass(tinds) = burned_volume(tinds) * FAE7_prop_density;
    
    pzz_hist(:, :, tinds) = pzz;
    pzz = pzz_next;
    
    im = imagesc(pzz);
    hold on
    im_overlay = imagesc(pzz, 'alphadata', pzz ~= 0);
    im_overlay.CData = im_overlay.CData;
    set(gca, 'colormap', copper)
    f.Color = [1,1,1]; 
    set(gca, 'visible', 'off')
    drawnow
    saveFrame(gca, writerObj)
end
close(writerObj)

mass_flow_rate = diff(burned_mass) ./ diff(tspan);
mean_flow_rate = mean(mass_flow_rate);
mass_flow_rate_at_t = @(t) interp1(tspan(1:end-1), mass_flow_rate, t);
volumetric_loading = 1 - burned_area(1)/SPRINT_stage1_combustor_area;
mass_flow_rate_error = (mean_flow_rate - SPRINT_stage_1_m_dot) / SPRINT_stage_1_m_dot * 100
total_burn_time_to_empty = (burned_mass(end) - STAGE1_time * SPRINT_stage_1_m_dot) / mean_flow_rate + 1.2


t=initTiles(2,2);
nexttile
plot(tspan, burned_mass, 'linewidth', 2.5)
texit("Burned Propellant Mass", "Time in [$$s$$]", "Mass in [$$kg$$]")

nexttile
plot(tspan(1:end-1), mass_flow_rate, 'linewidth', 2.5)
texit("Propellant Mass Flow Rate", "Time in [$$s$$]", "Mass Flow Rate in [$$\frac{kg}{s}$$]")
ylim([0, Inf])

nexttile
plot(tspan, burned_volume, 'linewidth', 2.5)
texit("Open Combustor Volume", "Time in [$$s$$]", "Volume [$$m^3$$]")

nexttile
cla
set(gca, 'colormap', gray)
step_alpha = 0.4;
plot_intermediate = [80 60 40 20];
for i = plot_intermediate
    im = imagesc(pzz_hist(:,:,i) + 3*i);
    im.AlphaData = pzz_hist(:,:,i) == 0;
    im.CDataMapping = 'direct';
    hold on
end
im_init = imagesc(pzz_hist(:,:,1));
im_init.AlphaData = pzz_hist(:,:,1) == 0;

im_overlay = imagesc(pzz_hist(:,:,end), 'alphadata', pzz_hist(:,:,end) ~= 0);
im_overlay.CData = -im_overlay.CData;
texit("Grain Burn Progress", "", "")
grid off
shg

function nvals = pix_neighbors(pzz, i,j)
    if i > 1 &&i < size(pzz,1)
        if j > 1 && j < size(pzz,2)
            nvals = pzz(i-1:i+1, j-1:j+1);
            return
        end
    end
    nvals = NaN;
end
