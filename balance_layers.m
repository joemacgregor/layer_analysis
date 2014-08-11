% BALANCE_LAYERS Map balance velocity history across Greenland.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/04/14

clear

age_min                     = 3.5e3;
age_max                     = 9e3;
decim                       = 5; % decimation, indices but also effectively km
filt_type                   = 'exp'; % triang or exp
scale_type                  = 'thick'; % ind or thick
ind_filt                    = 5; % width of filter, indices, used only when scale_type=ind, should be an odd integer
std_exp                     = 4; % decay length of filtering exponential, indices, used only when filt_type=exp and scale_type=ind
thick_filt_ref              = 5; % number of ice thicknesses over which to run filter, used only when scale_type=thick
plotting                    = false;

% miscellaneous variables
load mat/grl_coast num_coast x_coast y_coast
load mat/grl_basin num_basin x_basin y_basin
load mat/mask_D1 mask_D1 x_D1_ord y_D1_ord
load mat/gimp_1km elev_surf_gimp
load mat/core_int name_core name_core_short num_core x_core y_core

% x/y limits
[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);

% surface speeds
load mat/speed_greenland_2007-2010 x_grd y_grd speed_x speed_y
[x_speed, y_speed]          = deal(x_grd, y_grd);
clear x_grd y_grd

% trim speed grid
[speed_x, speed_y]          = deal(speed_x(((y_speed(:, 1) >= y_min) & (y_speed(:, 1) <= y_max)), ((x_speed(1, :) >= x_min) & (x_speed(1, :) <= x_max))), speed_y(((y_speed(:, 1) >= y_min) & (y_speed(:, 1) <= y_max)), ((x_speed(1, :) >= x_min) & (x_speed(1, :) <= x_max))));
[x_speed, y_speed]          = deal(x_speed(((y_speed(:, 1) >= y_min) & (y_speed(:, 1) <= y_max)), ((x_speed(1, :) >= x_min) & (x_speed(1, :) <= x_max))), y_speed(((y_speed(:, 1) >= y_min) & (y_speed(:, 1) <= y_max)), ((x_speed(1, :) >= x_min) & (x_speed(1, :) <= x_max))));

% clean up speeds a bit to deal with NaNs in the south
[ind_clean_x, ind_clean_y]  = deal(find((x_speed(1, :) >= -100) & (x_speed(1, :) <= 300)), find((y_speed(:, 1) >= -2750) & (y_speed(:, 1) <= -2000)));
for ii = ind_clean_y'
    if ~any(isnan(speed_x(ii, ind_clean_x)))
        continue
    end
    speed_x(ii, ind_clean_x(isnan(speed_x(ii, ind_clean_x)))) ...
                            = interp1(x_speed(ii, ind_clean_x(~isnan(speed_x(ii, ind_clean_x)))), speed_x(ii, ind_clean_x(~isnan(speed_x(ii, ind_clean_x)))), x_speed(ii, ind_clean_x(isnan(speed_x(ii, ind_clean_x)))));
end
for ii = ind_clean_x
    if ~any(isnan(speed_y(ind_clean_y, ii)))
        continue
    end
    speed_y(ind_clean_y(isnan(speed_y(ind_clean_y, ii))), ii) ...
                            = interp1(y_speed(ind_clean_y(~isnan(speed_y(ind_clean_y, ii))), ii), speed_y(ind_clean_y(~isnan(speed_y(ind_clean_y, ii))), ii), y_speed(ind_clean_y(isnan(speed_y(ind_clean_y, ii))), ii));
end

speed                       = sqrt((speed_x .^ 2) + (speed_y .^ 2));

% load Bamber et al. (2013) data
load mat/greenland_bed_v3 mask_gris mask_land x y
[mask_gris, mask_land]      = deal(mask_gris(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))), mask_land(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))));
clear x y
mask_combo                  = double(mask_land);
mask_combo(mask_gris)       = 2;

% load Morlighem et al. (2014)
load mat/greenland_mc_bed_1km thick

% load CISM data
load mat/greenland_cism accum x y
accum_RACMO2                = accum(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max)));
clear accum x y

% load gridded isochrones
% load mat/age_grd2 age_iso depth_iso2 x_grd y_grd
load mat/age_grd2_krige age_iso depth_iso2 x_grd y_grd
[num_y, num_x]              = size(x_grd);

% valid isochrone ages
ind_age_iso                 = find((age_iso >= age_min) & (age_iso <= age_max));
num_iso                     = length(ind_age_iso);

% trim gridded isochrone set
age_iso                     = age_iso(ind_age_iso);
depth_iso                   = depth_iso2(:, :, ind_age_iso);
clear depth_iso2

% load and sort strain-rate modeling
% [accum, strain_rate, thick_shear] ...
[strain_rate, thick_shear]  = deal(NaN(size(x_grd, 1), size(x_grd, 2), num_iso));
% for ii = 1:num_iso
%     tmp                     = load(['mat/strain_all_' num2str(1e-3 * age_iso(ii)) 'ka.mat'], 'accum_shallow_grd', 'strain_rate_grd', 'thick_shear_grd');
%     [accum(:, :, ii), strain_rate(:, :, ii), thick_shear(:, :, ii)] ...
%                             = deal(tmp.accum_shallow_grd, tmp.strain_rate_grd, tmp.thick_shear_grd);
% end
load mat/strain_all_krige accum_shallow_grd strain_rate_grd thick_shear_grd
% load mat/strain_all_9ka strain_rate_grd thick_shear_grd
[accum, strain_rate(:, :, end), thick_shear(:, :, end)] ...
                            = deal(accum_shallow_grd, strain_rate_grd(:, :, end), thick_shear_grd);
clear accum_shallow_grd strain_rate_grd thick_shear_grd

% load Marcott et al. (2013) temperature reconstruction
load mat/marcott_2013_temp age_marcott temp_global temp_global_uncert temp_nh temp_nh_uncert

disp('Trimming and decimating non-model data...')

% decimate regular grids
length_grd                  = 1e3 * decim; % grid size, m
ind_x                       = find(~mod(x_grd(1, :), decim), 1):decim:find(~mod(x_grd(1, :), decim), 1, 'last');
ind_y                       = find(~mod(y_grd(:, 1), decim), 1):decim:find(~mod(y_grd(:, 1), decim), 1, 'last');
[x_filt, y_filt]            = deal(x_grd(ind_y, ind_x), y_grd(ind_y, ind_x));

% grid size
[num_decim_y, num_decim_x]  = deal(length(ind_y), length(ind_x));

% % decimate all standard grids
[accum_RACMO2_decim, elev_surf_gimp_decim, mask_gris_decim, mask_combo_decim, thick_decim, x_decim, y_decim, mask_D1_decim] ...
                            = deal(accum_RACMO2(ind_y, ind_x), elev_surf_gimp(ind_y, ind_x), mask_gris(ind_y, ind_x), mask_combo(ind_y, ind_x), thick(ind_y, ind_x), x_grd(ind_y, ind_x), y_grd(ind_y, ind_x), mask_D1(ind_y, ind_x));

% decimate speed grid
% ind_x_speed                 = find(~mod(x_speed(1, :), decim), 1):(2 * decim):find(~mod(x_speed(1, :), decim), 1, 'last');
% ind_y_speed                 = find(~mod(y_speed(:, 1), decim), 1):(2 * decim):find(~mod(y_speed(:, 1), decim), 1, 'last');
[speed_x_decim, speed_y_decim, x_speed_decim, y_speed_decim] ...
                            = deal(speed_x(1:2:end, 1:2:end), speed_y(1:2:end, 1:2:end), x_speed(1:2:end, 1:2:end), y_speed(1:2:end, 1:2:end));
speed_decim                 = sqrt((speed_x_decim .^ 2) + (speed_y_decim .^ 2));

disp('Building smoothing filter...')

% build smoothing filter
switch filt_type
    case 'triang' % triangle/chapeau
        switch scale_type
            case 'ind' % fixed number of indices
                filt_decim  = repmat(triang(ind_filt)', ind_filt, 1) .* repmat(triang(ind_filt), 1, ind_filt);
                filt_decim  = filt_decim ./ sum(filt_decim(:));
            case 'thick' % variable thickness filter
                % only need number of filter options depending on ice thickness range
                filt_decim  = cell(1, max([1 round(1e-3 * max(thick(:)) * thick_filt_ref)]));
                for ii = 1:length(filt_decim)
                    if iseven(ii)
                        jj  = ii + 1;
                    else
                        jj  = ii;
                    end
                    filt_decim{ii} ...
                            = repmat(triang(jj)', jj, 1) .* repmat(triang(jj), 1, jj);
                    filt_decim{ii} ...
                            = filt_decim{ii} ./ sum(filt_decim{ii}(:));
                end
        end
    case 'exp' % exponential decay
        switch scale_type
            case 'ind'
                filt_decim  = repmat(exp([-(((ind_filt - 1) / 2):-1:1) 0 -(1:((ind_filt - 1) / 2))] ./ std_exp), ind_filt, 1) .* repmat(exp([-(((ind_filt - 1) / 2):-1:1) 0 -(1:((ind_filt - 1) / 2))] ./ std_exp)', 1, ind_filt);
                filt_decim  = filt_decim ./ sum(filt_decim(:));
            case 'thick'
                filt_decim  = cell(1, max([1 round(1e-3 * max(thick(:)) * thick_filt_ref)]));
                for ii = 1:length(filt_decim)
                    filt_decim{ii} ...
                            = repmat(exp([-(ii:-1:1) 0 -(1:ii)] ./ thick_filt_ref), ((2 * ii) + 1), 1) .* repmat(exp([-(ii:-1:1) 0 -(1:ii)] ./ thick_filt_ref)', 1, ((2 * ii) + 1));
                    filt_decim{ii} ...
                            = filt_decim{ii} ./ sum(filt_decim{ii}(:));
                end
        end
end

disp('Filtering non-model data...')

% filter critical decimated grids
vars                        = {'elev_surf_gimp' 'speed_x_decim' 'speed_y_decim' 'thick'};
vars_filt                   = {'elev_surf_gimp_filt' 'speed_x_filt' 'speed_y_filt' 'thick_filt'};
num_var                     = length(vars);

switch scale_type
    case 'ind'
        for ii = 1:num_var
            eval([vars_filt{ii} ' = filter2(filt_decim, ' vars{ii} ');'])
        end
        filt_decim_ref      = [];
    case 'thick'
        
        for ii = 1:num_var
            eval([vars_filt{ii} ' = NaN(num_decim_y, num_decim_x);'])
        end
        
        % reference set of ice thicknesses
        thick_ref           = (1:length(filt_decim)) ./ (1e-3 * thick_filt_ref);
        % closest matching set of ice thicknesses
        filt_decim_ref      = interp1(thick_ref, 1:length(filt_decim), thick(ind_y, ind_x), 'nearest', 'extrap');
        
        for jj = 1:num_decim_x
            if ~mod(jj, 50)
                disp(jj)
            end
            for ii = 1:num_decim_y
                
                % current filter matrix and half-size
                filt_decim_curr ...
                            = filt_decim{filt_decim_ref(ii, jj)};
                num_filt_curr ...
                            = (size(filt_decim_curr, 1) - 1) / 2;
                
                % indices within original matrix
                ind_tmp     = -num_filt_curr:num_filt_curr;
                ii_tmp      = repmat((ind_y(ii) + ind_tmp)', ((2 * num_filt_curr) + 1), 1);
                jj_tmp      = ind_x(jj) + ind_tmp(ones(((2 * num_filt_curr) + 1), 1), :);
                jj_tmp      = jj_tmp(:);
                ind_good    = find((ii_tmp > 0) & (jj_tmp > 0) & (ii_tmp <= num_y) & (jj_tmp <= num_x));
                ind_curr    = sub2ind([num_y num_x], ii_tmp(ind_good), jj_tmp(ind_good));
                
                % filter each variable
                for kk = 1:num_var
                    eval([vars_filt{kk} '(ii, jj) = nansum(filt_decim_curr(ind_good) .* ' vars{kk} '(ind_curr)) * ((numel(filt_decim_curr) ^ 2) / (length(ind_good) * length(find(~isnan(' vars{kk} '(ind_curr))))));'])
                end
            end
        end
end

% filtered modern surface speed
speed_filt                  = sqrt((speed_x_filt .^ 2) + (speed_y_filt .^ 2));

% apply GrIS mask
[accum_RACMO2_decim(~mask_gris_decim), elev_surf_gimp_filt(~mask_gris_decim)] ...
                            = deal(NaN);

% surface elevation hollow filling
elev_surf_hollow            = elev_surf_gimp_filt;
for ii = 2:(num_decim_y - 1)
    for jj = 2:(num_decim_x - 1)
        curr_min            = min([elev_surf_hollow((ii - 1):(ii + 1), (jj - 1)); elev_surf_hollow([(ii - 1) (ii + 1)], jj); elev_surf_hollow((ii - 1):(ii + 1), (jj + 1))]);
        if (elev_surf_hollow(ii, jj) < curr_min)
            elev_surf_hollow(ii, jj) ...
                            = curr_min;
        end
    end
end

% surface-velocity and elevation-gradient flow azimuths
az_speed                    = atan2(speed_y_filt, speed_x_filt);
[elev_grad_x, elev_grad_y]  = gradient(elev_surf_hollow, length_grd, length_grd);
az_elev                     = atan2(-elev_grad_y, -elev_grad_x);

% smooth azimuths, extract useful elements
[az_sin_cat, az_cos_cat]    = deal(NaN(num_decim_y, num_decim_x, 2));
az_sin_cat(:, :, 1)         = sin(az_speed);
az_cos_cat(:, :, 1)         = cos(az_speed);
az_sin_cat(:, :, 2)         = sin(az_elev);
az_cos_cat(:, :, 2)         = cos(az_elev);

% weight filter exponentially using reference speed
speed_az_decay              = 25;
wt_az_elev                  = exp(-speed_filt ./ speed_az_decay);
wt_az_speed                 = 1 - wt_az_elev;
wt_az_elev(isnan(az_speed)) = 1;
wt_az_speed(isnan(az_elev)) = 1;
az_mean                     = atan2(((az_sin_cat(:, :, 1) .* wt_az_speed) + (az_sin_cat(:, :, 2) .* wt_az_elev)), ((az_cos_cat(:, :, 1) .* wt_az_speed) + (az_cos_cat(:, :, 2) .* wt_az_elev)));

% % old way just simple average
% az_mean                     = atan2(nanmean(az_sin_cat, 3), nanmean(az_cos_cat, 3));

az_sin                      = sin(az_mean);
az_cos                      = cos(az_mean);
az_sin_sign                 = sign(az_sin);
az_cos_sign                 = sign(az_cos);
az_sin_abs_rel              = abs(az_sin) ./ (abs(az_sin) + abs(az_cos));
az_cos_abs_rel              = abs(az_cos) ./ (abs(az_sin) + abs(az_cos));

% reproject speeds
[speed_x_filt, speed_y_filt]= deal((speed_filt .* az_cos), (speed_filt .* az_sin));
speed_filt                  = sqrt((speed_x_filt .^ 2) + (speed_y_filt .^ 2));

% sort grid indices in terms of descending elevation
[~, ind_sort]               = sort(elev_surf_hollow(:), 1, 'descend');
ind_sort                    = ind_sort(~isnan(ind_sort));

disp('Decimating and filtering strain-rate models, then doing flux balance...')
%%
[accum_filt, depth_filt, shape_filt, speed_bal, speed_bal_surf, strain_rate_filt, thick_shear_filt, vel_vert_filt] ...
                            = deal(NaN(num_decim_y, num_decim_x, num_iso));
for ii = num_iso
    disp([num2str(1e-3 * age_iso(ii)) ' ka...'])
    [speed_bal(:, :, ii), speed_bal_surf(:, :, ii), accum_filt(:, :, ii), strain_rate_filt(:, :, ii), thick_shear_filt(:, :, ii), shape_filt(:, :, ii), depth_filt(:, :, ii), vel_vert_filt(:, :, ii)] ...
                            = bal_vel(age_iso(ii), x_grd, y_grd, squeeze(accum(:, :, ii)), squeeze(depth_iso(:, :, ii)), squeeze(strain_rate(:, :, end)), squeeze(thick_shear(:, :, end)), mask_D1_decim, mask_gris_decim, thick_decim, az_sin_abs_rel, az_cos_abs_rel, az_sin_sign, az_cos_sign, ...
                                      ind_sort, length_grd, decim, true, scale_type, filt_decim, filt_decim_ref);
end

% [speed_bal_full, speed_bal_full_surf] ...
%                             = bal_vel(NaN, x_grd, y_grd, accum_RACMO2_decim, thick_decim, [], squeeze(thick_shear(:, :, end)), mask_D1_decim, mask_gris_decim, thick_decim, az_sin_abs_rel, az_cos_abs_rel, az_sin_sign, az_cos_sign, ind_sort, length_grd, decim, true, scale_type, filt_decim);

% % convert long-term averages to short-term history
% [accum_hist, speed_bal_surf_hist] ...
%                             = deal(NaN(num_decim_y, num_decim_x, (num_iso - 1)));
% for ii = (num_iso - 1):-1:1
%     speed_bal_surf_hist(:, :, ii) ...
%                             = (age_iso(ii + 1) / diff(age_iso([ii (ii + 1)]))) .* squeeze(speed_bal_surf(:, :, (ii + 1)) - ((age_iso(ii) / age_iso(ii + 1)) .* speed_bal_surf(:, :, ii)));
%     accum_hist(:, :, ii)    = (age_iso(ii + 1) / diff(age_iso([ii (ii + 1)]))) .* squeeze(accum_filt(:, :, (ii + 1)) - ((age_iso(ii) / age_iso(ii + 1)) .* accum_filt(:, :, ii)));
% end

% differences with modern speeds and accumulation rates
speed_change                = speed_bal_surf - speed_filt(:, :, ones(1, 1, num_iso));
speed_change_rel            = 100 .* (speed_change ./ speed_filt(:, :, ones(1, 1, num_iso)));
% speed_change_full           = speed_bal_full_surf - speed_filt;
% speed_change_rel_full       = 100 .* (speed_change_full ./ speed_filt);
% speed_hist_change           = speed_bal_surf_hist - speed_filt(:, :, ones(1, 1, (num_iso - 1)));
% speed_hist_change_rel       = 100 .* (speed_hist_change ./ speed_filt(:, :, ones(1, 1, (num_iso - 1))));
% accum_change                = accum_hist - accum_RACMO2_decim(:, :, ones(1, 1, (num_iso - 1)));

accum_diff                  = squeeze(accum_filt(:, :, end)) - accum_RACMO2_decim;
accum_diff_rel              = accum_diff ./ accum_RACMO2_decim;


%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
    
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
    
%% ISOCHRONE DEPTH, ACCUMULATION RATE, STRAIN RATE AND VERTICAL VELOCITY

    figure('position', [200 200 1600 600], 'name', [num2str(1e-3 * age_max) ' ka'])
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)])
    ranges                  = [200 2200; 0 50; 0 30; 0 15];
    incs                    = [100 2.5 1.5 0.75];
    plots                   = {'squeeze(depth_filt(:, :, end))' '1e2 .* squeeze(accum_filt(:, :, end))' '1e5 .* squeeze(strain_rate_filt(:, :, end))' '1e2 .* squeeze(vel_vert_filt(:, :, end))'};
    titles                  = {['$z_{' num2str(1e-3 * age_max) '\,\mathrm{ka}}$'] ['$\dot{b}_{' num2str(1e-3 * age_max) '\,\mathrm{ka}}$'] ['$\dot{\epsilon}_{' num2str(1e-3 * age_max) '\,\mathrm{ka}}$'] ['$w_{' num2str(1e-3 * age_max) '\,\mathrm{ka}}$']};
    spaces                  = [8 8 8 9];
    units                   = {'m' 'cm a^{-1}' '10^{-5} a^{-1}' 'cm a^{-1}'};
    letters                 = 'A':'D';
    [ax, cb1, cb2]          = deal(zeros(1, 4));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    for ii = 1:4
        ax(ii)              = subplot('position', [(0.005 + (0.25 * (ii - 1))) 0.005 0.245 0.95]);
        hold on
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        plot_tmp            = eval(plots{ii});
        plot_tmp            = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, plot_tmp, 'nearest', 'extrap');
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        [pp, pm, pb]        = deal(zeros(1, length(x_paral)), zeros(1, length(x_merid)), zeros(1, length(num_basin)));
        for jj = 1:length(x_paral)
            pp(jj)          = plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:length(x_merid)
            pm(jj)          = plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:num_basin
            pb(jj)          = plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2);
        end
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 2)
        caxis([1 23])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [])
        if (ii == 1)
            text(-530, -1045, 'NS', 'color', 'k', 'fontsize', 16, 'rotation', 20, 'fontweight', 'bold')
            text(-355, -1345, 'CC', 'color', 'k', 'fontsize', 16, 'fontweight', 'bold')
            text(-25, -2630, 'DYE-3', 'color', 'k', 'fontsize', 16, 'fontweight', 'bold')
            for jj = 1:2
                plot(x_core(jj), y_core(jj), 'k^', 'markerfacecolor', 'w', 'markersize', 12)
            end
            text(900, -525, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 20)            
        else
            text(850, -525, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 20)
        end
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-550, -800, [letters(ii) repmat(' ', 1, spaces(ii))], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        if (ii == 2)
            text(-450, -800, titles{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'interpreter', 'latex')
        else
            text(-450, -810, titles{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'interpreter', 'latex')
        end
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        tick_str            = cell(1, 21);
        for jj = 2:2:21
            tick_str{jj}    = '';
        end
        for jj = 1:2:21
            tick_str{jj}    = num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
        end
        if (ii == 1)
            tick_str{1}     = ['<' tick_str{1}];
        end
        tick_str{end}       = ['>' tick_str{end}];
        cb1(ii)             = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', ranges(ii, :), 'ytick', ranges(ii, 1):incs(ii):ranges(ii, 2), 'yticklabel', tick_str);
        set(cb1(ii), 'visible', 'off')
        cb2(ii)             = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2(ii), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str, 'position', get(cb1(ii), 'position'));
    end
    linkaxes(ax)
    
%% BALANCE-VELOCITY COMPARISON (FIG 2)

    figure('position', [200 200 1600 600], 'name', [num2str(1e-3 * age_max) ' ka'])
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)])
    ranges                  = [0 2; 0 2; -20 20; -100 100];
    incs                    = [0.1  0.1  2 10];
    plots                   = {'log10(speed_bal_surf(:, :, end))' 'log10(speed_filt)' 'speed_change(:, :, end)' 'speed_change_rel(:, :, end)'};
    titles                  = {['$u_s^{' num2str(1e-3 * age_max) '\,\mathrm{ka}}$'] '$u_s^{2007-2010}$' '$\Delta u_s$' '${\Delta}u_s/u_s^{2007-2010}$'};
    units                   = {'m a^{-1}' 'm a^{-1}' 'm a^{-1}' '%'};
    unitsx                  = [0 20 20 50];
    letters                 = 'A':'D';
    spaces                  = [8 17 8 25];
    [ax, cb, cb2]           = deal(zeros(1, 4));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([x_min x_max], [y_min y_max], 5, 10);
    for ii = 1:4
        ax(ii)              = subplot('position', [(0.005 + (0.25 * (ii - 1))) 0.005 0.245 0.95]);
        hold on
        plot_tmp            = eval(plots{ii});
        if (ii ~= 2)
            plot_tmp(~mask_D1_decim) ...
                            = NaN;
        end
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        plot_tmp            = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, plot_tmp, 'nearest', 'extrap');
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 2)
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        caxis([1 23])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
%         if (ii > 1)
            text((825 + unitsx(ii)), -525, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 20)
%         end
        text(-550, -800, [letters(ii) repmat(' ', 1, spaces(ii))], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        text(-450, -810, titles{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'interpreter', 'latex')
        tick_str        = cell(1, 21);
        if (ii > 2)
            for jj = 1:2:21
                tick_str{jj}= num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
            end
            for jj = 2:2:21
                tick_str{jj}= '';
            end
        else
            for jj = 1:21
                tick_str{jj}= '';
            end
            tick_str{1} = num2str(10 ^ (ranges(ii, 1)));
            tick_str{11}= '10';
            tick_str{end} ...
                = num2str(10 ^ (ranges(ii, 2)));
        end
        tick_str{1}     = ['<' tick_str{1}];
        if (ii >= 2)
            tick_str{end} ...
                = ['>' tick_str{end}];
        end
        cb(ii)          = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', ranges(ii, :), 'ytick', ranges(ii, 1):incs(ii):ranges(ii, 2), 'yticklabel', tick_str);
        set(cb(ii), 'visible', 'off')
        cb2(ii)         = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2(ii), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str, 'position', get(cb(ii), 'position'));
%         if (ii == 1)
%             set(cb2(ii), 'visible', 'off')
%         end
        if (any(ii == [1 2]))
            freezeColors
            if (ii == 2)
                colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; redblue(20)])
                set(cb2(ii), 'position', get(cb(ii), 'position'), 'ylim', [4 24], 'ytick', 4:24)
            end
        end
    end
    linkaxes(ax)
    
%% SHEAR-LAYER THICKNESS AND SHAPE FACTOR (FIG. S3)

    figure('position', [200 200 880 640])
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)])
    ranges                  = [0 2000; 0.9 1];
    incs                    = [100     0.005];
    plots                   = {'squeeze(thick_shear_filt(:, :, end))' 'squeeze(shape_filt(:, :, end))'};
    titles                  = {'$h$' '$\phi(z_{9\,\mathrm{ka}})$'};
    spaces                  = [4 13];
    units                   = {'m' ''};
    letters                 = 'A':'B';
    [ax, cb1, cb2]          = deal(zeros(1, 2));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    for ii = 1:2
        ax(ii)              = subplot('position', [(0.005 + (0.49 * (ii - 1))) 0.03 0.48 0.90]);
        hold on
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        plot_tmp            = eval(plots{ii});
        plot_tmp(~mask_D1_decim) ...
                            = NaN;
        plot_tmp            = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, plot_tmp, 'nearest', 'extrap');
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        [pp, pm, pb]        = deal(zeros(1, length(x_paral)), zeros(1, length(x_merid)), zeros(1, length(num_basin)));
        for jj = 1:length(x_paral)
            pp(jj)          = plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:length(x_merid)
            pm(jj)          = plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:num_basin
            pb(jj)          = plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2);
        end
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 1)
        caxis([1 23])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xticklabel', {}, 'yticklabel', {})
        if ~isempty(units{ii})
            text(900, -550, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 20)
        end
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-550, -800, [letters(ii) repmat(' ', 1, spaces(ii))], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        text(-450, -810, titles{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'interpreter', 'latex')
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        tick_str            = cell(1, 21);
        for jj = 2:2:21
            tick_str{jj}= '';
        end
        for jj = 1:2:21
            tick_str{jj}= num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
        end
        tick_str{1}         = ['<' tick_str{1}];
        tick_str{end}       = ['>' tick_str{end}];
        cb1(ii)              = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', ranges(ii, :), 'ytick', ranges(ii, 1):incs(ii):ranges(ii, 2), 'yticklabel', tick_str);
        set(cb1(ii), 'visible', 'off')
        cb2(ii)             = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2(ii), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str, 'position', get(cb1(ii), 'position'));
    end
    linkaxes(ax)
   
%% BALANCE-VELOCITY AND CHANGE ONLY

    figure('position', [200 200 880 640])
    colormap(jet(20))
    ranges                  = [0 2; -20 20];
    incs                    = [0.1  2];
    plots                   = {'log10(speed_bal_surf)' 'speed_change'};
    units                   = {'m a^{-1}' 'm a^{-1}'};
    unitsx                  = [0 20];
    letters                 = 'a':'b';
    spaces                  = [11 8];
    [ax, cb]                = deal(zeros(1, 2));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([x_min x_max], [y_min y_max], 5, 10);
    for ii = 1:2
        ax(ii)              = subplot('position', [(0.005 + (0.49 * (ii - 1))) 0.03 0.48 0.90]);
        hold on
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        if (ii > 1)
            range_tmp([1 end]) ...
                            = [-Inf Inf];
        else
            range_tmp([1 end]) ...
                            = [-1 3];
        end
        plot_tmp            = eval(plots{ii});
        plot_tmp(~mask_D1_filt) ...
                            = NaN;
        c                   = pcolor(x_grd_filt, y_grd_filt, plot_tmp);
        shading interp
        uistack(c, 'bottom')
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 1)
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        caxis(ranges(ii, :))
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xticklabel', {}, 'yticklabel', {})
        text((825 + unitsx(ii)), -550, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 18)
        text(-550, -800, ['(' letters(ii) ')'], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')
        tick_str        = cell(1, 21);
        if (ii > 1)
            for jj = 1:2:21
                tick_str{jj} ...
                        = num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
            end
            for jj = 2:2:21
                tick_str{jj} ...
                        = '';
            end
        else
            for jj = 1:21
                tick_str{jj}= '';
            end
            tick_str{1} = num2str(10 ^ (ranges(ii, 1)));
            tick_str{11}= '10';
            tick_str{end} ...
                        = num2str(10 ^ (ranges(ii, 2)));
        end
        tick_str{1}     = ['<' tick_str{1}];
        tick_str{end}   = ['>' tick_str{end}];
        cb(ii)          = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', ranges(ii, :), 'ytick', ranges(ii, 1):incs(ii):ranges(ii, 2), 'yticklabel', tick_str);
        if (ii == 1)
            freezeColors
            cb_pos      = get(cb(ii), 'position');
            cb2         = cbfreeze(cb(ii));
            colormap(redblue(20))
            set(cb2, 'position', cb_pos, 'ylim', [0 2], 'ytick', ranges(ii, 1):incs(ii):ranges(ii, 2))
        end
    end
    linkaxes(ax)
    
%% BALANCE-VELOCITY COMPARISON (FULL)

    figure('position', [200 200 1600 800], 'name', 'full-thickness balance velocity')
    colormap(jet(20))
    ranges                  = [0 2; 0 2; -20 20; -75 75];
    incs                    = [0.1  0.1  2       7.5];
    plots                   = {'log10(speed_bal_full)' 'log10(speed_filt)' 'speed_change_full' 'speed_change_rel_full'};
    titles                  = {'$u_s^H$' '$u_s^{2007-2010}$' '$\Delta u_s^H$' '${\Delta}u_s^H/u_s^H$'};
    units                   = {'m a^{-1}' 'm a^{-1}' 'm a^{-1}' '%'};
    unitsx                  = [0 20 20 70];
    letters                 = 'A':'D';
    spaces                  = [11 17 8 19];
    [ax, cb]                = deal(zeros(1, 4));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([x_min x_max], [y_min y_max], 5, 10);
    for ii = 1:4
        ax(ii)              = subplot('position', [(0.01 + (0.255 * (ii - 1))) 0.04 0.20 0.85]);
        hold on
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        if (ii > 2)
            range_tmp([1 end]) ...
                            = [-Inf Inf];
        else
            range_tmp([1 end]) ...
                            = [-1 3];
        end
        c                   = pcolor(x_decim, y_decim, eval(plots{ii}));
        shading interp
        uistack(c, 'bottom')
        if (ii ~= 2)
            plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 1)
        end
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        caxis(ranges(ii, :))
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xticklabel', {}, 'yticklabel', {})
        if (ii > 1)
            text((825 + unitsx(ii)), -550, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 18)
        end
        text(-550, -800, [letters(ii) repmat(' ', 1, spaces(ii))], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        text(-450, -810, titles{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'interpreter', 'latex')
        if (ii > 1)
            tick_str        = cell(1, 21);
            if (ii > 2)
                for jj = 1:2:21
                    tick_str{jj}= num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
                end
                for jj = 2:2:21
                    tick_str{jj}= '';
                end
            else
                for jj = 1:21
                    tick_str{jj}= '';
                end
                tick_str{1} = num2str(10 ^ (ranges(ii, 1)));
                tick_str{11}= '10';
                tick_str{end} ...
                            = num2str(10 ^ (ranges(ii, 2)));
            end
            tick_str{1}     = ['<' tick_str{1}];
            tick_str{end}   = ['>' tick_str{end}];
            cb(ii)          = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', ranges(ii, :), 'ytick', ranges(ii, 1):incs(ii):ranges(ii, 2), 'yticklabel', tick_str);
            set(gca, 'position', [(0.01 + (0.255 * (ii - 1))) 0.04 0.20 0.85])
            set(gca, 'position', (get(gca, 'position') + [-0.04 0 0 0]))
        end
        if (any(ii == [1 2]))
            freezeColors
            if (ii == 2)
                cb_pos      = get(cb(ii), 'position');
                cb2         = cbfreeze(cb(ii));
                colormap(redblue(20))
                set(cb2, 'position', cb_pos, 'ylim', [0 2], 'ytick', ranges(ii, 1):incs(ii):ranges(ii, 2))
            end
        end
    end
    linkaxes(ax)
        
%% SINGLE GRID PLOTS

    ranges                  = [0 40; 0 2000; -0.1 0.1];
    incs                    = [2 100 0.01; ];
    plots                   = {'1e5 .* strain_rate_filt' 'thick_shear_decim' 'melt_bed_decim' 'log10(speed_change ./ 9e3)'};
    units                   = {'10^{-5} a^{-1}' 'm' 'm a^{-1}'};
    [cb1, cb2]              = deal(zeros(1, 3));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    for ii = 4
        figure('position', [200 200 600 800])
        if (ii < 3)
            colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)])
        else
            colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; redblue(20)])
        end
        subplot('position', [0.005 0.03 0.9 0.90]);
        hold on
        plot_tmp            = eval(plots{ii});
        if (ii ~= 2)
            plot_tmp(~mask_D1_decim) ...
                            = NaN;
        end
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        plot_tmp            = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, plot_tmp, 'nearest', 'extrap');
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
        imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp)
        
        [pp, pm, pb]        = deal(zeros(1, length(x_paral)), zeros(1, length(x_merid)), zeros(1, length(num_basin)));
        for jj = 1:length(x_paral)
            pp(jj)          = plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:length(x_merid)
            pm(jj)          = plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:num_basin
            pb(jj)          = plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2);
        end
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 1)
        caxis([1 24])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xticklabel', {}, 'yticklabel', {})
        text(850, -550, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 20)
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)        
        tick_str            = cell(1, 21);
        for jj = 2:2:21
            tick_str{jj}    = '';
        end
        for jj = 1:2:21
            tick_str{jj}    = num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
        end
        if (ii > 1)
            tick_str{1}     = ['<' tick_str{1}];
        end
        tick_str{end}       = ['>' tick_str{end}];
        cb1(ii)              = colorbar('fontsize', 20, 'location', 'eastoutside');
        set(cb1(ii), 'visible', 'off')
        cb2(ii)             = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2(ii), 'position', get(cb1(ii), 'position'), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str)
    end
    
%% ACCUMULATION-RATE, BALANCE-VELOCITY AND NORTHERN HEMISPHERE TEMPERATURE HISTORY

    fgh                     = figure('position', [10 10 1280 720], 'color', 'w');
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20); redblue(20)])
    ranges_ref              = [0 60; 0 2];
    ranges                  = [-10 10; -20 20];
    incs_ref                = [3 0.1];
    incs                    = [1 2];
    plots_ref               = {'1e2 .* accum_RACMO2_decim' 'log10(speed_filt)'};
    plots                   = {'1e2 .* accum_change' 'speed_hist_change'};
    units                   = {'cm a^{-1}' 'm a^{-1}'};
    titles                  = {'A: Accumulation rate' 'B: Surface speed' {'C: Northern hemisphere'; '     temperature anomaly'}};
    [ax, p_ref]             = deal(zeros(1, 3));
    [cb, cb2, cb2_ref, pcr] = deal(zeros(1, 2));
    [tick_str, tick_str_ref]= deal(cell(1, 2));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([x_min x_max], [y_min y_max], 5, 10);
    plot_ind                = cell(2, (num_iso - 1));
    for ii = 1:2
        ax(ii)              = subplot('position', [(0.02 + (0.33 * (ii - 1))) 0.04 0.30 0.84]);
        hold on
        for jj = (num_iso - 1):-1:1
            plot_ind{ii, jj}= squeeze(eval([plots{ii} '(:, :, jj)']));
            plot_ind{ii, jj}(~mask_D1_decim) ...
                            = NaN;
            range_tmp       = ranges(ii, 1):incs(ii):ranges(ii, 2);
            plot_ind{ii, jj}= 23 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, plot_ind{ii, jj}, 'nearest', 'extrap');
            plot_ind{ii, jj}(isnan(plot_ind{ii, jj})) ...
                            = mask_combo_decim(isnan(plot_ind{ii, jj})) + 1;
        end
        plot_tmp_ref        = eval(plots_ref{ii});
        plot_tmp_ref(~mask_gris_decim) ...
                            = NaN;
        range_tmp           = ranges_ref(ii, 1):incs_ref(ii):ranges_ref(ii, 2);
        plot_tmp_ref        = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, plot_tmp_ref, 'nearest', 'extrap');
        plot_tmp_ref(isnan(plot_tmp_ref)) ...
                            = mask_combo_decim(isnan(plot_tmp_ref)) + 1;
        pcr(ii)             = imagesc(x_decim(1, :), y_decim(:, 1), plot_tmp_ref);
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 2)
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        caxis([1 43])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xticklabel', {}, 'yticklabel', {}, 'xtick', [], 'ytick', [])
        title(titles{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'normal')
        text(850, -525, ['(' units{ii} ')'], 'fontsize', 20, 'color', 'k')
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        [tick_str{ii}, tick_str_ref{ii}] ...
                            = deal(cell(1, 21));
        for jj = 1:21
            [tick_str{ii}{jj}, tick_str_ref{ii}{jj}] ...
                            = deal('');
        end
        for jj = 1:2:21
            tick_str{ii}{jj}= num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
            if (ii == 1)
                tick_str_ref{ii}{jj} ...
                            = num2str(ranges_ref(ii, 1) + (incs_ref(ii) * (jj - 1)));
            end
        end
        if (ii == 2)
            tick_str_ref{ii}{1} ...
                            = ['<' num2str(10 ^ (ranges_ref(ii, 1)))];
            tick_str_ref{ii}{11} ...
                            = '10';
            tick_str_ref{ii}{end} ...
                            = num2str(10 ^ (ranges_ref(ii, 2)));
        end
        tick_str{ii}{1}     = ['<' tick_str{ii}{1}];
        tick_str{ii}{end}   = ['>' tick_str{ii}{end}];
        tick_str_ref{ii}{end} ...
                            = ['>' tick_str_ref{ii}{end}];
        cb(ii)              = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', ranges_ref(ii, :), 'ytick', ranges_ref(ii, 1):incs(ii):ranges_ref(ii, 2), 'yticklabel', tick_str{ii});
        set(cb(ii), 'visible', 'off')
        cb2_ref(ii)         = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2_ref(ii), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str_ref{ii}, 'position', get(cb(ii), 'position'));
        set(ax(ii), 'userdata', [])
    end
    p_title                 = annotation('textbox', [0.21 0.93 0.4 0.05], 'string', 'Modern values', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none', 'horizontalalignment', 'center');
    p_age                   = annotation('textbox', [0.775 0.14 0.15 0.04], 'string', '', 'fontsize', 20, 'fontweight', 'bold', 'backgroundcolor', 'w', 'horizontalalignment', 'center', 'verticalalignment', 'baseline');
    p_ref                   = annotation('textbox', [0.70 0.06 0.3 0.05], 'string', {'A: Ettema et al. [2009; Geophys. Res. Lett.]'; 'B: Joughin et al. [2010; J. Glaciol.]'; 'C: Marcott et al. [2013; Science]'}, 'fontsize', 18, 'linestyle', 'none');
    ax(3)                   = subplot('position', [0.71 0.28 0.27 0.6]);
    hold on
    ind_age_curr            = find(~isnan(temp_nh));
    pmf                     = fill((1e-3 .* [age_marcott(ind_age_curr); flipud(age_marcott(ind_age_curr))]), [(temp_nh(ind_age_curr) - temp_nh_uncert(ind_age_curr)); flipud(temp_nh(ind_age_curr) + temp_nh_uncert(ind_age_curr))], [1 0.75 1], 'edgecolor', 'none');
    plot((1e-3 .* age_marcott(ind_age_curr)), temp_nh(ind_age_curr), 'color', [0.5 0.5 0.5], 'linewidth', 3)
    axis([(1e-3 .* age_marcott([1 end])') -2 1.5])
    set(gca, 'fontsize', 20, 'xdir', 'reverse')
    xlabel('Age (ka B.P.)')
    ylabel('(\circC)')
    title(titles{3}, 'color', 'k', 'fontsize', 20, 'fontweight', 'normal')
    box on
    grid on
    
%% ACCUMULATION/FLOW HISTORY ANIMATION

    figure(fgh)
    rect                    = get(fgh, 'position');
    mp4                     = VideoWriter('fig/gris_hist.mp4', 'MPEG-4');
    mp4.FrameRate           = 1;
    mp4.Quality             = 100;
    open(mp4)
    writeVideo(mp4, getframe(fgh, [0 0 rect(3:4)]));
    writeVideo(mp4, getframe(fgh, [0 0 rect(3:4)]));
    writeVideo(mp4, getframe(fgh, [0 0 rect(3:4)]));
    writeVideo(mp4, getframe(fgh, [0 0 rect(3:4)]));
    pause(1)
    pc                      = zeros(1, 2);
    for ii = (num_iso - 1):-1:1
        if (ii < (num_iso - 1))
            delete([pc p_fill_curr p_age_curr])
        else
            set(p_title, 'string', 'Differences from modern values')
            delete([pcr p_ref cb2_ref])
            for jj = 1:2
                ax_pos      = get(ax(jj), 'position');
                axes(ax(jj))
                cb2(jj)     = cbarf([24 44], 24:44, 'vertical', 'linear');
                set(ax(jj), 'position', ax_pos)
                set(cb2(jj), 'fontsize', 20, 'ytick', 24:44, 'yticklabel', tick_str{jj}, 'position', get(cb(jj), 'position'))
            end
        end
        for jj = 1:2
            axes(ax(jj))
            pc(jj)          = imagesc(x_decim(1, :), y_decim(:, 1), plot_ind{jj, ii});
            uistack(pc(jj), 'bottom')
        end
        set(p_age, 'string', [sprintf('%1.2f', (1e-3 * age_iso(ii + 1))) ' - ' sprintf('%1.2f', (1e-3 * age_iso(ii))) ' ka B.P.'])
        axes(ax(3))
        ind_age_curr        = find((age_marcott <= age_iso(ii + 1)) & (age_marcott >= age_iso(ii)) & ~isnan(temp_nh));
        p_fill_curr         = fill((1e-3 .* [age_marcott(ind_age_curr); flipud(age_marcott(ind_age_curr))]), [(temp_nh(ind_age_curr) - temp_nh_uncert(ind_age_curr)); flipud(temp_nh(ind_age_curr) + temp_nh_uncert(ind_age_curr))], [1 0.25 1], 'edgecolor', 'none');
        p_age_curr          = plot((1e-3 .* age_marcott(ind_age_curr)), temp_nh(ind_age_curr), 'k', 'linewidth', 3);
        pause(1)
        writeVideo(mp4, getframe(fgh, [0 0 rect(3:4)]));
        writeVideo(mp4, getframe(fgh, [0 0 rect(3:4)]));
    end
    close(mp4)

%% ACCUMULATION-RATE COMPARISON (FIG. S4)

    figure('position', [200 200 1600 600], 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)])
    ranges                  = [0 60; 0 60; -10 10; -50 50];
    incs                    = [3 3 1 5];
    plots                   = {'1e2 .* squeeze(accum_filt(:, :, end))' '1e2 .* accum_RACMO2_decim' '1e2 .* accum_diff' '1e2 .* accum_diff_rel'};
    titles                  = {'$\dot{b}_{9\,\mathrm{ka}}$' '$\dot{b}_{1958-2007}$' '$\Delta \dot{b}$' '$\Delta \dot{b}/\dot{b}_{1958-2007}$'};
    units                   = {'(cm a^{-1})' '(cm a^{-1})' '(cm a^{-1})' '(%)'};
    spaces                  = [8 17 7 23];
    unitsx                  = [0 0 25 60];
    letters                 = 'A':'D';
    [ax, cb, cb2]           = deal(zeros(1, 4));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([x_min x_max], [y_min y_max], 5, 10);
    for ii = 1:4
        ax(ii)              = subplot('position', [(0.01 + (0.25 * (ii - 1))) 0.005 0.245 0.95]);
        hold on
        plot_tmp            = eval(plots{ii});
        if (ii ~= 2)
            plot_tmp(~mask_D1_decim) ...
                            = NaN;
        end
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        plot_tmp            = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, plot_tmp, 'nearest', 'extrap');
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo_decim(isnan(plot_tmp)) + 1;
        imagesc(x_grd(1, :), y_grd(:, 1), plot_tmp)
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), '--', 'linewidth', 0.5, 'color', [0.5 0.5 0.5])
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 2)
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        caxis([1 23])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [], 'xticklabel', {}, 'yticklabel', {})
        text((825 + unitsx(ii)), -525, units{ii}, 'color', 'k', 'fontsize', 20)
        text(-607, -3215, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(483, -2771, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(717, -2897, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        text(-550, -800, [letters(ii) repmat(' ', 1, spaces(ii))], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        text(-450, -810, titles{ii}, 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'interpreter', 'latex')
        tick_str            = cell(1, 21);
        for jj = 1:2:21
            tick_str{jj}    = num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
        end
        for jj = 2:2:21
            tick_str{jj}    = '';
        end
        if (ii > 2)
            tick_str{1}     = ['<' tick_str{1}];
        end
        tick_str{end}       = ['>' tick_str{end}];
        cb(ii)              = colorbar('fontsize', 20, 'location', 'eastoutside');
        set(cb(ii), 'visible', 'off')
        cb2(ii)             = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2(ii), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str, 'position', get(cb(ii), 'position'));
        if (any(ii == [1 2]))
            freezeColors
            if (ii == 2)
                colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; redblue(20)])
            end
        end
    end
    linkaxes(ax)
    
%%
end