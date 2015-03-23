function atten_layers(hfcp_model, do_set, year_good, dist_decim, depth_range, thick_range, num_layer_min, az_diff_max, exp_geo, slope_max, slope_loss_max, envelope_bound, frac_test, conf_uncert)
% ATTEN_LAYERS Radar attenuation-rate and inferred temperature from dated radar reflections.
% 
% ATTEN_LAYERS(HFCP_MODEL,DO_SET,YEAR_GOOD,DIST_DECIM,DEPTH_RANGE,THICK_RANGE,NUM_LAYER_MIN,AZ_DIFF_MAX,EXP_GEO,SLOPE_MAX,SLOPE_LOSS_MAX,ENVELOPE_BOUND,FRAC_TEST,CONF_UNCERT)
% calculates depth-averaged radar-attenuation rates and equivalent
% temperatures using traced radiostratigraphy, with inputs as follows:
% 
% HFCP_MODEL: Ice-conductivity model (see HFCONDPROP for options).
% DO_SET: 1x5 logical vector determining whether to only analyze core
% intersections, to ignore transect portions with large azimuth changes, to
% correct echo intensities for slope, to use a percentile envelope instead
% of the arithmetic mean and to use spatially varying chemistry.
% YEAR_GOOD: Vector of indices of campaigns to analyze.
% DIST_DECIM: Distance over which to decimate results (km).
% DEPTH_RANGE: Minimum/maximum depth between which to include reflections
% in analysis (m).
% THICK_RANGE: Minimum/maximum thickness of layer depth range relative to
% measured ice thickness.
% NUM_LAYER_MIN: Minimum number of layers for analysis.
% AZ_DIFF_MAX: Maximum change in along-track azimuth between fixed
% intervals (deg).
% EXP_GEO: Geometric spreading exponent.
% SLOPE_MAX: Maximum Doppler slope (deg).
% SLOPE_LOSS_MAX: Maximum correction for slope-induced power loss (dB)
% ENVELOPE_BOUND: Envelope bounds (%).
% FRAC_TEST: Range around best-fit model parameter about which to test
% (fractions of model parameter).
% CONF_UNCERT: Confidence bound for uncertainty calculations, e.g., 0.95 is
% 95%, 0.68 ~ 1 sigma.
% 
% Joe MacGregor (UTIG)
% Last updated: 03/23/15

if (nargin ~= 14)
    error('atten_layers:nargin', 'Incorrect number of arguments (must be 14).')
end
if ~ischar(hfcp_model)
    error('atten_layers:hfcp_model', 'HFCP_MODEL is not a string.')
elseif ~any(strcmp(hfcp_model, {'M07-std' 'M07-adj' 'M07-adj-NH4' 'M07-wt' 'W97' 'W97-adj'}))
    error('atten_layers:hfcp_model_type', 'HFCP_MODEL must be either ''M07-std'', ''M07-adj'', ''M07-adj-NH4'', ''M07-wt'', ''W97'', ''W97-adj''.')
end
if (~islogical(do_set) || ~isvector(do_set) || (length(do_set) ~= 5))
    error('atten_layers:do_set', 'DO_SET must be a 5-element logical vector.')
end
if (~isnumeric(year_good) || ~isvector(year_good) || any(mod(year_good, 1)))
    error('atten_layers:year_good', 'YEAR_GOOD must be a vector of integers.')
end
if (~isnumeric(dist_decim) || ~isscalar(dist_decim) || (dist_decim <= 0))
    error('atten_layers:dist_decim', 'DIST_DECIM must be a numeric positive scalar.')
end
if (~isnumeric(depth_range) || ~isvector(depth_range) || (length(depth_range) ~= 2) || any(depth_range < 0))
    error('atten_layers:depth_range', 'DEPTH_RANGE must be a 2-element numeric positive vector.')
end
if (~isnumeric(thick_range) || ~isvector(thick_range) || (length(thick_range) ~= 2) || any(thick_range < 0))
    error('atten_layers:thick_range', 'THICK_RANGE must be a 2-element numeric positive vector.')
end
if (~isnumeric(num_layer_min) || ~isscalar(num_layer_min) || mod(num_layer_min, 1))
    error('atten_layers:num_layer_min', 'NUM_LAYER_MIN must be a scalar integer.')
end
if (~isnumeric(az_diff_max) || ~isscalar(az_diff_max) || (az_diff_max <= 0))
    error('atten_layers:az_diff_max', 'AZ_DIFF_MAX must be a numeric positive scalar.')
end
if (~isnumeric(exp_geo) || ~isscalar(exp_geo) || (exp_geo <= 0))
    error('atten_layers:exp_geo', 'EXP_GEO must be a numeric positive scalar.')
end
if (~isnumeric(slope_max) || ~isscalar(slope_max) || (slope_max <= 0))
    error('atten_layers:slope_max', 'SLOPE_MAX must be a numeric positive scalar.')
end
if (~isnumeric(slope_loss_max) || ~isscalar(slope_loss_max) || (slope_loss_max <= 0))
    error('atten_layers:slope_loss_max', 'SLOPE_LOSS_MAX must be a numeric positive scalar.')
end
if (~isnumeric(envelope_bound) || ~isvector(envelope_bound) || (length(envelope_bound) ~= 2) || any(envelope_bound < 0) || any(envelope_bound > 100))
    error('atten_layers:envelope_bound', 'ENVELOPE_BOUND must be a 2-element numeric positive vector with values between 0 and 100.')
end
if (~isnumeric(frac_test) || ~isvector(frac_test) || any(frac_test < 0))
    error('atten_layers:frac_test', 'FRAC_TEST must be a numeric positive vector.')
end
if (~isnumeric(conf_uncert) || ~isscalar(conf_uncert) || (conf_uncert <= 0) || (conf_uncert >= 1))
    error('atten_layers:conf_uncert', 'CONF_UNCERT must be a numeric positive scalar between 0 and 1.')
end
if nargout
    error('atten_layers:nargout', 'ATTEN_LAYERS has no outputs.')
end

% switches and reassignments
[do_core, do_az, do_slope, do_envelope, do_chem] ...
                            = deal(do_set(1), do_set(2), do_set(3), do_set(4), do_set(5));
[depth_min, depth_max]      = deal(depth_range(1), depth_range(2));
[thick_rel_min, thick_rel_max] ...
                            = deal(thick_range(1), thick_range(2));
hfcp_struct                 = hfcondprop(hfcp_model);

load mat/xy_all name_trans name_year num_trans num_year
load mat/merge_all dist depth_smooth elev_air_gimp elev_smooth_gimp elev_surf_gimp int_smooth lat lon num_layer num_subtrans thick twtt_smooth twtt_surf x_pk y_pk
load mat/date_all age
load mat/core_int_temp int_core_merge
load mat/grip_chem age_chem Cl H NH4
load mat/post_pulse_comp_gain_corr gain_corr twtt_gain_corr

if do_slope
    load mat/doppler_win doppler_bw doppler_win
end

permitt_ice                 = 3.15; % CReSIS standard real part of complex relative permittivity
letters                     = 'a':'z';
twtt_combine_wf_2013        = [9e-6 1e-5]; % waveform combination traveltimes for 2013 P3
ind_trans_switch_2013       = 17; % transect number at which waveform combination traveltime switches to second value

load mat/greenland_cism temp_surf x y
[temp_surf_grd, x_cism, y_cism] ...
                            = deal((temp_surf + 273.15), x, y);
clear temp_surf x y

opt                         = optimset('fminsearch');
opt.Display                 = 'off'; % avoid display of having reach maximum number of iterations
opt.FunValCheck             = 'on'; % sometimes strain-rate models get complex/NaN due to, e.g., logs of a negative argument; this switch stops that

[H_mean, Cl_mean, NH4_mean] = deal(mean(H), mean(Cl), mean(NH4));

if license('checkout', 'distrib_computing_toolbox')
    pool_check              = gcp('nocreate');
    if isempty(pool_check)
        try
            pool            = parpool('local', 4);
        catch
            pool            = parpool('local');
        end
    end
    parallel_check          = true;
else
    parallel_check          = false;
end

% initialize numerous variables
atten_core                  = [int_core_merge(:, 1:4) NaN(size(int_core_merge, 1), 8)];

[atten_rate, atten_rate_uncert, az_trans, depth_bound_decim, depth_min_atten, depth_max_atten, ind_decim, ind_decim_mid, ind_layer_int, int_smooth_corr, int_smooth_corr_uncert, mse_atten, num_decim, res_iso, slope_layer, temp_iso, ...
 temp_iso_uncert, thick_atten, thick_decim, Cl_trans, H_trans, NH4_trans] ...
                            = deal(cell(1, num_year));
for ii = year_good
    [atten_rate{ii}, atten_rate_uncert{ii}, az_trans{ii}, depth_bound_decim{ii}, depth_min_atten{ii}, depth_max_atten{ii}, ind_decim{ii}, ind_decim_mid{ii}, ind_layer_int{ii}, int_smooth_corr{ii}, int_smooth_corr_uncert{ii}, mse_atten{ii}, num_decim{ii}, ...
     res_iso{ii}, slope_layer{ii}, temp_iso{ii}, temp_iso_uncert{ii}, thick_atten{ii}, thick_decim{ii}, Cl_trans{ii}, H_trans{ii}, NH4_trans{ii}] ...
                            = deal(cell(1, num_trans(ii)));
    for jj = 1:num_trans(ii)
        [atten_rate{ii}{jj}, atten_rate_uncert{ii}{jj}, az_trans{ii}{jj}, depth_bound_decim{ii}{jj}, depth_min_atten{ii}{jj}, depth_max_atten{ii}{jj}, ind_decim{ii}{jj}, ind_decim_mid{ii}{jj}, ind_layer_int{ii}{jj}, int_smooth_corr{ii}{jj}, ...
         int_smooth_corr_uncert{ii}{jj}, mse_atten{ii}{jj}, res_iso{ii}{jj}, slope_layer{ii}{jj}, temp_iso{ii}{jj}, temp_iso_uncert{ii}{jj}, thick_atten{ii}{jj}, thick_decim{ii}{jj}, Cl_trans{ii}{jj}, H_trans{ii}{jj}, NH4_trans{ii}{jj}] ...
                            = deal(cell(1, num_subtrans{ii}(jj)));
        num_decim{ii}{jj}   = NaN(1, num_subtrans{ii}(jj));
        for kk = 1:num_subtrans{ii}(jj)
            dist_int        = dist{ii}{jj}{kk}(1):dist_decim:dist{ii}{jj}{kk}(end);
            dist_int(end)   = dist{ii}{jj}{kk}(end);
            ind_decim{ii}{jj}{kk} ...
                            = unique(interp1(unique(dist{ii}{jj}{kk}), 1:length(unique(dist{ii}{jj}{kk})), dist_int, 'nearest', 'extrap'));
            ind_decim_mid{ii}{jj}{kk} ...
                            = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));
            num_decim{ii}{jj}(kk) ...
                            = length(ind_decim{ii}{jj}{kk}) - 1;
            [int_smooth_corr{ii}{jj}{kk}, int_smooth_corr_uncert{ii}{jj}{kk}, slope_layer{ii}{jj}{kk}] ...
                            = deal(NaN(num_layer{ii}{jj}(kk), num_decim{ii}{jj}(kk)));
            [atten_rate{ii}{jj}{kk}, atten_rate_uncert{ii}{jj}{kk}, az_trans{ii}{jj}{kk}, depth_min_atten{ii}{jj}{kk}, depth_max_atten{ii}{jj}{kk},  mse_atten{ii}{jj}{kk}, res_iso{ii}{jj}{kk}, ...
             temp_iso{ii}{jj}{kk}, temp_iso_uncert{ii}{jj}{kk}, thick_atten{ii}{jj}{kk}, thick_decim{ii}{jj}{kk}] ...
                            = deal(NaN(1, num_decim{ii}{jj}(kk)));
            [depth_bound_decim{ii}{jj}{kk}, ind_layer_int{ii}{jj}{kk}, Cl_trans{ii}{jj}{kk}, H_trans{ii}{jj}{kk}, NH4_trans{ii}{jj}{kk}] ...
                            = deal(cell(1, num_decim{ii}{jj}(kk)));
        end
    end
end

if parallel_check
    pctRunOnAll warning('off', 'MATLAB:lscov:RankDefDesignMat')
else
    warning('off', 'MATLAB:lscov:RankDefDesignMat')
end

for ii = year_good
    
    disp(['Campaign: ' name_year{ii} '...'])
    
    for jj = 1:num_trans(ii)
        
        if do_slope
            angle_doppler_curr ...
                            = linspace(-(doppler_bw{ii}(jj) / 2), (doppler_bw{ii}(jj) / 2), 1e3);
            switch doppler_win{ii}{jj}
                case 'blackman'
                    loss_doppler_curr ...
                            = 20 .* log10(blackman(1e3));
                case 'hanning'
                    loss_doppler_curr ...
                            = 20 .* log10(hanning(1e3));
                case 'tukeywin'
                    loss_doppler_curr ...
                            = 20 .* log10(tukeywin(1e3, 0.2));
            end
        end
        
        for kk = 1:num_subtrans{ii}(jj)
            
            if (do_core && isempty(find(ismember(int_core_merge(:, 1:3), [ii jj kk], 'rows'), 1))) % require a core intersection
                continue
            end
            if (num_layer{ii}{jj}(kk) < num_layer_min) % require at least four layers
                continue
            else
                disp([name_trans{ii}{jj} letters(kk) '...'])
            end
            
            % along-track azimuth
            if license('checkout', 'map_toolbox')
                az_trans{ii}{jj}{kk} ...
                            = azimuth(lat{ii}{jj}{kk}(ind_decim{ii}{jj}{kk}(1:(end - 1))), lon{ii}{jj}{kk}(ind_decim{ii}{jj}{kk}(1:(end - 1))), lat{ii}{jj}{kk}(ind_decim{ii}{jj}{kk}(2:end)), lon{ii}{jj}{kk}(ind_decim{ii}{jj}{kk}(2:end)), wgs84Ellipsoid);
            else
                az_trans{ii}{jj}{kk} ...
                            = atan2d(diff(y_pk{ii}{jj}{kk}(ind_decim{ii}{jj}{kk})), diff(x_pk{ii}{jj}{kk}(ind_decim{ii}{jj}{kk})));
            end
            
            % correct 2013 echo intensities for waveform combination error
            if strcmp(name_year{ii}, '2013_p3')
                
                % initial gain correction for observed traveltimes
                gain_corr_curr ...
                            = cell(1, 2);
                for mm = 1:2
                    gain_corr_curr{mm} ...
                            = interp1(twtt_gain_corr{mm}, gain_corr{mm}, twtt_smooth{ii}{jj}{kk});
                end
                
                % waveform combination traveltime, based on transect index (changes towards end of 2013)
                twtt_combine_wf_curr ...
                            = twtt_combine_wf_2013(1);
                if (jj >= ind_trans_switch_2013)
                    twtt_combine_wf_curr ...
                            = twtt_combine_wf_2013(2);
                end
                
                % combine two waveforms' gain corrections as appropriate
                gain_corr_combine ...
                            = gain_corr_curr{1}; % start out with first waveform's correction
                gain_corr_combine(twtt_smooth{ii}{jj}{kk} < twtt_gain_corr{1}(1)) ...
                            = gain_corr_curr{1}(1);
                gain_corr_combine(twtt_smooth{ii}{jj}{kk} >= (repmat(twtt_surf{ii}{jj}{kk}, num_layer{ii}{jj}(kk), 1) + twtt_combine_wf_curr)) ...
                            = gain_corr_curr{2}(twtt_smooth{ii}{jj}{kk} >= (repmat(twtt_surf{ii}{jj}{kk}, num_layer{ii}{jj}(kk), 1) + twtt_combine_wf_curr));
                gain_corr_combine(twtt_smooth{ii}{jj}{kk} > twtt_gain_corr{2}(end)) ...
                            = gain_corr_curr{2}(end);
                
                % add gain correction to current sub-transect's echo intensities
                int_smooth{ii}{jj}{kk} ...
                            = int_smooth{ii}{jj}{kk} + gain_corr_combine;
            end
            
            % geometric spreading loss correction
            int_smooth_corr_full ...
                            = int_smooth{ii}{jj}{kk} + ((1e1 * exp_geo) .* log10(repmat((elev_air_gimp{ii}{jj}{kk} - elev_surf_gimp{ii}{jj}{kk}), num_layer{ii}{jj}(kk), 1) + (depth_smooth{ii}{jj}{kk} ./ sqrt(permitt_ice))));
            
            depth_smooth_decim ...
                            = NaN(num_layer{ii}{jj}(kk), num_decim{ii}{jj}(kk));
            
            if ~isempty(find(ismember(int_core_merge(:, 1:3), [ii jj kk], 'rows'), 1))
                ind_core_full_curr ...
                            = find(ismember(int_core_merge(:, 1:3), [ii jj kk], 'rows')); % row indices in core intersection list that match current transect
                if (do_core && isempty(ind_core_full_curr))
                    continue
                else
                    ind_core_decim_curr ...
                            = interp1(ind_decim_mid{ii}{jj}{kk}, 1:num_decim{ii}{jj}(kk), int_core_merge(ind_core_full_curr, 5), 'nearest')';
                end
            else
                [ind_core_full_curr, ind_core_decim_curr] ...
                            = deal([]);
            end
            
            % along-track modern surface temperature
            temp_surf_curr  = interp2(x_cism, y_cism, temp_surf_grd, x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}), y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}), 'linear', NaN);
            
            % initialize variables for extraction within each interval
            [depth_bound_decim{ii}{jj}{kk}, depth_inc_curr, loss_curr, loss_uncert_curr, H_curr, Cl_curr, NH4_curr] ...
                            = deal(cell(1, num_decim{ii}{jj}(kk)));
            ind_good        = false(1, num_decim{ii}{jj}(kk)); % logical test for each interval
            num_layer_curr  = zeros(1, num_decim{ii}{jj}(kk));
            
            if do_core
                ind_decim_curr ...
                            = unique([ind_core_decim_curr (ind_core_decim_curr - 1) (ind_core_decim_curr + 1)]);
                ind_decim_curr ...
                            = ind_decim_curr((ind_decim_curr > 0) & (ind_decim_curr <= num_decim{ii}{jj}(kk)));
            else
                ind_decim_curr ...
                            = 1:num_decim{ii}{jj}(kk);
            end
            
            for ll = ind_decim_curr
                
                ind_curr    = ind_decim{ii}{jj}{kk}(ll):ind_decim{ii}{jj}{kk}(ll + 1);
                
                thick_decim{ii}{jj}{kk}(ll) ...
                            = nanmean(thick{ii}{jj}{kk}(ind_curr));
                if isnan(thick_decim{ii}{jj}{kk}(ll))
                    continue
                end
                
                % smoothed layer depths
                depth_smooth_decim(:, ll) ...
                            = nanmean(depth_smooth{ii}{jj}{kk}(:, ind_curr), 2);
                
                % layer slope at fixed intervals
                for mm = 1:num_layer{ii}{jj}(kk)
                    if (length(find(~isnan(elev_smooth_gimp{ii}{jj}{kk}(mm, ind_curr)))) > 3)
                        slope_layer_tmp ...
                            = polyfit(dist{ii}{jj}{kk}(ind_curr(~isnan(elev_smooth_gimp{ii}{jj}{kk}(mm, ind_curr)))), elev_smooth_gimp{ii}{jj}{kk}(mm, ind_curr(~isnan(elev_smooth_gimp{ii}{jj}{kk}(mm, ind_curr)))), 1);
                        slope_layer{ii}{jj}{kk}(mm, ll) ...
                            = slope_layer_tmp(1);
                    end
                end
                
                if do_envelope % echo intensity at fixed intervals and standard deviation within percentile bounds
                    int_smooth_corr_full_curr ...
                            = 1e1 .^ (1e-1 .* int_smooth_corr_full(:, ind_curr));
                    int_smooth_prctile ...
                            = prctile(int_smooth_corr_full_curr, envelope_bound, 2); % envelope bounds
                    for mm = find(~isnan(int_smooth_prctile(:, 1)))'
                        int_smooth_set ...
                            = int_smooth_corr_full_curr(mm, (int_smooth_corr_full_curr(mm, :) >= int_smooth_prctile(mm, 1)) & (int_smooth_corr_full_curr(mm, :) <= int_smooth_prctile(mm, 2)));
                        int_smooth_corr{ii}{jj}{kk}(mm, ll) ...
                            = 1e1 .* log10(mean(int_smooth_set));
                        int_smooth_corr_uncert{ii}{jj}{kk}(mm, ll) ...
                            = nanstd((1e1 .* log10(int_smooth_set)), 0, 2);
                    end
                else % straight average
                    int_smooth_corr{ii}{jj}{kk}(:, ll) ...
                            = 1e1 .* log10(nanmean((1e1 .^ (1e-1 .* int_smooth_corr_full(:, ind_curr))), 2));
                    int_smooth_corr_uncert{ii}{jj}{kk}(:, ll) ...
                            = nanstd(int_smooth_corr_full(:, ind_curr), 0, 2);
                end
                
                % correct for power loss due to reflection slope
                if do_slope
                    % surface slope, m/km
                    slope_surf_curr ...
                            = polyfit(dist{ii}{jj}{kk}(ind_curr(~isnan(elev_surf_gimp{ii}{jj}{kk}(ind_curr)))), elev_surf_gimp{ii}{jj}{kk}(ind_curr(~isnan(elev_surf_gimp{ii}{jj}{kk}(ind_curr)))), 1);
                    slope_surf_curr ...
                            = slope_surf_curr(1);
                    % Doppler layer slopes, deg
                    slope_layer_doppler ...
                            = asind((sqrt(permitt_ice) * sind(atand(1e-3 .* slope_layer{ii}{jj}{kk}(:, ll))))) - atand(1e-3 * slope_surf_curr);%- atand(1e-3 * slope_air_curr);
                    for mm = find(~isnan(int_smooth_corr{ii}{jj}{kk}(:, ll)))'
                        if (abs(interp1(angle_doppler_curr, loss_doppler_curr, slope_layer_doppler(mm), 'linear', 'extrap')) <= slope_loss_max) % don't let slope correction get too big
                            int_smooth_corr{ii}{jj}{kk}(mm, ll) ...
                                = int_smooth_corr{ii}{jj}{kk}(mm, ll) - interp1(angle_doppler_curr, loss_doppler_curr, slope_layer_doppler(mm), 'linear', 'extrap');
                        else
                            int_smooth_corr{ii}{jj}{kk}(mm, ll) ...
                                = NaN;
                        end
                    end
                end
                
                % layers that meet thresholds
                if do_slope
                    ind_layer_curr ...
                            = find(~isnan(depth_smooth_decim(:, ll)) & ~isnan(int_smooth_corr{ii}{jj}{kk}(:, ll)) & ~isnan(int_smooth_corr_uncert{ii}{jj}{kk}(:, ll)) & (depth_smooth_decim(:, ll) >= depth_min) & ...
                                   (depth_smooth_decim(:, ll) <= min([depth_max (thick_rel_max * thick_decim{ii}{jj}{kk}(ll))])) & (abs(slope_layer_doppler) <= slope_max));
                else
                    ind_layer_curr ...
                            = find(~isnan(depth_smooth_decim(:, ll)) & ~isnan(int_smooth_corr{ii}{jj}{kk}(:, ll)) & ~isnan(int_smooth_corr_uncert{ii}{jj}{kk}(:, ll)) & (depth_smooth_decim(:, ll) >= depth_min) & ...
                                   (depth_smooth_decim(:, ll) <= min([depth_max (thick_rel_max * thick_decim{ii}{jj}{kk}(ll))])));
                end
                
                % require a minimum number of useful layers
                if (length(ind_layer_curr) < num_layer_min)
                    continue
                else
                    num_layer_curr(ll) ...
                            = length(ind_layer_curr) - 1;
                end
                
                % breakout current depths / echo intensities
                [depth_bound_decim{ii}{jj}{kk}{ll}, ind_layer_int{ii}{jj}{kk}{ll}] ...
                            = sort(depth_smooth_decim(ind_layer_curr, ll));
                ind_layer_int{ii}{jj}{kk}{ll} ...
                            = ind_layer_curr(ind_layer_int{ii}{jj}{kk}{ll});
                [depth_min_atten{ii}{jj}{kk}(ll), depth_max_atten{ii}{jj}{kk}(ll), thick_atten{ii}{jj}{kk}(ll)] ...
                            = deal(min(depth_bound_decim{ii}{jj}{kk}{ll}), max(depth_bound_decim{ii}{jj}{kk}{ll}), range(depth_bound_decim{ii}{jj}{kk}{ll}));
                
                % require minimum ice thickness
                if (thick_atten{ii}{jj}{kk}(ll) < (thick_rel_min * thick_decim{ii}{jj}{kk}(ll)))
                    continue
                end
                
                depth_inc_curr{ll} ...
                            = diff(depth_bound_decim{ii}{jj}{kk}{ll});
                [loss_curr{ll}, loss_uncert_curr{ll}] ...
                            = deal((int_smooth_corr{ii}{jj}{kk}(ind_layer_int{ii}{jj}{kk}{ll}(2:end), ll) - int_smooth_corr{ii}{jj}{kk}(ind_layer_int{ii}{jj}{kk}{ll}(1), ll)), ...
                                   sqrt((int_smooth_corr_uncert{ii}{jj}{kk}(ind_layer_int{ii}{jj}{kk}{ll}(2:end), ll) .^ 2) + (int_smooth_corr_uncert{ii}{jj}{kk}(ind_layer_int{ii}{jj}{kk}{ll}(1), ll) ^ 2)));
                
                % initialize impurity-concentration profiles as mean values
                [H_curr{ll}, Cl_curr{ll}, NH4_curr{ll}] ...
                            = deal(H_mean(ones(num_layer_curr(ll), 1)), Cl_mean(ones(num_layer_curr(ll), 1)), NH4_mean(ones(num_layer_curr(ll), 1)));
                
                % extract current ages
                age_curr_all= age{ii}{jj}{kk}(ind_layer_int{ii}{jj}{kk}{ll});
                age_curr    = age_curr_all(2:end);
                
                % rescale impurity concentration profiles if any dated layers
                if (do_chem && ~isempty(find(~isnan(age_curr_all), 1)))
                    
                    % dated vertical layers within current interval and boundary conditions
                    age_bound ...
                            = [0; age_curr_all(~isnan(age_curr_all)); Inf];
                    
                    % evaluate each vertical layer
                    for mm = 1:(length(age_bound) - 1)
                        ind_curr ...
                            = find(age_curr > age_bound(mm) & (age_curr <= age_bound(mm + 1))); % indices of dated layers within current layer
                        ind_chem_curr  ...
                            = find((age_chem > age_bound(mm)) & (age_chem <= age_bound(mm + 1))); % indices of chemistry values within current layer
                        if (~isempty(ind_curr) && ~isempty(ind_chem_curr))
                            [H_curr{ll}(ind_curr), Cl_curr{ll}(ind_curr), NH4_curr{ll}(ind_curr)] ...
                                = deal(mean(H(ind_chem_curr)), mean(Cl(ind_chem_curr)), mean(NH4(ind_chem_curr))); % assign mean value of chemistry within this layer
                        end
                    end
                end
                
                % good interval
                ind_good(ll)= true;
                
            end
            
            % limit by azimuth change
            if do_az
                ind_good(abs(diff(az_trans{ii}{jj}{kk})) > az_diff_max) ...
                            = false;
            end
            
            % skip sub-transect if no good traces
            if isempty(find(ind_good, 1))
                continue
            end
            
            % speed up core intersections
            if do_core
                ind_good_new= false(1, num_decim{ii}{jj}(kk));
                ind_good_new(ind_decim_curr) ...
                            = true;
                ind_good_new(~ind_good) ...
                            = false;
                ind_good    = ind_good_new;
            end
                       
            % slice key variables
            [depth_inc_curr, loss_curr, loss_uncert_curr, num_layer_curr, temp_surf_curr, H_curr, Cl_curr, NH4_curr] ...
                            = deal(depth_inc_curr(ind_good), loss_curr(ind_good), loss_uncert_curr(ind_good), num_layer_curr(ind_good), temp_surf_curr(ind_good), H_curr(ind_good), Cl_curr(ind_good), NH4_curr(ind_good));
            ind_good        = find(ind_good);
            
            % determine best-fit attenuation rate and equivalent temperature
            [atten_rate_curr, atten_rate_uncert_curr, mse_atten_curr, res_iso_curr, temp_iso_curr, temp_iso_uncert_curr] ...
                            = deal(NaN(1, length(ind_good)));
            
            if parallel_check
                parfor ll = 1:length(ind_good)
                    [atten_rate_tmp, atten_rate_uncert_tmp, mse_atten_curr(ll)] ...
                            = lscov([ones(num_layer_curr(ll), 1) cumsum(depth_inc_curr{ll})], loss_curr{ll}, (1 ./ (loss_uncert_curr{ll} .^ 2)));
                    [atten_rate_curr(ll), atten_rate_uncert_curr(ll)] ...
                            = deal((-0.5 * atten_rate_tmp(2)), (0.5 * atten_rate_uncert_tmp(2))); %#ok<PFOUS>
                    if (isnan(atten_rate_curr(ll)) || (atten_rate_curr(ll) < 0))
                        continue
                    end
                    try
                        [temp_iso_curr(ll), res_iso_curr(ll)] ...
                            = fminsearch(@atten_fit, temp_surf_curr(ll), opt, atten_rate_curr(ll), atten_rate_uncert_curr(ll), depth_inc_curr{ll}, H_curr{ll}, Cl_curr{ll}, NH4_curr{ll}, hfcp_struct);
                        if ~isnan(temp_iso_curr(ll))
                            tmp = atten_bound(frac_test, conf_uncert, atten_rate_curr(ll), atten_rate_uncert_curr(ll), depth_inc_curr{ll}, H_curr{ll}, Cl_curr{ll}, NH4_curr{ll}, hfcp_struct, temp_iso_curr(ll), res_iso_curr(ll));
                            temp_iso_uncert_curr(ll) ...
                                = nanmean([(temp_iso_curr(ll) - tmp(1)) (tmp(2) - temp_iso_curr(ll))]);
                        end
                    catch
                        disp([num2str([ii jj kk ll]) ' failed iso...'])
                    end
                end
            else
                for ll = 1:length(ind_good)
                    [atten_rate_tmp, atten_rate_uncert_tmp, mse_atten_curr(ll)] ...
                            = lscov([ones(num_layer_curr(ll), 1) cumsum(depth_inc_curr{ll})], loss_curr{ll}, (1 ./ (loss_uncert_curr{ll} .^ 2)));
                    [atten_rate_curr(ll), atten_rate_uncert_curr(ll)] ...
                            = deal((-0.5 * atten_rate_tmp(2)), (0.5 * atten_rate_uncert_tmp(2)));
                    if (isnan(atten_rate_curr(ll)) || (atten_rate_curr(ll) < 0))
                        continue
                    end
                    try
                        [temp_iso_curr(ll), res_iso_curr(ll)] ...
                            = fminsearch(@atten_fit, temp_surf_curr(ll), opt, atten_rate_curr(ll), atten_rate_uncert_curr(ll), depth_inc_curr{ll}, H_curr{ll}, Cl_curr{ll}, NH4_curr{ll}, hfcp_struct);
                        if ~isnan(temp_iso_curr(ll))
                            tmp = atten_bound(frac_test, conf_uncert, atten_rate_curr(ll), atten_rate_uncert_curr(ll), depth_inc_curr{ll}, H_curr{ll}, Cl_curr{ll}, NH4_curr{ll}, hfcp_struct, temp_iso_curr(ll), res_iso_curr(ll));
                            temp_iso_uncert_curr(ll) ...
                                = nanmean([(temp_iso_curr(ll) - tmp(1)) (tmp(2) - temp_iso_curr(ll))]);
                        end
                    catch
                        disp([num2str([ii jj kk ll]) ' failed iso...'])
                    end
                end
            end
            
            [atten_rate_curr(atten_rate_curr <= 0), atten_rate_uncert_curr(atten_rate_curr <= 0)] ...
                                = deal(NaN);
            
            % redistribute sliced attenuation-rate and best-fit temperature
            [atten_rate{ii}{jj}{kk}(ind_good), atten_rate_uncert{ii}{jj}{kk}(ind_good), mse_atten{ii}{jj}{kk}(ind_good), res_iso{ii}{jj}{kk}(ind_good), ...
             temp_iso{ii}{jj}{kk}(ind_good), temp_iso_uncert{ii}{jj}{kk}(ind_good), Cl_trans{ii}{jj}{kk}(ind_good), H_trans{ii}{jj}{kk}(ind_good), NH4_trans{ii}{jj}{kk}(ind_good)] ...
                                = deal((1e3 .* atten_rate_curr), (1e3 .* atten_rate_uncert_curr), mse_atten_curr, res_iso_curr, (temp_iso_curr - 273.15), temp_iso_uncert_curr, Cl_curr, H_curr, NH4_curr);
            
            % record attenuation rate and temperature if interval intersects an ice core
            if isempty(ind_core_decim_curr)
                continue
            end
            for mm = ind_core_decim_curr % loop through each decimated point
                if ~isnan(atten_rate{ii}{jj}{kk}(mm))
                    atten_core(ind_core_full_curr(mm == ind_core_decim_curr), 5:11) ...
                        = [mm atten_rate{ii}{jj}{kk}(mm) atten_rate_uncert{ii}{jj}{kk}(mm) temp_iso{ii}{jj}{kk}(mm) temp_iso_uncert{ii}{jj}{kk}(mm) depth_bound_decim{ii}{jj}{kk}{mm}([1 end])'];
                elseif ~isnan(atten_rate{ii}{jj}{kk}(mm - 1))
                    atten_core(ind_core_full_curr(mm == ind_core_decim_curr), 5:11) ...
                        = [(mm - 1) atten_rate{ii}{jj}{kk}(mm - 1) atten_rate_uncert{ii}{jj}{kk}(mm - 1) temp_iso{ii}{jj}{kk}(mm - 1) temp_iso_uncert{ii}{jj}{kk}(mm - 1) depth_bound_decim{ii}{jj}{kk}{mm - 1}([1 end])'];
                elseif ~isnan(atten_rate{ii}{jj}{kk}(mm + 1))
                    atten_core(ind_core_full_curr(mm == ind_core_decim_curr), 5:11) ...
                        = [(mm + 1) atten_rate{ii}{jj}{kk}(mm + 1) atten_rate_uncert{ii}{jj}{kk}(mm + 1) temp_iso{ii}{jj}{kk}(mm + 1) temp_iso_uncert{ii}{jj}{kk}(mm + 1) depth_bound_decim{ii}{jj}{kk}{mm + 1}([1 end])'];
                end
            end
        end
    end
end

if parallel_check
    pctRunOnAll warning('on', 'MATLAB:lscov:RankDefDesignMat')
else
    warning('on', 'MATLAB:lscov:RankDefDesignMat')
end

disp('Saving attenuation-rate modeling...')
if do_core
    if (exp_geo == 2)
        save(['mat/atten_core_' hfcp_model], '-v7.3', 'atten_core', 'atten_rate', 'atten_rate_uncert', 'az_trans', 'az_diff_max', 'depth_bound_decim', 'depth_min', 'depth_min_atten', 'depth_max', 'depth_max_atten', 'hfcp_model', 'ind_decim', 'ind_decim_mid', 'ind_layer_int', 'int_smooth_corr', ...
                                                      'int_smooth_corr_uncert', 'num_layer_min', 'res_iso', 'slope_layer', 'temp_iso', 'temp_iso_uncert', 'thick_atten', 'thick_decim', 'thick_rel_min', 'thick_rel_max', 'Cl_trans', 'H_trans', 'NH4_trans')
        disp(['Saved attenuation-rate modeling as mat/atten_core_' hfcp_model '.mat.'])
    else
        save(['mat/atten_core_' hfcp_model '_' num2str(exp_geo)], '-v7.3', 'atten_core', 'atten_rate', 'atten_rate_uncert', 'az_trans', 'az_diff_max', 'depth_bound_decim', 'depth_min', 'depth_min_atten', 'depth_max', 'depth_max_atten', 'hfcp_model', 'ind_decim', 'ind_decim_mid', ...
                                                                           'ind_layer_int', 'int_smooth_corr', 'int_smooth_corr_uncert', 'num_decim', 'num_layer_min', 'res_iso', 'slope_layer', 'temp_iso', 'temp_iso_uncert', 'thick_atten', 'thick_decim', 'thick_rel_min', 'thick_rel_max', ...
                                                                           'Cl_trans', 'H_trans', 'NH4_trans')
    end
else
    if do_chem
        if (exp_geo == 2)
            save mat/atten_all -v7.3 atten_core atten_rate atten_rate_uncert az_trans az_diff_max depth_bound_decim depth_min depth_min_atten depth_max depth_max_atten hfcp_model ind_decim ind_decim_mid ind_layer_int int_smooth_corr int_smooth_corr_uncert num_decim num_layer_min res_iso ...
                                     slope_layer temp_iso temp_iso_uncert thick_atten thick_decim thick_rel_min thick_rel_max Cl_trans H_trans NH4_trans
        else
            save(['mat/atten_core_' hfcp_model '_' num2str(exp_geo)], '-v7.3', 'atten_core', 'atten_rate', 'atten_rate_uncert', 'az_trans', 'az_diff_max', 'depth_bound_decim', 'depth_min', 'depth_min_atten', 'depth_max', 'depth_max_atten', 'hfcp_model', 'ind_decim', 'ind_decim_mid', ...
                                                                               'ind_layer_int', 'int_smooth_corr', 'int_smooth_corr_uncert', 'num_decim', 'num_layer_min', 'res_iso', 'slope_layer', 'temp_iso', 'temp_iso_uncert', 'thick_atten', 'thick_decim', 'thick_rel_min', 'thick_rel_max', ...
                                                                               'Cl_trans', 'H_trans', 'NH4_trans')
        end
    else
        save mat/atten_all_nochem -v7.3 temp_iso
    end
    disp('Saved attenuation-rate modeling as mat/atten_all.mat.')
end

if parallel_check
    delete(pool)
end