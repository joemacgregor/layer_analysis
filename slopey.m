% SLOPEY Analyze layer and bed slopes
% 
% Joe MacGregor (UTIG)
% Last updated: 08/07/13

do_load                     = true;

% load necessary transect data
if do_load
    clear all
    load mat/xy_all name_year name_trans num_trans num_year
    load mat/merge_all dist elev_bed_gimp elev_surf_gimp elev_smooth_gimp ind_fence num_fence num_layer
    do_save                 = true;
end

letters                     = 'a':'z';

dist_int_decim              = 5; % distance interval over which to calculate layer/bed slope
num_slope_min               = 5; % minimum number of non-NaN values over which to calculate layer/bed slope
num_layer_predict           = 5; % minimum number of layers whose slopes have been calculated with which to predict bed slope
depth_diff_max              = 1000; % maximum vertical distance between deepest layer and bed
ord_poly_layer              = 2; % order of polynomial fit between layer elevation and slope

% initializations
[ind_decim5, num_decim5, poly_slope_layer, res_slope_bed, slope_bed, slope_bed_predict, slope_layer] ...
                            = deal(cell(1, num_year));

warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')

% load all merged picks files from all years/campaigns that have been analyzed through FENCEGUI
for ii = 1:num_year
    
    disp(['Campaign: ' name_year{ii} '...'])
    
    [ind_decim5{ii}, num_decim5{ii}, poly_slope_layer{ii}, res_slope_bed{ii}, slope_bed{ii}, slope_bed_predict{ii}, slope_layer{ii}] ...
                            = deal(cell(1, num_trans(ii)));
    
    for jj = find(num_fence{ii})
        
        [ind_decim5{ii}{jj}, poly_slope_layer{ii}{jj}, res_slope_bed{ii}{jj}, slope_bed{ii}{jj}, slope_bed_predict{ii}{jj}, slope_layer{ii}{jj}] ...
                            = deal(cell(1, length(ind_fence{ii}{jj})));
        num_decim5{ii}{jj}  = zeros(1, length(ind_fence{ii}{jj}));
        
        for kk = find(ind_fence{ii}{jj})
            
            if (num_layer{ii}{jj}(kk) > 5)
                disp(['Transect: ' name_trans{ii}{jj} letters(kk) '...'])
            else
                continue
            end
            
            % break up data into intervals equal to dist_int_decim in kilometers
            dist_int        = dist{ii}{jj}{kk}(1):dist_int_decim:dist{ii}{jj}{kk}(end);
            dist_int(end)   = dist{ii}{jj}{kk}(end);
            ind_decim5{ii}{jj}{kk} ...
                            = unique(interp1(unique(dist{ii}{jj}{kk}), 1:length(unique(dist{ii}{jj}{kk})), dist_int, 'nearest', 'extrap'));
            num_decim5{ii}{jj}(kk) ...
                            = length(ind_decim5{ii}{jj}{kk}) - 1;
            
            poly_slope_layer{ii}{jj}{kk} ...
                            = NaN(3, num_decim5{ii}{jj}(kk));
            slope_layer{ii}{jj}{kk} ...
                            = NaN(num_layer{ii}{jj}(kk), num_decim5{ii}{jj}(kk));
            [slope_bed{ii}{jj}{kk}, slope_bed_predict{ii}{jj}{kk}] ...
                            = deal(NaN(1, num_decim5{ii}{jj}(kk)));
            
            % evaluate each interval
            for ll = 1:num_decim5{ii}{jj}(kk)
                
                ind_curr    = ind_decim5{ii}{jj}{kk}(ll):ind_decim5{ii}{jj}{kk}(ll + 1);
                dist_curr   = dist{ii}{jj}{kk}(ind_curr);
                elev_smooth_decim ...
                            = NaN(num_layer{ii}{jj}(kk), 1);
                
                for mm = 1:num_layer{ii}{jj}(kk)
                    
                    elev_layer_curr ...
                            = elev_smooth_gimp{ii}{jj}{kk}(mm, ind_curr);
                    
                    % require a minimum number of values required for layer slope measurement
                    if (length(find(~isnan(elev_layer_curr))) < num_slope_min)
                        continue
                    end
                    
                    % measure layer slope (m/km) for all available layers within the current horizontal interval
                    slope_layer_tmp ...
                            = polyfit(dist_curr(~isnan(elev_layer_curr)), elev_layer_curr(~isnan(elev_layer_curr)), 1);
                    slope_layer{ii}{jj}{kk}(mm, ll) ...
                            = slope_layer_tmp(1);
                    elev_smooth_decim(mm) ...
                            = nanmean(elev_layer_curr);
                end
                
                % again require a minimum number of values to measure bed slope (m/km)
                if (length(find(~isnan(elev_bed_gimp{ii}{jj}{kk}(ind_curr)))) < num_slope_min)
                    continue
                end
                
                elev_bed_curr ...
                            = elev_bed_gimp{ii}{jj}{kk}(ind_curr);
                slope_bed_tmp ...
                            = polyfit(dist_curr(~isnan(elev_bed_gimp{ii}{jj}{kk}(ind_curr))), elev_bed_curr(~isnan(elev_bed_curr)), 1);
                slope_bed{ii}{jj}{kk}(ll) ...
                            = slope_bed_tmp(1);
                elev_bed_decim ...
                            = nanmean(elev_bed_curr);
                
                % require a minimum number of layers for bed slope prediction and a layer that reaches a minimum normalized depth
                slope_layer_curr ...
                            = slope_layer{ii}{jj}{kk}(:, ll);
                if ((length(find(~isnan(slope_layer_curr))) < num_layer_predict) || ((nanmin(elev_smooth_decim) - elev_bed_decim) > depth_diff_max))
                    continue
                end
                
                % polynomial fit between layer elevation and slope
                poly_slope_layer{ii}{jj}{kk}(:, ll) ...
                            = polyfit(elev_smooth_decim(~isnan(slope_layer_curr)), slope_layer_curr(~isnan(slope_layer_curr)), ord_poly_layer)';
                
                % layer-predicted bed slope based on layer elevation-slope polynomial fit and mean bed elevation over the interval, m/km
                slope_bed_predict{ii}{jj}{kk}(ll) ...
                            = polyval(poly_slope_layer{ii}{jj}{kk}(:, ll), elev_bed_decim);
            end
            
            % residual between measured and layer-predicted bed slope, m/km
            res_slope_bed{ii}{jj}{kk} ...
                            = slope_bed{ii}{jj}{kk} - slope_bed_predict{ii}{jj}{kk};
        end
    end
end

warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale')

if do_save
    disp('Saving slope data...')
    save('mat/slope_all', '-v7.3', 'ind_decim5', 'num_decim5', 'poly_slope_layer', 'res_slope_bed', 'slope_bed', 'slope_bed_predict', 'slope_layer')
    disp('Saved slope data in mat/slope_all.mat.')
end