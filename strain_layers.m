function strain_layers(age_max, radar_type, do_strain, do_fix, do_uncert, do_bad, do_smooth, do_grd, do_save, num_smooth, strain_rate_ref, strain_rate_surf)
% STRAIN_LAYERS 1-D ice-flow modeling using dated layers.
% 
% STRAIN_LAYERS(AGE_MAX,RADAR_TYPE,DO_STRAIN,DO_FIX,DO_UNCERT,DO_BAD,DO_SMOOTH,DO_GRD,DO_SAVE,NUM_SMOOTH,STRAIN_RATE_REF,SPREAD_RATE)
% calculates the best-fit parameters for several 1-D ice-flow models using
% dated layers younger than or equal to AGE_MAX (units of years).
% RADAR_TYPE is a string that determines which radar dataset to analyze,
% e.g., 'deep' or 'accum'. DO_STRAIN is a 7-element logical vector that
% determines which models to run (true). The order is: Nye,
% Dansgaard-Johnsen, Nye+melt, Dansgaard-Johnsen+melt,
% Dansgaard-Johnsen+lateral strain, Dansgaard-Johnsen+surface
% strain,shallow strain. DO_FIX, DO_UNCERT, DO_BAD, DO_SMOOTH, DO_GRD and
% DO_SAVE are logical scalars that determine whether the vertical strain
% rate is fixed in the shallow strain model, uncertainties are calculated,
% whether bad traces are NaN'd out, whether the along-transect values are
% smoothed, whether grids are calculated and whether the results are saved,
% respectively. If DO_SMOOTH=TRUE, then the number of samples NUM_SMOOTH
% over which to smooth the along-transect ice-flow model parameters must
% also be given. If DO_FIX=TRUE, then a previously modeled vertical strain
% rate is used instead (STRAIN_RATE_REF). STRAIN_RATE_SURF is a structure
% that contains the longitudinal and lateral strain rates at the surface
% used by DJ_LAT_FIT/DJ_SURF_FIT and their x/y grids.
% 
% See also DJ_FIT, DJ_LAT_FIT, DJ_MELT_FIT, NYE_FIT, NYE_MELT_FIT,
% SHALLOW_STRAIN_FIT, SMOOTH_LOWESS and STRAIN_UNCERT.
% 
% Joe MacGregor (UTIG), Mark Fahnestock (UAF)
% Last updated: 08/20/15

if ~exist('dj_fit', 'file')
    error('strain_layers:djfit', 'Function DJ_FIT is not available within this user''s path.')
end
if ~exist('dj_lat_fit', 'file')
    error('strain_layers:djlatfit', 'Function DJ_LAT_FIT is not available within this user''s path.')
end
if ~exist('dj_melt_fit', 'file')
    error('strain_layers:djmeltfit', 'Function DJ_MELT_FIT is not available within this user''s path.')
end
if ~exist('dj_surf_fit', 'file')
    error('strain_layers:djsurffit', 'Function DJ_SURF_FIT is not available within this user''s path.')
end
if ~exist('nye_fit', 'file')
    error('strain_layers:nyefit', 'Function NYE_FIT is not available within this user''s path.')
end
if ~exist('nye_melt_fit', 'file')
    error('strain_layers:nyemeltfit', 'Function NYE_MELT_FIT is not available within this user''s path.')
end
if ~exist('shallow_strain_fit', 'file')
    error('strain_layers:shallowstrainfit', 'Function SHALLOW_STRAIN_FIT is not available within this user''s path.')
end
if ~exist('shallow_strain_fix', 'file')
    error('strain_layers:shallowstrainfix', 'Function SHALLOW_STRAIN_FIX is not available within this user''s path.')
end
if ~exist('strain_uncert', 'file')
    error('strain_layers:strainuncert', 'Function STRAIN_UNCERT is not available within this user''s path.')
end
if ~exist('smooth_lowess', 'file')
    error('strain_layers:smoothlowess', 'Function SMOOTH_LOWESS is not available within this user''s path.')
end
if ~any(nargin == 9:12)
    error('strain_layers:nargin', ['Number of arguments (' num2str(nargin) ') is not 9, 10, 11 or 12.'])
end
if (~ischar(radar_type) || ~any(strcmp(radar_type, {'deep' 'accum'})))
    error('strain_layers:radar_type', 'RADAR_TYPE is not a string equal to ''deep'' or ''accum''.')
end
if (~islogical(do_strain) || ~isvector(do_strain) || (length(do_strain) ~= 7) || ~any(do_strain))
    error('strain_layers:dostrain', 'DO_STRAIN is not a logical 7-element vector that contains at least one true value.')
end
if (~islogical(do_fix) || ~isscalar(do_fix))
    error('strain_layers:dofix', 'DO_FIX is not a logical scalar.')
end
if (~islogical(do_uncert) || ~isscalar(do_uncert))
    error('strain_layers:douncert', 'DO_UNCERT is not a logical scalar.')
end
if (~islogical(do_bad) || ~isscalar(do_bad))
    error('strain_layers:dobad', 'DO_BAD is not a logical scalar.')
end
if (~islogical(do_smooth) || ~isscalar(do_smooth))
    error('strain_layers:dosmooth', 'DO_SMOOTH is not a logical scalar.')
end
if (~islogical(do_grd) || ~isscalar(do_grd))
    error('strain_layers:dogrd', 'DO_GRD is not a logical scalar.')
end
if (~islogical(do_save) || ~isscalar(do_save))
    error('strain_layers:dosave', 'DO_SAVE is not a logical scalar.')
end
if (do_smooth && (nargin < 10))
    error('strain_layers:narginsmooth', 'Number of arguments must be at least 10 if DO_SMOOTH=TRUE.')
end
if do_smooth
    if (~isnumeric(num_smooth) || ~isscalar(num_smooth) || mod(num_smooth, 1))
        error('strain_layers:numsmooth', 'NUM_SMOOTH is not a numeric scalar integer.')
    end
end
if (do_fix && (nargin ~= 11))
    error('strain_layers:narginfix', 'Number of arguments must be 11 if DO_FIX=TRUE.')
end
if do_fix
    if ~iscell(strain_rate_ref)
        error('strain_layers:strainrateref', 'STRAIN_RATE_REF is not a cell.')
    end
end
if (any(do_strain(5:6)) && (nargin ~= 12))
    error('strain_layers:latnargin', 'ANY(DO_STRAIN(5:6))=TRUE but STRAIN_RATE_SURF is not an argument.')
end
if (nargin == 12)
    if ~isstruct(strain_rate_surf)
        error('strain_layers:strainratesurf', 'STRAIN_RATE_SURF is not a structure.')
    end
end
if nargout
    error('strain_layers:nargout', 'STRAIN_LAYERS has no outputs.')
end

disp('Starting strain-rate modeling...')

switch radar_type
    case 'accum'
        load mat/xy_all_accum name_year name_trans num_year num_trans
        load mat/merge_all_accum depth_smooth ind_decim ind_decim_mid num_decim num_subtrans thick_decim x_pk y_pk
        load mat/date_all_accum age age_uncert
    case 'deep'
        load mat/xy_all name_year name_trans num_year num_trans
        load mat/merge_all depth_smooth ind_decim ind_decim_mid num_decim num_subtrans thick_decim x_pk y_pk
        load mat/date_all age age_uncert
        % temporary fix for incomplete 2014 P3 tracing
        if (num_year > 20)
            num_year        = 20;
        end
end

% initial strain-rate model parameters
thick_shear_def             = 200; % initial basal shear layer thickness, m ("little h")
melt_rate_def               = 0; % initial basal melt rate, m/a
frac_slide_def              = 0; % initial sliding fraction (dimensionless)
frac_test                   = [0.01 0.025 0.05 0.1 0.25 0.5 1 5 10]; % range around best-fit model parameter about which to test
conf_uncert                 = 0.95; % confidence bound for uncertainty calculations, e.g., 0.95 is 95%

% sub-transect letters
letters                     = 'a':'z';

% set names of strain-rate model parameters
strain_param_def            = {'accum_nye_start' 'ind_strain_fail' 'num_layer_strain'};
strain_param                = {{'accum_nye' 'res_nye'} {'accum_dj' 'res_dj' 'shape_fact' 'thick_shear'} {'accum_nye_melt' 'melt_bed' 'res_nye_melt'} {'accum_dj_melt' 'frac_slide' 'melt_bed_dj' 'res_dj_melt' 'shape_fact_melt' 'thick_shear_melt'} ...
                               {'accum_dj_lat' 'res_dj_lat' 'shape_fact_lat' 'thick_shear_lat'} {'res_dj_surf' 'shape_fact_surf' 'thick_shear_surf'} {'accum_shallow' 'strain_rate' 'res_shallow'}};
if do_fix
    strain_param{7}         = strain_param{7}([1 3]);
end
num_param_def               = length(strain_param_def);
num_param                   = zeros(1, length(strain_param));
for ii = 1:length(strain_param)
    num_param(ii)           = length(strain_param{ii});
end

% initialize model parameters
for ii = 1:num_param_def
    eval([strain_param_def{ii} ' = cell(1, num_year);']);
end
for ii = find(do_strain)
    for jj = 1:num_param(ii)
        eval([strain_param{ii}{jj} ' = cell(1, num_year);']);
        if do_uncert
            eval([strain_param{ii}{jj} '_uncert = cell(1, num_year);']);
        end
    end
end
for ii = 1:num_year
    for jj = 1:num_param_def
        eval([strain_param_def{jj} '{ii} = cell(1, num_trans(ii));']);
    end
    for jj = find(do_strain)
        for kk = 1:num_param(jj)
            eval([strain_param{jj}{kk} '{ii} = cell(1, num_trans(ii));']);
            if do_uncert
                eval([strain_param{jj}{kk} '_uncert{ii} = cell(1, num_trans(ii));']);
            end
        end
    end
    for jj = 1:num_trans(ii)
        for kk = 1:num_param_def
            eval([strain_param_def{kk} '{ii}{jj} = cell(1, num_subtrans{ii}(jj));']);
        end
        for kk = find(do_strain)
            for ll = 1:num_param(kk)
                eval([strain_param{kk}{ll} '{ii}{jj} = cell(1, num_subtrans{ii}(jj));']);
                if do_uncert
                    eval([strain_param{kk}{ll} '_uncert{ii}{jj} = cell(1, num_subtrans{ii}(jj));']);
                end
            end
        end
        for kk = 1:num_subtrans{ii}(jj)
            for ll = 1:num_param_def
                eval([strain_param_def{ll} '{ii}{jj}{kk} = NaN(1, num_decim{ii}{jj}(kk));']);
            end
            for ll = find(do_strain)
                for mm = 1:num_param(ll)
                    eval([strain_param{ll}{mm} '{ii}{jj}{kk} = NaN(1, num_decim{ii}{jj}(kk));']);
                    if do_uncert
                        eval([strain_param{ll}{mm} '_uncert{ii}{jj}{kk} = NaN(2, num_decim{ii}{jj}(kk));']);
                    end
                end
            end
        end
    end
end

if any(do_strain(5:6))
    try
        [x_grd_strain, y_grd_strain, lon_strain_rate_grd, lat_strain_rate_grd] ...
                            = deal(strain_rate_surf.x_grd_strain, strain_rate_surf.y_grd_strain, strain_rate_surf.lon_strain_rate_grd, strain_rate_surf.lat_strain_rate_grd);
    catch me
        disp(['does not contain the required variables. Error: ' me.message])
        return
    end
end

% fminsearch options setup
opt                         = optimset('fminsearch');
opt.Display                 = 'off'; % avoid display of having reach maximum number of iterations
opt.FunValCheck             = 'on'; % sometimes strain-rate models get complex/NaN due to, e.g., logs of a negative argument; this switch stops that

% if license('checkout', 'distrib_computing_toolbox')
%     pool_check              = gcp('nocreate');
%     if isempty(pool_check)
%         try
%             pool            = parpool('local', 4);
%         catch
%             pool            = parpool('local');
%         end
%     end
%     parallel_check          = true;
% else
    parallel_check          = false;
% end

for ii = 1:num_year
    
    disp(['Campaign: ' name_year{ii} '...'])
    
    for jj = 1:num_trans(ii)
        
        for kk = 1:num_subtrans{ii}(jj)
            
            if isempty(age{ii}{jj}{kk})
                continue
            end
            
            % layers that are dated and whose age is less than age_max
            ind_layer_good  = find(~isnan(age{ii}{jj}{kk}) & (age{ii}{jj}{kk} < age_max));
            
            if (length(ind_layer_good) < 4) % insufficient layer ages
                continue
            else
                disp([name_trans{ii}{jj} letters(kk) '...'])
            end
            
            depth_smooth_decim ...
                            = NaN(length(ind_layer_good), num_decim{ii}{jj}(kk));
            for ll = 1:num_decim{ii}{jj}(kk)
                depth_smooth_decim(:, ll) ...
                            = mean(depth_smooth{ii}{jj}{kk}(ind_layer_good, ind_decim{ii}{jj}{kk}(ll):ind_decim{ii}{jj}{kk}(ll + 1)), 2, 'omitnan');
            end
            
            % set uncertainty of zero-uncertainty layers to mean relative age uncertainty of other layers
            age_uncert_trim = age_uncert{ii}{jj}{kk}(ind_layer_good);
            age_uncert_trim(~age_uncert_trim) ...
                            = age{ii}{jj}{kk}(ind_layer_good(~age_uncert_trim)) .* mean(age_uncert_trim(age_uncert_trim > 0) ./ age{ii}{jj}{kk}(ind_layer_good(age_uncert_trim > 0)));
            age_trim        = age{ii}{jj}{kk}(ind_layer_good);
            
            if parallel_check
                
                % iterate at traces with at least four dated layers and an ice thickness
                curr_ind    = find((sum(~isnan(depth_smooth_decim)) > 3) & ~isnan(thick_decim{ii}{jj}{kk}));
                thick_slice = thick_decim{ii}{jj}{kk}(curr_ind);
                
                % dated layers at current trace
                [curr_age, curr_age_uncert, depth_smooth_cell] ...
                            = deal(cell(1, length(curr_ind)));
                for ll = 1:length(curr_ind)
                    curr_layer ...
                            = find(~isnan(depth_smooth_decim(:, curr_ind(ll))));
                    curr_age{ll} ...
                            = age_trim(curr_layer);
                    curr_age_uncert{ll} ...
                            = age_uncert_trim(curr_layer);
                    depth_smooth_cell{ll} ...
                            = depth_smooth_decim(curr_layer, curr_ind(ll));
                end
                
                % initialize slices
                for ll = 1:num_param_def
                    eval([strain_param_def{ll} '_slice = NaN(1, length(curr_ind));']);
                end
                for ll = find(do_strain)
                    for mm = 1:num_param(ll)
                        eval([strain_param{ll}{mm} '_slice = NaN(1, length(curr_ind));']);
                        if do_uncert
                            eval([strain_param{ll}{mm} '_uncert_slice = NaN(2, length(curr_ind));']);
                        end
                    end
                end
                
                if do_fix
                    strain_rate_slice ...
                            = strain_rate_ref{ii}{jj}{kk}(curr_ind);
                end
                
                if any(do_strain(5:6))
                    [lon_strain_rate_slice, lat_strain_rate_slice] ...
                            = deal(interp2(x_grd_strain, y_grd_strain, lon_strain_rate_grd, x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}(curr_ind)), y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}(curr_ind)), '*linear', NaN), ...
                                   interp2(x_grd_strain, y_grd_strain, lat_strain_rate_grd, x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}(curr_ind)), y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}(curr_ind)), '*linear', NaN));
                end
                
                parfor ll = 1:length(curr_ind)
                    
                    % number layers used in current model
                    num_layer_strain_slice(ll) ...
                            = length(curr_age{ll}); %#ok<PFOUS,NASGU>
                    
                    % LLA accumulation rate, i.e., Nye model
                    [accum_nye_start_slice(ll), accum_ref] ...
                            = deal(mean(-log(1 - (depth_smooth_cell{ll} ./ thick_slice(ll))) .* (thick_slice(ll) ./ curr_age{ll}))); %#ok<PFOUS>
                    
                    if do_strain(1) %#ok<PFBNS>
                        try % standard Nye model (iterative)
                            [accum_nye_slice(ll), res_nye_slice(ll)] ...
                                = fminsearch(@nye_fit, accum_nye_start_slice(ll), opt, thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}); %#ok<PFOUS>
                            if do_uncert
                                accum_nye_uncert_slice(:, ll) ...
                                    = strain_uncert(frac_test, conf_uncert, 'nye', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, accum_nye_slice(ll), res_nye_slice(ll)); %#ok<PFOUS,NASGU>
                            end
                            accum_ref ...
                                = accum_nye_slice(ll);
                        catch
                            ind_strain_fail_slice(ll) ...
                                = 1; %#ok<PFOUS>
                        end
                    end
                    
                    if do_strain(2)
                        try % Dansgaard-Johnsen model (original)
                            [tmp1, res_dj_slice(ll)] ...
                                = fminsearch(@dj_fit, [accum_ref thick_shear_def], opt, thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}); %#ok<PFOUS>
                            [accum_dj_slice(ll), thick_shear_slice(ll)] ...
                                = deal(tmp1(1), tmp1(2)); %#ok<NASGU,PFOUS>
                            if do_uncert
                                [accum_dj_uncert_slice(:, ll), thick_shear_uncert_slice(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, tmp1, res_dj_slice(ll)); %#ok<PFOUS,NASGU>
                            end
                            accum_ref ...
                                = accum_dj_slice(ll);
                        catch
                            ind_strain_fail_slice(ll) ...
                                = 2;
                        end
                    end
                    
                    if do_strain(3)
                        try % Nye + melt model
                            [tmp2, res_nye_melt_slice(ll)] ...
                                = fminsearch(@nye_melt_fit, [accum_ref melt_rate_def], opt, thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}); %#ok<PFOUS>
                            [accum_nye_melt_slice(ll), melt_bed_slice(ll)] ...
                                = deal(tmp2(1), tmp2(2)); %#ok<PFOUS>
                            if do_uncert
                                [accum_nye_melt_uncert_slice(:, ll), melt_bed_uncert_slice(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'nye_melt', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, tmp2, res_nye_melt_slice(ll)); %#ok<PFOUS,NASGU>
                            end
                            accum_ref ...
                                = accum_nye_melt_slice(ll);
                        catch
                            ind_strain_fail_slice(ll) ...
                                = 3;
                        end
                    end
                    
                    if do_strain(4)
                        try % Dansgaard-Johnsen model + basal melting
                            [tmp3, res_dj_melt_slice(ll)] ...
                                = fminsearch(@dj_melt_fit, [accum_ref thick_shear_def melt_bed_slice(ll) frac_slide_def], opt, thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}); %#ok<PFOUS>
                            [accum_dj_melt_slice(ll), thick_shear_melt_slice(ll), melt_bed_dj_slice(ll), frac_slide_slice(ll)] ...
                                = deal(tmp3(1), tmp3(2), tmp3(3), tmp3(4)); %#ok<PFOUS,NASGU>
                            if do_uncert
                                [accum_dj_melt_uncert_slice(:, ll), thick_shear_melt_uncert_slice(:, ll), melt_bed_dj_uncert_slice(:, ll), frac_slide_uncert_slice(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj_melt', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, tmp3, res_dj_melt_slice(ll)); %#ok<PFOUS,NASGU>
                            end
                            accum_ref ...
                                = accum_dj_melt_slice(ll);
                        catch
                            ind_strain_fail_slice(ll) ...
                                = 4;
                        end
                    end
                    
                    if do_strain(5)
                        try % Dansgaard-Johnsen model + lateral strain
                            [tmp4, res_dj_lat_slice(ll)] ...
                                = fminsearch(@dj_lat_fit, [accum_ref thick_shear_def], opt, thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, lon_strain_rate_slice(ll), lat_strain_rate_slice(ll)); %#ok<PFOUS>
                            [accum_dj_lat_slice(ll), thick_shear_lat_slice(ll)] ...
                                = deal(tmp4(1), tmp4(2)); %#ok<PFOUS,NASGU>
                            if do_uncert
                                [accum_dj_lat_uncert_slice(:, ll), thick_shear_lat_uncert_slice(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj_lat', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, tmp4, res_dj_lat_slice(ll), lon_strain_rate_slice(ll), lat_strain_rate_slice(ll)); %#ok<PFOUS,NASGU>
                            end
                            accum_ref ...
                                = accum_dj_lat_slice(ll);
                        catch
                            ind_strain_fail_slice(ll) ...
                                = 5;
                        end
                    end

                    if do_strain(6)
                        try % Dansgaard-Johnsen model + surface strain
                            [thick_shear_surf_slice(ll), res_dj_surf_slice(ll)] ...
                                = fminsearch(@dj_surf_fit, thick_shear_def, opt, thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, lon_strain_rate_slice(ll), lat_strain_rate_slice(ll)); %#ok<PFOUS>
                            if do_uncert
                                thick_shear_surf_uncert_slice(:, ll) ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj_surf', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, thick_shear_surf_slice(ll), res_dj_surf_slice(ll), lon_strain_rate_slice(ll), lat_strain_rate_slice(ll)); %#ok<PFOUS,NASGU>
                            end
                        catch
                            ind_strain_fail_slice(ll) ...
                                = 6;
                        end
                    end
                    
                    if do_strain(7)
                        try % shallow strain model
                            if do_fix
                                [accum_shallow_slice(ll), res_shallow_slice(ll)] ...
                                    = fminsearch(@shallow_strain_fix, accum_ref, opt, depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, strain_rate_slice(ll)); %#ok<PFOUS>
                                if do_uncert
                                    accum_shallow_uncert_slice(:, ll) ...
                                        = strain_uncert(frac_test, conf_uncert, 'shallow_strain', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, accum_shallow_slice(ll), res_shallow_slice(ll)); %#ok<PFOUS>
                                end
                            else
                                [tmp5, res_shallow_slice(ll)] ...
                                    = fminsearch(@shallow_strain_fit, [accum_ref (accum_ref / thick_slice(ll))], opt, depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll});
                                [accum_shallow_slice(ll), strain_rate_slice(ll)] ...
                                    = deal(tmp5(1), tmp5(2));
                                if do_uncert
                                    [accum_shallow_uncert_slice(:, ll), strain_rate_uncert_slice(:, ll)] ...
                                        = strain_uncert(frac_test, conf_uncert, 'shallow_strain', thick_slice(ll), depth_smooth_cell{ll}, curr_age{ll}, curr_age_uncert{ll}, tmp5, res_shallow_slice(ll)); %#ok<PFOUS,NASGU>
                                end
                            end
                        catch
                            ind_strain_fail_slice(ll) ...
                                = 7;
                        end
                    end
                end
                
                % reassign slices
                for ll = 1:num_param_def
                    eval([strain_param_def{ll} '{ii}{jj}{kk}(curr_ind) = ' strain_param_def{ll} '_slice;']);
                end
                for ll = find(do_strain)
                    for mm = 1:num_param(ll)
                        eval([strain_param{ll}{mm} '{ii}{jj}{kk}(curr_ind) = ' strain_param{ll}{mm} '_slice;']);
                        if do_uncert
                            eval([strain_param{ll}{mm} '_uncert{ii}{jj}{kk}(:, curr_ind) = ' strain_param{ll}{mm} '_uncert_slice;']);
                        end
                    end
                end
                
            else
                
                ind_strain_fail{ii}{jj}{kk} ...
                            = NaN(2, num_decim{ii}{jj}(kk)); %#ok<AGROW>
                
                % current along-track longitudinal and lateral strain rates
                if any(do_strain(5:6))
                    [lon_strain_rate_curr, lat_strain_rate_curr] ...
                            = deal(interp2(x_grd_strain, y_grd_strain, lon_strain_rate_grd, x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}), y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}), '*linear', NaN), ...
                                   interp2(x_grd_strain, y_grd_strain, lat_strain_rate_grd, x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}), y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk}), '*linear', NaN));
                end
                
                % iterate at traces with at least two dated layers and an ice thickness
                for ll = find((sum(~isnan(depth_smooth_decim)) > 3) & ~isnan(thick_decim{ii}{jj}{kk}))
                    
                    % dated layers at current trace
                    curr_layer ...
                            = find(~isnan(depth_smooth_decim(:, ll)));
                    num_layer_strain{ii}{jj}{kk}(ll) ...
                            = length(curr_layer); %#ok<AGROW,NASGU>
                    
                    % initial Nye accumulation rate
                    [accum_nye_start{ii}{jj}{kk}(ll), accum_ref] ...
                            = deal(mean(-log(1 - (depth_smooth_decim(curr_layer, ll) ./ thick_decim{ii}{jj}{kk}(ll))) .* (thick_decim{ii}{jj}{kk}(ll) ./ age_trim(curr_layer)))); %#ok<AGROW>
                    
                    if do_strain(1)
                        try % standard Nye model (iterative)
                            [accum_nye{ii}{jj}{kk}(ll), res_nye{ii}{jj}{kk}(ll)] ...
                                = fminsearch(@nye_fit, accum_nye_start{ii}{jj}{kk}(ll), opt, thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer)); %#ok<AGROW>
                            if do_uncert
                                accum_nye_uncert{ii}{jj}{kk}(:, ll) ...
                                    = strain_uncert(frac_test, conf_uncert, 'nye', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), accum_nye{ii}{jj}{kk}(ll), res_nye{ii}{jj}{kk}(ll)); %#ok<AGROW,NASGU>
                            end
                            accum_ref ...
                                = accum_nye{ii}{jj}{kk}(ll);
                        catch
                            ind_strain_fail{ii}{jj}{kk}(:, ll) ...
                                = [ll; 1];
                        end
                    end
                    
                    if do_strain(2)
                        try % Dansgaard-Johnsen model (original)
                            [tmp6, res_dj{ii}{jj}{kk}(ll)] ...
                                = fminsearch(@dj_fit, [accum_ref thick_shear_def], opt, thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer)); %#ok<AGROW>
                            [accum_dj{ii}{jj}{kk}(ll), thick_shear{ii}{jj}{kk}(ll)] ...
                                = deal(tmp6(1), tmp6(2)); %#ok<AGROW>
                            if do_uncert
                                [accum_dj_uncert{ii}{jj}{kk}(:, ll), thick_shear_uncert{ii}{jj}{kk}(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), tmp6, res_dj{ii}{jj}{kk}(ll)); %#ok<AGROW,ASGLU>
                            end
                            accum_ref ...
                                = accum_dj{ii}{jj}{kk}(ll);
                        catch
                            ind_strain_fail{ii}{jj}{kk}(:, ll) ...
                                = [ll; 2];
                        end
                    end
                    
                    if do_strain(3)
                        try % Nye + melt model
                            [tmp7, res_nye_melt{ii}{jj}{kk}(ll)] ...
                                = fminsearch(@nye_melt_fit, [accum_ref melt_rate_def], opt, thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer)); %#ok<AGROW>
                            [accum_nye_melt{ii}{jj}{kk}(ll), melt_bed{ii}{jj}{kk}(ll)] ...
                                = deal(tmp7(1), tmp7(2)); %#ok<AGROW>
                            if do_uncert
                                [accum_nye_melt_uncert{ii}{jj}{kk}(:, ll), melt_bed_uncert{ii}{jj}{kk}(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'nye_melt', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), tmp7, res_nye_melt{ii}{jj}{kk}(ll)); %#ok<AGROW,NASGU>
                            end
                            accum_ref ...
                                = accum_nye_melt{ii}{jj}{kk}(ll);
                        catch
                            ind_strain_fail{ii}{jj}{kk}(:, ll) ...
                                = [ll; 3];
                        end
                    end
                    
                    if do_strain(4)
                        try % Dansgaard-Johnsen model + basal melting
                            [tmp8, res_dj_melt{ii}{jj}{kk}(ll)] ...
                                = fminsearch(@dj_melt_fit, [accum_ref thick_shear_def melt_bed{ii}{jj}{kk}(ll) frac_slide_def], opt, thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer)); %#ok<AGROW>
                            [accum_dj_melt{ii}{jj}{kk}(ll), thick_shear_melt{ii}{jj}{kk}(ll), melt_bed_dj{ii}{jj}{kk}(ll), frac_slide{ii}{jj}{kk}(ll)] ...
                                = deal(tmp8(1), tmp8(2), tmp8(3), tmp8(4)); %#ok<NASGU,AGROW>
                            if do_uncert
                                [accum_dj_melt_uncert{ii}{jj}{kk}(:, ll), thick_shear_melt_uncert{ii}{jj}{kk}(:, ll), melt_bed_dj_uncert{ii}{jj}{kk}(:, ll), frac_slide_uncert{ii}{jj}{kk}(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj_melt', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), tmp8, res_dj_melt{ii}{jj}{kk}(ll)); %#ok<NASGU,AGROW,ASGLU>
                            end
                            accum_ref ...
                                = accum_dj_melt{ii}{jj}{kk}(ll);
                        catch
                            ind_strain_fail{ii}{jj}{kk}(:, ll) ...
                                = [ll; 4];
                        end
                    end
                    
                    if do_strain(5)
                        try % Dansgaard-Johnsen model + lateral strain
                            [tmp9, res_dj_lat{ii}{jj}{kk}(ll)] ...
                                = fminsearch(@dj_lat_fit, [accum_ref thick_shear_def], opt, thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), lon_strain_rate_curr(ll), lat_strain_rate_curr(ll)); %#ok<AGROW>
                            [accum_dj_lat{ii}{jj}{kk}(ll), thick_shear_lat{ii}{jj}{kk}(ll)] ...
                                = deal(tmp9(1), tmp9(2)); %#ok<AGROW>
                            if do_uncert
                                [accum_dj_lat_uncert{ii}{jj}{kk}(:, ll), thick_shear_lat_uncert{ii}{jj}{kk}(:, ll)] ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj_lat', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), tmp9, res_dj_lat{ii}{jj}{kk}(ll), lon_strain_rate_curr(ll), lat_strain_rate_curr(ll)); ...
                                      %#ok<AGROW,ASGLU>
                            end
                            accum_ref ...
                                = accum_dj_lat{ii}{jj}{kk}(ll);
                        catch
                            ind_strain_fail{ii}{jj}{kk}(:, ll) ...
                                = [ll; 5];
                        end
                    end
                    
                    if do_strain(6)
                        try % Dansgaard-Johnsen model + surface strain
                            [thick_shear_surf{ii}{jj}{kk}(ll), res_dj_surf{ii}{jj}{kk}(ll)] ...
                                = fminsearch(@dj_surf_fit, thick_shear_def, opt, thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), lon_strain_rate_curr(ll), lat_strain_rate_curr(ll)); %#ok<AGROW>
                            if do_uncert
                                thick_shear_surf_uncert{ii}{jj}{kk}(:, ll) ...
                                    = strain_uncert(frac_test, conf_uncert, 'dj_surf', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), thick_shear_surf{ii}{jj}{kk}(ll), res_dj_surf{ii}{jj}{kk}(ll), lon_strain_rate_curr(ll), ...
                                                    lat_strain_rate_curr(ll)); %#ok<AGROW>
                            end
                        catch
                            ind_strain_fail{ii}{jj}{kk}(:, ll) ...
                                = [ll; 6];
                        end
                    end
                    
                    if do_strain(7)
                        try % shallow strain layer model
                            if do_fix
                                [accum_shallow{ii}{jj}{kk}(ll), res_shallow{ii}{jj}{kk}(ll)] ...
                                    = fminsearch(@shallow_strain_fix, accum_ref, opt, depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), strain_rate_ref{ii}{jj}{kk}(ll)); %#ok<AGROW>
                                if do_uncert
                                    accum_shallow_uncert{ii}{jj}{kk}(:, ll) ...
                                        = strain_uncert(frac_test, conf_uncert, 'shallow_strain', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), accum_shallow{ii}{jj}{kk}(ll), res_shallow{ii}{jj}{kk}(ll)); %#ok<AGROW>
                                end
                            else
                                [tmp10, res_shallow{ii}{jj}{kk}(ll)] ...
                                    = fminsearch(@shallow_strain_fit, [accum_ref (accum_ref / thick_decim{ii}{jj}{kk}(ll))], opt, depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer)); %#ok<AGROW>
                                [accum_shallow{ii}{jj}{kk}(ll), strain_rate_surf{ii}{jj}{kk}(ll)] ...
                                    = deal(tmp10(1), tmp10(2)); %#ok<AGROW>
                                if do_uncert
                                    [accum_shallow_uncert{ii}{jj}{kk}(:, ll), strain_rate_uncert{ii}{jj}{kk}(:, ll)] ...
                                        = strain_uncert(frac_test, conf_uncert, 'shallow_strain', thick_decim{ii}{jj}{kk}(ll), depth_smooth_decim(curr_layer, ll), age_trim(curr_layer), age_uncert_trim(curr_layer), tmp10, res_shallow{ii}{jj}{kk}(ll)); %#ok<AGROW,NASGU>
                                end
                            end
                        catch
                            ind_strain_fail{ii}{jj}{kk}(:, ll) ...
                                = [ll; 7];
                        end
                    end
                end
                ind_strain_fail{ii}{jj}{kk} ...
                            = ind_strain_fail{ii}{jj}{kk}(:, ~isnan(ind_strain_fail{ii}{jj}{kk}(1, :))); %#ok<AGROW>
            end
            
            if do_strain(2)
                shape_fact{ii}{jj}{kk} ...
                            = 1 - (thick_shear{ii}{jj}{kk} ./ (2 .* thick_decim{ii}{jj}{kk})); %#ok<AGROW,NASGU>
                if do_uncert
                    shape_fact_uncert{ii}{jj}{kk} ...
                            = 1 - (thick_shear_uncert{ii}{jj}{kk} ./ (2 .* repmat(thick_decim{ii}{jj}{kk}, 2, 1))); %#ok<AGROW,NASGU>
                end
            end
            if do_strain(4)
                shape_fact_melt{ii}{jj}{kk} ...
                            = 1 - (thick_shear_melt{ii}{jj}{kk} ./ (2 .* thick_decim{ii}{jj}{kk})); %#ok<AGROW,NASGU>
                if do_uncert
                    shape_fact_melt_uncert{ii}{jj}{kk} ...
                            = 1 - (thick_shear_melt_uncert{ii}{jj}{kk} ./ (2 .* repmat(thick_decim{ii}{jj}{kk}, 2, 1))); %#ok<AGROW,NASGU>
                end
            end
            if do_strain(5)
                shape_fact_lat{ii}{jj}{kk} ...
                            = 1 - (thick_shear_lat{ii}{jj}{kk} ./ (2 .* thick_decim{ii}{jj}{kk})); %#ok<AGROW,NASGU>
                if do_uncert
                    shape_fact_lat_uncert{ii}{jj}{kk} ...
                            = 1 - (thick_shear_lat_uncert{ii}{jj}{kk} ./ (2 .* repmat(thick_decim{ii}{jj}{kk}, 2, 1))); %#ok<AGROW,NASGU>
                end
            end
            if do_strain(6)
                shape_fact_surf{ii}{jj}{kk} ...
                            = 1 - (thick_shear_surf{ii}{jj}{kk} ./ (2 .* thick_decim{ii}{jj}{kk})); %#ok<AGROW,NASGU>
                if do_uncert
                    shape_fact_surf_uncert{ii}{jj}{kk} ...
                            = 1 - (thick_shear_surf_uncert{ii}{jj}{kk} ./ (2 .* repmat(thick_decim{ii}{jj}{kk}, 2, 1))); %#ok<AGROW,NASGU>
                end
            end
        end
    end
end

% NaN out poor strain-rate model results
if (strcmp(radar_type, 'deep') && do_bad)
    disp('Erasing results from known bad transect portions...')
    trans2nan               = [1  2   2 220  260;
                               4  7   5 30   130;
                               4  7   6 75   81;
                               4  7   7 0    180;
                               5  1   6 2715 2730;
                               5  9   4 2292 2296;
                               6  2   2 3181 3212;
                               6  3   1 2885 2890;
                               6  7   4 2428 2432;
                               6  8   1 5043 5044;
                               6  8   3 4219 4240;
                               6  11  4 3840 3861;
                               6  11  4 3885 4500;
                               9  1   1 1810 1820;
                               9  4   1 2280 2360;
                               11 6   1 1140 1261;
                               12 8   1 2700 2710;
                               17 104 1 2341 2345;
                               17 130 1 2110 2200;
                               18 3   2 830  865;
                               18 5   2 494  508;
                               19 19  1 1268 1330;
                               19 24  4 3750 3756;
                               19 27  1 3500 3668;
                               19 30  1 3750 3760;
                               19 30  1 3852 3856;
                               19 47  1 1010 1025;
                               19 51  2 3115 3130;
                               19 56  1 700  707;
                               20 4   2 2980 2985;
                               20 6   3 2200 2210;
                               20 7   1 20   45;
                               20 10  5 2295 2296;
                               20 13  2 1376 1444;
                               20 14  1 500  600;
                               20 15  2 1576 1579];
    load mat/merge_all dist
    for ii = 1:size(trans2nan, 1)
        ind2nan             = find((dist{trans2nan(ii, 1)}{trans2nan(ii, 2)}{trans2nan(ii, 3)}(ind_decim_mid{trans2nan(ii, 1)}{trans2nan(ii, 2)}{trans2nan(ii, 3)}) >= trans2nan(ii, 4)) & ...
                                   (dist{trans2nan(ii, 1)}{trans2nan(ii, 2)}{trans2nan(ii, 3)}(ind_decim_mid{trans2nan(ii, 1)}{trans2nan(ii, 2)}{trans2nan(ii, 3)}) <= trans2nan(ii, 5))); %#ok<NASGU>
        for jj = find(do_strain)
            for kk = 1:num_param(jj)
                eval([strain_param{jj}{kk} '{trans2nan(ii, 1)}{trans2nan(ii, 2)}{trans2nan(ii, 3)}(ind2nan) = NaN;'])
                if do_uncert
                    eval([strain_param{jj}{kk} '_uncert{trans2nan(ii, 1)}{trans2nan(ii, 2)}{trans2nan(ii, 3)}(:, ind2nan) = NaN;'])
                end
            end
        end
    end
    disp('Done erasing bad regions.')
end

% smooth along-transect values
if do_smooth
    disp('Smoothing strain-rate model parameters along-transect...')
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            for kk = 1:num_subtrans{ii}(jj)
                for ll = find(do_strain)
                    for mm = 1:num_param(ll)
                        if any(~isnan(eval([strain_param{ll}{mm} '{ii}{jj}{kk}'])))
                            eval([strain_param{ll}{mm} '{ii}{jj}{kk} = smooth_lowess(' strain_param_def{ll}{mm} '{ii}{jj}{kk}, num_smooth)'';'])
                        end
                        if do_uncert
                            if any(~isnan(eval([strain_param{ll}{mm} '_uncert{ii}{jj}{kk}(:)'])))
                                eval([strain_param{ll}{mm} '_uncert{ii}{jj}{kk}(1, :) = smooth_lowess(' strain_param_def{ll}{mm} '_uncert{ii}{jj}{kk}(1, :), num_smooth)'';'])
                                eval([strain_param{ll}{mm} '_uncert{ii}{jj}{kk}(2, :) = smooth_lowess(' strain_param_def{ll}{mm} '_uncert{ii}{jj}{kk}(2, :), num_smooth)'';'])
                            end
                        end
                    end
                end
            end
        end
    end
    disp('Done smoothing strain-rate model parameters.')
end

if do_grd
    
    load mat/greenland_bed_v3 x y
    [x_grd, y_grd]          = deal(x, y);
    clear x y
    [x_min, x_max, y_min, y_max] ...
                            = deal(-632, 846, -3344, -670);
    [x_tmp, y_tmp]          = deal(x_grd(((y_grd(:, 1) >= y_min) & (y_grd(:, 1) <= y_max)), ((x_grd(1, :) >= x_min) & (x_grd(1, :) <= x_max))), y_grd(((y_grd(:, 1) >= y_min) & (y_grd(:, 1) <= y_max)), ((x_grd(1, :) >= x_min) & (x_grd(1, :) <= x_max))));
    [x_grd, y_grd]          = deal(x_tmp, y_tmp); %#ok,ASGLU>
    
    for ii = find(do_strain)
        for jj = 1:num_param(ii)
            disp(['Gridding ' strain_param{ii}{jj} '...'])
            [x_all, y_all, param_all] ...
                            = deal([]);
            if do_uncert
                name_uncert = {'lo' 'hi'};
                param_all_uncert ...
                            = cell(1, 2);
            end
            for kk = 1:num_year
                for ll = 1:num_trans(kk)
                    for mm = 1:num_subtrans{kk}(ll)
                        if isempty(eval([strain_param{ii}{jj} '{kk}{ll}{mm}']))
                            continue
                        end
                        x_all ...
                            = [x_all; x_pk{kk}{ll}{mm}(ind_decim_mid{kk}{ll}{mm})']; %#ok<AGROW>
                        y_all ...
                            = [y_all; y_pk{kk}{ll}{mm}(ind_decim_mid{kk}{ll}{mm})']; %#ok<AGROW>
                        param_all ...
                            = [param_all; eval([strain_param{ii}{jj} '{kk}{ll}{mm}'''])]; %#ok<AGROW>
                        if do_uncert
                            for nn = 1:2
                                param_all_uncert{nn} ...
                                    = [param_all_uncert{nn}; eval([strain_param{ii}{jj} '_uncert{kk}{ll}{mm}(nn, :)'''])];
                            end
                        end
                    end
                end
            end
            [x_all_ref, y_all_ref] ...
                            = deal(x_all, y_all);
            ind_good        = find(~isnan(param_all));
            [x_all, y_all, param_all] ...
                            = deal(x_all_ref(ind_good), y_all_ref(ind_good), param_all(ind_good));
            [xy_all, ind_unique] ...
                            = unique([x_all y_all], 'rows');
            [x_all, y_all, param_all] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), param_all(ind_unique));
            try
                param_interpolant ...
                            = scatteredInterpolant(x_all, y_all, param_all, 'natural', 'none'); %#ok<NASGU>
                eval([strain_param{ii}{jj} '_grd = param_interpolant(x_grd, y_grd);'])
            catch
                disp([strain_param{ii}{jj} ' gridding failed...'])
            end
            if do_uncert
                for kk = 1:2
                    ind_good= find(~isnan(param_all_uncert{kk}));
                    [x_all_tmp, y_all_tmp, param_all_uncert{kk}] ...
                            = deal(x_all_ref(ind_good), y_all_ref(ind_good), param_all_uncert{kk}(ind_good));
                    [xy_all_tmp, ind_unique] ...
                            = unique([x_all_tmp y_all_tmp], 'rows');
                    [x_all_tmp, y_all_tmp, param_all_uncert{kk}] ...
                            = deal(xy_all_tmp(:, 1), xy_all_tmp(:, 2), param_all_uncert{kk}(ind_unique));
                    try
                        param_interpolant ...
                            = scatteredInterpolant(x_all_tmp, y_all_tmp, param_all_uncert{kk}, 'natural', 'none'); %#ok<NASGU>
                        eval([strain_param{ii}{jj} '_' name_uncert{kk} '_grd = param_interpolant(x_grd, y_grd);'])
                    catch
                        disp([strain_param{ii}{jj} '_' name_uncert{kk} ' gridding failed...'])
                    end
                end
            end
        end
    end
end

if do_save
    disp('Saving strain-rate modeling...')
    var2save                = '''age_max'',';
    for ii = 1:num_param_def
        var2save            = [var2save '''' strain_param_def{ii} ''',']; %#ok<AGROW>
    end
    for ii = find(do_strain)
        for jj = 1:num_param(ii)
            var2save        = [var2save '''' strain_param{ii}{jj} ''',']; %#ok<AGROW>
            if do_uncert
                var2save    = [var2save '''' strain_param{ii}{jj} '_uncert'',']; %#ok<AGROW>
            end
            if do_grd
                var2save    = [var2save '''' strain_param{ii}{jj} '_grd'',']; %#ok<AGROW>
                if do_uncert
                    var2save= [var2save '''' strain_param{ii}{jj} '_lo_grd'',''' strain_param{ii}{jj} '_hi_grd'',']; %#ok<AGROW>
                end
            end
        end
    end
    switch radar_type
        case 'accum'
            name_save       = ['mat/strain_all_' num2str(1e-3 * age_max) 'ka_accum.mat'];
        case 'deep'
            name_save       = ['mat/strain_all_' num2str(1e-3 * age_max) 'ka.mat'];
    end
    eval(['save(''' name_save ''', ''-v7.3'', ' var2save(1:(end - 1)) ')'])
    disp(['Saved strain-rate modeling in ' name_save '.'])
end

if parallel_check
    delete(pool)
end