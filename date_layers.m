% DATE_LAYERS Date layers using ice-core depth/age scales and 1D/2D interpolation/extrapolation or quasi-Nye dating of overlapping dated layers.
% 
% Joe MacGregor (UTIG)
% Last updated: 03/24/15

clear

% switches
radar_type                  = 'accum';

switch radar_type
    case 'accum'
        do_snr              = false;
        do_date             = true;
        do_age_check        = true;
        do_match            = [true true];
        do_interp           = true;
        interp_type         = 'quasi Nye';
        do_save             = true;
        do_grd1             = false;
        do_grd2             = false;
        do_nye_norm         = false;
    case 'deep'
        do_snr              = true;
        do_date             = true;
        do_age_check        = true;
        do_match            = [true true];
        do_interp           = true;
        interp_type         = 'quasi Nye';
        do_save             = true;
        do_grd1             = true;
        do_grd2             = true;
        do_nye_norm         = true;
end

% variables of interest
switch radar_type
    case 'accum'
        depth_uncert        = 5; % depth uncertainty representing unknown firn correction, m
        age_uncert_rel_max  = 0.1; % maximum relative age uncertainty
        dist_int_max        = 0.25; % km, +/- range to extract core-intersecting layer depths
        thick_diff_max      = 0.05; % fraction of ice thickness within which overlapping layer must lie
        layer_diff_max      = 1; % fraction of layer thickness within which overlapping layer must lie
        num_depth_norm      = 0; % number of normalized depths
        age_iso             = 1e3 .* [0.25 0.5 1 2]';
        num_age_iso         = length(age_iso); % number of isochrones to calculate
        thick_diff_max_iso  = 0.05; % fraction of ice thickness within which overlapping layer must lie (for age_norm1)
        layer_diff_max_iso  = 1; % fraction of layer thickness within which overlapping layer must lie (for age_norm1 and depth_iso1)
        age_diff_max_iso    = 5e2; % age range within which layer must exist (for depth_iso1)
        snr_ref             = 10 ^ (5 / 10); % linear ratio (5 is dB), reference SNR value if unknown
        depth_shallow_min   = 0.02; % minimum fraction of ice thickness to attempt shallow fitting
        depth_shallow_max   = 0.05; % maximum fraction of ice thickness to attempt shallow fitting
        num_date_loop_max   = 10; % maximum number of dating loops
        age_max             = 1e4; % maximum permitted reflector age
        year_ord            = [1 2 3]; % order in which to analyze years/campaigns based on their overall data quality
        num_gauss           = ones(1, 3); % number of Gaussians to use to fit Nye differences of normalized ages
        core_avail          = true(1, 6); % cores to include in dating, following order in core names
    case 'deep'
        depth_uncert        = 10; % depth uncertainty representing unknown firn correction, m
        age_uncert_rel_max  = 0.25; % maximum relative age uncertainty
        dist_int_max        = 0.25; % km, +/- range to extract core-intersecting layer depths
        thick_diff_max      = 0.2; % fraction of ice thickness within which overlapping layer must lie
        layer_diff_max      = 1; % fraction of layer thickness within which overlapping layer must lie
        num_depth_norm      = 25; % number of normalized depths
        age_iso             = 1e3 .* [9 11.7 29 57 115]';
        num_age_iso         = length(age_iso); % number of isochrones to calculate
        thick_diff_max_iso  = 0.2; % fraction of ice thickness within which overlapping layer must lie (for age_norm1)
        layer_diff_max_iso  = 1; % fraction of layer thickness within which overlapping layer must lie (for age_norm1 and depth_iso1)
        age_diff_max_iso    = 5e3; % age range within which layer must exist (for depth_iso1)
        snr_ref             = 10 ^ (5 / 10); % linear ratio (5 is dB), reference SNR value if unknown
        depth_shallow_min   = 0.03; % minimum fraction of ice thickness to attempt shallow fitting
        depth_shallow_max   = 0.2; % maximum fraction of ice thickness to attempt shallow fitting
        num_date_loop_max   = 10; % maximum number of dating loops
        age_max             = 1.5e5; % maximum permitted reflector age, a
        year_ord            = [17 19 20 6 7 8 4 9 5 2 3 12 1 13 15 16 18 11 10 14]; % order in which to analyze years/campaigns based on their overall data quality
        num_gauss           = [ones(1, 15) (2 .* ones(1, 8)) 1]; % number of Gaussians to use to fit Nye differences of normalized ages
        core_avail          = true(1, 6); % cores to include in dating, following order in core names
end

if strcmp(interp_type, 'quasi Nye')
    tol                     = 1e-2; % residual tolerance in m
    iter_max                = 100; % maximum number of iterations in search of best-fit strain rate
else
    [tol, iter_max]         = deal([]);
end

% load transect data and core intersection data
switch radar_type
    case 'accum'
        load mat/xy_all_accum num_year num_trans name_year name_trans
        load mat/core_int_accum int_core_merge num_core name_core_short
        load mat/id_layer_master_accum_stable id_layer_master_mat
        load mat/merge_all_accum depth_smooth dist ind_decim ind_decim_mid num_decim num_layer num_subtrans num_trace_tot thick thick_decim x_pk y_pk
        load mat/range_resolution_accum range_resolution
    case 'deep'
        load mat/xy_all num_year num_trans name_year name_trans
        load mat/core_int int_core_merge num_core name_core_short
        load mat/id_layer_master_stable id_layer_master_mat
        load mat/merge_all depth_smooth dist ind_decim ind_decim_mid num_decim num_layer num_subtrans num_trace_tot thick thick_decim x_pk y_pk
        load mat/range_resolution range_resolution
end
thick_trans                 = thick; % temporary switcheroo on name for along-transect ice thickness
load mat/greenland_mc_bed_1km thick x_grd y_grd % Mathieu's mass conservation bed
thick_grd                   = thick;
thick                       = thick_trans;
clear thick_trans

% do not permit matches between campaigns at far end of period from each other, i.e., ICORDS cannot be matched with MCoRDSv2 (rules established for MacGregor et al. [2015])
year_match                  = true(num_year);
switch radar_type
    case 'deep'
        year_match(1:4, 9:end) ...
                            = false;
        year_match(9:end, 1:4) ...
                            = false;
end

% start core pool if possible
if license('checkout', 'distrib_computing_toolbox')
    pool_check              = gcp('nocreate');
    if isempty(pool_check)
        try
            pool            = parpool('local', 4); % start 4 cores (no need to get crazy here)
        catch
            pool            = parpool('local');
        end
    end
    num_pool                = pool.NumWorkers; % number of workers the pool (will be less than 4 if less than that available)
    parallel_check          = true; % flag for parallelization
else
    num_pool                = 0;
    parallel_check          = false;
end

% load core depth-age scales
core                        = cell(1, num_core);
for ii = 1:num_core
    core{ii}                = load(['mat/depth_age/' name_core_short{ii} '_depth_age']);
    if isnan(core{ii}.age(end))
        [core{ii}.age, core{ii}.age_uncert, core{ii}.depth] ...
                            = deal(core{ii}.age(1:(end - 1)), core{ii}.age_uncert(1:(end - 1)), core{ii}.depth(1:(end - 1)));
    end
end

% index of NorthGRIP core, which is the backup core for age uncertainty
ind_ngrip                   = find(strcmp(name_core_short, 'ngrip'));

% dummy letters for sub-transect names
letters                     = 'a':'z';

% Greenland-centered, GIMP-projected 1-km grid limits
[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);

% load CISM accumulation data (i.e., RACMO2 from Ettema et al. [2009])
if (do_date || do_grd1 || (do_grd2 && do_nye_norm))
    load mat/greenland_cism accum x y
    [accum, x_cism, y_cism] = deal(double(accum), x, y);
    clear x y
    accum_grd               = accum(((y_cism(:, 1) >= y_min) & (y_cism(:, 1) <= y_max)), ((x_cism(1, :) >= x_min) & (x_cism(1, :) <= x_max)));
    clear accum x_cism y_cism
end

% layer signal-to-noise ratio
if do_snr
    switch radar_type
        case 'accum'
            load mat/snr_all_accum snr_all
        case 'deep'
            load mat/snr_all snr_all
    end
    % loop through transects to check for unavailable SNRs and assign them if necessary
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            if isempty(snr_all{ii}{jj})
                continue
            end
            for kk = 1:num_subtrans{ii}(jj)
                if isempty(snr_all{ii}{jj}{kk})
                    continue
                end
                snr_all{ii}{jj}{kk}(snr_all{ii}{jj}{kk} <= 0) ...
                            = NaN;
                snr_all{ii}{jj}{kk}(~isnan(snr_all{ii}{jj}{kk})) ...
                            = 10 .^ (snr_all{ii}{jj}{kk}(~isnan(snr_all{ii}{jj}{kk})) ./ 10);
            end
        end
    end
else % assign all SNRs to the reference value (snr_ref) otherwise
    snr_all                 = cell(1, num_year);
    for ii = 1:num_year
        snr_all{ii}         = deal(cell(1, num_trans(ii)));
        for jj = 1:num_trans(ii)
            snr_all{ii}{jj} = cell(1, num_subtrans{ii}(jj));
            for kk = 1:num_subtrans{ii}(jj)
                snr_all{ii}{jj}{kk} ...
                            = snr_ref .* ones(num_layer{ii}{jj}(kk), num_core);
            end
        end
    end
end

%%
if do_date
%%
    % initialize age-related cells
    date_counter            = 1;
    [age, age_core, age_match, age_ord, age_n, age_range, age_type, age_uncert] ...
                            = deal(cell(1, num_year));
    for ii = 1:num_year
        [age{ii}, age_core{ii}, age_match{ii}, age_ord{ii}, age_n{ii}, age_range{ii}, age_type{ii}, age_uncert{ii}] ...
                            = deal(cell(1, num_trans(ii)));
        for jj = 1:num_trans(ii)
            [age{ii}{jj}, age_core{ii}{jj}, age_match{ii}{jj}, age_ord{ii}{jj}, age_n{ii}{jj}, age_range{ii}{jj}, age_type{ii}{jj}, age_uncert{ii}{jj}] ...
                            = deal(cell(1, num_subtrans{ii}(jj)));
            for kk = 1:num_subtrans{ii}(jj)
                age_core{ii}{jj}{kk} ...
                            = NaN(num_layer{ii}{jj}(kk), num_core);
                [age{ii}{jj}{kk}, age_match{ii}{jj}{kk}, age_ord{ii}{jj}{kk}, age_n{ii}{jj}{kk}, age_range{ii}{jj}{kk}, age_type{ii}{jj}{kk}, age_uncert{ii}{jj}{kk}] ...
                            = deal(NaN(num_layer{ii}{jj}(kk), 1));
            end
        end
    end
    
    % clean up matching list by putting earlier year first
    for ii = 1:size(id_layer_master_mat, 1)
        if ((id_layer_master_mat(ii, 1) > id_layer_master_mat(ii, 5)) || ((id_layer_master_mat(ii, 1) == id_layer_master_mat(ii, 5)) && (id_layer_master_mat(ii, 2) > id_layer_master_mat(ii, 6))) || ...
           ((id_layer_master_mat(ii, 1) == id_layer_master_mat(ii, 5)) && (id_layer_master_mat(ii, 2) == id_layer_master_mat(ii, 6)) && (id_layer_master_mat(ii, 3) > id_layer_master_mat(ii, 7))) || ...
           ((id_layer_master_mat(ii, 1) == id_layer_master_mat(ii, 5)) && (id_layer_master_mat(ii, 2) == id_layer_master_mat(ii, 6)) && (id_layer_master_mat(ii, 3) == id_layer_master_mat(ii, 7)) && (id_layer_master_mat(ii, 4) > id_layer_master_mat(ii, 8))))
            id_layer_master_mat(ii, :) ...
                            = id_layer_master_mat(ii, [5:8 1:4]);
        end
    end
    id_layer_master_mat     = unique(id_layer_master_mat, 'rows');
    
    disp('Assigning master IDs within campaigns...')
    
    % intial master id
    master_id_curr          = 1;
    layer_bin               = {};
    
    % assign master ids within individual campaign first
    for ii = year_ord
        
        % all matches within current campaign
        match_curr          = id_layer_master_mat(ismember(id_layer_master_mat(:, [1 5]), [ii ii], 'rows'), :);
        match_curr_cell     = match_curr;
        match_curr_cell(~match_curr_cell) ...
                            = 1;
        
        % assign lowest possible master ID to match pair
        for jj = 1:size(match_curr, 1)
            
            % if no master ID assigned to either layer then make a new one
            if all(isnan([age_match{match_curr(jj, 1)}{match_curr(jj, 2)}{match_curr_cell(jj, 3)}(match_curr(jj, 4)) age_match{match_curr(jj, 5)}{match_curr(jj, 6)}{match_curr_cell(jj, 7)}(match_curr(jj, 8))]))
                
                % make new bin
                layer_bin{master_id_curr} ...
                            = [match_curr(jj, 1:4); match_curr(jj, 5:8)]; %#ok<SAGROW>
                
                % assign layer pair same ID
                [age_match{match_curr(jj, 1)}{match_curr(jj, 2)}{match_curr_cell(jj, 3)}(match_curr(jj, 4)), age_match{match_curr(jj, 5)}{match_curr(jj, 6)}{match_curr_cell(jj, 7)}(match_curr(jj, 8))] ...
                            = deal(master_id_curr);
                
                % increment ID
                master_id_curr ...
                            = master_id_curr + 1;
                
            else
                
                % current matching IDs
                age_match_curr ...
                            = sort([age_match{match_curr(jj, 1)}{match_curr(jj, 2)}{match_curr_cell(jj, 3)}(match_curr(jj, 4)) age_match{match_curr(jj, 5)}{match_curr(jj, 6)}{match_curr_cell(jj, 7)}(match_curr(jj, 8))]);
                
                % add matching layers to lower ID bin
                layer_bin{age_match_curr(1)} ...
                            = unique([layer_bin{age_match_curr(1)}; match_curr(jj, 1:4); match_curr(jj, 5:8)], 'rows'); %#ok<SAGROW>
                
                % add upper ID bin's layers to lower ID's bin
                if (~isnan(age_match_curr(2)) && (age_match_curr(1) ~= age_match_curr(2)))
                    layer_bin{age_match_curr(1)} ...
                            = unique([layer_bin{age_match_curr(1)}; layer_bin{age_match_curr(2)}], 'rows'); %#ok<SAGROW>
                    layer_bin{age_match_curr(2)} ...
                            = []; %#ok<SAGROW>
                end
                
                % assign all layers in lower bin to same ID
                for kk = 1:size(layer_bin{age_match_curr(1)}, 1)
                    age_match{layer_bin{age_match_curr(1)}(kk, 1)}{layer_bin{age_match_curr(1)}(kk, 2)}{max([1 layer_bin{age_match_curr(1)}(kk, 3)])}(layer_bin{age_match_curr(1)}(kk, 4)) ...
                            = age_match_curr(1);
                end
            end
        end
    end
    
    % bin size
    num_bin                 = zeros(1, length(layer_bin));
    for ii = 1:length(layer_bin)
        num_bin(ii)         = size(layer_bin{ii}, 1);
    end
    
    % bin descending order
    [~, ind_bin_ord]        = sort(num_bin, 'descend');
    
    % adjust age_match based on bin totals
    for ii = 1:length(layer_bin)
        for jj = 1:num_bin(ind_bin_ord(ii))
            age_match{layer_bin{ind_bin_ord(ii)}(jj, 1)}{layer_bin{ind_bin_ord(ii)}(jj, 2)}{max([1 layer_bin{ind_bin_ord(ii)}(jj, 3)])}(layer_bin{ind_bin_ord(ii)}(jj, 4)) ...
                            = ii;
        end
    end
    
    % trim bins
    layer_bin               = layer_bin(ind_bin_ord);
    layer_bin               = layer_bin(logical(num_bin(ind_bin_ord)));
    
    % reset new master id
    master_id_curr          = length(layer_bin) + 1;
    
    disp('Assigning master IDs between different campaigns...')
    
    % assign master ids for matches between different campaigns
    for ii = year_ord(1:(end - 1))
        
        disp(name_year{ii})
        
        for jj = year_ord((find(ii == year_ord) + 1):end)
            
            % skip if direct matching between campaign pair is not allowed
            if ~year_match(ii, jj)
                continue
            end
            
            % all matches for current campaign pair
            match_curr      = id_layer_master_mat((ismember(id_layer_master_mat(:, [1 5]), [ii jj], 'rows') | ismember(id_layer_master_mat(:, [1 5]), [jj ii], 'rows')), 1:8);
            match_curr_cell = match_curr;
            match_curr_cell(~match_curr_cell) ...
                            = 1;
            
            % assign lowest possible id to pair
            for kk = 1:size(match_curr, 1)
                if all(isnan([age_match{match_curr(kk, 1)}{match_curr(kk, 2)}{match_curr_cell(kk, 3)}(match_curr(kk, 4)) age_match{match_curr(kk, 5)}{match_curr(kk, 6)}{match_curr_cell(kk, 7)}(match_curr(kk, 8))]))
                    layer_bin{master_id_curr} ...
                            = [match_curr(kk, 1:4); match_curr(kk, 5:8)];
                    [age_match{match_curr(kk, 1)}{match_curr(kk, 2)}{match_curr_cell(kk, 3)}(match_curr(kk, 4)), age_match{match_curr(kk, 5)}{match_curr(kk, 6)}{match_curr_cell(kk, 7)}(match_curr(kk, 8))] ...
                            = deal(master_id_curr);
                    master_id_curr ...
                            = master_id_curr + 1;
                else
                    age_match_curr ...
                            = sort([age_match{match_curr(kk, 1)}{match_curr(kk, 2)}{match_curr_cell(kk, 3)}(match_curr(kk, 4)) age_match{match_curr(kk, 5)}{match_curr(kk, 6)}{match_curr_cell(kk, 7)}(match_curr(kk, 8))]);
                    layer_bin{age_match_curr(1)} ...
                            = unique([layer_bin{age_match_curr(1)}; match_curr(kk, 1:4); match_curr(kk, 5:8)], 'rows');
                    if (~isnan(age_match_curr(2)) && (age_match_curr(1) ~= age_match_curr(2)))
                        layer_bin{age_match_curr(1)} ...
                            = unique([layer_bin{age_match_curr(1)}; layer_bin{age_match_curr(2)}], 'rows');
                        layer_bin{age_match_curr(2)} ...
                            = [];
                    end
                    for ll = 1:size(layer_bin{age_match_curr(1)}, 1)
                        age_match{layer_bin{age_match_curr(1)}(ll, 1)}{layer_bin{age_match_curr(1)}(ll, 2)}{max([1 layer_bin{age_match_curr(1)}(ll, 3)])}(layer_bin{age_match_curr(1)}(ll, 4)) ...
                            = age_match_curr(1);
                    end
                end
            end
        end
    end
    
    % bin size
    num_bin                 = zeros(1, length(layer_bin));
    for ii = 1:length(layer_bin)
        num_bin(ii)         = size(layer_bin{ii}, 1);
    end
    
    % bin descending order
    [~, ind_bin_ord]        = sort(num_bin, 'descend');
    
    % adjust age_match based on bin totals
    for ii = 1:length(layer_bin)
        for jj = 1:num_bin(ind_bin_ord(ii))
            age_match{layer_bin{ind_bin_ord(ii)}(jj, 1)}{layer_bin{ind_bin_ord(ii)}(jj, 2)}{max([1 layer_bin{ind_bin_ord(ii)}(jj, 3)])}(layer_bin{ind_bin_ord(ii)}(jj, 4)) ...
                            = ii;
        end
    end
    
    % trim bins
    layer_bin               = layer_bin(ind_bin_ord);
    layer_bin               = layer_bin(logical(num_bin(ind_bin_ord)));
    num_bin                 = num_bin(ind_bin_ord);
    num_bin                 = num_bin(logical(num_bin)); % number of total bins
    
    % assign master ids to matches following current scheme
    id_layer_master_mat     = [id_layer_master_mat NaN(size(id_layer_master_mat, 1), 1)];
    for ii = 1:size(id_layer_master_mat, 1)
        id_curr             = id_layer_master_mat(ii, 1:8);
        id_curr(~id_curr)   = 1;
        id_layer_master_mat(ii, end) ...
                            = min([age_match{id_curr(1)}{id_curr(2)}{id_curr(3)}(id_curr(4)) age_match{id_curr(5)}{id_curr(6)}{id_curr(7)}(id_curr(8))]);
    end
    
    % remove matches with no master ID (i.e., retroactively not permitted matches)
    id_layer_master_mat     = id_layer_master_mat(~isnan(id_layer_master_mat(:, end)), :);
    
    % initialize master age list
    num_master              = length(layer_bin);
    age_master              = NaN(num_master, (num_core + 6));
    
%% 
    
    disp('Dating core-intersecting layers in core-intersecting transects...')
    
    % loop through years/campaigns
    for ii = year_ord
        
        disp(['Campaign: ' name_year{ii} '...'])
        
        % loop through merged picks files
        for jj = 1:num_trans(ii)
            
            % loop through all merged picks file for current transect
            for kk = 1:num_subtrans{ii}(jj)
                
                % skip if no core intersection or if no core intersection with available depth-age scale
                if isempty(ismember(int_core_merge(:, 1:3), [ii jj kk], 'rows'))
                    continue
                elseif isempty(find(core_avail(int_core_merge(ismember(int_core_merge(:, 1:3), [ii jj kk], 'rows'), 4)), 1))
                    continue
                end
                
                disp([name_trans{ii}{jj} letters(kk) '...'])
                
                % indices in int_core_merge of intersections with current transect
                ind_int_core= int_core_merge((ismember(int_core_merge(:, 1:3), [ii jj kk], 'rows') & core_avail(int_core_merge(:, 4))'), 4:5);
                
                % initialize uncertainty variables
                [age_uncert_depth, age_uncert_interp, age_uncert_radar] ...
                            = deal(NaN(num_layer{ii}{jj}(kk), num_core));
                
                % loop through good core intersections
                for ll = 1:size(ind_int_core, 1)
                    
                    disp([name_core_short{ind_int_core(ll, 1)} '...'])
                    
                    % unique values of distance vector
                    [dist_curr, ind_sort] ...
                            = unique(dist{ii}{jj}{kk});
                    
                    % smooth near-intersection layer depths if intersection index is unique
                    if ~isempty(find((ind_sort == ind_int_core(ll, 2)), 1))
                        depth_curr ...
                            = depth_smooth{ii}{jj}{kk}(:, ind_sort);
                        depth_curr ...
                            = mean(depth_curr(:, interp1(dist_curr, 1:length(dist_curr), (dist_curr(find((ind_sort == ind_int_core(ll, 2)), 1)) - dist_int_max), 'nearest', 'extrap'):...
                                                 interp1(dist_curr, 1:length(dist_curr), (dist_curr(find((ind_sort == ind_int_core(ll, 2)), 1)) + dist_int_max), 'nearest', 'extrap')), 2, 'omitnan');
                    else
                        depth_curr ...
                            = depth_smooth{ii}{jj}{kk}(:, ind_int_core(ll, 2));
                    end
                    
                    for mm = find(~isnan(depth_curr))' % loop through each of the merged picks file layers, only interpolate age if layer present at core intersection
                        
                        age_core_curr ...
                            = interp1(core{ind_int_core(ll, 1)}.depth(~isnan(core{ind_int_core(ll, 1)}.age)), core{ind_int_core(ll, 1)}.age(~isnan(core{ind_int_core(ll, 1)}.age)), depth_curr(mm), 'linear', 'extrap'); % interpolated age at core intersection
                        % adjust if transect intersects the same ice core multiple times
                        if isnan(age_core{ii}{jj}{kk}(mm, ind_int_core(ll, 1)))
                            age_core{ii}{jj}{kk}(mm, ind_int_core(ll, 1)) ...
                                = age_core_curr;
                        else
                            age_core{ii}{jj}{kk}(mm, ind_int_core(ll, 1)) ...
                                = mean([age_core_curr age_core{ii}{jj}{kk}(mm, ind_int_core(ll, 1))], 'omitnan');
                        end
                        
                        % get signal-to-noise ratio
                        if isnan(snr_all{ii}{jj}{kk}(mm, ind_int_core(ll, 1)))
                            snr_curr ...
                                = mean(snr_all{ii}{jj}{kk}(:, ind_int_core(ll, 1)), 'omitnan');
                            if isnan(snr_curr)
                                snr_curr ...
                                    = snr_ref;
                            end
                        else
                            snr_curr ...
                                = snr_all{ii}{jj}{kk}(mm, ind_int_core(ll, 1));
                        end
                        
                        % radar-induced uncertainty (range resolution / sqrt(snr))
                        age_uncert_radar_curr ...
                            = 0.5 * sum(abs(interp1(core{ind_int_core(ll, 1)}.depth(~isnan(core{ind_int_core(ll, 1)}.age)), core{ind_int_core(ll, 1)}.age(~isnan(core{ind_int_core(ll, 1)}.age)), ...
                                                    (depth_curr(mm) + ([range_resolution(ii) -range_resolution(ii)] ./ sqrt(snr_curr))), 'linear', 'extrap') - repmat(age_core{ii}{jj}{kk}(mm, ind_int_core(ll, 1)), 1, 2)));
                        if isnan(age_uncert_radar(mm, ind_int_core(ll, 1)))
                            age_uncert_radar(mm, ind_int_core(ll, 1)) ...
                                = age_uncert_radar_curr;
                        else
                            age_uncert_radar(mm, ind_int_core(ll, 1)) ...
                                = mean([age_uncert_radar_curr age_uncert_radar(mm, ind_int_core(ll, 1))], 'omitnan');
                        end
                        
                        % depth-induced uncertainty
                        age_uncert_depth_curr ...
                            = 0.5 * sum(abs(interp1(core{ind_int_core(ll, 1)}.depth(~isnan(core{ind_int_core(ll, 1)}.age)), core{ind_int_core(ll, 1)}.age(~isnan(core{ind_int_core(ll, 1)}.age)), (depth_curr(mm) + [depth_uncert -depth_uncert]), 'linear', 'extrap') - ...
                                            repmat(age_core{ii}{jj}{kk}(mm, ind_int_core(ll, 1)), 1, 2)));
                        if isnan(age_uncert_depth(mm, ind_int_core(ll, 1)))
                            age_uncert_depth(mm, ind_int_core(ll, 1)) ...
                                = age_uncert_depth_curr;
                        else
                            age_uncert_depth(mm, ind_int_core(ll, 1)) ...
                                = mean([age_uncert_depth_curr age_uncert_depth(mm, ind_int_core(ll, 1))], 'omitnan');
                        end
                        
                    end
                    
                    % core-reported age uncertainty
                    age_uncert_interp_curr ...
                            = interp1(core{ind_int_core(ll, 1)}.depth(~isnan(core{ind_int_core(ll, 1)}.age_uncert)), core{ind_int_core(ll, 1)}.age_uncert(~isnan(core{ind_int_core(ll, 1)}.age_uncert)), depth_curr(~isnan(depth_curr)), 'linear', 'extrap');
                    if isnan(age_uncert_radar(mm, ind_int_core(ll, 1)))
                        age_uncert_interp(~isnan(depth_curr), ind_int_core(ll, 1)) ...
                            = age_uncert_interp_curr;
                    else
                        age_uncert_interp(~isnan(depth_curr), ind_int_core(ll, 1)) ...
                            = mean([age_uncert_interp_curr age_uncert_interp(~isnan(depth_curr), ind_int_core(ll, 1))], 2, 'omitnan');
                    end
                    
                    % order in which layers were dated
                    age_ord{ii}{jj}{kk}(~isnan(depth_curr)) ...
                            = date_counter;
                    date_counter ...
                            = date_counter + 1;
                end
                
                % assign core-intersecting ages, uncertainties, range, number of layers used (1) and type (0 for layers at closest trace in core-intersecting transects)
                age{ii}{jj}{kk} ...
                            = mean(age_core{ii}{jj}{kk}, 2, 'omitnan');
                age_uncert{ii}{jj}{kk} ...
                            = sqrt(mean((age_uncert_interp .^ 2), 2, 'omitnan') + mean((age_uncert_radar .^ 2), 2, 'omitnan') + mean((age_uncert_depth .^ 2), 2, 'omitnan')); % RSS of core age uncertainty, range resolution and depth uncertainty
                age_n{ii}{jj}{kk}(~isnan(age{ii}{jj}{kk})) ...
                            = 1;
                age_range{ii}{jj}{kk} ...
                            = range(age_core{ii}{jj}{kk}, 2);
                age_type{ii}{jj}{kk}(~isnan(age{ii}{jj}{kk})) ...
                            = 0;
                disp([num2str(length(find(~isnan(age{ii}{jj}{kk})))) '/' num2str(num_layer{ii}{jj}(kk)) ' layers dated.'])
            end
        end
    end
    
%%
    
    if do_match(1)
        
        disp('Assigning ages to matched core-intersecting layers...')
        
        [age_master, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter] ...
                            = match_age(1:num_master, layer_bin, num_core, age_master, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter, true, 0, [0.5 1], do_age_check, name_trans, letters, depth_smooth);
        
        disp('...done assigning ages to matched core-intersecting layers.')
        
    end
    
%%  
    
    % determine how many layers are currently dated and which transects
    [num_layer_dated_all, num_layer_dated_old] ...
                            = deal(0);
    age_old                 = cell(1, num_year);
    for ii = 1:num_year
        age_old{ii}         = cell(1, num_trans(ii));
        for jj = 1:num_trans(ii)
            age_old{ii}{jj} = cell(1, num_subtrans{ii}(jj));
            for kk = 1:num_subtrans{ii}(jj)
                age_old{ii}{jj}{kk} ...
                            = NaN(num_layer{ii}{jj}(kk), 1);
                num_layer_dated_all ...
                            = num_layer_dated_all + length(find(~isnan(age{ii}{jj}{kk})));
            end
        end
    end
    
    % initialize dating loop number
    num_date_loop           = 1;
        
    % keep going until dating does not improve anymore
    while ((num_layer_dated_all > num_layer_dated_old) && (num_date_loop <= num_date_loop_max))
        
        % don't bother if dating interpolation not requested
        if ~do_interp
            break
        end
        
        disp(['Dating loop #' num2str(num_date_loop) '...'])
        
        % loop through years/campaigns
        for ii = year_ord
            
            if (ii == year_ord(1))
                disp('Dating overlapping layers in all transects...')
            end
            
            disp(['Campaign: ' name_year{ii} '...'])
            
            % loop through transects
            for jj = 1:num_trans(ii)
                
                % loop through all sub-transects for current transect
                for kk = 1:num_subtrans{ii}(jj)
                    
                    % move onto to next sub-transect if all or no layers were dated, or not enough dates, or no change from previous loop
                    if (~any(isnan(age{ii}{jj}{kk})) || all(isnan(age{ii}{jj}{kk})) || (length(find(~isnan(age{ii}{jj}{kk}))) < 2) || isempty(setdiff(find(~isnan(age{ii}{jj}{kk})), find(~isnan(age_old{ii}{jj}{kk})))))
                        continue
                    else
                        if (num_subtrans{ii}(jj) == 1)
                            curr_name ...
                                = name_trans{ii}{jj};
                        else
                            curr_name ...
                                = [name_trans{ii}{jj} letters(kk)];
                        end
                        disp([curr_name '...'])
                    end
                    
                    % Nye strain rate
                    if strcmp(interp_type, 'quasi Nye')
                        accum_curr ...
                                = interp2(x_grd, y_grd, accum_grd, x_pk{ii}{jj}{kk}, y_pk{ii}{jj}{kk}, 'linear', NaN);
                        strain_rate_curr ...
                                = accum_curr ./ thick{ii}{jj}{kk}; % start with Nye strain rate as initial guess
                        ind_nothick ...
                                = find(isnan(strain_rate_curr) | isinf(strain_rate_curr));
                        if ~isempty(ind_nothick)
                            strain_rate_curr(ind_nothick) ...
                                = accum_curr(ind_nothick) ./ interp2(x_grd, y_grd, thick_grd, x_pk{ii}{jj}{kk}(ind_nothick), y_pk{ii}{jj}{kk}(ind_nothick), 'linear', NaN);
                        end
                    else
                        strain_rate_curr ...
                                = [];
                    end
                    
                    % remember ages before attempting dating this time
                    age_old{ii}{jj}{kk} ...
                            = age{ii}{jj}{kk};
                    
                    % initial list of undated reflectors
                    [ind_layer_undated, ind_layer_undated_ref] ...
                            = deal(find(isnan(age{ii}{jj}{kk}))');
                    
                    % dummy list of undated reflectors to start first while loop
                    ind_layer_undated_old ...
                            = 1:num_layer{ii}{jj}(kk);
                    
                    disp([num2str(length(ind_layer_undated)) '/' num2str(num_layer{ii}{jj}(kk)) ' undated...'])
                    
                    % repeat dating until no more progress is made
                    while (length(ind_layer_undated) < length(ind_layer_undated_old))
                        
                        ind_layer_blacklist ...
                            = [];
                        
                        % progressively assign ages to layers in core-intersecting transects that did not intersect an ice core
                        while any(isnan(age{ii}{jj}{kk}(setdiff(ind_layer_undated, ind_layer_blacklist))))
                            [age{ii}{jj}{kk}, age_ord{ii}{jj}{kk}, age_n{ii}{jj}{kk}, age_range{ii}{jj}{kk}, age_type{ii}{jj}{kk}, age_uncert{ii}{jj}{kk}, date_counter, ind_layer_blacklist] ...
                                = date_interp(depth_smooth{ii}{jj}{kk}, thick{ii}{jj}{kk}, age{ii}{jj}{kk}, age_ord{ii}{jj}{kk}, age_n{ii}{jj}{kk}, age_range{ii}{jj}{kk}, age_type{ii}{jj}{kk}, age_uncert{ii}{jj}{kk}, num_layer{ii}{jj}(kk), ind_layer_undated, ind_layer_blacklist, ...
                                              thick_diff_max, layer_diff_max, curr_name, interp_type, strain_rate_curr, parallel_check, num_pool, tol, iter_max, age_max, do_age_check, age_uncert_rel_max, date_counter, (num_date_loop + 1));
                        end
                        
                        % updated undated layer set
                        ind_layer_undated_old ...
                            = ind_layer_undated;
                        ind_layer_undated ...
                            = find(isnan(age{ii}{jj}{kk}))';
                        
                    end
                    
                    % shallow quasi-Nye using surface as a last resort, if possible
                    for ll = find(isnan(age{ii}{jj}{kk}))'
                        [age{ii}{jj}{kk}, age_ord{ii}{jj}{kk}, age_n{ii}{jj}{kk}, age_range{ii}{jj}{kk}, age_type{ii}{jj}{kk}, age_uncert{ii}{jj}{kk}, date_counter] ...
                            = date_interp_shallow(depth_smooth{ii}{jj}{kk}, thick{ii}{jj}{kk}, age{ii}{jj}{kk}, age_ord{ii}{jj}{kk}, age_n{ii}{jj}{kk}, age_range{ii}{jj}{kk}, age_type{ii}{jj}{kk}, age_uncert{ii}{jj}{kk}, interp_type, date_counter, ll, age_uncert_rel_max, depth_shallow_min, ...
                                                  depth_shallow_max, age_max, do_age_check, curr_name, (num_date_loop + 1.25));
                    end
                    
                    % display number of additional layers that were dated
                    disp([num2str(length(ind_layer_undated_ref) - length(ind_layer_undated)) ' dated...'])
                end
            end
        end
        
        disp('...done dating overlapping layers.')
        
        if do_match(2)
            
            disp('Assigning ages to matched layers...')
            
            [age_master, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter] ...
                            = match_age(find(isnan(age_master(:, (num_core + 1))))', layer_bin, num_core, age_master, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter, false, (num_date_loop + 1.25), (num_date_loop + [1.5 1.75]), do_age_check, name_trans, letters, ...
                                        depth_smooth);
            
            disp('...done assigning ages to matched layers.')
            
        end
        
        % reassign old number of dated layers
        num_layer_dated_old = num_layer_dated_all;
        
        % new number of dated layers
        num_layer_dated_all = 0;
        for ii = 1:num_year
            for jj = 1:num_trans(ii)
                for kk = 1:num_subtrans{ii}(jj)
                    num_layer_dated_all ...
                            = num_layer_dated_all + length(find(~isnan(age{ii}{jj}{kk})));
                end
            end
        end
        
        disp([num2str(num_layer_dated_all - num_layer_dated_old) ' dated layers added in this dating loop...'])
        if ((num_layer_dated_all - num_layer_dated_old) <= 0)
            disp('Dating loop shutting down...')
        end
        
        % increment dating loop counter
        num_date_loop       = num_date_loop + 1;
        
    end
    
%%    
    if do_save
        disp('Saving layer ages...')
        switch radar_type
            case 'accum'
                save('mat/date_all_accum', '-v7.3', 'age', 'age_core', 'age_master', 'age_match', 'age_ord', 'age_n', 'age_range', 'age_type', 'age_uncert', 'age_uncert_rel_max', 'dist_int_max', 'layer_bin')
                disp('Layer ages saved as mat/date_all_accum.mat.')
            case 'deep'
                save('mat/date_all', '-v7.3', 'age', 'age_core', 'age_master', 'age_match', 'age_ord', 'age_n', 'age_range', 'age_type', 'age_uncert', 'age_uncert_rel_max', 'dist_int_max', 'layer_bin')
%                 save('mat/date_all_ngrip_qnye', '-v7.3', 'age', 'age_core', 'age_uncert', 'age_type', 'dist_int_max')
                disp('Layer ages saved as mat/date_all.mat.')
        end
    end
    
%%
elseif (do_grd1 || do_grd2)
%%
    disp('Loading layer ages...')
    switch radar_type
        case 'accum'
            load mat/date_all_accum age age_uncert
            disp('Loaded layer ages from mat/date_all_accum.mat')
        case 'deep'
            load mat/date_all age age_uncert
            disp('Loaded layer ages from mat/date_all.mat')
    end
%%
end

%%
if do_grd1
    
    disp('1-D inter/extrapolation of age at normalized depths and depth at specified ages...')
    
    % decimation
    if logical(num_depth_norm)
        depth_norm          = ((1 / num_depth_norm):(1 / num_depth_norm):1)'; % normalized/relative depth vector
    else
        depth_norm          = 0;
    end
    
    % initialize along-track variables
    [age_diff1, age_norm1, age_uncert_norm1, depth_iso1, depth_norm_grd1, depth_uncert_iso1, dist_grd1] ...
                            = deal(cell(1, num_year));
    
    for ii = 1:num_year
        
        [age_diff1{ii}, age_norm1{ii}, age_uncert_norm1{ii}, depth_iso1{ii}, depth_norm_grd1{ii}, depth_uncert_iso1{ii}, dist_grd1{ii}] ...
                            = deal(cell(1, num_trans(ii)));
        
        disp(['Campaign: ' name_year{ii} '...'])
        
        % loop through merged picks file that have been analyzed through FENCEGUI
        for jj = 1:num_trans(ii)
            
            [age_diff1{ii}{jj}, age_norm1{ii}{jj}, age_uncert_norm1{ii}{jj}, depth_iso1{ii}{jj}, depth_norm_grd1{ii}{jj}, depth_uncert_iso1{ii}{jj}, dist_grd1{ii}{jj}] ...
                            = deal(cell(1, num_subtrans{ii}(jj)));
            
            % loop through all merged picks file for current transect
            for kk = 1:num_subtrans{ii}(jj)
                
                % dated layer indices
                ind_layer_dated ...
                            = find(~isnan(age{ii}{jj}{kk}) & ~isnan(age_uncert{ii}{jj}{kk}));
                
                % not enough dated layers in this transect
                if (length(ind_layer_dated) < 2)
                    continue
                else
                    disp(['Transect: ' name_trans{ii}{jj} letters(kk) '...'])
                end
                
                % decimated ice thickness and layer depths
                dist_decim  = NaN(1, num_decim{ii}{jj}(kk));
                depth_smooth_decim ...
                            = NaN(length(ind_layer_dated), num_decim{ii}{jj}(kk));
                for ll = 1:num_decim{ii}{jj}(kk)
                    depth_smooth_decim(:, ll) ...
                            = mean(depth_smooth{ii}{jj}{kk}(ind_layer_dated, ind_decim{ii}{jj}{kk}(ll):ind_decim{ii}{jj}{kk}(ll + 1)), 2, 'omitnan');
                    dist_decim(ll) ...
                            = mean(dist{ii}{jj}{kk}(ind_decim{ii}{jj}{kk}(ll):ind_decim{ii}{jj}{kk}(ll + 1)), 'omitnan');
                end
                
                % decimated Nye strain rate
                if strcmp(interp_type, 'quasi Nye')
                    accum_curr ...
                            = interp2(x_grd, y_grd, accum_grd, x_pk{ii}{jj}{kk}, y_pk{ii}{jj}{kk}, 'linear', NaN);
                    accum_decim ...
                            = NaN(1, num_decim{ii}{jj}(kk));
                    for ll = 1:num_decim{ii}{jj}(kk)
                        accum_decim(ll) ...
                            = mean(accum_curr(ind_decim{ii}{jj}{kk}(ll):ind_decim{ii}{jj}{kk}(ll + 1)), 'omitnan');
                    end
                    strain_rate_decim ...
                            = accum_decim ./ thick_decim{ii}{jj}{kk}; % start with Nye strain rate as initial guess
                    strain_rate_decim(isnan(strain_rate_decim) | isinf(strain_rate_decim)) ...
                            = mean(strain_rate_decim(~isinf(strain_rate_decim)), 'omitnan');
                end
                
                % thickness-normalized layer depths
                depth_smooth_norm ...
                            = depth_smooth_decim ./ repmat(thick_decim{ii}{jj}{kk}, length(ind_layer_dated), 1);
                
                % transect-specific grid dimensions
                [dist_grd1{ii}{jj}{kk}, depth_norm_grd1{ii}{jj}{kk}] ...
                            = meshgrid(dist_decim, depth_norm);
                
                % dated layers and uncertainties
                [age_dated, age_uncert_dated] ...
                            = deal(age{ii}{jj}{kk}(ind_layer_dated), age_uncert{ii}{jj}{kk}(ind_layer_dated)); % only use dated layers
                
                [age_norm1{ii}{jj}{kk}, age_uncert_norm1{ii}{jj}{kk}] ...
                            = deal(NaN(num_depth_norm, num_decim{ii}{jj}(kk)));
                [depth_iso1{ii}{jj}{kk}, depth_uncert_iso1{ii}{jj}{kk}] ...
                            = deal(NaN(num_age_iso, num_decim{ii}{jj}(kk)));
                
                % do age and depth uncertainties first
                for ll = find(~isnan(thick_decim{ii}{jj}{kk}))
                    if ~thick_decim{ii}{jj}{kk}(ll)
                        continue
                    end
                    if (length(find(~isnan(depth_smooth_norm(:, ll)))) < 2)
                        continue
                    end
                    [tmp1, tmp2] ...
                            = unique(depth_smooth_norm(~isnan(depth_smooth_norm(:, ll)), ll));
                    age_uncert_curr ...
                            = age_uncert_dated(~isnan(depth_smooth_norm(:, ll)));
                    age_uncert_norm1{ii}{jj}{kk}(:, ll) ...
                            = interp1(tmp1, age_uncert_curr(tmp2), depth_norm, 'linear', 'extrap');
                    [tmp1, tmp2] ...
                            = unique(age_dated);
                    depth_uncert_curr ...
                            = mean(diff(core{ind_ngrip}.depth(interp1(core{ind_ngrip}.age, 1:length(core{ind_ngrip}.age), [(age_dated - age_uncert_dated) age_dated (age_dated + age_uncert_dated)], 'nearest', 'extrap')), 1, 2), 2, 'omitnan');
                    depth_uncert_iso1{ii}{jj}{kk}(:, ll) ...
                            = interp1(tmp1, depth_uncert_curr(tmp2), age_iso, 'linear', 'extrap');
               end
                
                % loop through each decimated trace and calculate 1D vertical age profile at uniform horizontal intervals in along-transect grid
                switch interp_type
                    
                    case 'linear'
                        
                        for ll = find(~isnan(thick_decim{ii}{jj}{kk}))
                            
                            % skip if zero thickness or not enough dated layers
                            if ~thick_decim{ii}{jj}{kk}(ll)
                                continue
                            end
                            if (length(find(~isnan(depth_smooth_norm(:, ll)))) < 2)
                                continue
                            end
                            
                            [tmp1, tmp2] ...
                                = unique(depth_smooth_norm(~isnan(depth_smooth_norm(:, ll)), ll));
                            
                            % linearly interpolate age at normalized depths
                            if logical(num_depth_norm)
                                age_curr ...
                                    = age_dated(~isnan(depth_smooth_norm(:, ll)));
                                age_norm1{ii}{jj}{kk}(:, ll) ...
                                    = interp1(tmp1, age_curr(tmp2), depth_norm, 'linear', 'extrap');
                            end
                            
                            % linearly interpolate depth at specified ages
                            [tmp1, tmp2] ...
                                = unique(age_dated(~isnan(depth_smooth_decim(:, ll))));
                            if (length(tmp1) > 1)
                                depth_smooth_decim_curr ...
                                    = depth_smooth_decim(~isnan(depth_smooth_decim(:, ll)), ll);
                                depth_iso1{ii}{jj}{kk}(:, ll) ...
                                    = interp1(tmp1, depth_smooth_decim_curr(tmp2), age_iso, 'linear', 'extrap');
                            end
                            
                            depth_iso1{ii}{jj}{kk}((depth_iso1{ii}{jj}{kk}(:, ll) > thick_decim{ii}{jj}{kk}(ll)), ll) ...
                                    = NaN;
                        end
                        
                    case 'quasi Nye'
                        
                        for ll = find(~isnan(thick_decim{ii}{jj}{kk}))
                            
                            if (length(find(~isnan(depth_smooth_norm(:, ll)))) < 2)
                                continue
                            end
                            if ~thick_decim{ii}{jj}{kk}(ll)
                                continue
                            end
                            
                            % current normalized layer depths
                            [depth_curr, thick_diff_max_curr, age_curr, depth_smooth_decim_curr] ...
                                = deal((depth_norm .* thick_decim{ii}{jj}{kk}(ll)), (thick_diff_max_iso * thick_decim{ii}{jj}{kk}(ll)), age_dated(~isnan(depth_smooth_decim(:, ll))), depth_smooth_decim(~isnan(depth_smooth_decim(:, ll)), ll));
                            [depth_smooth_decim_curr, ind_depth_ord] ...
                                = sort(depth_smooth_decim_curr);
                            age_curr ...
                                = age_curr(ind_depth_ord);
                            
                            for mm = 1:num_depth_norm
                                
                                [age_bound, depth_bound] ...
                                    = deal([]);
                                
                                [ind_top, ind_bot] ...
                                    = deal(find((depth_smooth_decim_curr < depth_curr(mm)), 2, 'last'), find((depth_smooth_decim_curr > depth_curr(mm)), 2));
                                
                                if (isempty(ind_top) && isempty(ind_bot))
                                    continue
                                end
                                
                                if (isempty(ind_top) && (length(ind_bot) == 2))
                                    if all((depth_smooth_decim_curr(ind_bot(1)) - depth_curr(mm)) <= [(layer_diff_max_iso * diff(depth_smooth_decim_curr(ind_bot))) thick_diff_max_curr])
                                        [age_bound, depth_bound] ...
                                            = deal(age_curr(ind_bot), depth_smooth_decim_curr(ind_bot));
                                    end
                                elseif (isempty(ind_bot) && (length(ind_top) == 2))
                                    if all((depth_curr(mm) - depth_smooth_decim_curr(ind_top(2))) <= [(layer_diff_max_iso * diff(depth_smooth_decim_curr(ind_top))) thick_diff_max_curr])
                                        [age_bound, depth_bound] ...
                                            = deal(age_curr(ind_top), depth_smooth_decim_curr(ind_top));
                                    end
                                elseif (~isempty(ind_top) && ~isempty(ind_bot))
                                    [age_bound, depth_bound] ...
                                        = deal(age_curr([ind_top(end); ind_bot(1)]), depth_smooth_decim_curr([ind_top(end); ind_bot(1)]));
                                end
                                
                                % quasi-Nye age
                                if (~isempty(depth_bound) && ~isempty(age_bound) && (diff(depth_bound) > 0) && (diff(age_bound) > 0))
                                    age_norm1{ii}{jj}{kk}(mm, ll) ...
                                        = date_quasi_nye(depth_bound, age_bound, (-diff(log(1 - (depth_bound ./ thick_decim{ii}{jj}{kk}(ll)))) / diff(age_bound)), depth_curr(mm), tol, iter_max, 'age');
                                    % sanity check
                                    if (((depth_curr(mm) > depth_bound(2)) && (age_norm1{ii}{jj}{kk}(mm, ll) < age_bound(2))) || ((depth_curr(mm) < depth_bound(1)) && (age_norm1{ii}{jj}{kk}(mm, ll) > age_bound(1))) || ...
                                        (((depth_curr(mm) > depth_bound(1)) && (depth_curr(mm) < depth_bound(2))) && ((age_norm1{ii}{jj}{kk}(mm, ll) < age_bound(1)) || (age_norm1{ii}{jj}{kk}(mm, ll) > age_bound(2)))))
                                        age_norm1{ii}{jj}{kk}(mm, ll) ...
                                            = NaN;
                                    end
                                end
                            end
                            
                            for mm = 1:num_age_iso
                                
                                [age_bound, depth_bound] ...
                                    = deal([]);
                                
                                [ind_top, ind_bot] ...
                                    = deal(find((age_curr < age_iso(mm)), 2, 'last'), find((age_curr > age_iso(mm)), 2));
                                
                                if (isempty(ind_top) && isempty(ind_bot))
                                    continue
                                end
                                
                                if (isempty(ind_top) && (length(ind_bot) == 2))
                                    if all((age_curr(ind_bot(1)) - age_iso(mm)) <= [(layer_diff_max_iso * diff(age_curr(ind_bot))) age_diff_max_iso])
                                        [age_bound, depth_bound] ...
                                            = deal(age_curr(ind_bot), depth_smooth_decim_curr(ind_bot));
                                    end
                                elseif (isempty(ind_bot) && (length(ind_top) == 2))
                                    if all((age_iso(mm) - age_curr(ind_top(2))) <= [(layer_diff_max_iso * diff(age_curr(ind_top))) age_diff_max_iso])
                                        [age_bound, depth_bound] ...
                                            = deal(age_curr(ind_top), depth_smooth_decim_curr(ind_top));
                                    end
                                elseif (~isempty(ind_top) && ~isempty(ind_bot))
                                    [age_bound, depth_bound] ...
                                        = deal(age_curr([ind_top(end); ind_bot(1)]), depth_smooth_decim_curr([ind_top(end); ind_bot(1)]));
                                end
                                
                                % quasi-Nye isochrone depth
                                if (~isempty(depth_bound) && ~isempty(age_bound) && (diff(depth_bound) > 0) && (diff(age_bound) > 0))
                                    depth_iso1{ii}{jj}{kk}(mm, ll) ...
                                        = date_quasi_nye(depth_bound, age_bound, (-diff(log(1 - (depth_bound ./ thick_decim{ii}{jj}{kk}(ll)))) / diff(age_bound)), age_iso(mm), tol, iter_max, 'depth');
                                    if (((age_iso(mm) > age_bound(2)) && (depth_iso1{ii}{jj}{kk}(mm, ll) < depth_bound(2))) || ((age_iso(mm) < age_bound(1)) && (depth_iso1{ii}{jj}{kk}(mm, ll) > depth_bound(1))) || ...
                                        (((age_iso(mm) > age_bound(1)) && (age_iso(mm) < age_bound(2))) && ((depth_iso1{ii}{jj}{kk}(mm, ll) < depth_bound(1)) || (depth_iso1{ii}{jj}{kk}(mm, ll) > depth_bound(2)))))
                                        depth_iso1{ii}{jj}{kk}(mm, ll) ...
                                            = NaN;
                                    end
                                end
                            end
                            depth_iso1{ii}{jj}{kk}((depth_iso1{ii}{jj}{kk}(:, ll) > thick_decim{ii}{jj}{kk}(ll)), ll) ...
                                = NaN;
                        end
                end
                
                % NaN out misbehaving values
                age_norm1{ii}{jj}{kk}(isinf(age_norm1{ii}{jj}{kk}) | (age_norm1{ii}{jj}{kk} < 0) | (age_norm1{ii}{jj}{kk} > age_max)) ...
                            = NaN;
                age_uncert_norm1{ii}{jj}{kk}(isinf(age_uncert_norm1{ii}{jj}{kk}) | (age_uncert_norm1{ii}{jj}{kk} < 0)) ...
                            = NaN;
                depth_iso1{ii}{jj}{kk}(isinf(depth_iso1{ii}{jj}{kk}) | (depth_iso1{ii}{jj}{kk} < 0)) ...
                            = NaN;
                depth_uncert_iso1{ii}{jj}{kk}(depth_uncert_iso1{ii}{jj}{kk} < 0) ...
                            = NaN;
                
                % reference Nye age field
                if logical(num_depth_norm)
                    age_nye1= [-(1 ./ strain_rate_decim(ones((num_depth_norm - 1), 1), :)) .* log(1 - depth_norm(1:(end - 1), ones(1, num_decim{ii}{jj}(kk)))); NaN(1, num_decim{ii}{jj}(kk))];
                    age_diff1{ii}{jj}{kk} ...
                            = (age_norm1{ii}{jj}{kk} - age_nye1) ./ age_nye1;
                else
                    age_diff1{ii}{jj}{kk} ...
                            = NaN(0, num_decim{ii}{jj}(kk));
                end
            end
        end
    end
    
    disp('...done 1-D gridding layer ages.')
%%
    if do_save
        disp('Saving 1-D age grids...')
        switch radar_type
            case 'accum'
                save('mat/age_grd1_accum', '-v7.3', 'age_diff1', 'age_norm1', 'age_uncert_norm1', 'age_iso', 'depth_norm', 'depth_iso1', 'depth_norm', 'depth_norm_grd1', 'depth_uncert_iso1', 'dist_grd1', 'num_age_iso', 'num_depth_norm')
                disp('Saved 1-D age grids as mat/age_grd1_accum.mat.')
            case 'deep'
                save('mat/age_grd1', '-v7.3', 'age_diff1', 'age_norm1', 'age_uncert_norm1', 'age_iso', 'depth_norm', 'depth_iso1', 'depth_norm', 'depth_norm_grd1', 'depth_uncert_iso1', 'dist_grd1', 'num_age_iso', 'num_depth_norm')
                disp('Saved 1-D age grids as mat/age_grd1.mat.')
        end
    end
%%
elseif do_grd2
%%
    disp('Loading 1-D age gridding...')
    switch radar_type
        case 'accum'
            load mat/age_grd1_accum age_norm1 age_uncert_norm1 depth_iso1 depth_norm depth_uncert_iso1 num_age_iso num_depth_norm
            if do_nye_norm
                load mat/age_grd1_accum age_diff1
            end
            disp('Loaded 1-D age gridding from mat/age_grd1_accum.mat')
        case 'deep'
            load mat/age_grd1 age_norm1 age_uncert_norm1 depth_iso1 depth_norm depth_uncert_iso1 num_age_iso num_depth_norm
            if do_nye_norm
                load mat/age_grd1 age_diff1
            end
            disp('Loaded 1-D age gridding from mat/age_grd1.mat')
    end
%%
end

if do_grd2
    
    disp('2-D gridding of normalized ages, isochrone depths and their uncertainties...')
    
    % prepare cells for ages/depths to be gridded
    [age_norm2, age_uncert_norm2] ...
                            = deal(NaN(size(x_grd, 1), size(x_grd, 2), num_depth_norm));
    [depth_iso2, depth_uncert_iso2] ...
                            = deal(NaN(size(x_grd, 1), size(x_grd, 2), num_age_iso));
    [thick_cat, x_all_age_cat, x_all_age_uncert_cat, x_all_depth_cat, x_all_depth_uncert_cat, y_all_age_cat, y_all_age_uncert_cat, y_all_depth_cat, y_all_depth_uncert_cat] ...
                            = deal([]);
    [age_all, age_all_min, age_all_max, age_uncert_all, x_all_age, x_all_age_min, x_all_age_max, x_all_age_uncert, y_all_age, y_all_age_min, y_all_age_max, y_all_age_uncert] ...
                            = deal(cell(1, num_depth_norm));
    if do_nye_norm
        age_diff_all        = cell(1, num_depth_norm);
    end
    [depth_all, depth_all_min, depth_all_max, depth_uncert_all, x_all_depth, x_all_depth_min, x_all_depth_max, x_all_depth_uncert, y_all_depth, y_all_depth_min, y_all_depth_max, y_all_depth_uncert] ...
                            = deal(cell(1, num_age_iso));
    
    % loop through each 1-D set and add its values to the mix
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            for kk = 1:num_subtrans{ii}(jj)
                if ~isempty(age_norm1{ii}{jj}{kk})
                    [x_all_age_cat, y_all_age_cat] ...
                            = deal([x_all_age_cat; x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})'], [y_all_age_cat; y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})']);
                    for ll = 1:num_depth_norm
                        age_all{ll} ...
                            = [age_all{ll}; age_norm1{ii}{jj}{kk}(ll, :)'];
                        if do_nye_norm
                            age_diff_all{ll} ...
                                = [age_diff_all{ll}; age_diff1{ii}{jj}{kk}(ll, :)'];
                        end
                        [age_min_tmp, age_max_tmp] ...
                            = deal(age_norm1{ii}{jj}{kk}(ll, :));
                        for mm = find(isnan(age_norm1{ii}{jj}{kk}(ll, :)))
                            if ~isempty(find(age_norm1{ii}{jj}{kk}(1:(ll - 1), mm), 1, 'last'))
                                age_min_tmp(ll) ...
                                    = age_norm1{ii}{jj}{kk}(find(age_norm1{ii}{jj}{kk}(1:(ll - 1), mm), 1, 'last'), mm);
                            end
                            if ~isempty(find(age_norm1{ii}{jj}{kk}((ll + 1):end, mm), 1))
                                age_max_tmp(ll) ...
                                    = age_norm1{ii}{jj}{kk}(find(age_norm1{ii}{jj}{kk}((ll + 1):end, mm), 1), mm);
                            end
                        end
                        [age_all_min{ll}, age_all_max{ll}] ...
                            = deal([age_all_min{ll}; age_min_tmp'], [age_all_max{ll}; age_max_tmp']);
                    end
                end
                if ~isempty(age_uncert_norm1{ii}{jj}{kk})
                    [x_all_age_uncert_cat, y_all_age_uncert_cat] ...
                            = deal([x_all_age_uncert_cat; x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})'], [y_all_age_uncert_cat; y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})']);
                    for ll = 1:num_depth_norm
                        age_uncert_all{ll} ...
                            = [age_uncert_all{ll}; age_uncert_norm1{ii}{jj}{kk}(ll, :)'];
                    end
                end
                if ~isempty(depth_iso1{ii}{jj}{kk})
                    [x_all_depth_cat, y_all_depth_cat, thick_cat] ...
                            = deal([x_all_depth_cat; x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})'], [y_all_depth_cat; y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})'], [thick_cat; double(thick_decim{ii}{jj}{kk}')]);
                    for ll = 1:num_age_iso
                        depth_all{ll} ...
                            = [depth_all{ll}; depth_iso1{ii}{jj}{kk}(ll, :)'];
                        [depth_min_tmp, depth_max_tmp] ...
                                = deal(depth_iso1{ii}{jj}{kk}(ll, :));
                        for mm = find(isnan(depth_iso1{ii}{jj}{kk}(ll, :)))
                            if ~isempty(find(depth_iso1{ii}{jj}{kk}(1:(ll - 1), mm), 1, 'last'))
                                depth_min_tmp(ll) ...
                                    = depth_iso1{ii}{jj}{kk}(find(depth_iso1{ii}{jj}{kk}(1:(ll - 1), mm), 1, 'last'), mm);
                            end
                            if ~isempty(find(depth_iso1{ii}{jj}{kk}((ll + 1):end, mm), 1))
                                depth_max_tmp(ll) ...
                                    = depth_iso1{ii}{jj}{kk}(find(depth_iso1{ii}{jj}{kk}((ll + 1):end, mm), 1), mm);
                            end
                        end
                        [depth_all_min{ll}, depth_all_max{ll}] ...
                            = deal([depth_all_min{ll}; depth_min_tmp'], [depth_all_max{ll}; depth_max_tmp']);
                    end
                end
                if ~isempty(depth_uncert_iso1{ii}{jj}{kk})
                    [x_all_depth_uncert_cat, y_all_depth_uncert_cat] ...
                            = deal([x_all_depth_uncert_cat; x_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})'], [y_all_depth_uncert_cat; y_pk{ii}{jj}{kk}(ind_decim_mid{ii}{jj}{kk})']);
                    for ll = 1:num_age_iso
                        depth_uncert_all{ll} ...
                            = [depth_uncert_all{ll}; depth_uncert_iso1{ii}{jj}{kk}(ll, :)'];
                    end
                end
            end
        end
    end
    
    % concatenate and remove the bad stuff prior to gridding
    for ii = 1:num_depth_norm
        [x_all_age{ii}, x_all_age_min{ii}, x_all_age_max{ii}] ...
                            = deal(x_all_age_cat);
        [y_all_age{ii}, y_all_age_min{ii}, y_all_age_max{ii}] ...
                            = deal(y_all_age_cat);
        ind_good            = find(~isnan(age_all{ii}));
        [x_all_age{ii}, y_all_age{ii}, age_all{ii}] ...
                            = deal(x_all_age{ii}(ind_good), y_all_age{ii}(ind_good), age_all{ii}(ind_good));
        if do_nye_norm
            age_diff_all{ii}= age_diff_all{ii}(ind_good);
        end
        [xy_all, ind_unique]= unique([x_all_age{ii} y_all_age{ii}], 'rows');
        [x_all_age{ii}, y_all_age{ii}, age_all{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all{ii}(ind_unique));
        if do_nye_norm
            age_diff_all{ii}= age_diff_all{ii}(ind_unique);
        end
        ind_good            = find(~isnan(age_all_min{ii}));
        [x_all_age_min{ii}, y_all_age_min{ii}, age_all_min{ii}] ...
                            = deal(x_all_age_min{ii}(ind_good), y_all_age_min{ii}(ind_good), age_all_min{ii}(ind_good));
        [xy_all, ind_unique]= unique([x_all_age_min{ii} y_all_age_min{ii}], 'rows');
        [x_all_age_min{ii}, y_all_age_min{ii}, age_all_min{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_min{ii}(ind_unique));
        ind_good            = find(~isnan(age_all_max{ii}));
        [x_all_age_max{ii}, y_all_age_max{ii}, age_all_max{ii}] ...
                            = deal(x_all_age_max{ii}(ind_good), y_all_age_max{ii}(ind_good), age_all_max{ii}(ind_good));
        [xy_all, ind_unique]= unique([x_all_age_max{ii} y_all_age_max{ii}], 'rows');
        [x_all_age_max{ii}, y_all_age_max{ii}, age_all_max{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_max{ii}(ind_unique));
        [x_all_age_uncert{ii}, y_all_age_uncert{ii}] ...
                            = deal(x_all_age_uncert_cat, y_all_age_uncert_cat);
        ind_good            = find(~isnan(age_uncert_all{ii}));
        [x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}] ...
                            = deal(x_all_age_uncert{ii}(ind_good), y_all_age_uncert{ii}(ind_good), age_uncert_all{ii}(ind_good));
        [xy_all, ind_unique]= unique([x_all_age_uncert{ii} y_all_age_uncert{ii}], 'rows');
        [x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_uncert_all{ii}(ind_unique));
    end
    
    for ii = 1:num_age_iso
        [x_all_depth{ii}, x_all_depth_min{ii}, x_all_depth_max{ii}] ...
                            = deal(x_all_depth_cat);
        [y_all_depth{ii}, y_all_depth_min{ii}, y_all_depth_max{ii}] ...
                            = deal(y_all_depth_cat);
        [thick_all_curr, thick_all_min_curr, thick_all_max_curr] ...
                            = deal(thick_cat);
        ind_good            = find(~isnan(depth_all{ii}) & ~isnan(thick_all_curr));
        [x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, thick_all_curr] ...
                            = deal(x_all_depth{ii}(ind_good), y_all_depth{ii}(ind_good), depth_all{ii}(ind_good), thick_all_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth{ii} y_all_depth{ii}], 'rows', 'stable');
        [x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, thick_all_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_all{ii}(ind_unique), thick_all_curr(ind_unique));
        depth_all{ii}       = depth_all{ii} ./ thick_all_curr;
        ind_good            = find(~isnan(depth_all_min{ii}) & ~isnan(thick_all_min_curr));
        [x_all_depth_min{ii}, y_all_depth_min{ii}, depth_all_min{ii}, thick_all_min_curr] ...
                            = deal(x_all_depth_min{ii}(ind_good), y_all_depth_min{ii}(ind_good), depth_all_min{ii}(ind_good), thick_all_min_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth_min{ii} y_all_depth_min{ii}], 'rows', 'stable');
        [x_all_depth_min{ii}, y_all_depth_min{ii}, depth_all_min{ii}, thick_all_min_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_all_min{ii}(ind_unique), thick_all_min_curr(ind_unique));
        depth_all_min{ii}   = depth_all_min{ii} ./ thick_all_min_curr;
        ind_good            = find(~isnan(depth_all_max{ii}) & ~isnan(thick_all_max_curr));
        [x_all_depth_max{ii}, y_all_depth_max{ii}, depth_all_max{ii}, thick_all_max_curr] ...
                            = deal(x_all_depth_max{ii}(ind_good), y_all_depth_max{ii}(ind_good), depth_all_max{ii}(ind_good), thick_all_max_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth_max{ii} y_all_depth_max{ii}], 'rows', 'stable');
        [x_all_depth_max{ii}, y_all_depth_max{ii}, depth_all_max{ii}, thick_all_max_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_all_max{ii}(ind_unique), thick_all_max_curr(ind_unique));
        depth_all_max{ii}   = depth_all_max{ii} ./ thick_all_max_curr;
        [x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, thick_uncert_curr] ...
                            = deal(x_all_depth_uncert_cat, y_all_depth_uncert_cat, thick_cat);
        ind_good            = find(~isnan(depth_uncert_all{ii}) & ~isnan(thick_uncert_curr));
        [x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, thick_uncert_curr] ...
                            = deal(x_all_depth_uncert{ii}(ind_good), y_all_depth_uncert{ii}(ind_good), depth_uncert_all{ii}(ind_good), thick_uncert_curr(ind_good));
        [xy_all, ind_unique]= unique([x_all_depth_uncert{ii} y_all_depth_uncert{ii}], 'rows', 'stable');
        [x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, thick_uncert_curr] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), depth_uncert_all{ii}(ind_unique), thick_uncert_curr(ind_unique));
        depth_uncert_all{ii}= depth_uncert_all{ii} ./ thick_uncert_curr;
    end
    
    clear *_cat
    
    % use Nye normalization and Gaussian fits to remove major anomalies from age_norm1 sets to be gridded
    if do_nye_norm
        for ii = 1:(num_depth_norm - 1)
            fit_curr        = fit((min(age_diff_all{ii}):0.1:max(age_diff_all{ii}))', hist(age_diff_all{ii}, min(age_diff_all{ii}):0.1:max(age_diff_all{ii}))', ['gauss' num2str(num_gauss(ii))]); % Gaussian fit to histogran of age differences from Nye model
            switch num_gauss(ii)
                case 1
                    age_diff_range ...
                            = norminv([0.05 0.99], fit_curr.b1, fit_curr.c1); % lower/upper bounds of acceptable Nye differences
                case 2
                    age_diff_range ...
                            = zeros(2);
                    age_diff_range(1, :) ...
                            = norminv([0.05 0.99], fit_curr.b1, fit_curr.c1);
                    age_diff_range(2, :) ...
                            = norminv([0.05 0.99], fit_curr.b2, fit_curr.c2);
                    age_diff_range ...
                            = [min(age_diff_range(:, 1)) max(age_diff_range(:, 2))]; % take min/max of the two fits' lower/upper bounds, respectively
            end
            ind_good        = find((age_diff_all{ii} >= age_diff_range(1)) & (age_diff_all{ii} <= age_diff_range(2))); % indices of ages with Nye differences within spec
            [x_all_age{ii}, y_all_age{ii}, age_all{ii}, age_diff_all{ii}] ...
                            = deal(x_all_age{ii}(ind_good), y_all_age{ii}(ind_good), age_all{ii}(ind_good), age_diff_all{ii}(ind_good));
            [xy_all, ~, ind_common] ...
                            = intersect([x_all_age{ii} y_all_age{ii}], [x_all_age_min{ii} y_all_age_min{ii}], 'rows', 'stable');
            [x_all_age_min{ii}, y_all_age_min{ii}, age_all_min{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_min{ii}(ind_common));
            [xy_all, ~, ind_common] ...
                            = intersect([x_all_age{ii} y_all_age{ii}], [x_all_age_max{ii} y_all_age_max{ii}], 'rows', 'stable');
            [x_all_age_max{ii}, y_all_age_max{ii}, age_all_max{ii}] ...
                            = deal(xy_all(:, 1), xy_all(:, 2), age_all_max{ii}(ind_common));
        end
    end
%%    
    if parallel_check
        disp('...fixed-depth layers...')
        parfor ii = 1:num_depth_norm
            disp([num2str(ii) '/' num2str(num_depth_norm) '...'])
            try
                age_interpolant ...
                            = scatteredInterpolant(x_all_age{ii}, y_all_age{ii}, age_all{ii}, 'natural', 'none');
                age_norm2(:, :, ii) ...
                            = age_interpolant(x_grd, y_grd);
                age_uncert_interpolant ...
                            = scatteredInterpolant(x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}, 'natural', 'none');
                age_uncert_norm2(:, :, ii) ...
                            = age_uncert_interpolant(x_grd, y_grd);
            catch
                disp([num2str(depth_norm(ii)) ' isopach failed...'])
            end
        end
        disp('...fixed-age layers...')
        parfor ii = 1:num_age_iso
            disp([num2str(ii) '/' num2str(num_age_iso) '...'])
            try
                depth_interpolant ...
                            = scatteredInterpolant(x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, 'natural', 'none');
                depth_iso2(:, :, ii) ...
                            = depth_interpolant(x_grd, y_grd) .* thick_grd;
                depth_uncert_interpolant ...
                            = scatteredInterpolant(x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, 'natural', 'none');
                depth_uncert_iso2(:, :, ii) ...
                            = depth_uncert_interpolant(x_grd, y_grd) .* thick_grd;
            catch
                disp([num2str(age_iso(ii)) '-ka isochrone failed...'])
            end
        end
    else
        disp('...fixed-depth layers...')
        for ii = 1:num_depth_norm
            disp([num2str(ii) '/' num2str(num_depth_norm) '...'])
            try
                age_interpolant ...
                            = scatteredInterpolant(x_all_age{ii}, y_all_age{ii}, age_all{ii}, 'natural', 'none');
                age_norm2(:, :, ii) ...
                            = age_interpolant(x_grd, y_grd);
                age_uncert_interpolant ...
                            = scatteredInterpolant(x_all_age_uncert{ii}, y_all_age_uncert{ii}, age_uncert_all{ii}, 'natural', 'none');
                age_uncert_norm2(:, :, ii) ...
                            = age_uncert_interpolant(x_grd, y_grd);
            catch
                disp([num2str(depth_norm(ii)) ' isopach failed...'])
            end
        end
        disp('...fixed-age layers...')
        for ii = 1:num_age_iso
            disp([num2str(ii) '/' num2str(num_age_iso) '...'])
            try
                depth_interpolant ...
                            = scatteredInterpolant(x_all_depth{ii}, y_all_depth{ii}, depth_all{ii}, 'natural', 'none');
                depth_iso2(:, :, ii) ...
                            = depth_interpolant(x_grd, y_grd) .* thick_grd;
                depth_uncert_interpolant ...
                            = scatteredInterpolant(x_all_depth_uncert{ii}, y_all_depth_uncert{ii}, depth_uncert_all{ii}, 'natural', 'none');
                depth_uncert_iso2(:, :, ii) ...
                            = depth_uncert_interpolant(x_grd, y_grd) .* thick_grd;
            catch
                disp([num2str(age_iso(ii)) '-ka isochrone failed...'])
            end
        end
    end
    
    % trim bad ages/depths
    age_norm2((age_norm2 < 0) | (age_norm2 > age_max)) ...
                            = NaN;
    age_uncert_norm2(age_uncert_norm2 < 0) ...
                            = NaN;
    depth_iso2(depth_iso2 < 0) ...
                            = NaN;
    depth_iso2(depth_iso2 > thick_grd(:, :, ones(1, 1, num_age_iso))) ...
                            = NaN;
    depth_uncert_iso2(depth_uncert_iso2 < 0) ...
                            = NaN;
    
    disp('...done 2-D gridding age.')
%%
    if do_save
        disp('Saving 2-D age grids...')
        switch radar_type
            case 'accum'
                save('mat/age_grd2_accum', '-v7.3', 'age_norm2', 'age_uncert_norm2', 'age_iso', 'depth_iso2', 'depth_norm', 'depth_uncert_iso2', 'num_age_iso', 'num_depth_norm', 'x_grd', 'y_grd')
                disp('Saved 2-D age grids as mat/age_grd2_accum.mat.')
            case 'deep'
                save('mat/age_grd2', '-v7.3', 'age_norm2', 'age_uncert_norm2', 'age_iso', 'depth_iso2', 'depth_norm', 'depth_uncert_iso2', 'num_age_iso', 'num_depth_norm', 'x_grd', 'y_grd')
                disp('Saved 2-D age grids as mat/age_grd2.mat.')
        end
    end
%%
end

if parallel_check
    delete(pool)
end