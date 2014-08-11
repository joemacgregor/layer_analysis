% MERGEINT Determine intersections of merged radar transects.
% 
%   MERGEINT determines the positions of intersections between each
%   merged transect and all other merged transects.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/08/14

clear

radar_type                  = 'deep';
letters                     = 'a':'z';
dir_save                    = 'mat/';
do_xy                       = true;
do_int                      = false;
do_fence                    = false;
dist_density                = 2.5; % maximum intersection density in km
angle_threshold             = 2.5; % minimum intersection angle in degrees
decim                       = 25;
dist_int_max                = 0.25; % km, +/- range to extract local layer depths
depth_diff_2012             = -15; % m, 2012 P3 offset
plotting                    = false;

switch radar_type
    case 'accum'
        load mat/xy_all_accum name_trans name_year num_trans num_year
        load mat/range_resolution_accum range_resolution
    case 'deep'
        load mat/xy_all name_trans name_year num_trans num_year
        load mat/range_resolution range_resolution
end

%% concatenate x/y positions for each campaign

if do_xy
    
    disp('Extracting x/y positions from merged picks files...')
    
    [depth_smooth, dist, dist_diff, ind_merge, num_merge, x, y] ...
                            = deal(cell(1, num_year));
    
    for ii = 1:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        [depth_smooth{ii}, dist{ii}, dist_diff{ii}, ind_merge{ii}, x{ii}, y{ii}] ...
                            = deal(cell(1, num_trans(ii)));
        num_merge{ii}       = deal(zeros(1, num_trans(ii)));
        
        for jj = 1:num_trans(ii)
            
            disp([name_trans{ii}{jj} ' (' num2str(jj) ' / ' num2str(num_trans(ii)) ')...'])
            
            if exist([name_year{ii} '/merge/' name_trans{ii}{jj} '_pk_merge.mat'], 'file')
                disp(['Transect: ' name_trans{ii}{jj} '...'])
                pk          = load([name_year{ii} '/merge/' name_trans{ii}{jj} '_pk_merge.mat']); % load a transect's merge file (no a/b/c/etc, just a single file)
                pk          = pk.pk;
                [depth_smooth{ii}{jj}{1}, dist{ii}{jj}{1}, x{ii}{jj}{1}, y{ii}{jj}{1}] ...
                            = deal(pk.depth_smooth, pk.dist, pk.x, pk.y);
                dist_diff{ii}{jj}{1} ...
                            = diff(dist{ii}{jj}{1});
                num_merge{ii}(jj) ...
                            = 1;
                ind_merge{ii}{jj} ...
                            = 1;
            else
                for kk = 1:length(letters)
                    if ~exist([name_year{ii} '/merge/' name_trans{ii}{jj} letters(kk) '_pk_merge.mat'], 'file') % load merge file (a/b/c/etc)
                        continue
                    end
                    disp(['Transect: ' name_trans{ii}{jj} letters(kk) '...'])
                    pk      = load([name_year{ii} '/merge/' name_trans{ii}{jj} letters(kk) '_pk_merge.mat']);
                    pk      = pk.pk;
                    num_merge{ii}(jj) ...
                            = num_merge{ii}(jj) + 1;
                    [depth_smooth{ii}{jj}{kk}, dist{ii}{jj}{kk}, x{ii}{jj}{kk}, y{ii}{jj}{kk}] ...
                            = deal(pk.depth_smooth, pk.dist, pk.x, pk.y);
                    dist_diff{ii}{jj}{kk} ...
                            = diff(dist{ii}{jj}{kk});
                    ind_merge{ii}{jj} ...
                            = [ind_merge{ii}{jj} kk];
                end
            end
        end
    end
    
    switch radar_type
        case 'accum'
            save([dir_save 'merge_xy_all_accum'], '-v7.3', 'depth_smooth', 'dist', 'dist_diff', 'ind_merge', 'num_merge', 'x', 'y')
        case 'deep'
            save([dir_save 'merge_xy_all'], '-v7.3', 'depth_smooth', 'dist', 'dist_diff', 'ind_merge', 'num_merge', 'x', 'y')
    end
    disp(['Done saving key merge data in ' dir_save '.'])
    
else
    switch radar_type
        case 'accum'
            load([dir_save 'merge_xy_all_accum'], 'depth_smooth', 'dist', 'dist_diff', 'ind_merge', 'num_merge', 'x', 'y')
        case 'deep'
            load([dir_save 'merge_xy_all'], 'depth_smooth', 'dist', 'dist_diff', 'ind_merge', 'num_merge', 'x', 'y')
    end    
    disp(['Loaded key merge data from ' dir_save '.'])
end

%% Calculate intersections between each transect all others (all others in that campaign and all other years)

if do_int
    
    disp('Finding all intersections...')
    
    num_int                 = deal(cell(1, num_year));
    int_all                 = [];
    
    for ii = 1:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        num_int{ii}         = zeros(1, num_trans(ii));
        
        for jj = ii:num_year
            
            disp(['...compared with ' name_year{jj} '...'])
            
            for kk = find(num_merge{ii})
                
                ll_tmp      = find(num_merge{jj});
                if (ii == jj)
                    ll_tmp  = ll_tmp(ll_tmp >= kk);
                end
                
                if isempty(ll_tmp)
                    continue
                end
                
                disp(['...transect ' name_trans{ii}{kk} ' compared with...' ])
                
                for ll = ll_tmp
                    
                    disp(['...' name_trans{jj}{ll} '...'])
                    
                    tmp     = size(int_all, 1);
                    
                    for nn = ind_merge{ii}{kk}
                        
                        oo_tmp                              = ind_merge{jj}{ll};
                        if ((ii == jj) && (kk == ll))
                            oo_tmp                          = oo_tmp(oo_tmp > nn);
                        end
                        
                        if isempty(oo_tmp)
                            continue
                        end
                        
                        ind_tmp1                            = 1:decim:length(x{ii}{kk}{nn});
                        ind_tmp1(end)                       = length(x{ii}{kk}{nn});
                        
                        for oo = oo_tmp
                            
                            ind_tmp2                        = 1:decim:length(x{jj}{ll}{oo});
                            ind_tmp2(end)                   = length(x{jj}{ll}{oo});
                            
                            % intersection of two vectors using linear interpolation
                            [~, ~, ind_int1, ind_int2]      = intersecti(x{ii}{kk}{nn}(ind_tmp1), y{ii}{kk}{nn}(ind_tmp1), x{jj}{ll}{oo}(ind_tmp2), y{jj}{ll}{oo}(ind_tmp2));
                            
                            if isempty(ind_int1)
                                continue
                            end
                            
                            for mm = 1:length(ind_int1)
                                
                                if (ind_int1(mm) == 1)
                                    ind_tmp3                = 1:ind_tmp1(2);
                                elseif (ind_int1(mm) == length(ind_tmp1))
                                    ind_tmp3                = ind_tmp1(end - 1):ind_tmp1(end);
                                else
                                    ind_tmp3                = ind_tmp1(ind_int1(mm) - 1):ind_tmp1(ind_int1(mm) + 1);
                                end
                                
                                if (ind_int2(mm) == 1)
                                    ind_tmp4                = 1:ind_tmp2(2);
                                elseif (ind_int2(mm) == length(ind_tmp2))
                                    ind_tmp4                = ind_tmp2(end - 1):ind_tmp2(end);
                                else
                                    ind_tmp4                = ind_tmp2(ind_int2(mm) - 1):ind_tmp2(ind_int2(mm) + 1);
                                end
                                
                                % intersection of two vectors using linear interpolation
                                [~, ~, ind_int3, ind_int4]  = intersecti(x{ii}{kk}{nn}(ind_tmp3), y{ii}{kk}{nn}(ind_tmp3), x{jj}{ll}{oo}(ind_tmp4), y{jj}{ll}{oo}(ind_tmp4));
                                
                                if ~isempty(ind_int3) % only when there is an intersection will we make a matrix for ind_match_trans
                                    [ind_int3, ind_int4]    = deal((ind_int3 + ind_tmp3(1) - 1), (ind_int4 + ind_tmp4(1) - 1));
                                    if ~isempty(find((diff(dist{ii}{kk}{nn}(ind_int3)) < dist_density), 1)) % ensure density of consecutive intersections is not too great
                                        tmp1                = find((diff(dist{ii}{kk}{nn}(ind_int3)) >= dist_density));
                                        [ind_int3, ind_int4]= deal(ind_int3(tmp1), ind_int4(tmp1));
                                    end
                                    if ~isempty(find((diff(dist{jj}{ll}{oo}(ind_int4)) < dist_density), 1))
                                        tmp1                = find(diff(dist{jj}{ll}{oo}(ind_int4)) >= dist_density);
                                        [ind_int3, ind_int4]= deal(ind_int3(tmp1), ind_int4(tmp1));
                                    end
                                    if ~isempty(ind_int3)
                                        disp(['...' letters(nn) ' ' letters(oo) '...'])
                                        int_all = [int_all; ...
                                                   repmat([ii kk nn], length(ind_int3), 1) ind_int3 NaN(length(ind_int3), 1) repmat([jj ll oo], length(ind_int3), 1) ind_int4 NaN(length(ind_int3), 1)]; %#ok<AGROW>
                                    end
                                end
                            end
                        end
                    end
                    
                    ind2remove      = [];
                    
                    % remove intersections due to transect jumps or small intersection angles or more than one intersection per km
                    for nn = (tmp + 1):size(int_all, 1)
                        
                        ind_tmp1    = (int_all(nn, 4) - 1):(int_all(nn, 4) + 1);
                        ind_tmp1    = ind_tmp1((ind_tmp1 > 0) & (ind_tmp1 <= length(dist_diff{ii}{kk}{int_all(nn, 3)})));
                        
                        if any(dist_diff{ii}{kk}{int_all(nn, 3)}(ind_tmp1) > (5 * median(dist_diff{ii}{kk}{int_all(nn, 3)})))
                            ind2remove ...
                                    = [ind2remove nn]; %#ok<AGROW>
                        end
                        
                        ind_tmp2    = (int_all(nn, 9) - 1):(int_all(nn, 9) + 1);
                        ind_tmp2    = ind_tmp2((ind_tmp2 > 0) & (ind_tmp2 <= length(dist_diff{jj}{ll}{int_all(nn, 8)})));
                        
                        if any(dist_diff{jj}{ll}{int_all(nn, 8)}(ind_tmp2) > (5 * median(dist_diff{jj}{ll}{int_all(nn, 8)})))
                            ind2remove ...
                                    = [ind2remove nn]; %#ok<AGROW>
                        end
                        
                        int_all(nn, [5 10]) ...
                                    = [atan2d(diff(y{ii}{kk}{int_all(nn, 3)}(ind_tmp1([1 end]))), diff(x{ii}{kk}{int_all(nn, 3)}(ind_tmp1([1 end])))) ...
                                       atan2d(diff(y{jj}{ll}{int_all(nn, 8)}(ind_tmp2([1 end]))), diff(x{jj}{ll}{int_all(nn, 8)}(ind_tmp2([1 end]))))]; %#ok<SAGROW>
                        
                        % transect intersection angle must be greater than angle_threshold
                        if (abs(diff(int_all(nn, [5 10]))) < angle_threshold)
                            ind2remove ...
                                    = [ind2remove nn]; %#ok<AGROW>
                        end
                    end
                    
                    ind2remove      = unique(ind2remove);
                    
                    if ~isempty(ind2remove)
                        int_all     = int_all(setdiff(1:size(int_all, 1), ind2remove), :); % the transect intersection where data actually exist
                    end
                    
                    num_int{ii}(jj) = length((tmp + 1):size(int_all, 1));
                    
                end
            end
        end
    end
    
    for ii = find(int_all(:, 3) == 1)'
        if exist([name_year{int_all(ii, 1)} '/merge/' name_trans{int_all(ii, 1)}{int_all(ii, 2)} '_pk_merge.mat'], 'file')
            int_all(ii, 3)  = 0; %#ok<SAGROW>
        end
    end
    for ii = find(int_all(:, 8) == 1)'
        if exist([name_year{int_all(ii, 6)} '/merge/' name_trans{int_all(ii, 6)}{int_all(ii, 7)} '_pk_merge.mat'], 'file')
            int_all(ii, 8)  = 0; %#ok<SAGROW>
        end
    end
    
    switch radar_type
        case 'accum'
            save([dir_save 'merge_int_all_accum'], '-v7.3', 'int_all')
        case 'deep'
            save([dir_save 'merge_int_all'], '-v7.3', 'int_all')
    end
    disp(['Done calculating all intersections and saved in ' dir_save '.'])
    
else
    switch radar_type
        case 'accum'
            load([dir_save 'merge_int_all_accum'], 'int_all')
        case 'deep'
            load([dir_save 'merge_int_all'], 'int_all')
    end
    disp(['Loaded all intersections from ' dir_save '.'])
end

%% Populate fence log

if do_fence
    
%     disp('Preparing fencing spreadsheet...')
%     
%     fence                       = cell(size(int_all, 1), 4);
%     for ii = 1:size(int_all, 1)
%         fence{ii, 1}            = name_year{int_all(ii, 1)};
%         if ~int_all(ii, 3)
%             fence{ii, 2}        = [name_trans{int_all(ii, 1)}{int_all(ii, 2)}];
%         else
%             fence{ii, 2}        = [name_trans{int_all(ii, 1)}{int_all(ii, 2)} letters(int_all(ii, 3))];
%         end
%         fence{ii, 3}            = name_year{int_all(ii, 6)};
%         if ~int_all(ii, 8)
%             fence{ii, 4}        = [name_trans{int_all(ii, 6)}{int_all(ii, 7)}];
%         else
%             fence{ii, 4}        = [name_trans{int_all(ii, 6)}{int_all(ii, 7)} letters(int_all(ii, 8))];
%         end
%     end
%     
%     [~, ind_unique]             = unique(int_all(:, [1 2 3 6 7 8]), 'rows');
%     fence                       = fence(ind_unique, :);
    
    disp('Auto-populating fenced list...')
    
    % initialize master layer intersection lists
    id_layer_master_mat     = [];
    [id_layer_master_cell, ind_layer_auto] ...
                            = deal(cell(1, num_year));
    for ii = 1:num_year
        [id_layer_master_cell{ii}, ind_layer_auto{ii}] ...
                            = deal(cell(1, num_trans(ii)));
        for jj = 1:num_trans(ii)
            [id_layer_master_cell{ii}{jj}, ind_layer_auto{ii}{jj}] ...
                            = deal(cell(1, max(ind_merge{ii}{jj})));
        end
    end
    
    % correct intersection list for 0/1 subtransect ambiguity
    int_all_fix             = int_all;
    int_all_fix(~int_all_fix(:, 3), 3) ...
                            = 1;
    int_all_fix(~int_all_fix(:, 8), 8) ...
                            = 1;
    
    % evaluate each transect intersection
    for ii = 1:size(int_all_fix, 1)
        
        disp([num2str(ii) '/' num2str(size(int_all, 1)) '...'])
        
        % extract unique current distances and depths for each transect
        [dist_curr1, ind_sort1] ...
                            = unique(dist{int_all_fix(ii, 1)}{int_all_fix(ii, 2)}{int_all_fix(ii, 3)});
        if (length(dist_curr1) < length(dist{int_all_fix(ii, 1)}{int_all_fix(ii, 2)}{int_all_fix(ii, 3)}))
            if int_all(ii, 3)
                disp([name_trans{int_all(ii, 1)}{int_all(ii, 2)} letters(int_all(ii, 3)) ' has a funky distance vector (' ...
                      num2str(length(dist{int_all_fix(ii, 1)}{int_all_fix(ii, 2)}{int_all_fix(ii, 3)})) ' total vs. ' num2str(length(dist_curr1)) ' unique).'])
            else
                disp([name_trans{int_all(ii, 1)}{int_all(ii, 2)} ' has a funky distance vector (' ...
                      num2str(length(dist{int_all_fix(ii, 1)}{int_all_fix(ii, 2)}{int_all_fix(ii, 3)})) ' total vs. ' num2str(length(dist_curr1)) ' unique).'])
            end
        end
        if ~any(ind_sort1 == int_all_fix(ii, 4)) % skip if original intersection index no longer exists
            continue
        else
            depth_smooth_curr1 ...
                            = depth_smooth{int_all_fix(ii, 1)}{int_all_fix(ii, 2)}{int_all_fix(ii, 3)}(:, ind_sort1);
            ind_sort1       = find((ind_sort1 == int_all_fix(ii, 4)), 1);
        end
        [dist_curr2, ind_sort2] ...
                            = unique(dist{int_all_fix(ii, 6)}{int_all_fix(ii, 7)}{int_all_fix(ii, 8)});
        if (length(dist_curr2) < length(dist{int_all_fix(ii, 6)}{int_all_fix(ii, 7)}{int_all_fix(ii, 8)}))
            if int_all(ii, 8)
                disp([name_trans{int_all(ii, 6)}{int_all(ii, 7)} letters(int_all(ii, 8)) ' has a funky distance vector (' ...
                      num2str(length(dist{int_all_fix(ii, 6)}{int_all_fix(ii, 7)}{int_all_fix(ii, 8)})) ' total vs. ' num2str(length(dist_curr2)) ' unique).'])
            else
                disp([name_trans{int_all(ii, 6)}{int_all(ii, 7)} ' has a funky distance vector (' ...
                      num2str(length(dist{int_all_fix(ii, 6)}{int_all_fix(ii, 7)}{int_all_fix(ii, 8)})) ' total vs. ' num2str(length(dist_curr2)) ' unique).'])
            end
        end
        if ~any(ind_sort2 == int_all_fix(ii, 9))
            continue
        else
            depth_smooth_curr2 ...
                            = depth_smooth{int_all_fix(ii, 6)}{int_all_fix(ii, 7)}{int_all_fix(ii, 8)}(:, ind_sort2);
            ind_sort2       = find((ind_sort2 == int_all_fix(ii, 9)), 1);
        end
        
        % near-intersection layer depths
        ind_curr1           = interp1(dist_curr1, 1:length(dist_curr1), (dist_curr1(ind_sort1) - dist_int_max), 'nearest', 'extrap'):...
                              interp1(dist_curr1, 1:length(dist_curr1), (dist_curr1(ind_sort1) + dist_int_max), 'nearest', 'extrap');
        depth_smooth_curr1  = nanmean(depth_smooth_curr1(:, ind_curr1), 2);
        
        ind_curr2           = interp1(dist_curr2, 1:length(dist_curr2), (dist_curr2(ind_sort2) - dist_int_max), 'nearest', 'extrap'):...
                              interp1(dist_curr2, 1:length(dist_curr2), (dist_curr2(ind_sort2) + dist_int_max), 'nearest', 'extrap');
        depth_smooth_curr2  = nanmean(depth_smooth_curr2(:, ind_curr2), 2);
        
        % number of layer matches
        kk                  = 0;
        
        % skip if no layers to match
        if (isempty(find(~isnan(depth_smooth_curr1), 1)) || isempty(find(~isnan(depth_smooth_curr2), 1)))
            continue
        end
        
        % adjust 2012 data
        if (int_all_fix(ii, 1) == num_year)
            depth_smooth_curr1 ...
                            = depth_smooth_curr1 + depth_diff_2012;
        end
        if (int_all_fix(ii, 6) == num_year)
            depth_smooth_curr2 ...
                            = depth_smooth_curr2 + depth_diff_2012;
        end
        
        % set threshold depth difference to minimum range resolution of layer pair's respective campaigns
        depth_diff_curr     = min(range_resolution(int_all_fix(ii, [1 6])));
        
        for jj = find(~isnan(depth_smooth_curr1))'
            
            % only keep going if a layer pair meets the current maximum range resolution
            if (nanmin(abs(depth_smooth_curr1(jj) - depth_smooth_curr2)) > depth_diff_curr)
                continue
            end
            
            % layer index in intersecting (second) transect of the pair that matches the master (first) transect layer depth
            [~, curr_layer] = nanmin(abs(depth_smooth_curr1(jj) - depth_smooth_curr2));
            curr_layer      = curr_layer(1);
            
            % skip if layer pair has already been matched
            if (ii > 1)
                if ~isempty(find(ismember(id_layer_master_mat(:, 1:8), [int_all(ii, 1:3) jj int_all(ii, 6:8) curr_layer], 'rows'), 1))
                    continue
                elseif ~isempty(find(ismember(id_layer_master_mat(:, 1:8), [int_all(ii, 6:8) curr_layer int_all(ii, 1:3) jj], 'rows'), 1))
                    continue
                end
            end
            
            % add to master list
            id_layer_master_mat ...
                            = [id_layer_master_mat; int_all(ii, 1:3) jj int_all(ii, 6:8) curr_layer NaN]; %#ok<AGROW>
            kk              = kk + 1; % increment number of matches
            
            % assign master ID; only assign new master ID if neither of the layer pair have already been matched
            curr_ind        = find(ismember(id_layer_master_mat(1:(end - 1), 1:4), [int_all(ii, 1:3) jj], 'rows') | ismember(id_layer_master_mat(1:(end - 1), 5:8), [int_all(ii, 1:3) jj], 'rows') | ...
                                   ismember(id_layer_master_mat(1:(end - 1), 1:4), [int_all(ii, 6:8) curr_layer], 'rows') | ismember(id_layer_master_mat(1:(end - 1), 5:8), [int_all(ii, 6:8) curr_layer], 'rows'));
            if ~isempty(curr_ind)
                id_layer_master_mat([curr_ind; end], end) ...
                            = min(id_layer_master_mat(curr_ind, end));
            else
                id_layer_master_mat(end, end) ...
                            = max([0; id_layer_master_mat(:, end)]) + 1;
            end
        end
        
        if kk
            disp([num2str(kk) ' matches.'])
        end
    end
    
    % recheck for repeated layers
    for ii = 1:size(id_layer_master_mat, 1)
        curr_ind            = find(ismember(id_layer_master_mat(:, 1:4), id_layer_master_mat(ii, 1:4), 'rows') | ismember(id_layer_master_mat(:, 5:8), id_layer_master_mat(ii, 1:4), 'rows') | ...
                                   ismember(id_layer_master_mat(:, 1:4), id_layer_master_mat(ii, 5:8), 'rows') | ismember(id_layer_master_mat(:, 5:8), id_layer_master_mat(ii, 5:8), 'rows'));
        id_layer_master_mat(curr_ind, end) ...
                            = min(id_layer_master_mat(curr_ind, end)); %#ok<SAGROW>
    end
    
    % add layer matches to cell lists
    for ii = 1:size(id_layer_master_mat, 1)
        curr_id         = id_layer_master_mat(ii, :);
        curr_id(~curr_id) ...
                        = 1;
        id_layer_master_cell{curr_id(1)}{curr_id(2)}{curr_id(3)} ...
                        = [id_layer_master_cell{curr_id(1)}{curr_id(2)}{curr_id(3)}; id_layer_master_mat(ii, 4:end)];
        id_layer_master_cell{curr_id(5)}{curr_id(6)}{curr_id(7)} ...
                        = [id_layer_master_cell{curr_id(5)}{curr_id(6)}{curr_id(7)}; id_layer_master_mat(ii, [8 1:4 9])];
        ind_layer_auto{curr_id(1)}{curr_id(2)}{curr_id(3)} ...
                        = [ind_layer_auto{curr_id(1)}{curr_id(2)}{curr_id(3)}; id_layer_master_mat(ii, 4:end)];
        ind_layer_auto{curr_id(5)}{curr_id(6)}{curr_id(7)} ...
                        = [ind_layer_auto{curr_id(5)}{curr_id(6)}{curr_id(7)}; id_layer_master_mat(ii, [8 1:4 9])];
    end
    
    switch radar_type
        case 'accum'
            save([dir_save 'id_layer_master_auto_accum'], '-v7.3', 'id_layer_master_cell', 'id_layer_master_mat', 'ind_layer_auto')
        case 'deep'
            save([dir_save 'id_layer_master_auto'], '-v7.3', 'id_layer_master_cell', 'id_layer_master_mat', 'ind_layer_auto')
    end    
end