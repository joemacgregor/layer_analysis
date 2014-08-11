% FENCE_CHECK Evaluate fenced intersections of merged radar transects.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/07/14

clear

do_check                    = false;
do_save                     = false;
plotting                    = false;

letters                     = 'a':'z';
dist_int_max                = 0.25; % km, +/- range to extract local layer depths
range_factor                = 10;
depth_diff_2012             = -15; % m, 2012 P3 offset

load mat/xy_all name_trans num_year

%%
if do_check
%%  
    disp('Evaluating fenced list...')
    
    load mat/range_resolution range_resolution
    load mat/merge_all depth_smooth dist
    load mat/id_layer_master id_layer_master_mat
    load mat/merge_int_all int_all
    
    int_all                     = int_all(:, [1:4 6:9]);
    
    int_bad                     = [];
    
    % evaluate each existing transect intersection
    for ii = 1:size(id_layer_master_mat, 1)
        
        % current transect/layer ID pair
        [id1, id2]          = deal(id_layer_master_mat(ii, 1:4), id_layer_master_mat(ii, 5:8));
        if ~id1(3)
            id1(3)          = 1;
        end
        if ~id2(3)
            id2(3)          = 1;
        end
        
        % intersection indices in merged files
        int_set             = [int_all((ismember(int_all(:, 1:3), id_layer_master_mat(ii, 1:3), 'rows') & ismember(int_all(:, 5:7), id_layer_master_mat(ii, 5:7), 'rows')), [4 8]); ...
                               int_all((ismember(int_all(:, 5:7), id_layer_master_mat(ii, 1:3), 'rows') & ismember(int_all(:, 1:3), id_layer_master_mat(ii, 5:7), 'rows')), [8 4])];
        
        % extract unique current distances and depths for each transect's layer
        [dist_curr1, ind_sort1] ...
                            = unique(dist{id1(1)}{id1(2)}{id1(3)});
        depth_smooth1       = double(depth_smooth{id1(1)}{id1(2)}{id1(3)}(id1(4), ind_sort1));
        
        % adjust 2012 data
        if (id1(1) == num_year)
            depth_smooth1   = depth_smooth1 + depth_diff_2012;
        end
        
        [dist_curr2, ind_sort2] ... 
                            = unique(dist{id2(1)}{id2(2)}{id2(3)});
        depth_smooth2       = double(depth_smooth{id2(1)}{id2(2)}{id2(3)}(id2(4), ind_sort2));
        if (id2(1) == num_year)
            depth_smooth2   = depth_smooth2 + depth_diff_2012;
        end
        
        % set threshold depth difference to be proportional to range resolution of layer pair's respective campaigns
        depth_diff_curr     = range_factor * mean(range_resolution([id1(1) id2(1)]));
        
        % for this layer pair, loop through each intersection of its transects and sanity-check their depths there
        for jj = 1:size(int_set, 1)
            
            ind_start1      = find((ind_sort1 == int_set(jj, 1)), 1);
            ind_start2      = find((ind_sort2 == int_set(jj, 2)), 1);
            
            % near-intersection layer depths
            ind_curr1       = interp1(dist_curr1, 1:length(dist_curr1), (dist_curr1(ind_start1) - dist_int_max), 'nearest', 'extrap'):interp1(dist_curr1, 1:length(dist_curr1), (dist_curr1(ind_start1) + dist_int_max), 'nearest', 'extrap');
            depth_smooth_int1 ...
                            = nanmean(depth_smooth1(:, ind_curr1));
            
            ind_curr2       = interp1(dist_curr2, 1:length(dist_curr2), (dist_curr2(ind_start2) - dist_int_max), 'nearest', 'extrap'):interp1(dist_curr2, 1:length(dist_curr2), (dist_curr2(ind_start2) + dist_int_max), 'nearest', 'extrap');
            depth_smooth_int2 ...
                            = nanmean(depth_smooth2(:, ind_curr2));
            
            if (isnan(depth_smooth_int1) || isnan(depth_smooth_int2))
                continue
            end
            
            if (abs(depth_smooth_int1 - depth_smooth_int2) > depth_diff_curr)
                int_bad     = [int_bad; id_layer_master_mat(ii, 1:4) int_set(jj, 1) dist_curr1(ind_start1) id_layer_master_mat(ii, 5:8) int_set(jj, 2) dist_curr2(ind_start2) (depth_smooth_int1 - depth_smooth_int2)]; %#ok<AGROW>
                disp([num2str(ii) '/' num2str(size(id_layer_master_mat, 1))])
            end
        end
    end
    
    if do_save
        save mat/int_bad -v7.3 int_bad
    end
    
else
    load mat/int_bad int_bad
end
%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%%
    [~, ind_sort]           = sort(int_bad(:, 13));
    int_bad_sort            = int_bad(ind_sort, :);
    name_trans_bad          = cell(size(int_bad_sort, 1), 2);
    for ii = 1:size(int_bad_sort, 1)
        if int_bad_sort(ii, 3)
            name_trans_bad{ii, 1} ...
                            = [name_trans{int_bad_sort(ii, 1)}{int_bad_sort(ii, 2)} letters(int_bad_sort(ii, 3))];
        else
            name_trans_bad{ii, 1} ...
                            = name_trans{int_bad_sort(ii, 1)}{int_bad_sort(ii, 2)};
        end
        if int_bad_sort(ii, 9)
            name_trans_bad{ii, 2} ...
                            = [name_trans{int_bad_sort(ii, 7)}{int_bad_sort(ii, 8)} letters(int_bad_sort(ii, 9))];
        else
            name_trans_bad{ii, 2} ...
                            = name_trans{int_bad_sort(ii, 7)}{int_bad_sort(ii, 8)};
        end
    end
    figure('name', 'Bad intersection table')
    uitable('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'fontsize', 18, ...
            'columnname', {'Name (1)' 'Year(1)' 'Transect (1)' 'Sub-transect (1)' 'Layer # (1)' 'Distance (km) (1)' ...
                           'Name (2)' 'Year (2)' 'Transect (2)' 'Sub-transect (2)' 'Layer # (2)', 'Distance (km) (2)' 'Depth difference (m)'}, ...
            'data', [name_trans_bad(:, 1) mat2cell(int_bad_sort(:, [1:4 6]), ones(size(int_bad_sort, 1), 1), ones(1, 5)) name_trans_bad(:, 2) mat2cell(int_bad_sort(:, [7:10 12 13]), ones(size(int_bad_sort, 1), 1), ...
                     ones(1, 6))])
end