% FENCE_REL_CHECK Evaluate fenced intersections of merged radar transects.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/07/14

clear

do_check                    = false;
do_save                     = false;
plotting                    = false;

letters                     = 'a':'z';
dist_int_max                = 0.25; % km, +/- range to extract local layer depths
depth_diff_rel_max          = 0.20;
depth_diff_2012             = -15; % m, 2012 P3 offset
dist_check_max              = 100;

load mat/xy_all name_trans num_year

%%
if do_check
%%  
    disp('Evaluating fenced list...')
    
    load mat/merge_all depth_smooth dist thick
    load mat/id_layer_master id_layer_master_mat
    load mat/merge_int_all int_all
    
    int_all                     = int_all(:, [1:4 6:9]);
    
    int_bad_rel                 = [];
    
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
        
        % for each layer pair, loop through each intersection of its transects and sanity-check their depths there
        for jj = 1:size(int_set, 1)
            
            [depth_smooth_int1_rel, depth_smooth_int2_rel] ...
                            = deal(NaN);
            
            ind_start1      = find((ind_sort1 == int_set(jj, 1)), 1);
            ind_start2      = find((ind_sort2 == int_set(jj, 2)), 1);
            
            % near-intersection layer depths
            ind_curr1       = interp1(dist_curr1, 1:length(dist_curr1), (dist_curr1(ind_start1) - dist_int_max), 'nearest', 'extrap'):...
                              interp1(dist_curr1, 1:length(dist_curr1), (dist_curr1(ind_start1) + dist_int_max), 'nearest', 'extrap');
            depth_smooth_int1 ...
                            = nanmean(depth_smooth1(ind_curr1));
            
            ind_curr2       = interp1(dist_curr2, 1:length(dist_curr2), (dist_curr2(ind_start2) - dist_int_max), 'nearest', 'extrap'):...
                              interp1(dist_curr2, 1:length(dist_curr2), (dist_curr2(ind_start2) + dist_int_max), 'nearest', 'extrap');
            depth_smooth_int2 ...
                            = nanmean(depth_smooth2(ind_curr2));
            
            if isnan(depth_smooth_int1)
                ind_tmp1    = find(~isnan(depth_smooth1));
                ind_bound1  = NaN(1, 2);
                if ~isempty(find((ind_tmp1 < ind_start1), 1, 'last'))
                    ind_bound1(1) ...
                            = ind_tmp1(find((ind_tmp1 < ind_start1), 1, 'last'));
                end
                if ~isempty(find((ind_tmp1 > ind_start1), 1))
                    ind_bound1(2) ...
                            = ind_tmp1(find((ind_tmp1 > ind_start1), 1));
                end
                if (~isnan(ind_bound1(1)) && isnan(ind_bound1(2)))
                    [dist_diff1, ind_min1] ...
                            = deal(abs(dist_curr1(ind_bound1(1)) - dist_curr1(ind_start1)), 1);
                elseif (isnan(ind_bound1(1)) && ~isnan(ind_bound1(2)))
                    [dist_diff1, ind_min1] ...
                            = deal(abs(dist_curr1(ind_bound1(2)) - dist_curr1(ind_start1)), 2);
                elseif all(isnan(ind_bound1))
                    continue
                else
                    [dist_diff1, ind_min1] ...
                            = nanmin(abs(dist_curr1(ind_bound1) - dist_curr1(ind_start1)));
                end
                if (dist_diff1 <= dist_check_max)
                    [depth_smooth_int1, thick1] ...
                            = deal(depth_smooth1(ind_bound1(ind_min1)), thick{id1(1)}{id1(2)}{id1(3)}(ind_bound1(ind_min1)));
                end
            else
                thick1      = nanmean(thick{id1(1)}{id1(2)}{id1(3)}(ind_curr1));
            end
            
            if isnan(depth_smooth_int2)
                ind_tmp2    = find(~isnan(depth_smooth2));
                ind_bound2  = NaN(1, 2);
                if ~isempty(find((ind_tmp2 < ind_start2), 1, 'last'))
                    ind_bound2(1) ...
                            = ind_tmp2(find((ind_tmp2 < ind_start2), 1, 'last'));
                end
                if ~isempty(find((ind_tmp2 > ind_start2), 1))
                    ind_bound2(2) ...
                            = ind_tmp2(find((ind_tmp2 > ind_start2), 1));
                end
                if (~isnan(ind_bound2(1)) && isnan(ind_bound2(2)))
                    [dist_diff2, ind_min2] ...
                            = deal(abs(dist_curr2(ind_bound2(1)) - dist_curr2(ind_start2)), 1);
                elseif (isnan(ind_bound2(1)) && ~isnan(ind_bound2(2)))
                    [dist_diff2, ind_min2] ...
                            = deal(abs(dist_curr2(ind_bound2(2)) - dist_curr2(ind_start2)), 2);
                elseif all(isnan(ind_bound2))
                    continue
                else
                    [dist_diff2, ind_min2] ...
                            = nanmin(abs(dist_curr2(ind_bound2) - dist_curr2(ind_start2)));
                end
                if (dist_diff2 <= dist_check_max)
                    [depth_smooth_int2, thick2] ...
                            = deal(depth_smooth2(ind_bound2(ind_min2)), thick{id2(1)}{id2(2)}{id2(3)}(ind_bound2(ind_min2)));
                end
            else
                thick2      = nanmean(thick{id2(1)}{id2(2)}{id2(3)}(ind_curr2));
            end
            
            if any(isnan([depth_smooth_int1 depth_smooth_int2]))
                continue
            end
            
            [depth_smooth_int1_rel, depth_smooth_int2_rel] ...
                            = deal((depth_smooth_int1 ./ thick1), (depth_smooth_int2 ./ thick2));
            
            if (any(isinf([depth_smooth_int1_rel depth_smooth_int2_rel])) || any(isnan([depth_smooth_int1_rel depth_smooth_int2_rel])))
                continue
            end
            
            if (abs(depth_smooth_int1_rel - depth_smooth_int2_rel) > depth_diff_rel_max)
                int_bad_rel     = [int_bad_rel; id_layer_master_mat(ii, 1:8) (depth_smooth_int1_rel - depth_smooth_int2_rel)]; %#ok<AGROW>
                disp([num2str(ii) '/' num2str(size(id_layer_master_mat, 1))])
            end
        end
    end
    
    if do_save
        save mat/int_bad_rel -v7.3 int_bad_rel
    end
    
else
    load mat/int_bad_rel int_bad_rel
end
%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%%
    [~, int_bad_rel_sort]   = unique(int_bad_rel(:, 1:8), 'rows');
    int_bad_rel_sort        = sortrows(int_bad_rel(int_bad_rel_sort, :), [1 2 3 5 6 7 4 8 9]);
%     [~, ind_sort]           = sort(abs(int_bad_rel_sort(:, end)), 'descend');
%     int_bad_rel_sort        = int_bad_rel_sort(ind_sort, :);
    name_trans_bad          = cell(size(int_bad_rel_sort, 1), 2);
    for ii = 1:size(int_bad_rel_sort, 1)
        if int_bad_rel_sort(ii, 3)
            name_trans_bad{ii, 1} ...
                            = [name_trans{int_bad_rel_sort(ii, 1)}{int_bad_rel_sort(ii, 2)} letters(int_bad_rel_sort(ii, 3))];
        else
            name_trans_bad{ii, 1} ...
                            = name_trans{int_bad_rel_sort(ii, 1)}{int_bad_rel_sort(ii, 2)};
        end
        if int_bad_rel_sort(ii, 7)
            name_trans_bad{ii, 2} ...
                            = [name_trans{int_bad_rel_sort(ii, 5)}{int_bad_rel_sort(ii, 6)} letters(int_bad_rel_sort(ii, 7))];
        else
            name_trans_bad{ii, 2} ...
                            = name_trans{int_bad_rel_sort(ii, 5)}{int_bad_rel_sort(ii, 6)};
        end
    end
    figure('name', 'Bad relative depth intersection table')
    uitable('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'fontsize', 18, ...
            'columnname', {'Name (1)' 'Year(1)' 'Transect (1)' 'Sub-transect (1)' 'Layer # (1)' 'Name (2)' 'Year (2)' 'Transect (2)' 'Sub-transect (2)' 'Layer # (2)', 'Relative depth difference (%)'}, ...
            'data', [name_trans_bad(:, 1) mat2cell(int_bad_rel_sort(:, 1:4), ones(size(int_bad_rel_sort, 1), 1), ones(1, 4)) name_trans_bad(:, 2) mat2cell([int_bad_rel_sort(:, 5:8) (100 .* int_bad_rel_sort(:, 9))], ones(size(int_bad_rel_sort, 1), 1), ones(1, 5))]);
end