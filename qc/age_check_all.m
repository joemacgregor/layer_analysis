% AGE_CHECK Check for age overturning.
% 
% Joe MacGregor (UTIG)
% Last updated: 06/17/14

clear

do_age                      = true;
do_save                     = true;
plotting                    = false;

letters                     = 'a':'z';

load mat/xy_all name_trans num_trans num_year

%%
if do_age
%%  
    disp('Evaluating each transect''s layers for jumps...')
    
    load mat/merge_all depth_smooth ind_fence num_fence num_layer num_trace_tot
    load mat/date_all age age_match age_ord
    
    layer_bad               = [];
%%    
    % evaluate each existing layer intersection
    for ii = 1:num_year
        
        for jj = find(num_fence{ii})
            
            for kk = find(ind_fence{ii}{jj})
                
                if (num_layer{ii}{jj}(kk) < 2)
                    continue
                end
                
                disp([name_trans{ii}{jj} letters(kk)])
                
                % check for and report overturned ages
                ind_layer_overturn ...
                            = [];
                
                % sort depths, get ages at those sorted depths, difference at each trace, and then find the bad ones
                [depth_sort, ind_sort_depth] ...
                            = sort(depth_smooth{ii}{jj}{kk});
                age_sort    = age{ii}{jj}{kk}(ind_sort_depth);
                age_sort(isnan(depth_sort)) ...
                            = NaN;
                age_diff    = diff(age_sort);
                [ind_age_bad, ind_trace_bad] ...
                            = find(age_diff <= 0);
                
                % only necessary if bad ages found
                if ~isempty(ind_age_bad)
                    ind_layer_overturn ...
                            = NaN(length(ind_age_bad), 5);
                    for ll = 1:length(ind_age_bad)
                        ind_layer_overturn(ll, 1:2) ...
                            = ind_sort_depth([ind_age_bad(ll) (ind_age_bad(ll) + 1)], ind_trace_bad(ll));
                    end
                    ind_layer_overturn(:, 3) ...
                            = 1e-3 .* age_diff(sub2ind([(num_layer{ii}{jj}(kk) - 1) num_trace_tot{ii}{jj}(kk)], ind_age_bad, ind_trace_bad));
                    ind_layer_overturn(:, 4:5) ...
                            = age_ord{ii}{jj}{kk}(ind_layer_overturn(:, 1:2));
                    ind_layer_overturn ...
                            = unique(ind_layer_overturn, 'rows');
                end
                
                age_match_overturn ...
                            = NaN(size(ind_layer_overturn, 1), 2);
                for ll = 1:size(ind_layer_overturn, 1)
                    age_match_overturn(ll, :) ...
                            = age_match{ii}{jj}{kk}(ind_layer_overturn(ll, 1:2))';
                end
                
                if ~isempty(ind_layer_overturn)
                    disp('Layer pairs overturned or equal...')
                    disp(num2str(ind_layer_overturn))
                    layer_bad ...
                            = [layer_bad; repmat([ii jj kk], size(ind_layer_overturn, 1), 1) ind_layer_overturn age_match_overturn]; %#ok<AGROW>
                end
            end
        end
    end
%%    
    if do_save
        save mat/age_bad -v7.3 layer_bad
    end
    
else
    load mat/age_bad layer_bad
end
%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%%
    layer_bad_sort          = sortrows(layer_bad, [6 1 2 3 7 8]);
    name_trans_bad          = cell(size(layer_bad_sort, 1), 2);
    for ii = 1:size(layer_bad_sort, 1)
        if layer_bad_sort(ii, 3)
            name_trans_bad{ii, 1} ...
                            = [name_trans{layer_bad_sort(ii, 1)}{layer_bad_sort(ii, 2)} letters(layer_bad_sort(ii, 3))];
        else
            name_trans_bad{ii, 1} ...
                            = name_trans{layer_bad_sort(ii, 1)}{layer_bad_sort(ii, 2)};
        end
    end
    figure('name', 'Overturned age table')
    uitable('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'fontsize', 18, 'columnname', {'Name' 'Year' 'Transect' 'Sub-transect' 'Layer 1 #' 'Layer 2 #' 'Age difference (ka)' 'Date order 1 #' 'Date order 2 #' 'Master bin #1' 'Master bin #2'}, ...
            'data', [name_trans_bad(:, 1) mat2cell(layer_bad_sort, ones(size(layer_bad_sort, 1), 1), ones(1, 10))]);
end