function [age_master, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter] ...
                            = match_age(ind_layer_curr, layer_bin, num_core, age_master, age_core, age, age_uncert, age_range, age_n, age_type, age_ord, date_counter, do_core, age_type_max, age_type_assign, do_age_check, name_trans, letters, depth)
% MATCH_AGE Match ages of dated and undated layers.
% 
% To be used within DATE_LAYERS only.
% 
% Joe MacGregor
% Last updated: 05/07/14

% loop through master layers to consider
for ii = ind_layer_curr
    age_cat                 = [];
    id                      = layer_bin{ii};
    id(~id(:, 3), 3)        = 1;
    
    % loop through all dated/fenced matches to test whether they are dated; concatenate them if they are dated
    for jj = 1:size(id, 1)
        if (~isnan(age{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4))) && (age_type{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)) <= age_type_max))
            if do_core
                age_cat     = [age_cat; age_core{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4), :) age_uncert{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4))]; %#ok<AGROW>
            else
                age_cat     = [age_cat; age{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)) age_uncert{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4))]; %#ok<AGROW>
            end
        end
    end
    
    % pass to next iteration if no dated layers
    if isempty(age_cat)
        continue
    end
    
    % break out concatenated dated layers
    age_master_tmp          = NaN(1, (num_core + 6));
    if do_core
        age_master_tmp(1:num_core) ...
                            = nanmean(age_cat(:, 1:num_core), 1);
        age_master_tmp(num_core + 1) ...
                            = nanmean(age_master_tmp(1:num_core));
        age_master_tmp(num_core + 2) ...
                            = nanmean(age_cat(:, (num_core + 1)));                        
        age_master_tmp(num_core + 3) ...
                            = range(age_cat(:, (num_core + 1)));
    else
        age_master_tmp(num_core + 1) ...
                            = nanmean(age_cat(:, 1));
        age_master_tmp(num_core + 2) ...
                            = nanmean(age_cat(:, 2));        
        age_master_tmp(num_core + 3) ...
                            = range(age_cat(:, 1));
    end
    age_master_tmp(num_core + 4) ...
                            = size(age_cat, 1);
    age_master_tmp(num_core + 5) ...
                            = age_type_assign(1);
    age_master_tmp(num_core + 6) ...
                            = date_counter; % dating order increment
    
    % assign best age set to master layer
    age_master(ii, :)       = age_master_tmp;
    
    % reassign ages of fenced layers to that of master layer, depending on status of fenced layer and if age does not create overturning
    for jj = 1:size(id, 1)
        curr_name           = name_trans{id(jj, 1)}{id(jj, 2)};
        if layer_bin{ii}(jj, 3)
            curr_name       = [curr_name letters(id(jj, 3))]; %#ok<AGROW>
        end
        if (do_age_check && ~age_check(age{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}, id(jj, 4), age_master_tmp(num_core + 1), depth{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}, curr_name))
            continue
        end
        if isnan(age{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)))
            age_type{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)) ...
                            = age_type_assign(2);
        else
            age_type{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)) ...
                            = age_type_assign(1);
        end
        if do_core
            age_core{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4), :) ...
                            = age_master_tmp(1:num_core);
        end
        [age{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)), age_uncert{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)), age_range{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)), age_n{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4)), age_ord{id(jj, 1)}{id(jj, 2)}{id(jj, 3)}(id(jj, 4))] ...
                            = deal(age_master_tmp(num_core + 1), age_master_tmp(num_core + 2), age_master_tmp(num_core + 3), age_master_tmp(num_core + 4), age_master_tmp(num_core + 6));
    end
    
    date_counter            = date_counter + 1; % dating order increment
end