function age_test           = age_check(age, curr_layer, age_candidate, depth, curr_name)
% AGE_CHECK Test whether a reflector age would lead to overturning prior to assigning it.
% 
% To be used within DATE_LAYERS only.
% 
% Joe MacGregor
% Last updated: 05/02/14

% assume candidate age will pass
age_test                    = true;

% skip if only one layer or no dated layers
if ((length(age) == 1) || isempty(find(~isnan(age), 1)))
    return
end

% assign candidate age for testing
age(curr_layer)             = age_candidate;

% determine dated layers
ind_layer_good              = find(~isnan(age));

if (length(ind_layer_good) == 1)
    return
end

% index of current layer
ind_curr_layer              = find(ind_layer_good == curr_layer);

% trim age to dated layers
age                         = age(ind_layer_good);

% sort layers and dated layers by increasing depth only where current layer exists
[depth_sort, ind_depth_sort]= sort(depth(ind_layer_good, ~isnan(depth(curr_layer, :))), 1);

% sort ages and NaN out empty ones
age_sort                    = age(ind_depth_sort);
age_sort(isnan(depth_sort)) = NaN;

% difference between sorted ages column-wise
age_diff                    = diff(age_sort);

% review if negative/zero age differences found
if ~isempty(find((age_diff <= 0), 1))
    
    [i_bad, j_bad]          = find(age_diff <= 0); % row/column indices of overturns
    
    % check if current layer is part of the overturn
    if ~isempty(find((ind_depth_sort(sub2ind(size(ind_depth_sort), [i_bad; min([(i_bad + 1) (size(ind_depth_sort, 1) .* ones(length(j_bad), 1))], [], 2)], repmat(j_bad, 2, 1))) == ind_curr_layer), 1)) 
        disp(['!!! ' curr_name ' #' num2str(curr_layer) ' would have been age-overturned so not dated here.'])
        age_test            = false;
        return
    end
end