function [age, age_ord, age_n, age_range, age_type, age_uncert, date_counter, ind_layer_blacklist] ...
                            = date_interp(depth, thick, age, age_ord, age_n, age_range, age_type, age_uncert, num_layer, ind_layer_undated, ind_layer_blacklist, curr_name, interp_type, parallel_check, num_pool, date_counter, strain_rate, age_uncert_rel_max, thick_diff_max, layer_diff_max, tol, ...
                                          iter_max, age_max, do_age_check, age_type_assign)
% DATE_INTERP Date an undated layer.
%
% To be used within DATE_LAYERS only.
%
% Joe MacGregor
% Last updated: 08/18/14

% determine which traces are potentially usable because they have at least two dated layers
ind_trace_usable            = find(sum(~isnan(depth(~isnan(age), :)), 1) > 1);

% stop if insufficient usable traces for overlapping dating
if isempty(ind_trace_usable)
    ind_layer_blacklist     = find(isnan(age))';
    return
end

% loop through undated layers to find number of overlapping traces with dated layers
num_overlap                 = zeros(num_layer, 1);
for ii = intersect(find(isnan(age))', setdiff(ind_layer_undated, ind_layer_blacklist), 'stable')
    num_overlap(ii)         = length(find(~isnan(depth(ii, ind_trace_usable))));
end

% stop if there are no more undated layers that overlap with dated layers
if ~any(num_overlap > 1)
    ind_layer_blacklist     = find(isnan(age))';
    return
end

% layer with most overlap
[~, curr_layer]             = max(num_overlap);

% along-transect overlapping indices of layer with most overlap
ind_overlap                 = ind_trace_usable(~isnan(depth(curr_layer, ind_trace_usable)));
num_overlap_max             = length(ind_overlap);

% prepare for loop iteration by slicing all necessary variables
[age_slice, age_uncert_slice, depth_slice] ...
                            = deal(NaN(2, num_overlap_max));

% initially include all depth and then whittle them down to the bounding layers
[age_slice_tmp, age_uncert_slice_tmp, depth_slice_tmp, depth_curr, ind_overlap_bounded, thick_diff_max_curr] ...
                            = deal(age(~isnan(age)), age_uncert(~isnan(age_uncert)), depth(~isnan(age), ind_overlap), depth(curr_layer, ind_overlap), false(1, num_overlap_max), (thick(ind_overlap) .* thick_diff_max));

% sort dated layers by increasing depth
[depth_slice_tmp, ind_depth_ord] ...
                            = sort(depth_slice_tmp, 1);
[age_slice_tmp, age_uncert_slice_tmp] ...
                            = deal(age_slice_tmp(:, ones(1, num_overlap_max)), age_uncert_slice_tmp(:, ones(1, num_overlap_max)));
[age_slice_tmp, age_uncert_slice_tmp] ...
                            = deal(age_slice_tmp(ind_depth_ord), age_uncert_slice_tmp(ind_depth_ord));

% trim slices to upper/lower bounds
for ii = 1:num_overlap_max
    
    % upper/lower bounding layer indices
    [ind_top, ind_bot]      = deal(find((depth_slice_tmp(:, ii) < depth_curr(ii)), 2, 'last'), find((depth_slice_tmp(:, ii) > depth_curr(ii)), 2));
    
    % skip if no bounding dated reflectors
    if (isempty(ind_top) && isempty(ind_bot))
        continue
    end
        
    % bin ages and depths that survive criteria if unbounded above or below (not both)
    if (isempty(ind_top) && (length(ind_bot) == 2))
        if all((depth_slice_tmp(ind_bot(1), ii) - depth_curr(ii)) <= [(layer_diff_max * diff(depth_slice_tmp(ind_bot, ii))) thick_diff_max_curr(ii)])
            [age_slice(:, ii), age_uncert_slice(:, ii), depth_slice(:, ii), ind_overlap_bounded(ii)] ...
                            = deal(age_slice_tmp(ind_bot, ii), age_uncert_slice_tmp(ind_bot, ii), depth_slice_tmp(ind_bot, ii), true);
        end
    elseif (isempty(ind_bot) && (length(ind_top) == 2))
        if all((depth_curr(ii) - depth_slice_tmp(ind_top(2), ii)) <= [(layer_diff_max * diff(depth_slice_tmp(ind_top, ii))) thick_diff_max_curr(ii)])
            [age_slice(:, ii), age_uncert_slice(:, ii), depth_slice(:, ii), ind_overlap_bounded(ii)] ...
                            = deal(age_slice_tmp(ind_top, ii), age_uncert_slice_tmp(ind_top, ii), depth_slice_tmp(ind_top, ii), true);
        end
    elseif (~isempty(ind_top) && ~isempty(ind_bot))
        [age_slice(:, ii), age_uncert_slice(:, ii), depth_slice(:, ii), ind_overlap_bounded(ii)] ...
                            = deal(age_slice_tmp([ind_top(end); ind_bot(1)], ii), age_uncert_slice_tmp([ind_top(end); ind_bot(1)], ii), depth_slice_tmp([ind_top(end); ind_bot(1)], ii), true);
    end
end

% skip traces with bad age bounds
if ~isempty(find((diff(age_slice) <= 0), 1))
    ind_overlap_bounded(diff(age_slice) <= 0) ...
                            = false;
    disp([curr_name ' layer #' num2str(curr_layer) ' is sometimes badly age-bounded (' num2str(1e-3 * mean(diff(age_slice(:, (diff(age_slice) <= 0))))) ' ka mean top/bottom difference).'])
end

% remove unbounded or poorly dated traces
if isempty(find(ind_overlap_bounded, 1))
    ind_layer_blacklist     = [ind_layer_blacklist curr_layer];
    return
elseif any(~ind_overlap_bounded)
    [age_slice, age_uncert_slice, depth_slice, ind_overlap, depth_curr] ...
                            = deal(age_slice(:, ind_overlap_bounded), age_uncert_slice(:, ind_overlap_bounded), depth_slice(:, ind_overlap_bounded), ind_overlap(ind_overlap_bounded), depth_curr(ind_overlap_bounded));
    num_overlap_max         = length(ind_overlap);
end
age_overlap                 = NaN(1, num_overlap_max);

% date layer using either linear or quasi-Nye vertical strain rate method
switch interp_type
    
    case 'lin'
        
        % parallelized loop through each trace where current layer exists, linearly inter/extrapolating age
        if (parallel_check && (num_overlap_max > (3 * num_pool)))
            parfor ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = interp1(depth_slice(:, ii), age_slice(:, ii), depth_curr(ii), 'linear', 'extrap');
            end
        else
            % non-parallel loop doing same thing
            for ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = interp1(depth_slice(:, ii), age_slice(:, ii), depth_curr(ii), 'linear', 'extrap');
            end
        end
        
    case 'quasi Nye'
        
        % quasi-Nye dating using bounding or sufficiently near dated reflectors
        if (parallel_check && (num_overlap_max > (3 * num_pool)))
            strain_rate_bounded ...
                            = strain_rate(ind_overlap);
            parfor ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = date_strain(depth_slice(:, ii), age_slice(:, ii), strain_rate_bounded(ii), depth_curr(ii), tol, iter_max, 'age');
            end
        else
            for ii = 1:num_overlap_max
                age_overlap(ii) ...
                            = date_strain(depth_slice(:, ii), age_slice(:, ii), strain_rate(ind_overlap(ii)), depth_curr(ii), tol, iter_max, 'age');
            end
        end
end

% sanity checks
for ii = 1:num_overlap_max
    if ~isreal(age_overlap(ii))
        age_overlap(ii)     = NaN;
    end
end
age_overlap(isinf(age_overlap) | (age_overlap < 0) | (age_overlap > age_max) | ((depth_curr > depth_slice(2, :)) & (age_overlap < age_slice(2, :))) | ((depth_curr < depth_slice(1, :)) & (age_overlap > age_slice(1, :)))) ...
                            = NaN;

% weighted mean
age_overlap_mean            = nanmean(age_overlap, 2);

% check if new age causes overturning
if (do_age_check && ~age_check(age, curr_layer, age_overlap_mean, depth, curr_name))
    ind_layer_blacklist     = [ind_layer_blacklist curr_layer];
    return
end

% weighted uncertainty as the sum of uncertainties of both bounding layers
age_overlap_uncert          = sqrt(nanvar(age_overlap) + nanmean(nanmean((age_uncert_slice .^ 2), 1), 2));

% limit age assignment to well-constrained layers
if (isempty(find(~isnan(age_overlap), 1)) || ((age_overlap_uncert / age_overlap_mean) > age_uncert_rel_max) || (age_overlap_mean <= 0))
    ind_layer_blacklist     = [ind_layer_blacklist curr_layer];
    return
end

% breakout age information for newly dated layer
age(curr_layer)             = age_overlap_mean;
age_uncert(curr_layer)      = age_overlap_uncert;
age_ord(curr_layer)         = date_counter;
date_counter                = date_counter + 1;
age_n(curr_layer)           = 1;
age_range(curr_layer)       = range(age_overlap);
age_type(curr_layer)        = age_type_assign;