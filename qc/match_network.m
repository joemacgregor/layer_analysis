% MATCH_NETWORK Network graph analysis of layer matches.
% 
% Joe MacGregor (UTIG)
% Last updated: 04/07/14

clear

plotting                    = false;

load mat/xy_all num_trans num_year name_trans name_year
load mat/merge_all depth_smooth ind_fence num_fence num_layer
load mat/id_layer_master id_layer_master_mat

colors_def                  = [0    0       0.75;
                               0    0       1;
                               0    0.25    1;
                               0    0.50    1;
                               0    0.75    1;
                               0    1       1;
                               0.25 1       0.75;
                               0.50 1       0.50;
                               0.75 1       0.25;
                               1    1       0;
                               1    0.75    0;
                               1    0.50    0;
                               1    0.25    0;
                               1    0       0;
                               0.75 0       0];

% load mat/merge_all depth_smooth num_layer
% load mat/id_layer_master_auto id_layer_master_mat ind_fence num_fence

% order in which to analyze years/campaigns based on their overall data quality
year_ord                    = [17 19 6 7 8 4 9 5 2 3 12 1 13 15 16 18 11 10 14];

% do not permit matches between campaigns at far end of period from each other, i.e., ICORDS cannot be matched with any system after ICORDS v2
year_match                  = true(num_year);
year_match(1:4, 9:end)      = false;
year_match(9:end, 1:4)      = false;

age_match                   = cell(1, num_year);
for ii = 1:num_year
    if ~any(num_fence{ii})
        continue
    end
    age_match{ii}           = cell(1, num_trans(ii));
    for jj = find(num_fence{ii})
        age_match{ii}{jj}   = cell(1, length(ind_fence{ii}{jj}));
        for kk = find(ind_fence{ii}{jj})
            age_match{ii}{jj}{kk} ...
                            = NaN(num_layer{ii}{jj}(kk), 1);
        end
    end
end

% clean up matching list
for ii = 1:size(id_layer_master_mat, 1)
    if ((id_layer_master_mat(ii, 1) > id_layer_master_mat(ii, 5)) || ((id_layer_master_mat(ii, 1) == id_layer_master_mat(ii, 5)) && (id_layer_master_mat(ii, 2) > id_layer_master_mat(ii, 6))) || ...
       ((id_layer_master_mat(ii, 1) == id_layer_master_mat(ii, 5)) && (id_layer_master_mat(ii, 2) == id_layer_master_mat(ii, 6)) && (id_layer_master_mat(ii, 3) > id_layer_master_mat(ii, 7))) || ...
       ((id_layer_master_mat(ii, 1) == id_layer_master_mat(ii, 5)) && (id_layer_master_mat(ii, 2) == id_layer_master_mat(ii, 6)) && (id_layer_master_mat(ii, 3) == id_layer_master_mat(ii, 7)) && (id_layer_master_mat(ii, 4) > id_layer_master_mat(ii, 8))))
         id_layer_master_mat(ii, :) ...
                            = id_layer_master_mat(ii, [5:8 1:4]);
    end
end
id_layer_master_mat         = unique(id_layer_master_mat, 'rows');

%%

disp('Assigning master IDs within campaigns...')

% intial master id
master_id_curr              = 1;
layer_bin                   = {};

% assign master ids within individual campaign first
for ii = year_ord
    
    % all matches within current campaign
    match_curr              = id_layer_master_mat(ismember(id_layer_master_mat(:, [1 5]), [ii ii], 'rows'), :);
    match_curr_cell         = match_curr;
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
            master_id_curr  = master_id_curr + 1;
            
        else
            
            % current matching IDs
            age_match_curr  = sort([age_match{match_curr(jj, 1)}{match_curr(jj, 2)}{match_curr_cell(jj, 3)}(match_curr(jj, 4)) age_match{match_curr(jj, 5)}{match_curr(jj, 6)}{match_curr_cell(jj, 7)}(match_curr(jj, 8))]);
            
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
num_bin                     = zeros(1, length(layer_bin));
for ii = 1:length(layer_bin)
    num_bin(ii)             = size(layer_bin{ii}, 1);
end

% bin descending order
[~, ind_bin_ord]            = sort(num_bin, 'descend');

% adjust age_match based on bin totals
for ii = 1:length(layer_bin)
    for jj = 1:num_bin(ind_bin_ord(ii))
        age_match{layer_bin{ind_bin_ord(ii)}(jj, 1)}{layer_bin{ind_bin_ord(ii)}(jj, 2)}{max([1 layer_bin{ind_bin_ord(ii)}(jj, 3)])}(layer_bin{ind_bin_ord(ii)}(jj, 4)) ...
                            = ii;
    end
end

% trim bins
layer_bin                   = layer_bin(ind_bin_ord);
layer_bin                   = layer_bin(logical(num_bin(ind_bin_ord)));
num_bin                     = num_bin(ind_bin_ord);
num_bin                     = num_bin(logical(num_bin)); %#ok<NASGU>

% reset new master id
master_id_curr              = length(layer_bin) + 1;

%%

disp('Assigning master IDs between different campaigns...')

% assign master ids for matches between different campaigns
for ii = year_ord(1:(end - 1))
    
    disp([name_year{ii} ' ' num2str(find(year_ord == ii)) '/' num2str(num_year - 1)])
    
    for jj = year_ord((find(ii == year_ord) + 1):end)
        
        % skip if direct matching between campaign pair is not allowed
        if ~year_match(ii, jj)
            continue
        end
        
        match_curr          = id_layer_master_mat((ismember(id_layer_master_mat(:, [1 5]), [ii jj], 'rows') | ismember(id_layer_master_mat(:, [1 5]), [jj ii], 'rows')), 1:8);
        match_curr_cell     = match_curr;
        match_curr_cell(~match_curr_cell) ...
                            = 1;
        
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

num_bin                     = zeros(1, length(layer_bin));
for ii = 1:length(layer_bin)
    num_bin(ii)             = size(layer_bin{ii}, 1);
end

[~, ind_bin_ord]            = sort(num_bin, 'descend');

for ii = 1:length(layer_bin)
    for jj = 1:num_bin(ind_bin_ord(ii))
        age_match{layer_bin{ind_bin_ord(ii)}(jj, 1)}{layer_bin{ind_bin_ord(ii)}(jj, 2)}{max([1 layer_bin{ind_bin_ord(ii)}(jj, 3)])}(layer_bin{ind_bin_ord(ii)}(jj, 4)) ...
                            = ii;
    end
end

layer_bin                   = layer_bin(ind_bin_ord);
layer_bin                   = layer_bin(logical(num_bin(ind_bin_ord)));
num_bin                     = num_bin(ind_bin_ord);
num_bin                     = num_bin(logical(num_bin));

%%

disp('Evaluating within each bin for overlapping layers...')

% cell of logical vectors with overlap evaluation for each transect
layer_overlap               = cell(1, length(layer_bin));

% evaluate all bins
for ii = 1:length(layer_bin)
    
    % assume no layer overlap within current bin
    layer_overlap{ii}       = zeros(num_bin(ii), 1);
    mm                      = 1;
    
    % transects represented in current bin
    trans_unique            = unique(layer_bin{ii}(:, 1:3), 'rows');
    
    % evaluate all transects within bin
    for jj = 1:size(trans_unique, 1)
        
        % skip if only one layer from that transect within bin
        if (length(find(ismember(layer_bin{ii}(:, 1:3), trans_unique(jj, :), 'rows'))) == 1)
            continue
        end
        
        % indices within layer_bin that from current transect
        ind_bin_risk        = find(ismember(layer_bin{ii}(:, 1:3), trans_unique(jj, :), 'rows'));
        
        ind_bin_bad         = [];
        for kk = 1:(length(ind_bin_risk) - 1)
            ll              = kk + 1;
            while ((ll <= length(ind_bin_risk)) && ~any(ll == ind_bin_bad))
                diff_depth_smooth_curr ...
                            = diff(depth_smooth{trans_unique(jj, 1)}{trans_unique(jj, 2)}{max([1 trans_unique(jj, 3)])}(layer_bin{ii}(ind_bin_risk([kk ll]), 4), :));
                if ~isempty(find((~isnan(diff_depth_smooth_curr) & (diff_depth_smooth_curr ~= 0)), 1))
                    ind_bin_bad ...
                            = [ind_bin_bad; kk; ll]; %#ok<AGROW>
                    break
                end
                ll          = ll + 1;
            end
        end
        
        if isempty(ind_bin_bad)
            continue
        end
        
        layer_overlap{ii}(ind_bin_risk(ind_bin_bad)) ...
                            = mm;
        mm                  = mm + 1;
        if (mm > size(colors_def, 1))
            mm              = 1;
        end
    end
end

num_overlap                 = zeros(1, length(layer_bin));
for ii = 1:length(layer_bin)
    num_overlap(ii)         = length(find(layer_overlap{ii}));
end

%%

disp('Evaluating between each bin for overlapping layers...')

bin_overlap                 = [];

% loop through all bins except last
for ii = 1:(length(layer_bin) - 1)
    
    disp(ii)
    
    % transects represented in current bin1
    trans_unique1           = unique(layer_bin{ii}(:, 1:3), 'rows');
    
    % loop through all bins after current bin1
    for jj = (ii + 1):length(layer_bin)
        
        trans_unique2       = unique(layer_bin{jj}(:, 1:3), 'rows');
        
        sign_test           = NaN(size(layer_bin{ii}, 1), size(layer_bin{jj}, 1));
        
        % loop through all transects in bin1
        for kk = 1:size(trans_unique1, 1)
            
            % skip if transect is not in both bin1 and bin2
            if isempty(find(ismember(trans_unique2, trans_unique1(kk, :), 'rows'), 1))
                continue
            end
            
            % indices within each bin under consideration for each transect
            ind_bin_risk1   = find(ismember(layer_bin{ii}(:, 1:3), trans_unique1(kk, :), 'rows'));
            ind_bin_risk2   = find(ismember(layer_bin{jj}(:, 1:3), trans_unique1(kk, :), 'rows'));
            
            % loop through risky layers in bin1
            for ll = 1:length(ind_bin_risk1)
                
                % loop through risky layers in bin2
                for mm = 1:length(ind_bin_risk2)
                    
                    % sign of difference between each layer pair
                    sign_diff_depth_smooth ...
                            = sign(diff(depth_smooth{trans_unique1(kk, 1)}{trans_unique1(kk, 2)}{max([1 trans_unique1(kk, 3)])}([layer_bin{ii}(ind_bin_risk1(ll), 4) layer_bin{jj}(ind_bin_risk2(mm), 4)], :)));
                    
                    % skip if layer pair doesn't overlap
                    if isempty(find(~isnan(sign_diff_depth_smooth), 1))
                        continue
                    end
                    
                    % reduce to overlapping portions
                    sign_diff_depth_smooth ...
                            = sign_diff_depth_smooth(~isnan(sign_diff_depth_smooth));
                    
                    % reversal found in sign of difference of a layer pair
                    if (~isempty(find(diff(sign_diff_depth_smooth), 1)) && ~all(sign_diff_depth_smooth == 0))
                        bin_overlap ...
                            = [bin_overlap; ii ind_bin_risk1(ll) jj ind_bin_risk2(mm)]; %#ok<AGROW>
                    else
                        sign_test(ind_bin_risk1(ll), ind_bin_risk2(mm)) ...
                            = sign_diff_depth_smooth(1);
                    end
                end
            end
        end
        
        % most common sign difference
        mode_curr           = mode(sign_test(~isnan(sign_test)));
        
        % skip if no inconsistent sign differences
        if isempty(find(((sign_test ~= mode_curr) & ~isnan(sign_test)), 1))
            continue
        end
        
        % bin indices of errant signs
        [ind_bin_bad1, ind_bin_bad2] ...
                            = ind2sub([size(layer_bin{ii}, 1) size(layer_bin{jj}, 1)], find(((sign_test ~= mode_curr) & ~isnan(sign_test))));
        bin_overlap         = [bin_overlap; ii(ones(length(ind_bin_bad1), 1), 1) ind_bin_bad1 jj(ones(length(ind_bin_bad1), 1), 1) ind_bin_bad2]; %#ok<AGROW>
    end
end

bin_overlap                 = unique(bin_overlap, 'rows');


%%
bin_overlap_simple          = zeros(size(bin_overlap, 1), 8);
for ii = 1:size(bin_overlap, 1)
    bin_overlap_simple(ii, :) ...
                            = [layer_bin{bin_overlap(ii, 1)}(bin_overlap(ii, 2), :) layer_bin{bin_overlap(ii, 3)}(bin_overlap(ii, 4), :)];
end
% bin_overlap_simple          = sortrows(bin_overlap_simple, [1:3 5:7 4 8]);

%%

% unique layers that are matched
curr_bin                    = 1;
layer_match                 = unique(layer_bin{curr_bin}, 'rows');
num_layer_match             = size(layer_match, 1);

% simplify matching set to only those relevant for these layers
id_layer_master_simple      = NaN(size(id_layer_master_mat, 1), 2);
for ii = 1:num_layer_match
    if ~mod(ii, 100)
        disp(ii)
    end
    id_layer_master_simple(ismember(id_layer_master_mat(:, 1:4), layer_match(ii, :), 'rows'), 1) ...
                            = ii;
    id_layer_master_simple(ismember(id_layer_master_mat(:, 5:8), layer_match(ii, :), 'rows'), 2) ...
                            = ii;
end
id_layer_master_simple      = id_layer_master_simple((~isnan(id_layer_master_simple(:, 1)) & ~isnan(id_layer_master_simple(:, 2))), :);

% sparse matching network and graph
layer_match_net             = sparse(id_layer_master_simple(:, 1), id_layer_master_simple(:, 2), ones(size(id_layer_master_simple, 1), 1), num_layer_match, num_layer_match);
layer_match_str             = cell(num_layer_match, 1);

layer_match_graph           = biograph(layer_match_net, num2str((1:num_layer_match)'), 'layouttype', 'radial', 'showarrows', 'off', 'edgecallback', @do_break, 'nodeautosize', 'off');
% layer_match_graph           = biograph(layer_match_net, num2str((1:num_layer_match)'), 'layouttype', 'radial', 'showarrows', 'off', 'showtextinnodes', 'ID', 'edgecallback', @do_fence_snapshot);
for ii = 1:num_layer_match
    layer_match_graph.Nodes(ii).Label ...
                            = ' ';%num2str(layer_match(ii, :));
    layer_match_graph.Nodes(ii).FontSize ...
                            = 12;
    if ~isempty(find(layer_match_net(ii, :), 1))
        layer_match_graph.Nodes(ii).Size ...
                            = sqrt(length(find(layer_match_net(ii, :)))) .* [10 10];
    else
        layer_match_graph.Nodes(ii).Size ...
                            = sqrt(length(find(layer_match_net(:, ii)))) .* [10 10];
    end
    layer_match_graph.Nodes(ii).Description ...
                            = [num2str(layer_match(ii, :)) ' / ' idtrans(layer_match(ii, 1), layer_match(ii, 2), layer_match(ii, 3)) ' #' num2str(layer_match(ii, 4))];
    if layer_overlap{curr_bin}(ii)
        layer_match_graph.Nodes(ii).Color ...
                            = colors_def(layer_overlap{curr_bin}(ii), :);
        if any(layer_overlap{curr_bin}(ii) == [1 2 3])
            layer_match_graph.Nodes(ii).TextColor ...
                            = [1 1 1];
        end
    else
        layer_match_graph.Nodes(ii).Color ...
                            = [1 1 1];
%         layer_match_graph.Nodes(ii).Label...
%                             = ' ';
    end
end
for ii = 1:size(id_layer_master_simple, 1)
    layer_match_graph.Edges(ii).UserData ...
                            = [layer_match(id_layer_master_simple(ii, 1), :); layer_match(id_layer_master_simple(ii, 2), :)];
end

view(layer_match_graph)

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%%
end