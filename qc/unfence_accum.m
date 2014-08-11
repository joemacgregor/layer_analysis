function unfence_accum(year, trans, subtrans, layer)
% UNFENCE_ACCUM Unfence a particular AR transect or layer of any matches.
%
% Joe MacGregor (UTIG)
% Last updated: 04/17/14

load mat/id_layer_master_accum id_layer_master_mat id_layer_master_cell
load mat/xy_all_accum name_year name_trans

letters                     = 'a':'m';

if (nargin == 3)
    if ~subtrans
        disp(['Unfencing all matches with transect ' name_trans{year}{trans} '...'])
    else
        disp(['Unfencing all matches with transect ' name_trans{year}{trans} letters(subtrans) '...'])
    end
else
    if ~subtrans
        disp(['Unfencing all matches with layer #' num2str(layer) ' from transect ' name_trans{year}{trans} '...'])
    else
        disp(['Unfencing all matches with layer #' num2str(layer) ' from transect ' name_trans{year}{trans} letters(subtrans) '...'])
    end
end

if (nargin == 3)
    ind2unfence{1}          = find((id_layer_master_mat(:, 1) == year) & (id_layer_master_mat(:, 2) == trans) & (id_layer_master_mat(:, 3) == subtrans));
    ind2unfence{2}          = find((id_layer_master_mat(:, 5) == year) & (id_layer_master_mat(:, 6) == trans) & (id_layer_master_mat(:, 7) == subtrans));
else
    ind2unfence{1}          = find((id_layer_master_mat(:, 1) == year) & (id_layer_master_mat(:, 2) == trans) & (id_layer_master_mat(:, 3) == subtrans) & (id_layer_master_mat(:, 4) == layer));
    ind2unfence{2}          = find((id_layer_master_mat(:, 5) == year) & (id_layer_master_mat(:, 6) == trans) & (id_layer_master_mat(:, 7) == subtrans) & (id_layer_master_mat(:, 8) == layer));
end

if isempty([ind2unfence{1}; ind2unfence{2}])
    disp('No matches to unfence.')
    return
end

trans2unfence               = unique([id_layer_master_mat(ind2unfence{1}, 5:7); id_layer_master_mat(ind2unfence{2}, 1:3)], 'rows');

if (nargin == 3)
    if ~subtrans
        id_layer_master_cell{year}{trans} ...
                            = [];
    else
        id_layer_master_cell{year}{trans}{subtrans} ...
                            = [];
    end
    for ii = 1:size(trans2unfence, 1)
        if ~trans2unfence(ii, 3)
            id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} ...
                            = id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(setdiff(1:size(id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}, 1), ...
                              find((id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(:, 2) == year) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(:, 3) == trans) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(:, 4) == subtrans))), :);
        else
            id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)} ...
                            = id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}...
                              (setdiff(1:size(id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}, 1), ...
                              find((id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}(:, 2) == year) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}(:, 3) == trans) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}(:, 4) == subtrans))), :);
        end
    end
else
    if ~subtrans
        id_layer_master_cell{year}{trans} ...
                            = id_layer_master_cell{year}{trans}(setdiff(1:size(id_layer_master_cell{year}{trans}, 1), find(id_layer_master_cell{year}{trans}(:, 1) == layer)), :);
    else
        id_layer_master_cell{year}{trans}{subtrans} ...
                            = id_layer_master_cell{year}{trans}{subtrans}(setdiff(1:size(id_layer_master_cell{year}{trans}{subtrans}, 1), find(id_layer_master_cell{year}{trans}{subtrans}(:, 1) == layer)), :);
    end
    for ii = 1:size(trans2unfence, 1)
        if ~trans2unfence(ii, 3)
            id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} ...
                            = id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(setdiff(1:size(id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}, 1), ...
                              find((id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(:, 2) == year) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(:, 3) == trans) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(:, 4) == subtrans) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}(:, 5) == layer))), :);
        else
            id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)} ...
                            = id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}...
                              (setdiff(1:size(id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}, 1), ...
                              find((id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}(:, 2) == year) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}(:, 3) == trans) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}(:, 4) == subtrans) & ...
                                   (id_layer_master_cell{trans2unfence(ii, 1)}{trans2unfence(ii, 2)}{trans2unfence(ii, 3)}(:, 5) == layer))), :);
        end
    end
end

id_layer_master_mat         = id_layer_master_mat(setdiff(1:size(id_layer_master_mat, 1), [ind2unfence{1}; ind2unfence{2}]), :);

disp(['Unfenced ' num2str(length([ind2unfence{1}; ind2unfence{2}])) ' matches from id_layer_master.mat...']);

save mat/id_layer_master_accum -v7.3 id_layer_master_mat id_layer_master_cell

disp('Now unfencing transect''s merged picks...')

if ~subtrans
    load([name_year{year} '/merge/' name_trans{year}{trans} '_pk_merge'], 'pk')    
else
    load([name_year{year} '/merge/' name_trans{year}{trans} letters(subtrans) '_pk_merge'], 'pk')
end

if (nargin == 3)
    pk.ind_layer            = [];
else
    pk.ind_layer            = pk.ind_layer(setdiff(1:size(pk.ind_layer, 1), find(pk.ind_layer(:, 1) == layer)), :); %#ok<NODEF>
end

if ~subtrans
    save([name_year{year} '/merge/' name_trans{year}{trans} '_pk_merge'], '-v7.3', 'pk')
else
    save([name_year{year} '/merge/' name_trans{year}{trans} letters(subtrans) '_pk_merge'], '-v7.3', 'pk')
end

disp('Now unfencing merged picks files fenced with transect...')

for ii = 1:size(trans2unfence, 1)
    if ~trans2unfence(ii, 3)
        disp([name_trans{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} '...'])
        load([name_year{trans2unfence(ii, 1)} '/merge/' name_trans{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} '_pk_merge'], 'pk')
    else
        disp([name_trans{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} letters(trans2unfence(ii, 3)) '...'])
        load([name_year{trans2unfence(ii, 1)} '/merge/' name_trans{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} letters(trans2unfence(ii, 3)) '_pk_merge'], 'pk')
    end
    if (nargin == 3)
        pk.ind_layer        = pk.ind_layer(setdiff(1:size(pk.ind_layer, 1), find((pk.ind_layer(:, 2) == year) & (pk.ind_layer(:, 3) == trans) & (pk.ind_layer(:, 4) == subtrans))), :);
    else
        pk.ind_layer        = pk.ind_layer(setdiff(1:size(pk.ind_layer, 1), find((pk.ind_layer(:, 2) == year) & (pk.ind_layer(:, 3) == trans) & (pk.ind_layer(:, 4) == subtrans) & (pk.ind_layer(:, 5) == layer))), :);
    end
    if ~trans2unfence(ii, 3)
        save([name_year{trans2unfence(ii, 1)} '/merge/' name_trans{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} '_pk_merge'], '-v7.3', 'pk')
    else
        save([name_year{trans2unfence(ii, 1)} '/merge/' name_trans{trans2unfence(ii, 1)}{trans2unfence(ii, 2)} letters(trans2unfence(ii, 3)) '_pk_merge'], '-v7.3', 'pk')
    end
end

disp('DONE unfencing.')