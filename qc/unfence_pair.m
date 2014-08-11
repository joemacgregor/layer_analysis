function unfence_pair(year1, trans1, subtrans1, layer1, year2, trans2, subtrans2, layer2)
% UNFENCE_PAIR Unfence a particular transect/layer pair.
%
% Joe MacGregor (UTIG)
% Last updated: 03/30/14

if (nargin ~= 8)
    error('unfence_pair:nargin', 'Incorrect number of input arguments.')
end

load mat/id_layer_master id_layer_master_mat id_layer_master_cell
load mat/xy_all name_year name_trans

letters                     = 'a':'z';

ind2unfence                 = [find(ismember(id_layer_master_mat(:, 1:4), [year1 trans1 subtrans1 layer1], 'rows') & ismember(id_layer_master_mat(:, 5:8), [year2 trans2 subtrans2 layer2], 'rows')); ...
                               find(ismember(id_layer_master_mat(:, 5:8), [year1 trans1 subtrans1 layer1], 'rows') & ismember(id_layer_master_mat(:, 1:4), [year2 trans2 subtrans2 layer2], 'rows'))];

if isempty(ind2unfence)
    disp('No matches to unfence.')
    return
end

id_layer_master_mat         = id_layer_master_mat(setdiff(1:size(id_layer_master_mat, 1), ind2unfence), :);

if ~subtrans1
    id_layer_master_cell{year1}{trans1} ...
                            = id_layer_master_cell{year1}{trans1}(setdiff(1:size(id_layer_master_cell{year1}{trans1}, 1), ...
                                                                  find(ismember(id_layer_master_cell{year1}{trans1}(:, 1:5), [layer1 year2 trans2 subtrans2 layer2], 'rows'))), :);
else
    id_layer_master_cell{year1}{trans1}{subtrans1} ...
                            = id_layer_master_cell{year1}{trans1}{subtrans1}(setdiff(1:size(id_layer_master_cell{year1}{trans1}{subtrans1}, 1), ...
                                                                             find(ismember(id_layer_master_cell{year1}{trans1}{subtrans1}(:, 1:5), [layer1 year2 trans2 subtrans2 layer2], 'rows'))), :);
end
if ~subtrans2
    id_layer_master_cell{year2}{trans2} ...
                            = id_layer_master_cell{year2}{trans2}(setdiff(1:size(id_layer_master_cell{year2}{trans2}, 1), ...
                                                                  find(ismember(id_layer_master_cell{year2}{trans2}(:, 1:5), [layer2 year1 trans1 subtrans1 layer1], 'rows'))), :);
else
    id_layer_master_cell{year2}{trans2}{subtrans2} ...
                            = id_layer_master_cell{year2}{trans2}{subtrans2}(setdiff(1:size(id_layer_master_cell{year2}{trans2}{subtrans2}, 1), ...
                                                                             find(ismember(id_layer_master_cell{year2}{trans2}{subtrans2}(:, 1:5), [layer2 year1 trans1 subtrans1 layer1], 'rows'))), :);
end

disp(['Unfenced ' num2str(length(ind2unfence)) ' matches from id_layer_master.mat...']);

save mat/id_layer_master -v7.3 id_layer_master_mat id_layer_master_cell

disp('Now unfencing transect''s merged picks...')

if ~subtrans1
    load([name_year{year1} '/merge/' name_trans{year1}{trans1} '_pk_merge'], 'pk')    
else
    load([name_year{year1} '/merge/' name_trans{year1}{trans1} letters(subtrans1) '_pk_merge'], 'pk')
end

pk.ind_layer            = pk.ind_layer(setdiff(1:size(pk.ind_layer, 1), find(ismember(pk.ind_layer(:, 1:5), [layer1 year2 trans2 subtrans2 layer2], 'rows'))), :); %#ok<NODEF>

if ~subtrans1
    save([name_year{year1} '/merge/' name_trans{year1}{trans1} '_pk_merge'], '-v7.3', 'pk')
else
    save([name_year{year1} '/merge/' name_trans{year1}{trans1} letters(subtrans1) '_pk_merge'], '-v7.3', 'pk')
end

if ~subtrans2
    load([name_year{year2} '/merge/' name_trans{year2}{trans2} '_pk_merge'], 'pk')    
else
    load([name_year{year2} '/merge/' name_trans{year2}{trans2} letters(subtrans2) '_pk_merge'], 'pk')
end

pk.ind_layer            = pk.ind_layer(setdiff(1:size(pk.ind_layer, 1), find(ismember(pk.ind_layer(:, 1:5), [layer2 year1 trans1 subtrans1 layer1], 'rows'))), :);

if ~subtrans2
    save([name_year{year2} '/merge/' name_trans{year2}{trans2} '_pk_merge'], '-v7.3', 'pk')
else
    save([name_year{year2} '/merge/' name_trans{year2}{trans2} letters(subtrans2) '_pk_merge'], '-v7.3', 'pk')
end

disp('DONE unfencing.')