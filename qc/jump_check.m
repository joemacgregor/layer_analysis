% JUMP_CHECK Evaluate layer jumps in transects.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/07/14

clear

do_jump                     = false;
do_save                     = false;
plotting                    = false;

letters                     = 'a':'z';
dist_jump_max               = 100;

load mat/xy_all name_trans num_trans num_year

%%
if do_jump
%%  
    disp('Evaluating each transect''s layers for jumps...')
    
    load mat/merge_all depth_smooth dist ind_fence num_layer
    layer_bad               = [];
    
    % evaluate each existing layer intersection
    for ii = 1:num_year
        
        for jj = 1:num_trans(ii)
            
            for kk = 1:length(ind_fence{ii}{jj})
                
                for ll = 1:num_layer{ii}{jj}(kk)
                    
                    tmp1    = diff(~isnan(depth_smooth{ii}{jj}{kk}(ll, :)));
                    ind_cut = find(tmp1);
                    
                    if (length(ind_cut) == 1)
                        continue
                    end
                    
                    for mm = 1:(length(ind_cut) - 1)
                        if (isequal(tmp1(ind_cut([mm (mm + 1)])), [-1 1]) && (diff(dist{ii}{jj}{kk}(ind_cut([mm (mm + 1)]))) > dist_jump_max))
                            layer_bad ...
                                = [layer_bad; ii jj kk ll diff(dist{ii}{jj}{kk}(ind_cut([mm (mm + 1)])))]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
%%    
    if do_save
        save mat/layer_bad -v7.3 dist_jump_max layer_bad
    end
    
else
    load mat/layer_bad layer_bad
end
%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%%
    name_trans_bad          = cell(size(layer_bad, 1), 2);
    for ii = 1:size(layer_bad, 1)
        if layer_bad(ii, 3)
            name_trans_bad{ii, 1} ...
                            = [name_trans{layer_bad(ii, 1)}{layer_bad(ii, 2)} letters(layer_bad(ii, 3))];
        else
            name_trans_bad{ii, 1} ...
                            = name_trans{layer_bad(ii, 1)}{layer_bad(ii, 2)};
        end
    end
    figure('name', 'Jumping layer table')
    uitable('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'fontsize', 18, 'columnname', {'Name' 'Year' 'Transect' 'Sub-transect' 'Layer #' 'Distance jump (km)'}, 'data', [name_trans_bad(:, 1) mat2cell(layer_bad, ones(size(layer_bad, 1), 1), [1 1 1 1 1])])
end