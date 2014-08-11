% TRIM_LAYER Trim merged layer.
% 
% Joe MacGregor (UTIG)
% Last updated: 04/07/14

pk                          = load('/Volumes/icebridge/data/2011_p3/merge/20110426_10_pk_merge.mat');
pk                          = pk.pk;

%%

vars                        = {'depth' 'depth_smooth' 'elev' 'elev_gimp' 'elev_smooth' 'elev_smooth_gimp' 'ind_y' 'ind_y_flat_smooth' 'ind_y_smooth' 'int' 'int_smooth' 'twtt' 'twtt_ice' 'twtt_ice_smooth' 'twtt_smooth'};
layer2edit                  = 31;
range2edit                  = '3635:3640';

for ii = 1:length(vars)
    if isfield(pk, vars{ii})
        eval(['pk.' vars{ii} '(layer2edit, ' range2edit ') = NaN(1, length(' range2edit '));'])
    end
end

%%
save('/Volumes/icebridge/data/2011_p3/merge/20110426_10_pk_merge.mat', '-v7.3', 'pk')