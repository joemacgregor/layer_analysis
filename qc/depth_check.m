% DEPTH_CHECK Check depths for reversals of sign difference.
% 
% Joe MacGregor
% Last updated: 04/07/14

clear

load mat/xy_all num_trans num_year name_trans name_year
load mat/merge_all depth_smooth ind_fence num_fence num_layer

layer_overlap                 = [];

for ii = 1:num_year
    
    disp(name_year{ii})
        
    for jj = 1:num_trans(ii)
        
        disp(name_trans{ii}{jj})
        
        for kk = find(ind_fence{ii}{jj})
            
            for ll = 1:num_layer{ii}{jj}(kk)
                
                diff_curr   = nanmean(sign(depth_smooth{ii}{jj}{kk} - depth_smooth{ii}{jj}{kk}((ll .* ones(num_layer{ii}{jj}(kk), 1)), :)), 2);
                
                for mm = find(~isnan(diff_curr))'
                    if ~any(diff_curr(mm) == [-1 0 1])
                        layer_overlap ...
                                = [layer_overlap; ii jj kk ll mm]; %#ok<AGROW>
                        disp(layer_overlap(end, :))
                    end
                end
            end
        end
    end
end

save mat/layer_overlap -v7.3 layer_overlap