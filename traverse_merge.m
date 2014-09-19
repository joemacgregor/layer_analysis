% TRAVERSE_MERGE Extract merged layer data from melt/icebridge/data/.
% 
% Joe MacGregor (UTIG)
% Last updated: 09/19/14

clear

radar_type                  = 'deep';

% load transect data
switch radar_type
    case 'accum'
        load mat/xy_all_accum num_year num_trans name_year name_trans
    case 'deep'
        load mat/xy_all num_year num_trans name_year name_trans
end

% dummy letters
letters                     = 'a':'z';
dist_int_decim              = 1; % decimation index distance, km

% initializations
[depth_smooth, dist, elev_air_gimp, elev_bed_gimp, elev_smooth_gimp, elev_surf_gimp, file_block, ind_decim, ind_decim_mid, ind_fence, ind_layer, ind_trace_start, int_bed, int_smooth, int_surf, num_decim, num_fence, num_layer, num_trace_tot, thick, thick_decim, time, twtt_surf, x_pk, y_pk] ...
                            = deal(cell(1, num_year));

disp('Collecting all merged picks files...')

% load all merged picks files from all years/campaigns
for ii = 1:num_year
    
    disp(['Campaign: ' name_year{ii} '...'])
    
    [depth_smooth{ii}, dist{ii}, elev_air_gimp{ii}, elev_bed_gimp{ii}, elev_smooth_gimp{ii}, elev_surf_gimp{ii}, file_block{ii}, ind_decim{ii}, ind_decim_mid{ii}, ind_fence{ii}, ind_layer{ii}, ind_trace_start{ii}, int_bed{ii}, int_smooth{ii}, int_surf{ii}, num_decim{ii}, num_layer{ii}, ...
     num_trace_tot{ii}, thick{ii}, thick_decim{ii}, time{ii}, twtt_surf{ii}, x_pk{ii}, y_pk{ii}] ...
                            = deal(cell(1, num_trans(ii)));
    num_fence{ii}           = deal(zeros(1, num_trans(ii)));
    
    for jj = 1:num_trans(ii)
        
        if exist([name_year{ii} '/merge/' name_trans{ii}{jj} '_pk_merge.mat'], 'file')
            
            pk              = load([name_year{ii} '/merge/' name_trans{ii}{jj} '_pk_merge.mat']); % load a transect's merge file (no a/b/c/etc, just a single file)
            pk              = pk.pk;
            
            disp(['Transect: ' name_trans{ii}{jj} '...'])
            
            [depth_smooth{ii}{jj}{1}, dist{ii}{jj}{1}, elev_air_gimp{ii}{jj}{1}, elev_bed_gimp{ii}{jj}{1}, elev_smooth_gimp{ii}{jj}{1}, elev_surf_gimp{ii}{jj}{1}, file_block{ii}{jj}{1}, ind_trace_start{ii}{jj}{1}, int_bed{ii}{jj}{1}, int_smooth{ii}{jj}{1}, int_surf{ii}{jj}{1}, ...
             num_layer{ii}{jj}, num_trace_tot{ii}{jj}(1), time{ii}{jj}{1}, twtt_surf{ii}{jj}{1}, x_pk{ii}{jj}{1}, y_pk{ii}{jj}{1}] ...
                            = deal(pk.depth_smooth, pk.dist, pk.elev_air, pk.elev_bed_gimp, pk.elev_smooth_gimp, pk.elev_surf_gimp, pk.file_block, pk.ind_trace_start, pk.int_bed, pk.int_smooth, pk.int_surf, pk.num_layer, pk.num_trace_tot, pk.time, pk.twtt_surf, pk.x, pk.y);
            thick{ii}{jj}{1}= elev_surf_gimp{ii}{jj}{1} - elev_bed_gimp{ii}{jj}{1};
            thick{ii}{jj}{1}(isinf(thick{ii}{jj}{1})) ...
                            = NaN;
            % indices to be kept (decimation vector)
            dist_int        = dist{ii}{jj}{1}(1):dist_int_decim:dist{ii}{jj}{1}(end);
            dist_int(end)   = dist{ii}{jj}{1}(end);
            ind_decim{ii}{jj}{1} ...
                            = unique(interp1(unique(dist{ii}{jj}{1}), 1:length(unique(dist{ii}{jj}{1})), dist_int, 'nearest', 'extrap'));
            ind_decim_mid{ii}{jj}{1} ...
                            = round(ind_decim{ii}{jj}{1}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{1}) ./ 2));
            num_decim{ii}{jj}(1) ...
                            = length(ind_decim{ii}{jj}{1}) - 1;
            
            thick_decim{ii}{jj}{1} ...
                            = NaN(1, num_decim{ii}{jj}(1));
            for ll = 1:num_decim{ii}{jj}(1)
                thick_decim{ii}{jj}{1}(ll) ...
                            = nanmean(thick{ii}{jj}{1}(ind_decim{ii}{jj}{1}(ll):ind_decim{ii}{jj}{1}(ll + 1)));
            end
            
            if isfield(pk, 'ind_layer') % only true if merge file has been reviewed through FENCEGUI
                ind_layer{ii}{jj}{1} ...
                            = pk.ind_layer;
                if ~isempty(pk.ind_layer)
                    num_fence{ii}(jj) ...
                            = 1; % number of merge files for this transect
                    ind_fence{ii}{jj} ...
                            = true;
                else
                    ind_fence{ii}{jj} ...
                            = false;
                end
            end
            
        else
            
            % loop through possible sub-transects (letters a-z)
            for kk = 1:length(letters)
                
                if ~exist([name_year{ii} '/merge/' name_trans{ii}{jj} letters(kk) '_pk_merge.mat'], 'file') % load merge file (a/b/c/etc)
                    continue
                end
                
                pk          = load([name_year{ii} '/merge/' name_trans{ii}{jj} letters(kk) '_pk_merge.mat']);
                pk          = pk.pk;
                
                disp(['Transect: ' name_trans{ii}{jj} letters(kk) '...'])
                
                [depth_smooth{ii}{jj}{kk}, dist{ii}{jj}{kk}, elev_air_gimp{ii}{jj}{kk}, elev_bed_gimp{ii}{jj}{kk}, elev_smooth_gimp{ii}{jj}{kk}, elev_surf_gimp{ii}{jj}{kk}, file_block{ii}{jj}{kk}, ind_trace_start{ii}{jj}{kk}, int_bed{ii}{jj}{kk}, int_smooth{ii}{jj}{kk}, int_surf{ii}{jj}{kk}, ...
                 num_layer{ii}{jj}(kk), num_trace_tot{ii}{jj}(kk), time{ii}{jj}{kk}, twtt_surf{ii}{jj}{kk}, x_pk{ii}{jj}{kk}, y_pk{ii}{jj}{kk}] ...
                            = deal(pk.depth_smooth, pk.dist, pk.elev_air, pk.elev_bed_gimp, pk.elev_smooth_gimp, pk.elev_surf_gimp, pk.file_block, pk.ind_trace_start, pk.int_bed, pk.int_smooth, pk.int_surf, pk.num_layer, pk.num_trace_tot, pk.time, pk.twtt_surf, pk.x, pk.y);
                thick{ii}{jj}{kk} ...
                            = elev_surf_gimp{ii}{jj}{kk} - elev_bed_gimp{ii}{jj}{kk};
                thick{ii}{jj}{kk}(isinf(thick{ii}{jj}{kk})) ...
                            = NaN;
                % indices to be kept (decimation vector)
                dist_int    = dist{ii}{jj}{kk}(1):dist_int_decim:dist{ii}{jj}{kk}(end);
                dist_int(end) ...
                            = dist{ii}{jj}{kk}(end);
                ind_decim{ii}{jj}{kk} ...
                            = unique(interp1(unique(dist{ii}{jj}{kk}), 1:length(unique(dist{ii}{jj}{kk})), dist_int, 'nearest', 'extrap'));
                ind_decim_mid{ii}{jj}{kk} ...
                            = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));
                num_decim{ii}{jj}(kk) ...
                            = length(ind_decim{ii}{jj}{kk}) - 1;
                
                % ice thickness at fixed intervals
                thick_decim{ii}{jj}{kk} ...
                            = NaN(1, num_decim{ii}{jj}(kk));
                for ll = 1:num_decim{ii}{jj}(kk)
                    thick_decim{ii}{jj}{kk}(ll) ...
                            = nanmean(thick{ii}{jj}{kk}(ind_decim{ii}{jj}{kk}(ll):ind_decim{ii}{jj}{kk}(ll + 1)));
                end
                
                if isfield(pk, 'ind_layer') % only true if merge file has been reviewed through FENCEGUI
                    ind_layer{ii}{jj}{kk} ...
                            = pk.ind_layer;
                    if ~isempty(pk.ind_layer)
                        num_fence{ii}(jj) ...
                            = num_fence{ii}(jj) + 1;
                        ind_fence{ii}{jj}(kk) ...
                            = true;
                    else
                        ind_fence{ii}{jj}(kk) ...
                            = false;
                    end
                else
                    ind_fence{ii}{jj}(kk) ...
                            = false;
                end
            end
        end
    end
end

disp('Saving merged data...')
switch radar_type
    case 'accum'
        save('mat/merge_all_accum', '-v7.3', 'depth_smooth', 'dist', 'elev_air_gimp', 'elev_bed_gimp', 'elev_smooth_gimp', 'elev_surf_gimp', 'file_block', 'ind_decim', 'ind_decim_mid', 'ind_fence', 'ind_layer', 'ind_trace_start', 'int_bed', 'int_smooth', 'int_surf', 'num_decim', 'num_fence', ...
                                             'num_layer', 'num_trace_tot', 'thick', 'thick_decim', 'time', 'twtt_surf', 'x_pk', 'y_pk')
        disp('Saved merged data in mat/merge_all_accum.mat.')
    case 'deep'
        save('mat/merge_all', '-v7.3', 'depth_smooth', 'dist', 'elev_air_gimp', 'elev_bed_gimp', 'elev_smooth_gimp', 'elev_surf_gimp', 'file_block', 'ind_decim', 'ind_decim_mid', 'ind_fence', 'ind_layer', 'ind_trace_start', 'int_bed', 'int_smooth', 'int_surf', 'num_decim', 'num_fence', ...
                                       'num_layer', 'num_trace_tot', 'thick', 'thick_decim', 'time', 'twtt_surf', 'x_pk', 'y_pk')
        disp('Saved merged data in mat/merge_all.mat.')
end