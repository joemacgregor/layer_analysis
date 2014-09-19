% XY_INT_ACCUM Determine intersections of radar transects.
% 
%   XY_INT_ACCUM determines the positions of intersections between each
%   transect and itself, all other transects in the campaign, and finally
%   all other transects in all other campaigns. It first loads the
%   positions from original data frames and concatenates them for each
%   transect.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/31/14

clear

dir_save                    = 'mat/';
do_xy                       = true;
decim                       = 25;
plotting                    = false;

load([dir_save 'gimp_proj'], 'gimp_proj')

wgs84                       = wgs84Ellipsoid;

%% concatenate x/y positions for each campaign

if do_xy
       
    % extract names and ditch the ones we don't need
    name_year               = {'2011_accum' '2012_accum' '2013_accum'};
    num_year                = length(name_year);
    
    [lat, lon, dist, dist_diff, num_frame, name_trans, time, x, y] ...
                            = deal(cell(1, num_year));
    num_trans               = zeros(1, num_year);
        
    disp('Extracting x/y positions from data frames...')
    
    for ii = 1:num_year
        
        disp(['Campaign ' name_year{ii} ' (' num2str(ii) ' / ' num2str(num_year) ')...'])
        
        name_trans{ii}      = dir([name_year{ii} '/data/']);
        ind_dir             = false(1, length(name_trans{ii}));
        for jj = 1:length(name_trans{ii})
            if name_trans{ii}(jj).isdir
                ind_dir(jj) = true;
            end
        end
        
        name_trans{ii}      = {name_trans{ii}(ind_dir).name};
        name_trans{ii}      = name_trans{ii}(~strcmp(name_trans{ii}, '.'));
        name_trans{ii}      = name_trans{ii}(~strcmp(name_trans{ii}, '..'));
        
        num_trans(ii)       = length(name_trans{ii});
        num_frame{ii}       = zeros(1, num_trans(ii));
        [lat{ii}, lon{ii}, dist{ii}, dist_diff{ii}, time{ii}, x{ii}, y{ii}] ...
                            = deal(cell(1, num_trans(ii)));
        
        for jj = 1:num_trans(ii)
            
            disp([ name_trans{ii}{jj} ' (' num2str(jj) ' / ' num2str(num_trans(ii)) ')...'])
            
            dir_curr        = dir([name_year{ii} '/data/' name_trans{ii}{jj} '/*.mat']);
            dir_curr        = {dir_curr.name};
            num_frame{ii}(jj) ...
                            = length(dir_curr);
            
            for kk = 1:num_frame{ii}(jj)
                
                tmp         = load([name_year{ii} '/data/' name_trans{ii}{jj} '/' dir_curr{kk}], 'Latitude', 'Longitude', 'GPS_time');
                [lat_curr, lon_curr, time_curr] ...
                            = deal(tmp.Latitude, tmp.Longitude, tmp.GPS_time);
                if any(isnan(lat_curr))
                    [lon_curr, time_curr] ...
                            = deal(lon_curr(~isnan(lat_curr)), time_curr(~isnan(lat_curr)));
                    lat_curr= lat_curr(~isnan(lat_curr));
                end
                
                % ignore overlapping data
                if (kk > 1)
                    if any(time_curr < time_old(end))
                        ind_tmp ...
                            = find(time_curr > time_old(end));
                        if ~isempty(ind_tmp)
                            [lat_curr, lon_curr, time_curr] ...
                                = deal(lat_curr(ind_tmp), lon_curr(ind_tmp), time_curr(ind_tmp));
                        else
                            continue % no points in this frame that don't overlap with previous one
                        end
                    end
                end
                
                [x_curr, y_curr] ...
                            = projfwd(gimp_proj, lat_curr, lon_curr);
                [lat{ii}{jj}, lon{ii}{jj}, time{ii}{jj}, x{ii}{jj}, y{ii}{jj}] ...
                            = deal([lat{ii}{jj} lat_curr], [lon{ii}{jj} lon_curr], [time{ii}{jj} time_curr], [x{ii}{jj} (1e-3 .* x_curr)], [y{ii}{jj} (1e-3 .* y_curr)]);
                
                time_old    = time_curr;
                
            end
            
            dist{ii}{jj}    = 1e-3 .* cumsum([0 distance([lat{ii}{jj}(1:(end - 1))' lon{ii}{jj}(1:(end - 1))'], [lat{ii}{jj}(2:end)' lon{ii}{jj}(2:end)'], wgs84)']);
            dist_diff{ii}{jj} ...
                            = diff(dist{ii}{jj});
            
        end
    end
    
    save([dir_save 'xy_all_accum'], 'dist', 'dist_diff', 'lat', 'lon', 'num_frame', 'num_trans', 'num_year', 'name_year', 'name_trans', 'time', 'x', 'y')
    disp(['Done extracting x/y positions and saved in ' dir_save '.'])
    
else
    load([dir_save 'xy_all_accum'])
    disp(['Loaded x/y positions from ' dir_save '.'])
end

%%
if plotting
    
%% all transects and intersections
    figure
    hold on
    colors                  = colormap(jet(num_year));
    plot(x_gl, y_gl, 'k', 'linewidth', 1)
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            plot(x{ii}{jj}(1:decim:end), y{ii}{jj}(1:decim:end), '.', 'markersize', 12, 'color', colors(ii, :))
        end
    end
    set(gca, 'fontsize', 20)
    xlabel('Polar stereographic X (km)')
    ylabel('Polar stereographic Y (km)')
    axis equal tight
    grid on
    box on
    
%%
end