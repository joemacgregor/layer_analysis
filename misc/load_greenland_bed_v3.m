% LOAD_GREENLAND_BED_V3 Load new Greenland bed-elevation data and related variables.
% 
% Load Greenland bed-elevation data from Bamber et al. [2013, The Cyrosphere, 7, 499-510].
% Updated for v3 correction mentioned on Cryolist on 07/16/13.
% 
% Joe MacGregor (UTIG), joemac@ig.utexas.edu
% Last updated: 08/07/14

clear

do_load                     = false;
do_save                     = false;
plotting                    = false;

load mat/grl_coast num_coast x_coast y_coast

if do_load
    
    % NetCDF path/file
    file_nc                 = 'mat/Greenland_bedrock_topography_V3.nc';
    
    disp(['Loading ' file_nc '...'])
    
    % get variable names
    nc_info                 = ncinfo(file_nc);
    num_var                 = length(nc_info.Variables);
    load mat/gimp_proj gimp_proj
    
    % load variables as their original names
    var_long                = cell(1, num_var);
    for ii = 1:num_var
        var_long{ii}        = nc_info.Variables(ii).Name;
        if (ii == 1)
            polar_stereo    = ncread(file_nc, var_long{ii});
            
        else
            eval([var_long{ii} ' = ncread(file_nc, var_long{ii});'])
        end
    end
    
    % rename variables with shorter/simpler names
    var_short               = {'polar stereo' 'x' 'y' 'elev_bed' 'elev_surf' 'thick' 'elev_surf_uncert' 'elev_bed_uncert' 'mask_land' 'num_pt' 'elev_geoid' 'mask_elev_bed_change' 'mask_shelf' 'elev_bed_raw' 'thick_raw' 'mask_bathy'};
    for ii = 2:num_var
        eval([var_short{ii} ' = ' var_long{ii} ';']) % variable name switch
        eval([var_short{ii} '(' var_short{ii} ' == -9999) = NaN;']) % NaN out unknown values
        if any(strcmp(var_short{ii}, {'mask_land' 'mask_elev_bed_change' 'mask_shelf' 'mask_bathy'})) % make masks logical
            tmp             = false(size(eval(var_short{ii})));
            eval(['tmp(~isnan(' var_short{ii} ' )) = logical(' var_short{ii} '(~isnan(' var_short{ii} ')));'])
            eval([var_short{ii} ' = tmp;'])
        end
    end
    
    % rotate matrices 90 deg clockwise
    for ii = 4:num_var
        eval([var_short{ii} ' = flipud(rot90(' var_short{ii} ', 1));'])
    end
    
    % grid polar stereographic coordinates
    [x_bamber, y_bamber]    = deal(x, y);
    clear x y
    [x_bamber, y_bamber]    = meshgrid(x_bamber, y_bamber);
    
    % projection details
    y_false                 = str2double(nc_info.Variables(1).Attributes(2).Value);
    x_false                 = str2double(nc_info.Variables(1).Attributes(3).Value);
    lat_0                   = str2double(nc_info.Variables(1).Attributes(5).Value);
    lat_std                 = str2double(nc_info.Variables(1).Attributes(6).Value);
    lon_std                 = str2double(nc_info.Variables(1).Attributes(7).Value);
    
    wgs84                   = almanac('earth', 'wgs84', 'meters');
    [lat, lon]              = polarstereo_inv(x_bamber, y_bamber, wgs84(1), wgs84(2), lat_std, lon_std);
    [x_tmp, y_tmp]          = projfwd(gimp_proj, double(lat), double(lon)); % convert positions to GIMP projection
    
    % m to km
    [x, y]                  = deal(-1200e3:1e3:1300e3, (-3500e3:1e3:-500e3)'); % reasonable x/y range for GIMP projection
    [x, y]                  = meshgrid(x, y);
    
    disp('Interpolating gridded variables onto GIMP projection...')
    
    for ii = 4:num_var
        disp(var_short{ii})
        if any(strcmp(var_short{ii}, {'mask_land' 'mask_elev_bed_change' 'mask_shelf' 'mask_bathy'})) % make masks logical
            F               = scatteredInterpolant(x_tmp(:), y_tmp(:), double(eval([var_short{ii} '(:)'])), 'nearest', 'none');
            eval([var_short{ii} ' = F(x, y);'])
            tmp             = false(size(eval(var_short{ii})));
            tmp(~isnan(eval(var_short{ii}))) ...
                            = logical(round(eval([var_short{ii} '(~isnan(' var_short{ii} '))'])));
            eval([var_short{ii} ' = tmp;'])
        elseif strcmp(var_short{ii}, 'num_pt')
            F               = scatteredInterpolant(x_tmp(:), y_tmp(:), double(eval([var_short{ii} '(:)'])), 'nearest', 'none');
            eval([var_short{ii} ' = F(x, y);'])
            eval([var_short{ii} ' = round(' var_short{ii} ');'])
        else
            F               = scatteredInterpolant(x_tmp(:), y_tmp(:), double(eval([var_short{ii} '(:)'])), 'natural', 'none');
            eval([var_short{ii} ' = F(x, y);'])
        end
    end
    
    [x, y]                  = deal((1e-3 .* x), (1e-3 .* y)); % m to km
    
    mask_greenland          = false(size(x));
    mask_greenland(mask_land & ~isnan(elev_surf)) ...
                            = true;
    mask_gris               = mask_greenland;
    mask_gris(isnan(thick)) = false;
    
    % save output in MATLAB format
    if do_save
        disp('Saving mat/greenland_bed_v3.mat...')
        save mat/greenland_bed_v3 nc_info x y lat lon elev_bed elev_surf thick elev_bed_uncert elev_surf_uncert mask_land num_pt elev_geoid mask_elev_bed_change mask_shelf elev_bed_raw thick_raw mask_bathy ...
                                  mask_greenland mask_gris
        disp('Saved mat/greenland_bed_v3.mat.')
    end
    
else
    disp('Loading mat/greenland_bed_v3.mat...')
    load mat/greenland_bed_v3
    disp('Loaded mat/greenland_bed_v3.mat.')
end

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
%% plot a Greenland bed v3 variable by reassigning curr_plot to one of the var_short(4:end) strings
    curr_plot               = 'elev_bed';
    decim                   = 1; % decimate the big matrices a bit to speed up plotting
    figure('position', [10 10 1000 1000], 'name', curr_plot)
    hold on
    pcolor(x(1:decim:end, 1:decim:end), y(1:decim:end, 1:decim:end), double(eval([curr_plot '(1:decim:end, 1:decim:end)'])))
    shading interp
    caxis([min(eval([curr_plot '(:)'])) max(eval([curr_plot '(:)']))])
    axis xy image
    p_gl                    = zeros(1, num_coast);
    for ii = 1:num_coast
        p_gl(ii)            = plot3(x_coast{ii}, y_coast{ii}, zeros(1, length(x_coast{ii})), 'color', [0.5 0.5 0.5], 'linewidth', 2);
    end
    set(gca, 'fontsize', 20, 'layer', 'top')
    xlabel('Polar stereographic X (km) ')
    ylabel('Polar stereographic Y (km) ')
    title([curr_plot ' '], 'interpreter', 'none')
    colorbar('fontsize', 20)
    grid on
    box on
%%
end