% LOAD_GREENLAND_CISM Load Greenland CISM reference data.
% 
%   Downloaded from http://websrv.cs.umt.edu/isis/images/a/ab/Greenland1km.nc.
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
    file_nc                 = 'mat/Greenland1km.nc';
    
    disp(['Loading ' file_nc '...'])
    
    % get variable names
    nc_info                 = ncinfo(file_nc);
    num_var                 = length(nc_info.Variables);
    load mat/gimp_proj gimp_proj
    
    % load variables as their original names
    var_long                = cell(1, num_var);
    for ii = 1:num_var
        var_long{ii}        = nc_info.Variables(ii).Name;
        eval([var_long{ii} ' = ncread(file_nc, var_long{ii});'])
    end
    
    % rename variables with shorter/simpler names
    var_short               = {'accum' '' 'y' 'x' 'flux_geo' 'temp_surf' 'elev_bed' 'thick' ''};
    for ii = [1 3:(num_var - 1)]
        eval([var_short{ii} ' = ' var_long{ii} ';']) % variable name switch
    end
    
    % meters water equivalent to meters ice equivalent
    accum                   = accum ./ 0.917;
    
    thick(~thick)           = NaN;
    
    % rotate appropriate variables 90 deg counter-clockwise and flip
    for ii = [1 5:(num_var - 1)]
        eval([var_short{ii} ' = flipud(rot90(' var_short{ii} ', 1));'])
    end
    
    % grid CISM's polar stereographic coordinates
    [x_cism, y_cism]        = deal(x, y);
    clear x y
    [x_cism, y_cism]        = meshgrid(x_cism, y_cism);
    
    % CISM polar stereographic projection details
    y_false                 = str2double(nc_info.Variables(2).Attributes(5).Value);
    x_false                 = str2double(nc_info.Variables(2).Attributes(7).Value);
    lat_0                   = str2double(nc_info.Variables(2).Attributes(2).Value);
    lat_std                 = str2double(nc_info.Variables(2).Attributes(4).Value);
    lon_std                 = str2double(nc_info.Variables(2).Attributes(3).Value);
    
    % convert CISM x/y into GIMP projection
    wgs84                   = wgs84Ellipsoid('m');
    [lat, lon]              = polarstereo_inv(x_cism, y_cism, wgs84.SemimajorAxis, wgs84.Eccentricity, lat_std, lon_std); % invert 
    [x_tmp, y_tmp]          = projfwd(gimp_proj, double(lat), double(lon)); % convert geographic coordinates to GIMP projection
    
    % reasonable x/y range for GIMP projection
    [x, y]                  = deal(-1200e3:1e3:1300e3, (-3500e3:1e3:-500e3)');
    [x, y]                  = meshgrid(x, y);
    
    disp('Interpolating gridded variables onto GIMP projection...')
    
    for ii = [1 5:(num_var - 1)]
        disp(var_short{ii})
        F                   = scatteredInterpolant(x_tmp(:), y_tmp(:), double(eval([var_short{ii} '(:)'])), 'natural', 'none');
        eval([var_short{ii} ' = F(x, y);'])
    end
    
    [x, y]                  = deal((1e-3 .* x), (1e-3 .* y)); % m to km
    
    % save output in MATLAB format
    if do_save
        disp('Saving mat/greenland_cism.mat...')
        save mat/greenland_cism nc_info -v7.3 accum x y flux_geo temp_surf elev_bed thick
        disp('Saved mat/greenland_cism.mat.')
    end
    
else
    disp('Loading mat/greenland_cism.mat...')
    load mat/greenland_cism
    disp('Loaded mat/greenland_cism.mat.')
end

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
%% plot a Greenland CISM variable by reassigning curr_plot to one of the var_short([1 5:(num_var - 1)]) strings
    curr_plot               = 'accum';
    decim                   = 5; % decimate the big matrices a bit to speed up plotting
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