% LOAD_MC_BED Load Morlighem et al. (2014) Greenland mass-conserving bed-elevation model.
% 
% Load Greenland mass-conserving bed-elevation model developed and described by Morlighem et al. [2014, Nature Geosci., 7, 418-422, doi:10.1038/ngeo2167].
% 
% Joe MacGregor (UTIG), joemac@ig.utexas.edu
% Last updated: 06/18/14

clear

do_load                     = false;
do_save                     = false;
plotting                    = false;

if do_load

    % NetCDF path/file
    file_nc                 = '~/Documents/research/data/greenland/bedmap/MCdataset-2014-05-19.nc';
    
    disp(['Loading ' file_nc '...'])
    
    % get variable names
    nc_info                 = ncinfo(file_nc);
    num_var                 = length(nc_info.Variables);
    
    % load variables as their original names
    var_orig                = cell(1, num_var);
    for ii = 1:num_var
        var_orig{ii}        = nc_info.Variables(ii).Name;
        eval([var_orig{ii} ' = ncread(file_nc, var_orig{ii});'])
    end
    
    % rename variables with shorter/simpler names and clear original ones
    var_new                 = {'x' 'y' 'mask_gris' 'elev_surf' 'thick' 'elev_bed' 'elev_bed_uncert' 'source_grd'};
    for ii = 1:num_var
        eval([var_new{ii} ' = double(' var_orig{ii} ');']) % variable name change
        eval([var_new{ii} '(' var_new{ii} ' == -9999) = NaN;']) % NaN out unknown values
        if (ii > 2)
            eval(['clear ' var_orig{ii}])
        end
    end
    
    % rotate matrices 90 deg clockwise
    for ii = 3:num_var
        eval([var_new{ii} ' = rot90(' var_new{ii} ', 1);'])
    end
    
    % convert x/y to km, convert to grids
    [x, y]                  = meshgrid((1e-3 .* x)', flipud(1e-3 .* y));
    
    % save output in MATLAB format
    if do_save
        disp('Saving mat/greenland_mc_bed.mat...')
        save mat/greenland_mc_bed -v7.3 elev_bed elev_bed_uncert elev_surf mask_gris nc_infosource_grd thick x y
        disp('Saved mat/greenland_mc_bed.mat.')
    end
    
else
    disp('Loading mat/greenland_mc_bed.mat...')
    load mat/greenland_mc_bed
    disp('Loaded mat/greenland_mc_bed.mat.')
end

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
    
%%
    set(0, 'DefaultFigureWindowStyle', 'default')
    
%% plot a Greenland mass-conserving bed variable by reassigning curr_plot to one of the var_new(3:end) strings
    curr_plot               = 'mask_gris';
    decim                   = 10; % decimate matrices to speed up plotting
    figure('position', [10 10 1000 1000], 'name', curr_plot)
    hold on
    pcolor(x(1:decim:end, 1:decim:end), y(1:decim:end, 1:decim:end), double(eval([curr_plot '(1:decim:end, 1:decim:end)'])))
    shading interp
    caxis([min(eval([curr_plot '(:)'])) max(eval([curr_plot '(:)']))])
    axis xy image
    set(gca, 'fontsize', 20, 'layer', 'top')
    xlabel('Polar stereographic X (km) ')
    ylabel('Polar stereographic Y (km) ')
    title([curr_plot ' '], 'interpreter', 'none')
    colorbar('fontsize', 20)
    grid on
    box on
%%
end