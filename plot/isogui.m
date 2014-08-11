function isogui
% ISOGUI Interactive visualization of dated reflectors and gridded isochrones.
%   
%   ISOGUI loads a three-dimensional GUI for interactive investigation of
%   dated traced/gridded reflectors/isochrones.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/07/14

%% Intialize variables

load mat/xy_all num_year num_trans
load mat/merge_all elev_smooth_gimp ind_decim ind_fence num_fence num_layer x_pk y_pk
load mat/date_all age %age_best layer_unique num_unique
load mat/core_int x_core y_core name_core num_core
load mat/grl_coast num_coast x_coast y_coast
load mat/greenland_bed_v3 elev_bed elev_surf mask_greenland mask_gris x y
load mat/age_grd2 age_iso depth_iso2 x_grd y_grd
load mat/grl_basin x_basin y_basin num_basin

% age_best(:, 1:2)            = 1e-3 .* age_best(:, 1:2);

int_year                    = 1e3 .* (5:10:125);
num_int_year                = length(int_year);

[elev_bed_grd, elev_surf_grd] ...
                            = deal(elev_bed, elev_surf);
elev_bed_grd(~mask_greenland) ...
                            = NaN;
elev_surf_grd(~mask_greenland | ~mask_gris) ...
                            = NaN;
clear elev_bed elev_surf mask_greenland mask_gris

% decim                       = 2;
% [depth_iso2, elev_bed_grd, elev_surf_grd] ...
%                             = deal(depth_iso2(1:decim:end, 1:decim:end, :), elev_bed_grd(1:decim:end, 1:decim:end), elev_surf_grd(1:decim:end, 1:decim:end), x_grd(1:decim:end, 1:decim:end), ...
%                                    y_grd(1:decim:end, 1:decim:end));
[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);
[elev_bed_grd, elev_surf_grd] ...
                            = deal(elev_bed_grd(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))), elev_surf_grd(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))));
clear x y

% x/y defaults
[x_min_ref, x_max_ref]      = deal(-632, 846);
[y_min_ref, y_max_ref]      = deal(-3344, -670);
[z_min_ref, z_max_ref]      = deal(-1000, 3500);
[cb_min_ref, cb_max_ref]    = deal(0, 130);
[az_min, az_max]            = deal(-180, 180);
[el_min, el_max]            = deal(-90, 90);

[x_min, x_max]              = deal(x_min_ref, x_max_ref);
[y_min, y_max]              = deal(y_min_ref, y_max_ref);
[z_min, z_max]              = deal(z_min_ref, z_max_ref);
[cb_min, cb_max]            = deal(cb_min_ref, cb_max_ref);
aspect_ratio                = 10;

% allocate a bunch of variables
iso_str                     = cell(num_int_year, 1);
for ii = 1:length(age_iso)
    iso_str{ii}             = num2str(1e-3 .* age_iso(ii));
end
[az2, el2, ii, p_pk, tmp1]  = deal(0);
curr_iso                    = 1;
curr_dim                    = '3D';

%% draw first GUI

set(0, 'DefaultFigureWindowStyle', 'docked')
if ispc % windows switch
    igui                    = figure('toolbar', 'figure', 'name', 'ISOGUI', 'position', [1920 940 1 1], 'menubar', 'none', 'keypressfcn', @keypress, 'color', 'w');
    ax                      = subplot('position', [0.08 0.10 0.90 0.83]);
    size_font               = 14;
    width_slide             = 0.01;
else
    igui                    = figure('toolbar', 'figure', 'name', 'ISOGUI', 'position', [1864 1100 1 1], 'menubar', 'none', 'keypressfcn', @keypress, 'color', 'w');
    ax                      = subplot('position', [0.08 0.10 0.86 0.81]);
    size_font               = 18;
    width_slide             = 0.02;
end
hold on
axis([x_min_ref x_max_ref y_min_ref y_max_ref z_min_ref z_max_ref])
box off
grid on
[az3, el3]                  = view(3);
set(gca, 'fontsize', size_font, 'layer', 'top', 'dataaspectratio', [1 1 aspect_ratio])
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Elevation (m)')
h_pan                       = pan;
set(h_pan, 'actionpostcallback', @panzoom)
h_zoom                      = zoom;
set(h_zoom, 'actionpostcallback', @panzoom)
h_rot                       = rotate3d;
set(h_rot, 'actionpostcallback', @rot)
colors                      = colormap([copper(64); 0.75 0.75 1; jet(num_int_year)]);
p_bedgrd                    = surf(x_grd, y_grd, elev_bed_grd);
tmp1                        = elev_bed_grd;
tmp1                        = tmp1 - min(tmp1(:)) + 1;
tmp1                        = 1 + round(63 .* (tmp1 ./ (max(tmp1(:)) - min(tmp1(:)))));
tmp1                        = 64 - tmp1;
set(p_bedgrd, 'cdatamapping', 'direct', 'cdata', tmp1)
shading interp
alpha(p_bedgrd, 0.75)
shading interp
alim([0 1])
p_surfgrd                   = surf(x_grd, y_grd, elev_surf_grd);
shading interp
set(p_surfgrd, 'cdata', 65 .* ones(size(elev_surf_grd)), 'cdatamapping', 'direct')
shading interp
alpha(p_surfgrd, 0.25)
shading interp
p_pk                        = cell(1, num_year);
for ii = 1:num_year %#ok<*FXUP>
    if ~any(num_fence{ii})
        continue
    end
    p_pk{ii}                = cell(1, num_trans(ii));
    for jj = find(num_fence{ii})
        p_pk{ii}{jj}        = cell(1, length(ind_fence{ii}{jj}));
        for kk = find(ind_fence{ii}{jj})
            ind_decim_mid   = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));
            p_pk{ii}{jj}{kk}= zeros(1, num_layer{ii}{jj}(kk));
            for ll = 1:num_layer{ii}{jj}(kk)
                if (any(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid))) && ~isnan(age{ii}{jj}{kk}(ll))) % layer default gray
                    p_pk{ii}{jj}{kk}(ll) ...
                            = plot3(x_pk{ii}{jj}{kk}(ind_decim_mid(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid)))), ...
                                    y_pk{ii}{jj}{kk}(ind_decim_mid(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid)))), ...
                                    elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid)))), ...
                                    '.', 'color', colors((65 + interp1(int_year, 1:num_int_year, age{ii}{jj}{kk}(ll), 'nearest', 'extrap')), :), 'markersize', 6, 'tag', num2str([ii jj kk ll]), 'visible', 'off');
                elseif any(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid)))
                    p_pk{ii}{jj}{kk}(ll) ...
                            = plot3(x_pk{ii}{jj}{kk}(ind_decim_mid(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid)))), ...
                                    y_pk{ii}{jj}{kk}(ind_decim_mid(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid)))), ...
                                    elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid(~isnan(elev_smooth_gimp{ii}{jj}{kk}(ll, ind_decim_mid)))), ...
                                    '.', 'color', [0.5 0.5 0.5], 'markersize', 4, 'tag', num2str([ii jj kk ll]), 'visible', 'off');
                end
            end
        end
    end
end
p_iso                       = zeros(1, num_int_year);
for ii = 1:length(age_iso)
    p_iso(ii)               = surf(x_grd, y_grd, (elev_surf_grd - squeeze(depth_iso2(:, :, ii))), (ones(size(x_grd)) .* (65 + interp1(int_year, 1:num_int_year, age_iso(ii), 'nearest', 'extrap'))));
    set(p_iso(ii), 'cdatamapping', 'direct')
    alpha(p_iso(ii), 0.25)
    shading interp
    set(p_iso(ii), 'visible', 'off')
end
p_coast                     = zeros(1, num_coast);
for ii = 1:num_coast
    p_coast(ii)             = plot3(x_coast{ii}, y_coast{ii}, zeros(1, length(x_coast{ii})), 'color', [0.5 0.5 0.5], 'linewidth', 2);
end
p_basin                     = zeros(1, num_basin);
for ii = 1:num_basin
    p_basin(ii)             = plot3(x_basin{ii}, y_basin{ii}, interp2(x_grd, y_grd, elev_surf_grd, x_basin{ii}, y_basin{ii}), 'color', 'k', 'linewidth', 2);
end
[p_core, p_corename]        = deal(zeros(1, num_core));
for ii = 1:num_core
    p_core(ii)              = plot3(repmat(x_core(ii), 1, 2), repmat(y_core(ii), 1, 2), [z_min_ref z_max_ref], 'm', 'linewidth', 2);
    p_corename(ii)          = text(double(x_core(ii) + 1), double(y_core(ii) + 1), z_max_ref, name_core{ii}, 'color', 'm', 'fontsize', 16, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w');
end
set(gca, 'fontsize', size_font, 'layer', 'top', 'dataaspectratio', [1 1 10])
caxis([0 (65 + num_int_year)])
colorbar('fontsize', size_font, 'ylim', [65 (65 + num_int_year)], 'ytick', 65:(65 + num_int_year), 'yticklabel', 0:10:130, 'position', [0.955 0.09 0.0125 0.8150]);

% sliders
x_min_slide                 = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.11 0.04 0.2 0.02], 'callback', @slide_x_min, 'min', x_min_ref, 'max', x_max_ref, 'value', x_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
x_max_slide                 = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.11 0.01 0.2 0.02], 'callback', @slide_x_max, 'min', x_min_ref, 'max', x_max_ref, 'value', x_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
y_min_slide                 = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.665 0.04 0.2 0.02], 'callback', @slide_y_min, 'min', y_min_ref, 'max', y_max_ref, 'value', y_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
y_max_slide                 = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.665 0.01 0.2 0.02], 'callback', @slide_y_max, 'min', y_min_ref, 'max', y_max_ref, 'value', y_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
z_min_slide                 = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.07 width_slide 0.32], 'callback', @slide_z_min, 'min', z_min_ref, 'max', z_max_ref, 'value', z_min_ref, ...
                                               'sliderstep', [0.01 0.1]);
z_max_slide                 = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.005 0.52 width_slide 0.32], 'callback', @slide_z_max, 'min', z_min_ref, 'max', z_max_ref, 'value', z_max_ref, ...
                                               'sliderstep', [0.01 0.1]);
el_slide                    = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.035 0.07 width_slide 0.32], 'callback', @slide_el, 'min', el_min, 'max', el_max, 'value', el3, ...
                                               'sliderstep', [0.01 0.1]);
az_slide                    = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.035 0.52 width_slide 0.32], 'callback', @slide_az, 'min', az_min, 'max', az_max, 'value', az3, ...
                                               'sliderstep', [0.01 0.1]);
cb_min_slide                = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.92 0.07 width_slide 0.32], 'callback', @slide_cb_min, 'min', cb_min, 'max', cb_max, 'value', cb_min, ...
                                               'sliderstep', [0.01 0.1]);
cb_max_slide                = uicontrol(igui, 'style', 'slider', 'units', 'normalized', 'position', [0.92 0.52 width_slide 0.32], 'callback', @slide_cb_max, 'min', cb_min, 'max', cb_max, 'value', cb_max, ...
                                               'sliderstep', [0.01 0.1]);

% slider values
x_min_edit                  = annotation('textbox', [0.07 0.045 0.04 0.03], 'string', num2str(x_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
x_max_edit                  = annotation('textbox', [0.07 0.005 0.04 0.03], 'string', num2str(x_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
y_min_edit                  = annotation('textbox', [0.63 0.045 0.04 0.03], 'string', num2str(y_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
y_max_edit                  = annotation('textbox', [0.63 0.005 0.04 0.03], 'string', num2str(y_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_min_edit                  = annotation('textbox', [0.005 0.39 0.04 0.03], 'string', num2str(z_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
z_max_edit                  = annotation('textbox', [0.005 0.84 0.04 0.03], 'string', num2str(z_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_min_edit                 = annotation('textbox', [0.92 0.39 0.04 0.03], 'string', num2str(cb_min_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
cb_max_edit                 = annotation('textbox', [0.92 0.84 0.04 0.03], 'string', num2str(cb_max_ref), 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');

% push buttons
uicontrol(igui, 'style', 'pushbutton', 'string', 'test', 'units', 'normalized', 'position', [0.895 0.965 0.03 0.03], 'callback', @misctest, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset x/y/z', 'units', 'normalized', 'position', [0.94 0.925 0.055 0.03], 'callback', @reset_xyz, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.315 0.04 0.03 0.03], 'callback', @reset_x_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.315 0.005 0.03 0.03], 'callback', @reset_x_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.87 0.04 0.03 0.03], 'callback', @reset_y_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.87 0.005 0.03 0.03], 'callback', @reset_y_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.035 0.03 0.03], 'callback', @reset_z_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.005 0.48 0.03 0.03], 'callback', @reset_z_max, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.91 0.035 0.03 0.03], 'callback', @reset_cb_min, 'fontsize', size_font, 'foregroundcolor', 'r')
uicontrol(igui, 'style', 'pushbutton', 'string', 'reset', 'units', 'normalized', 'position', [0.91 0.48 0.03 0.03], 'callback', @reset_cb_max, 'fontsize', size_font, 'foregroundcolor', 'r')

% fixed text annotations
a(1)                        = annotation('textbox', [0.005 0.965 0.03 0.03], 'string', 'Cores', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(2)                        = annotation('textbox', [0.06 0.965 0.04 0.03], 'string', 'Coastline', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(3)                        = annotation('textbox', [0.22 0.965 0.06 0.03], 'string', 'Bed grid', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(4)                        = annotation('textbox', [0.28 0.965 0.08 0.03], 'string', 'Surface grid', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(5)                        = annotation('textbox', [0.88 0.925 0.05 0.03], 'string', 'gridlines', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(6)                        = annotation('textbox', [0.04 0.04 0.03 0.03], 'string', 'x_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(7)                        = annotation('textbox', [0.04 0.005 0.04 0.03], 'string', 'x_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(8)                        = annotation('textbox', [0.595 0.04 0.03 0.03], 'string', 'y_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(9)                        = annotation('textbox', [0.595 0.005 0.04 0.03], 'string', 'y_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(10)                       = annotation('textbox', [0.005 0.42 0.03 0.03], 'string', 'z_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(11)                       = annotation('textbox', [0.005 0.87 0.03 0.03], 'string', 'z_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(12)                       = annotation('textbox', [0.005 0.89 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(13)                       = annotation('textbox', [0.005 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(14)                       = annotation('textbox', [0.56 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(15)                       = annotation('textbox', [0.915 0.005 0.03 0.03], 'string', 'fix', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(16)                       = annotation('textbox', [0.825 0.965 0.04 0.03], 'string', 'aspect', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(17)                       = annotation('textbox', [0.35 0.965 0.10 0.03], 'string', 'Dated layers', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(18)                       = annotation('textbox', [0.43 0.965 0.12 0.03], 'string', 'Gridded isochrones', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
a(19)                       = annotation('textbox', [0.035 0.87 0.08 0.03], 'string', 'az', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(20)                       = annotation('textbox', [0.035 0.42 0.08 0.03], 'string', 'el', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(21)                       = annotation('textbox', [0.905 0.42 0.05 0.03], 'string', 'age_{min}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(22)                       = annotation('textbox', [0.905 0.87 0.05 0.03], 'string', 'age_{max}', 'fontsize', size_font, 'color', 'k', 'edgecolor', 'none');
a(23)                       = annotation('textbox', [0.12 0.965 0.11 0.03], 'string', 'Drainage basins', 'fontsize', size_font, 'color', 'b', 'edgecolor', 'none');
if ~ispc
    set(a, 'fontweight', 'bold')
end

dim_group                   = uibuttongroup('position', [0.93 0.965 0.06 0.03], 'selectionchangefcn', @choose_dim);
uicontrol(igui, 'style', 'text', 'parent', dim_group, 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', size_font)
dim_check(1)                = uicontrol(igui, 'style', 'radio', 'string', '2D', 'units', 'normalized', 'position', [0.01 0.1 0.45 0.8], 'parent', dim_group, 'fontsize', size_font, 'handlevisibility', 'off');
dim_check(2)                = uicontrol(igui, 'style', 'radio', 'string', '3D', 'units', 'normalized', 'position', [0.51 0.1 0.45 0.8], 'parent', dim_group, 'fontsize', size_font, 'handlevisibility', 'off');
set(dim_group, 'selectedobject', dim_check(2))

% value boxes
aspect_edit                 = uicontrol(igui, 'style', 'edit', 'string', num2str(aspect_ratio), 'units', 'normalized', 'position', [0.86 0.965 0.03 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                              'backgroundcolor', 'w', 'callback', @adj_aspect);
az_edit                     = uicontrol(igui, 'style', 'edit', 'string', num2str(az3), 'units', 'normalized', 'position', [0.035 0.84 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                              'backgroundcolor', 'w', 'callback', @adj_az);
el_edit                     = uicontrol(igui, 'style', 'edit', 'string', num2str(el3), 'units', 'normalized', 'position', [0.035 0.39 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                              'backgroundcolor', 'w', 'callback', @adj_el);
cb_min_edit                 = uicontrol(igui, 'style', 'edit', 'string', num2str(cb_min), 'units', 'normalized', 'position', [0.91 0.39 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                              'backgroundcolor', 'w', 'callback', @adj_cb_min);
cb_max_edit                 = uicontrol(igui, 'style', 'edit', 'string', num2str(cb_max), 'units', 'normalized', 'position', [0.91 0.84 0.04 0.03], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                              'backgroundcolor', 'w', 'callback', @adj_cb_max);

% menus
% layer_list                  = uicontrol(igui, 'style', 'popupmenu', 'string', 'N/A', 'value', 1, 'units', 'normalized', 'position', [0.35 0.955 0.10 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
%                                               'callback', @choose_layer);
iso_list                    = uicontrol(igui, 'style', 'popupmenu', 'string', iso_str, 'value', 1, 'units', 'normalized', 'position', [0.525 0.955 0.05 0.04], 'fontsize', size_font, 'foregroundcolor', 'k', ...
                                              'callback', @choose_iso);

% check boxes
xfix_check                  = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.03 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
yfix_check                  = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.57 0.005 0.01 0.03], 'fontsize', size_font, 'value', 0);
zfix_check                  = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.02 0.89 0.01 0.03], 'fontsize', size_font, 'value', 0);
grid_check                  = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.925 0.925 0.01 0.03], 'callback', @toggle_grid, 'fontsize', size_font, 'value', 1);
core_check                  = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.04 0.965 0.01 0.03], 'callback', @show_core, 'fontsize', size_font, 'value', 1);
coast_check                 = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.105 0.965 0.01 0.03], 'callback', @show_coast, 'fontsize', size_font, 'value', 1);
bedgrd_check                = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.265 0.965 0.01 0.03], 'callback', @show_bedgrd, 'fontsize', size_font, 'value', 1);
surfgrd_check               = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.335 0.965 0.01 0.03], 'callback', @show_surfgrd, 'fontsize', size_font, 'value', 1);
pk_check                    = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.415 0.965 0.01 0.03], 'callback', @show_pk, 'fontsize', size_font, 'value', 0);
iso_check                   = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.585 0.965 0.01 0.03], 'callback', @show_iso, 'fontsize', size_font, 'value', 0);
cbfix_check                 = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.93 0.01 0.01 0.03], 'fontsize', size_font, 'value', 0);
basin_check                 = uicontrol(igui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.20 0.965 0.01 0.03], 'callback', @show_basin, 'fontsize', size_font, 'value', 1);

%% Choose the current layer

    function choose_iso(source, eventdata) %#ok<*INUSD>
        curr_iso            = get(iso_list, 'value');
        for ii = 1:num_int_year
            alpha(p_iso(ii), 0.25)
        end
        if (logical(p_iso(curr_iso)) && ishandle(p_iso(curr_iso)))
            alpha(p_iso(curr_iso), 0.75)
        end
        if ~get(iso_check, 'value')
            set(iso_check, 'value', 1)
        end
        show_iso
    end

%% Update minimum x/y/z

    function slide_x_min(source, eventdata)
        if (get(x_min_slide, 'value') < x_max)
            if get(xfix_check, 'value')
                tmp1        = x_max - x_min;
            end
            x_min           = get(x_min_slide, 'value');
            if get(xfix_check, 'value')
                x_max       = x_min + tmp1;
                if (x_max > x_max_ref)
                    x_max   = x_max_ref;
                    x_min   = x_max - tmp1;
                    if (x_min < x_min_ref)
                        x_min = x_min_ref;
                    end
                    if (x_min < get(x_min_slide, 'min'))
                        set(x_min_slide, 'value', get(x_min_slide, 'min'))
                    else
                        set(x_min_slide, 'value', x_min)
                    end
                end
                if (x_max > get(x_max_slide, 'max'))
                    set(x_max_slide, 'value', get(x_max_slide, 'max'))
                else
                    set(x_max_slide, 'value', x_max)
                end
                set(x_max_edit, 'string', sprintf('%4.1f', x_max))
            end
            set(x_min_edit, 'string', sprintf('%4.1f', x_min))
            update_x_range
        else
            if (x_min < get(x_min_slide, 'min'))
                set(x_min_slide, 'value', get(x_min_slide, 'min'))
            else
                set(x_min_slide, 'value', x_min)
            end
        end
        set(x_min_slide, 'enable', 'off')
        drawnow
        set(x_min_slide, 'enable', 'on')
    end

    function slide_y_min(source, eventdata)
        if (get(y_min_slide, 'value') < y_max)
            if get(yfix_check, 'value')
                tmp1        = y_max - y_min;
            end
            y_min           = get(y_min_slide, 'value');
            if get(yfix_check, 'value')
                y_max       = y_min + tmp1;
                if (y_max > y_max_ref)
                    y_max   = y_max_ref;
                    y_min   = y_max - tmp1;
                    if (y_min < y_min_ref)
                        y_min = y_min_ref;
                    end
                    if (y_min < get(y_min_slide, 'min'))
                        set(y_min_slide, 'value', get(y_min_slide, 'min'))
                    else
                        set(y_min_slide, 'value', y_min)
                    end
                end
                if (y_max > get(y_max_slide, 'max'))
                    set(y_max_slide, 'value', get(y_max_slide, 'max'))
                else
                    set(y_max_slide, 'value', y_max)
                end
                set(y_max_edit, 'string', sprintf('%4.1f', y_max))
            end
            set(y_min_edit, 'string', sprintf('%4.1f', y_min))
            update_y_range
        else
            if (y_min < get(y_min_slide, 'min'))
                set(y_min_slide, 'value', get(y_min_slide, 'min'))
            else
                set(y_min_slide, 'value', y_min)
            end
        end
        set(y_min_slide, 'enable', 'off')
        drawnow
        set(y_min_slide, 'enable', 'on')
    end

    function slide_z_min(source, eventdata)
        if (get(z_min_slide, 'value') < z_max)
            if get(zfix_check, 'value')
                tmp1        = z_max - z_min;
            end
            z_min           = get(z_min_slide, 'value');
            if get(zfix_check, 'value')
                z_max       = z_min + tmp1;
                if (z_max > z_max_ref)
                    z_max   = z_max_ref;
                    z_min   = z_max - tmp1;
                    if (z_min < z_min_ref)
                        z_min = z_min_ref;
                    end
                    if (z_min < get(z_min_slide, 'min'))
                        set(z_min_slide, 'value', get(z_min_slide, 'min'))
                    else
                        set(z_min_slide, 'value', z_min)
                    end
                end
                set(z_max_edit, 'string', sprintf('%4.0f', z_max))
                if (z_max > get(z_max_slide, 'max'))
                    set(z_max_slide, 'value', get(z_max_slide, 'max'))
                else
                    set(z_max_slide, 'value', z_max)
                end
            end
            set(z_min_edit, 'string', sprintf('%4.0f', z_min))
            update_z_range
        else
            if (z_min < get(z_min_slide, 'min'))
                set(z_min_slide, 'value', get(z_min_slide, 'min'))
            else
                set(z_min_slide, 'value', z_min)
            end
        end
        set(z_min_slide, 'enable', 'off')
        drawnow
        set(z_min_slide, 'enable', 'on')
    end

    function slide_cb_min(source, eventdata)
        if (get(cb_min_slide, 'value') < cb_max)
            if get(cbfix_check, 'value')
                tmp1        = cb_max - cb_min;
            end
            cb_min          = get(cb_min_slide, 'value');
            if get(cbfix_check, 'value')
                cb_max      = cb_min + tmp1;
                if (cb_max > cb_max_ref)
                    cb_max  = cb_max_ref;
                    cb_min  = cb_max - tmp1;
                    if (cb_min < cb_min_ref)
                        cb_min = cb_min_ref;
                    end
                    if (cb_min < get(cb_min_slide, 'min'))
                        set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                    else
                        set(cb_min_slide, 'value', cb_min)
                    end
                end
                set(cb_max_edit, 'string', sprintf('%3.0f', cb_max))
                if (cb_max > get(cb_max_slide, 'max'))
                    set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                else
                    set(cb_max_slide, 'value', cb_max)
                end
            end
            set(cb_min_edit, 'string', sprintf('%3.0f', cb_min))
            update_cb_range
        else
            if (cb_min < get(cb_min_slide, 'min'))
                set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
            else
                set(cb_min_slide, 'value', cb_min)
            end
        end
        set(cb_min_slide, 'enable', 'off')
        drawnow
        set(cb_min_slide, 'enable', 'on')
    end

%% Update maximum x/y/z

    function slide_x_max(source, eventdata)
        if (get(x_max_slide, 'value') > x_min)
            if get(xfix_check, 'value')
                tmp1        = x_max - x_min;
            end
            x_max           = get(x_max_slide, 'value');
            if get(xfix_check, 'value')
                x_min       = x_max - tmp1;
                if (x_min < x_min_ref)
                    x_min   = x_min_ref;
                    x_max   = x_min + tmp1;
                    if (x_max > x_max_ref)
                        x_max = x_max_ref;
                    end
                    if (x_max > get(x_max_slide, 'max'))
                        set(x_max_slide, 'value', get(x_max_slide, 'max'))
                    else
                        set(x_max_slide, 'value', x_max)
                    end
                end
                if (x_min < get(x_min_slide, 'min'))
                    set(x_min_slide, 'value', get(x_min_slide, 'min'))
                else
                    set(x_min_slide, 'value', x_min)
                end
                set(x_min_edit, 'string', sprintf('%4.1f', x_min))
            end
            set(x_max_edit, 'string', sprintf('%4.1f', x_max))
            update_x_range
        else
            if (x_max > get(x_max_slide, 'max'))
                set(x_max_slide, 'value', get(x_max_slide, 'max'))
            else
                set(x_max_slide, 'value', x_max)
            end
        end
        set(x_max_slide, 'enable', 'off')
        drawnow
        set(x_max_slide, 'enable', 'on')
    end

    function slide_y_max(source, eventdata)
        if (get(y_max_slide, 'value') > y_min)
            if get(yfix_check, 'value')
                tmp1        = y_max - y_min;
            end
            y_max           = get(y_max_slide, 'value');
            if get(yfix_check, 'value')
                y_min       = y_max - tmp1;
                if (y_min < y_min_ref)
                    y_min   = y_min_ref;
                    y_max   = y_min + tmp1;
                    if (y_max > y_max_ref)
                        y_max = y_max_ref;
                    end
                    if (y_max > get(y_max_slide, 'max'))
                        set(y_max_slide, 'value', get(y_max_slide, 'max'))
                    else
                        set(y_max_slide, 'value', y_max)
                    end
                end
                if (y_min < get(y_min_slide, 'min'))
                    set(y_min_slide, 'value', get(y_min_slide, 'min'))
                else
                    set(y_min_slide, 'value', y_min)
                end
                set(y_min_edit, 'string', sprintf('%4.1f', y_min))
            end
            set(y_max_edit, 'string', sprintf('%4.1f', y_max))
            update_y_range
        else
            if (y_max > get(y_max_slide, 'max'))
                set(y_max_slide, 'value', get(y_max_slide, 'max'))
            else
                set(y_max_slide, 'value', y_max)
            end
        end
        set(y_max_slide, 'enable', 'off')
        drawnow
        set(y_max_slide, 'enable', 'on')
    end

    function slide_z_max(source, eventdata)
        if (get(z_max_slide, 'value') > z_min)
            if get(zfix_check, 'value')
                tmp1        = z_max - z_min;
            end
            z_max = get(z_max_slide, 'value');
            if get(zfix_check, 'value')
                z_min       = z_max - tmp1;
                if (z_min < z_min_ref)
                    z_min   = z_min_ref;
                    z_max   = z_min + tmp1;
                    if (z_max > get(z_max_slide, 'max'))
                        set(z_max_slide, 'value', get(z_max_slide, 'max'))
                    else
                        set(z_max_slide, 'value', z_max)
                    end
                end
                set(z_min_edit, 'string', sprintf('%4.0f', z_min))
                if (z_min < get(z_min_slide, 'min'))
                    set(z_min_slide, 'value', get(z_min_slide, 'min'))
                else
                    set(z_min_slide, 'value', z_min)
                end
            end
            set(z_max_edit, 'string', sprintf('%4.0f', z_max))
            update_z_range
        else
            if (z_max > get(z_max_slide, 'max'))
                set(z_max_slide, 'value', get(z_max_slide, 'max'))
            else
                set(z_max_slide, 'value', z_max)
            end
        end
        set(z_max_slide, 'enable', 'off')
        drawnow
        set(z_max_slide, 'enable', 'on')
    end

    function slide_cb_max(source, eventdata)
        if (get(cb_max_slide, 'value') > cb_min)
            if get(cbfix_check, 'value')
                tmp1        = cb_max - cb_min;
            end
            cb_max          = get(cb_max_slide, 'value');
            if get(cbfix_check, 'value')
                cb_min      = cb_max - tmp1;
                if (cb_min < cb_min_ref)
                    cb_min  = cb_min_ref;
                    cb_max  = cb_min + tmp1;
                    if (cb_max > get(cb_max_slide, 'max'))
                        set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
                    else
                        set(cb_max_slide, 'value', cb_max)
                    end
                end
                set(cb_min_edit, 'string', sprintf('%3.0f', cb_min))
                if (cb_min < get(cb_min_slide, 'min'))
                    set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
                else
                    set(cb_min_slide, 'value', cb_min)
                end
            end
            set(cb_max_edit, 'string', sprintf('%3.0f', cb_max))
            update_cb_range
        else
            if (cb_max > get(cb_max_slide, 'max'))
                set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
            else
                set(cb_max_slide, 'value', cb_max)
            end
        end
        set(cb_max_slide, 'enable', 'off')
        drawnow
        set(cb_max_slide, 'enable', 'on')
    end

%% Update azimuth/elevation

    function slide_az(source, eventdata)
        az3         = get(az_slide, 'value');
        set(az_edit, 'string', sprintf('%3.0f', az3))
        set(az_slide, 'enable', 'off')
        view(az3, el3)
        drawnow
        set(az_slide, 'enable', 'on')
    end

    function slide_el(source, eventdata)
        el3         = get(el_slide, 'value');
        set(el_edit, 'string', sprintf('%3.0f', el3))
        set(el_slide, 'enable', 'off')
        view(az3, el3)
        drawnow
        set(el_slide, 'enable', 'on')
    end

%% Reset minimum x/y/z/cb

    function reset_x_min(source, eventdata)
        if (x_min_ref < get(x_min_slide, 'min'))
            set(x_min_slide, 'value', get(x_min_slide, 'min'))
        else
            set(x_min_slide, 'value', x_min_ref)
        end
        set(x_min_edit, 'string', sprintf('%4.1f', x_min_ref))
        x_min               = x_min_ref;
        update_x_range
    end

    function reset_y_min(source, eventdata)
        if (y_min_ref < get(y_min_slide, 'min'))
            set(y_min_slide, 'value', get(y_min_slide, 'min'))
        else
            set(y_min_slide, 'value', y_min_ref)
        end
        set(y_min_edit, 'string', sprintf('%4.1f', y_min_ref))
        y_min               = y_min_ref;
        update_y_range
    end

    function reset_z_min(source, eventdata)
        if (z_min_ref < get(z_min_slide, 'min'))
            set(z_min_slide, 'value', get(z_min_slide, 'min'))
        else
            set(z_min_slide, 'value', z_min_ref)
        end
        set(z_min_edit, 'string', sprintf('%4.0f', z_min_ref))
        z_min               = z_min_ref;
        update_z_range
    end

    function reset_cb_min(source, eventdata)
        if (cb_min_ref < get(cb_min_slide, 'min'))
            set(cb_min_slide, 'value', get(cb_min_slide, 'min'))
        else
            set(cb_min_slide, 'value', cb_min_ref)
        end
        set(cb_min_edit, 'string', sprintf('%3.0f', cb_min_ref))
        cb_min              = cb_min_ref;
        update_cb_range
    end

%% Reset maximum x/y/z/cb

    function reset_x_max(source, eventdata)
        if (x_max_ref > get(x_max_slide, 'max'))
            set(x_max_slide, 'value', get(x_max_slide, 'max'))
        else
            set(x_max_slide, 'value', x_max_ref)
        end
        set(x_max_edit, 'string', sprintf('%4.1f', x_max_ref))
        x_max               = x_max_ref;
        update_x_range
    end

    function reset_y_max(source, eventdata)
        if (y_max_ref > get(y_max_slide, 'max'))
            set(y_max_slide, 'value', get(y_max_slide, 'max'))
        else
            set(y_max_slide, 'value', y_max_ref)
        end
        set(y_max_edit, 'string', sprintf('%4.1f', y_max_ref))
        y_max               = y_max_ref;
        update_y_range
    end

    function reset_z_max(source, eventdata)
        if (z_max_ref > get(z_max_slide, 'max'))
            set(z_max_slide, 'value', get(z_max_slide, 'max'))
        else
            set(z_max_slide, 'value', z_max_ref)
        end
        set(z_max_edit, 'string', sprintf('%4.0f', z_max_ref))
        z_max            = z_max_ref;
        update_z_range
    end

    function reset_cb_max(source, eventdata)
        if (cb_max_ref > get(cb_max_slide, 'max'))
            set(cb_max_slide, 'value', get(cb_max_slide, 'max'))
        else
            set(cb_max_slide, 'value', cb_max_ref)
        end
        set(cb_max_edit, 'string', sprintf('%3.0f', cb_max_ref))
        cb_max              = cb_max_ref;
        update_cb_range
    end

%% Reset all x/y/z

    function reset_xyz(source, eventdata)
        reset_x_min
        reset_x_max
        reset_y_min
        reset_y_max
        reset_z_min
        reset_z_max
        [az3, el3]= view(3);
        set(az_edit, 'string', num2str(az3))
        set(el_edit, 'string', num2str(el3))        
    end

%% Update x/y/z/cb range

    function update_x_range(source, eventdata)
        xlim([x_min x_max])
    end

    function update_y_range(source, eventdata)
        ylim([y_min y_max])
    end

    function update_z_range(source, eventdata)
        zlim([z_min z_max])
    end

    function update_cb_range(source, eventdata)
        for ii = 1:num_year %#ok<*FXUP>
            if ~any(num_fence{ii})
                continue
            end
            for jj = find(num_fence{ii})
                for kk = find(ind_fence{ii}{jj})
                    if (any(p_pk{ii}{jj}{kk}) && any(ishandle(p_pk{ii}{jj}{kk})))
                        if ~isempty(find((logical(p_pk{ii}{jj}{kk}) & ishandle(p_pk{ii}{jj}{kk})), 1))
                            set(p_pk{ii}{jj}{kk}(logical(p_pk{ii}{jj}{kk}) & ishandle(p_pk{ii}{jj}{kk})), 'visible', 'off')
                        end
                        if ~isempty(find((logical(p_pk{ii}{jj}{kk}) & ishandle(p_pk{ii}{jj}{kk}) & (age{ii}{jj}{kk}' >= (1e3 * cb_min)) & (age{ii}{jj}{kk}' <= (1e3 * cb_max))), 1))
                            set(p_pk{ii}{jj}{kk}(logical(p_pk{ii}{jj}{kk}) & ishandle(p_pk{ii}{jj}{kk}) & (age{ii}{jj}{kk}' >= (1e3 * cb_min)) & (age{ii}{jj}{kk}' <= (1e3 * cb_max))), 'visible', 'on')
                        end
                    end
                end
            end
        end
    end

%% Show dated layers

    function show_pk(source, eventdata)
        if get(pk_check, 'value')
            for ii = 1:num_year
                if ~any(num_fence{ii})
                    continue
                end
                for jj = find(num_fence{ii})
                    for kk = find(ind_fence{ii}{jj})
                        if (any(p_pk{ii}{jj}{kk}) && any(ishandle(p_pk{ii}{jj}{kk})))
                            set(p_pk{ii}{jj}{kk}(logical(p_pk{ii}{jj}{kk}) & ishandle(p_pk{ii}{jj}{kk})), 'visible', 'on')
                        end
                    end
                end
            end
        else
            for ii = 1:num_year
                if ~any(num_fence{ii})
                    continue
                end
                for jj = find(num_fence{ii})
                    for kk = find(ind_fence{ii}{jj})
                        if (any(p_pk{ii}{jj}{kk}) && any(ishandle(p_pk{ii}{jj}{kk})))
                            set(p_pk{ii}{jj}{kk}(logical(p_pk{ii}{jj}{kk}) & ishandle(p_pk{ii}{jj}{kk})), 'visible', 'off')
                        end
                    end
                end
            end
        end
    end

%% Show gridded isochones

    function show_iso(source, eventdata)
        if get(iso_check, 'value')
            if (any(p_iso) && any(ishandle(p_iso)))
                set(p_iso(logical(p_iso) & ishandle(p_iso)), 'visible', 'off')
            end
            if (logical(p_iso(curr_iso)) && ishandle(p_iso(curr_iso)))
                set(p_iso(curr_iso), 'visible', 'on')
                uistack(p_iso(curr_iso), 'top')
            end
        else
            if (any(p_iso) && any(ishandle(p_iso)))
                set(p_iso(logical(p_iso) & ishandle(p_iso)), 'visible', 'off')
            end
        end
    end

%% Show core locations

    function show_core(source, eventdata)
        if get(core_check, 'value')
            if (any(p_core) && any(ishandle(p_core)))
                set(p_core(logical(p_core) & ishandle(p_core)), 'visible', 'on')
                uistack(p_core(logical(p_core) & ishandle(p_core)), 'top')
            end
            if (any(p_corename) && any(ishandle(p_corename)))
                set(p_corename(logical(p_corename) & ishandle(p_corename)), 'visible', 'on')
                uistack(p_corename(logical(p_corename) & ishandle(p_corename)), 'top')
            end
        else
            if (any(p_core) && any(ishandle(p_core)))
                set(p_core(logical(p_core) & ishandle(p_core)), 'visible', 'off')
            end
            if (any(p_corename) && any(ishandle(p_corename)))
                set(p_corename(logical(p_corename) & ishandle(p_corename)), 'visible', 'off')
            end
        end
    end

%% Show coastline

    function show_coast(source, eventdata)
        if get(coast_check, 'value')
            if (any(p_coast) && any(ishandle(p_coast)))
                set(p_coast(logical(p_coast) & ishandle(p_coast)), 'visible', 'on')
                uistack(p_coast(logical(p_coast) & ishandle(p_coast)), 'top')
            end
        else
            if (any(p_coast) && any(ishandle(p_coast)))
                set(p_coast(logical(p_coast) & ishandle(p_coast)), 'visible', 'off')
            end
        end
    end

%% Show drainage basins

    function show_basin(source, eventdata)
        if get(basin_check, 'value')
            if (any(p_basin) && any(ishandle(p_basin)))
                set(p_basin(logical(p_basin) & ishandle(p_basin)), 'visible', 'on')
                uistack(p_basin(logical(p_basin) & ishandle(p_basin)), 'top')
            end
        else
            if (any(p_basin) && any(ishandle(p_basin)))
                set(p_basin(logical(p_basin) & ishandle(p_basin)), 'visible', 'off')
            end
        end
    end

%% Show bed grid

    function show_bedgrd(source, eventdata)
        if get(bedgrd_check, 'value')
            if (logical(p_bedgrd) && ishandle(p_bedgrd))
                set(p_bedgrd, 'visible', 'on')
                uistack(p_bedgrd, 'top')
            end
        else
            if (logical(p_bedgrd) && ishandle(p_bedgrd))
                set(p_bedgrd, 'visible', 'off')
            end
        end
    end

%% Show surface grid

    function show_surfgrd(source, eventdata)
        if get(surfgrd_check, 'value')
            if (logical(p_surfgrd) && ishandle(p_surfgrd))
                set(p_surfgrd, 'visible', 'on')
                uistack(p_surfgrd, 'top')
            end
        else
            if (logical(p_surfgrd) && ishandle(p_surfgrd))
                set(p_surfgrd, 'visible', 'off')
            end
        end
    end

%% Adjust 3D display aspect ratio / azimuth / elevation

    function adj_aspect(source, eventdata)
        aspect_ratio        = abs(round(str2double(get(aspect_edit, 'string'))));
        set(ax, 'dataaspectratio', [1 1 aspect_ratio])
        set(aspect_edit, 'string', num2str(aspect_ratio))
    end

    function adj_az(source, eventdata)
        if ((round(str2double(get(az_edit, 'string'))) >= az_min) && (round(str2double(get(el_edit, 'string'))) <= az_max))
            az3            = round(str2double(get(az_edit, 'string')));
            set(az_slide, 'value', az3)
            view(az3, el3)
            set(az_edit, 'string', num2str(az3))
        else
            set(az_edit, 'string', num2str(az3))
        end
    end

    function adj_el(source, eventdata)
        if ((round(str2double(get(el_edit, 'string'))) >= el_min) && (round(str2double(get(el_edit, 'string'))) <= el_max))
            el3        = round(str2double(get(el_edit, 'string')));
            set(el_slide, 'value', el3)
            view(az3, el3)
            set(el_edit, 'string', num2str(el3))
        else
            set(el_edit, 'string', num2str(el3))
        end
    end

    function adj_cb_min(source, eventdata)
        if ((round(str2double(get(cb_min_edit, 'string'))) >= cb_min_ref) && (round(str2double(get(cb_min_edit, 'string'))) <= cb_max_ref) && (round(str2double(get(cb_min_edit, 'string'))) < cb_max))
            cb_min          = round(str2double(get(cb_min_edit, 'string')));
            update_cb_range
            set(cb_min_slide, 'value', cb_min)
            set(cb_min_edit, 'string', num2str(cb_min))
        else
            set(cb_min_edit, 'string', num2str(cb_min))
        end
    end

    function adj_cb_max(source, eventdata)
        if ((round(str2double(get(cb_max_edit, 'string'))) >= cb_min_ref) && (round(str2double(get(cb_max_edit, 'string'))) <= cb_max_ref) && (round(str2double(get(cb_max_edit, 'string'))) > cb_min))
            cb_max          = round(str2double(get(cb_max_edit, 'string')));
            update_cb_range
            set(cb_max_slide, 'value', cb_max)
            set(cb_max_edit, 'string', num2str(cb_max))
        else
            set(cb_max_edit, 'string', num2str(cb_max))
        end
    end

%% Pan/zoom/rotate

    function panzoom(source, eventdata)
        tmp1                = get(gca, 'xlim');
        if (tmp1(1) < z_min_ref)
            reset_x_min
        else
            if (tmp1(1) < get(x_min_slide, 'min'))
                set(x_min_slide, 'value', get(x_min_slide, 'min'))
            else
                set(x_min_slide, 'value', tmp1(1))
            end
            set(x_min_edit, 'string', sprintf('%3.1f', tmp1(1)))
            x_min           = tmp1(1);
        end
        if (tmp1(2) > x_max_ref)
            reset_x_max
        else
            if (tmp1(2) > get(x_max_slide, 'max'))
                set(x_max_slide, 'value', get(x_max_slide, 'max'))
            else
                set(x_max_slide, 'value', tmp1(2))
            end
            set(x_max_edit, 'string', sprintf('%3.1f', tmp1(2)))
            x_max           = tmp1(2);
        end
        tmp1                = get(gca, 'ylim');
        if (tmp1(1) < y_min_ref)
            reset_y_min
        else
            if (tmp1(1) < get(y_min_slide, 'min'))
                set(y_min_slide, 'value', get(y_min_slide, 'min'))
            else
                set(y_min_slide, 'value', tmp1(1))
            end
            set(y_min_edit, 'string', sprintf('%4.0f', tmp1(1)))
            y_min           = tmp1(1);
        end
        if (tmp1(2) > y_max_ref)
            reset_y_max
        else
            if (tmp1(2) > get(y_max_slide, 'max'))
                set(y_max_slide, 'value', get(y_max_slide, 'max'))
            else
                set(y_max_slide, 'value', tmp1(2))
            end
            set(y_max_edit, 'string', sprintf('%4.0f', tmp1(2)))
            y_max           = tmp1(2);
        end
        tmp1                = get(gca, 'zlim');
        if (tmp1(1) < z_min_ref)
            reset_z_min
        else
            if (tmp1(1) < get(z_min_slide, 'min'))
                set(z_min_slide, 'value', get(z_min_slide, 'min'))
            else
                set(z_min_slide, 'value', tmp1(1))
            end
            set(z_min_edit, 'string', sprintf('%4.0f', tmp1(1)))
            z_min           = tmp1(1);
        end
        if (tmp1(2) > z_max_ref)
            reset_z_max
        else
            if (tmp1(2) > get(z_max_slide, 'max'))
                set(z_max_slide, 'value', get(z_max_slide, 'max'))
            else
                set(z_max_slide, 'value', tmp1(2))
            end
            set(z_max_edit, 'string', sprintf('%4.0f', tmp1(2)))
            z_max           = tmp1(2);
        end
    end

    function rot(source, eventdata)
        [az3, el3]= view;
        set(az_edit, 'string', num2str(az3))
        set(el_edit, 'string', num2str(el3))
    end

%% Toggle gridlines

    function toggle_grid(source, evenntdata)
        if get(grid_check, 'value')
            grid on
        else
            grid off
        end
    end

%% Switch viewing dimension

    function choose_dim(~, eventdata)
        curr_dim            = get(eventdata.NewValue, 'string');
        switch curr_dim
            case '2D'
                [az3, el3] ...
                            = view;
                if any([az2 el2])
                    view(az2, el2)
                else
                    view(2)
                end
            case '3D'
                [az2, el2] ...
                            = view;
                view(az3, el3)
        end
    end

%% Keyboard shortcuts for various functions

    function keypress(~, eventdata)
        switch eventdata.Key
            case '1'
                if get(core_check, 'value')
                    set(core_check, 'value', 0)
                else
                    set(core_check, 'value', 1)
                end
                show_core
            case '2'
                if get(coast_check, 'value')
                    set(coast_check, 'value', 0)
                else
                    set(coast_check, 'value', 1)
                end
                show_coast
            case '3'
                if get(basin_check, 'value')
                    set(basin_check, 'value', 0)
                else
                    set(basin_check, 'value', 1)
                end
                show_basin
            case '4'
                if get(bedgrd_check, 'value')
                    set(bedgrd_check, 'value', 0)
                else
                    set(bedgrd_check, 'value', 1)
                end
                show_bedgrd
            case '5'
                if get(surfgrd_check, 'value')
                    set(surfgrd_check, 'value', 0)
                else
                    set(surfgrd_check, 'value', 1)
                end
                show_surfgrd
            case 'c'
                if get(cbfix_check, 'value')
                    set(cbfix_check, 'value', 0)
                else
                    set(cbfix_check, 'value', 1)
                end
            case 'e'
                reset_xyz
            case 'g'
                if get(grid_check, 'value')
                    set(grid_check, 'value', 0)
                else
                    set(grid_check, 'value', 1)
                end
                toggle_grid
            case 'i'
                if get(iso_check, 'value')
                    set(iso_check, 'value', 0)
                else
                    set(iso_check, 'value', 1)
                end
                show_iso
            case 'p'
                if get(pk_check, 'value')
                    set(pk_check, 'value', 0)
                else
                    set(pk_check, 'value', 1)
                end
                show_pk
            case 't'
                misctest
            case 'x'
                if get(xfix_check, 'value')
                    set(xfix_check, 'value', 0)
                else
                    set(xfix_check, 'value', 1)
                end
            case 'y'
                if get(yfix_check, 'value')
                    set(yfix_check, 'value', 0)
                else
                    set(yfix_check, 'value', 1)
                end
            case 'z'
                if get(zfix_check, 'value')
                    set(zfix_check, 'value', 0)
                else
                    set(zfix_check, 'value', 1)
                end
            case 'spacebar'
                switch curr_dim
                    case '2D'
                        set(dim_group, 'selectedobject', dim_check(2))
                        curr_dim ...
                            = '3D';
                        [az2, el2] ...
                            = view;
                        view(az3, el3)
                        set(az_edit, 'string', az3)
                        set(el_edit, 'string', el3)                        
                    case '3D'
                        set(dim_group, 'selectedobject', dim_check(1))
                        curr_dim ...
                            = '2D';
                        [az3, el3] ...
                            = view;
                        if any([az2 el2])
                            view(az2, el2)
                        else
                            view(2)
                        end
                        set(az_edit, 'string', az2)
                        set(el_edit, 'string', el2)
                end
        end
    end

%% Test something

    function misctest(source, eventdata)
        
    end
%%
end