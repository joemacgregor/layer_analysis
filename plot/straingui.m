function straingui
% STRAINGUI GUI for plotting strain-rate model parameters.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/07/14

age_max                     = 9e3;

%   variable name       full name                               default min/max units/title
plots                       = ...
   {'accum_dj'          'accumulation rate (D-J)'               [0      0.5]    'm/a'           '$\dot{b}_{D-J}}$';
    'accum_dj_melt'     'accumulation rate (D-J + melt)'        [0      0.5]    'm/a'           '$\dot{b}_{D-J+melt}$';
    'accum_nye'         'accumulation rate (Nye)'               [0      0.5]    'm/a'           '$\dot{b}_{Nye}$';
    'accum_nye_melt'    'accumulation rate (Nye + melt)'        [0      0.5]    'm/a'           '$\dot{b}_{Nye+melt}$';
    'accum_nye_start'   'accumulation rate (Nye start)'         [0      0.5]    'm/a'           '$\dot{b}_{Nye start}$';
    'accum_shallow'     'accumulation rate (shallow strain)'    [0      0.5]    'm/a'           '$\dot{b}_{shallow}$';
    'melt_bed_dj'       'basal melt rate (D-J + melt)'          [-0.1   0.1]    'm/a'           '$\dot{m}_{D-J+melt}$';
    'melt_bed'          'basal melt rate (Nye + melt)'          [-0.1   0.1]    'm/a'           '$\dot{m}_{Nye+melt}$';
    'rate_elev_change'  'elevation change (Holocene)'           [-0.1   0]      'm/a'           '$dh/dt$';
    'flux_div'          'flux divergence'                       [-0.5   0.5]    'm/a'           '$\nabla \cdot (z u_s)$';
    'mass_con_diff'     'mass conservation difference'          [-0.5   0.5]    'm/a'           '${\Delta}q$';
    'mass_con_diff_rel' 'mass conservation difference (rel.)'   [-10    10]     '%'             '${\Delta}q/|u_s|$';
    'mass_con_diff_rhs' 'mass conservation (RHS)'               [0      0.25]   'm/a'           '$\dot{b} - w$';
    'num_layer_strain'  'number of layers used'                 [0      20]     ''              '$N$';
    'res_dj'            'residual (D-J)'                        [0      100]    ''              '$\chi^2_{D-J}$';
    'res_dj_melt'       'residual (D-J + melt)'                 [0      100]    ''              '$\chi^2_{D-J+melt}$';
    'res_nye'           'residual (Nye)'                        [0      100]    ''              '$\chi^2_{Nye}$';
    'res_nye_melt'      'residual (Nye + melt)'                 [0      100]    ''              '$\chi^2_{Nye+melt}$';
    'res_shallow'       'residual (shallow strain)'             [0      100]    ''              '$\chi^2_{shallow}$';
    'shape_fact'        'shape factor (D-J)'                    [0.75   1]      ''              '$f_{D-J}$';
    'shape_fact_melt'   'shape factor (D-J + melt)'             [0.75   1]      ''              '$f_{D-J+melt}$';
    'thick_shear'       'shear layer thickness (D-J)'           [0      2000]   'm'             '$h_{D-J}$';
    'thick_shear_melt'  'shear layer thickness (D-J + melt)'    [0      2000]   'm'             '$h_{D-J+melt}$';
    'frac_slide'        'sliding fraction (D-J + melt)'         [0      0.1]    ''              '$s_{D-J+melt}$';
    'strain_rate'       'shallow vertical strain rate'          [0      2.5e-4] '1/a'           '$\dot{\epsilon}$';
    'vel_vert'          'vertical velocity'                     [0      0.2]    'm/a'           '$w$';
    'slope_bed'         'local bed slope'                       [-50    50]     'm/km'          '$dz_b/dx$';
    'slope_bed_predict' 'local bed slope (predicted)'           [-50    50]     'm/km'          '$dz_b/dx$';
    'atten_rate'        'depth-averaged attenuation rate'       [5      30]     'dB km^{-1}'    'N_a';
    'temp_iso1'         'depth-averaged temperature'            [-25    -5]     '\circC'        'T_{mean}'
    'res_slope_bed'     'residual (bed slope)'                  [-50    50]     'm/km'          '$dz_b/dx$';
    'res_shallow / num_layer_strain' ...
                        'residual (shallow) / # layers'         [0      5]      ''              '$\chi^2_{shallow} / N$'};

load mat/xy_all num_year name_trans
load mat/merge_all dist ind_decim ind_fence num_fence x_pk y_pk
load mat/core_int x_core y_core num_core
load(['mat/strain_all_' num2str(1e-3 * age_max) 'ka'], 'accum_dj', 'accum_dj_melt', 'accum_nye', 'accum_nye_melt', 'accum_nye_start', 'accum_shallow', 'melt_bed', 'melt_bed_dj', 'frac_slide', 'num_layer_strain', 'res_dj', 'res_dj_melt', 'res_shallow', 'res_nye', 'res_nye_melt', 'shape_fact', ...
    'shape_fact_melt', 'strain_rate', 'thick_shear', 'thick_shear_melt')
% load(['mat/strain_all_' num2str(1e-3 * age_max) 'ka'], 'accum_dj', 'accum_dj_grd', 'accum_dj_bound', 'accum_dj_bound_lo_grd', 'accum_dj_bound_hi_grd', 'accum_dj_melt', 'accum_dj_melt_grd', 'accum_dj_melt_bound', ...
%     'accum_dj_melt_bound_lo_grd', 'accum_dj_melt_bound_hi_grd', 'accum_nye', 'accum_nye_grd', 'accum_nye_bound', 'accum_nye_bound_lo_grd', 'accum_nye_bound_hi_grd', 'accum_nye_melt', 'accum_nye_melt_grd', ...
%     'accum_nye_melt_bound', 'accum_nye_start', 'accum_nye_start_grd', 'accum_shallow', 'accum_shallow_grd', 'accum_shallow_bound', 'accum_shallow_bound_lo_grd', 'accum_shallow_bound_hi_grd', 'melt_bed', ...
%     'melt_bed_grd', 'melt_bed_bound', 'melt_bed_bound_lo_grd', 'melt_bed_bound_hi_grd', 'melt_bed_dj', 'melt_bed_dj_grd', 'melt_bed_dj_bound', 'melt_bed_dj_lo_grd', 'melt_bed_dj_hi_grd', 'frac_slide', ...
%     'frac_slide_grd', 'frac_slide_bound', 'frac_slide_bound_lo_grd', 'frac_slide_bound_hi_grd', 'num_layer_strain', 'res_dj', 'res_dj_melt', 'res_shallow', 'res_nye', 'res_nye_melt', 'shape_fact', 'shape_fact_grd', ...
%     'shape_fact_bound', 'shape_fact_bound_lo_grd', 'shape_fact_bound_hi_grd', 'shape_fact_melt', 'shape_fact_melt_grd', 'shape_fact_melt_bound', 'shape_fact_melt_bound_lo_grd', 'shape_fact_melt_bound_hi_grd', ...
%     'strain_rate', 'strain_rate_grd', 'strain_rate_bound', 'strain_rate_bound_lo_grd', 'strain_rate_bound_hi_grd', 'thick_shear', 'thick_shear_grd', 'thick_shear_bound', 'thick_shear_bound_lo_grd', ...
%     'thick_shear_bound_hi_grd', 'thick_shear_melt', 'thick_shear_melt_grd', 'thick_shear_melt_bound', 'thick_shear_melt_bound_lo_grd', 'thick_shear_melt_bound_hi_grd', 'x_grd', 'y_grd')
% load(['mat/flux_all_' num2str(1e-3 * age_max) 'ka'], 'flux_div', 'flux_div_grd', 'mass_con_diff', 'mass_con_diff_grd', 'mass_con_diff_rel', 'mass_con_diff_rel_grd', 'mass_con_rhs', 'mass_con_rhs_grd', 'rate_elev_change', 'rate_elev_change_grd', 'vel_vert', 'vel_vert_grd')
load('mat/age_grd1', 'age_iso', 'age_norm1', 'depth_iso1', 'depth_norm', 'num_depth_norm', 'num_age_iso')
% load('mat/atten_all', 'atten_rate' ,'temp_iso1')
load mat/grl_coast num_coast x_coast y_coast
load mat/grl_basin num_basin x_basin y_basin
load mat/xy_D1 x_D1 y_D1
% load('mat/slope_all', 'ind_decim5', 'res_slope_bed', 'slope_bed', 'slope_bed_predict')

num_normal                  = size(plots, 1);

plots_tmp1                  = cell(num_depth_norm, 5);
for ii = 1:num_depth_norm
    plots_tmp1{ii, 1}       = ['age_norm1_' num2str(ii)];
    plots_tmp1{ii, 2}       = ['Age at normalized depth (' num2str(ii) '/' num2str(num_depth_norm) ')'];
    plots_tmp1{ii, 3}       = [0 130];
    plots_tmp1{ii, 4}       = 'ka';
    plots_tmp1{ii, 5}       = '$\hat(z)$';
end

plots_tmp2                  = cell(num_age_iso, 5);
for ii = 1:num_age_iso
    plots_tmp2{ii, 1}       = ['depth_iso1_' num2str(ii)];
    plots_tmp2{ii, 2}       = ['Isochrone depth (' num2str(1e-3 * age_iso(ii)) ' ka)'];
    plots_tmp2{ii, 3}       = [0 3000];
    plots_tmp2{ii, 4}       = 'm';
    plots_tmp2{ii, 5}       = '$A$';
end
plots                       = [plots; plots_tmp1; plots_tmp2];

curr_clim                   = cell2mat(plots(:, 3)); % keep track of current plotting color ranges, start with assigned values as default

% plot range
[curr_year, curr_trans, curr_subtrans, dist_mid, dist_min, ii, ind_decim_mid, jj, kk, pdj, pdjm, pn, pnm, ps, sax, tmp, x_pk, y_pk] ...
                            = deal(0);
[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);
[p_basin, p_core]           = deal(zeros(num_core, 3));
p_map                       = cell(3, 2);
for ii = 1:num_year %#ok<*FXUP>
    for jj = find(num_fence{ii})
        for kk = 1:3
            for ll = 1:2
                p_map{kk, ll}{ii}{jj} ...
                            = zeros(1, length(ind_fence{ii}{jj}));
            end
        end
%         for kk = ind_fence{ii}{jj}
%             if ~isempty(depth_iso1{ii}{jj}{kk})
%                 depth_iso1{ii}{jj}{kk} ...
%                             = depth_iso1{ii}{jj}{kk}((age_iso == age_max), :); %#ok<AGROW>
%             end
%         end
    end
end
[ax, basin_check, curr_min, curr_max, map_type_group, map_subtype_group, max_box, min_box, plot_list, p_D1, p_grd, p_unit, trans_button] ...
                            = deal(zeros(1, 3));
curr_map                    = 1;
[curr_plot, curr_map_type, curr_map_subtype] ...
                            = deal(ones(1, 3));
map_type_check              = zeros(2, 3);
map_subtype_check           = zeros(3);

%%
set(0, 'DefaultFigureWindowStyle', 'docked')
sgui                        = figure('toolbar', 'figure', 'name', ['STRAIN GUI ' num2str(1e-3 * age_max) ' ka']);
for ii = 1%:3
    ax(ii)                  = axes('position', [(0.01 + (0.33 * (ii - 1))) 0.04 1 0.85]); %#ok<LAXES>
%     ax(ii)                  = axes('position', [(-0.25 + (0.33 * (ii - 1))) 0.04 0.92 0.85]); %#ok<LAXES>    
    hold on
    for jj = 1:num_coast
        plot(x_coast{jj}(1:3:end), y_coast{jj}(1:3:end), 'color', [0.5 0.5 0.5], 'linewidth', 1)
    end
    for jj = 1:num_basin
        p_basin(jj, ii)      = plot(x_basin{jj}(1:3:end), y_basin{jj}(1:3:end), 'k', 'linewidth', 2);
    end
    for jj = 1:num_core
        p_core(jj, ii)      = plot(x_core(jj), y_core(jj), 'k^', 'markerfacecolor', 'm', 'markersize', 12);
    end
    p_D1(ii)                = plot(x_D1(1:3:end), y_D1(1:3:end), 'm.', 'markersize', 6);
    colormap(jet)
    caxis([0 1])
    axis equal
%     axis([x_min x_max y_min y_max])
    ylim([y_min y_max])
    set(gca, 'fontsize', 20, 'layer', 'top')
    colorbar('fontsize', 20)
    p_unit(ii)              = annotation('textbox', [(0.01 + (0.33 * (ii - 1))) 0.915 0.06 0.03], 'string', '', 'color', 'k', 'fontsize', 20, 'edgecolor', 'none');
    min_box(ii)             = uicontrol(sgui, 'style', 'edit', 'string', 'min', 'units', 'normalized', 'position', [(0.005 + (0.33 * (ii - 1))) 0.955 0.035 0.04], 'callback', eval(['@do_min' num2str(ii)]), 'fontsize', 18);
    max_box(ii)             = uicontrol(sgui, 'style', 'edit', 'string', 'max', 'units', 'normalized', 'position', [(0.045 + (0.33 * (ii - 1))) 0.955 0.035 0.04], 'callback', eval(['@do_max' num2str(ii)]), 'fontsize', 18);
    uicontrol(sgui, 'style', 'pushbutton', 'string', 'default ', 'units', 'normalized', 'position', [(0.085 + (0.33 * (ii - 1))) 0.955 0.04 0.04], 'callback', eval(['@do_def' num2str(ii)]), 'fontsize', 18, 'foregroundcolor', 'k')
    trans_button(ii)        = uicontrol(sgui, 'style', 'pushbutton', 'string', 'transect', 'units', 'normalized', 'position', [(0.13 + (0.33 * (ii - 1))) 0.955 0.04 0.04], 'callback', eval(['@plot_trans' num2str(ii)]), 'fontsize', 18, 'foregroundcolor', 'b');
    uicontrol(sgui, 'style', 'pushbutton', 'string', 'pop', 'units', 'normalized', 'position', [(0.01 + (0.33 * (ii - 1))) 0.905 0.03 0.04], 'callback', eval(['@pop' num2str(ii)]), 'fontsize', 18, 'foregroundcolor', 'b')
    map_type_group(ii)      = uibuttongroup('position', [(0.175 + (0.33 * (ii - 1))) 0.96 0.07 0.03], 'selectionchangefcn', eval(['@choose_map_type' num2str(ii)]));
    uicontrol(sgui, 'style', 'text', 'parent', map_type_group(ii), 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', 18)
    map_type_check(1, ii)   = uicontrol(sgui, 'style', 'radio', 'string', 'lines', 'units', 'normalized', 'position', [0.01 0.1 0.49 0.8], 'parent', map_type_group(ii), 'fontsize', 18, 'handlevisibility', 'off');
    map_type_check(2, ii)   = uicontrol(sgui, 'style', 'radio', 'string', 'grid', 'units', 'normalized', 'position', [0.51 0.1 0.45 0.8], 'parent', map_type_group(ii), 'fontsize', 18, 'handlevisibility', 'off');
    set(map_type_group(ii), 'selectedobject', map_type_check(1, ii))
    map_subtype_group(ii)   = uibuttongroup('position', [(0.25 + (0.33 * (ii - 1))) 0.96 0.08 0.03], 'selectionchangefcn', eval(['@choose_map_subtype' num2str(ii)]));
    uicontrol(sgui, 'style', 'text', 'parent', map_subtype_group(ii), 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'fontsize', 18)
    map_subtype_check(1, ii)= uicontrol(sgui, 'style', 'radio', 'string', 'std', 'units', 'normalized', 'position', [0.01 0.1 0.36 0.8], 'parent', map_subtype_group(ii), 'fontsize', 18, 'handlevisibility', 'off');
    map_subtype_check(2, ii)= uicontrol(sgui, 'style', 'radio', 'string', 'lo', 'units', 'normalized', 'position', [0.38 0.1 0.28 0.8], 'parent', map_subtype_group(ii), 'fontsize', 18, 'handlevisibility', 'off');
    map_subtype_check(3, ii)= uicontrol(sgui, 'style', 'radio', 'string', 'hi', 'units', 'normalized', 'position', [0.67 0.1 0.28 0.8], 'parent', map_subtype_group(ii), 'fontsize', 18, 'handlevisibility', 'off');
    set(map_subtype_group(ii), 'selectedobject', map_subtype_check(1, ii))
    plot_list(ii)           = uicontrol(sgui, 'style', 'popupmenu', 'string', plots(:, 2), 'units', 'normalized', 'position', [(0.04 + (0.33 * (ii - 1))) 0.90 0.25 0.04], 'callback', eval(['@do_map' num2str(ii)]), 'fontsize', 18);
    annotation('textbox', [(0.295 + (0.33 * (ii - 1))) 0.915 0.05 0.03], 'string', 'basins', 'fontsize', 18, 'color', 'k', 'edgecolor', 'none');
    basin_check(ii)         = uicontrol(sgui, 'style', 'checkbox', 'units', 'normalized', 'position', [(0.33 + (0.33 * (ii - 1))) 0.915 0.01 0.03], 'callback', eval(['@show_basin' num2str(ii)]), 'fontsize', 18, 'value', 1);
    grid on
    box on
end
linkaxes(ax)

yr_check                    = zeros(1, num_year);
for ii = 1:num_year
    annotation('textbox', [0.9 (0.95 - ((ii - 1) * 0.05)) 0.01 0.03], 'string', num2str(ii), 'fontsize', 18, 'color', 'k', 'edgecolor', 'none');
    yr_check(ii)            = uicontrol(sgui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.95 (0.95 - ((ii - 1) * 0.05)) 0.01 0.03], 'callback', @do_yr_check, 'fontsize', 18, 'value', 1, 'tag', num2str(ii));
end

%% sub-functions

    function do_min1(source, eventdata) %#ok<DEFNU,*INUSD>
        curr_map            = 1;
        do_min
    end

    function do_min2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        do_min
    end

    function do_min3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        do_min
    end

    function do_min(source, eventdata)
        curr_min(curr_map)  = str2double(get(min_box(curr_map), 'string'));
        update_range
    end

    function do_max1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        do_max
    end

    function do_max2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        do_max
    end

    function do_max3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        do_max
    end

    function do_max(source, eventdata)
        curr_max(curr_map)  = str2double(get(max_box(curr_map), 'string'));
        update_range
    end

    function do_def1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        do_def
    end

    function do_def2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        do_def
    end

    function do_def3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        do_def
    end

    function do_def(source, eventdata) % revert current map to default caxis
        curr_min(curr_map)  = plots{curr_plot(curr_map), 3}(1);
        curr_max(curr_map)  = plots{curr_plot(curr_map), 3}(2);
        update_range
    end

    function update_range(source, eventdata)
        axes(ax(curr_map))
        set([min_box(curr_map) max_box(curr_map)], 'foregroundcolor', 'r')
        pause(0.1)
        caxis([curr_min(curr_map) curr_max(curr_map)]);
        curr_clim(curr_plot(curr_map), :) ...
                            = [curr_min(curr_map) curr_max(curr_map)]; % keep track for future maps
        set([min_box(curr_map) max_box(curr_map)], 'foregroundcolor', 'k')
        pause(0.1)
    end

    function do_map1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        do_map
    end

    function do_map2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        do_map
    end

    function do_map3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        do_map
    end

    function choose_map_type1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        switch get(map_type_group(curr_map), 'selectedobject')
            case map_type_check(1, curr_map)
                curr_map_type(curr_map) ...
                            = 1;
            case map_type_check(2, curr_map)
                curr_map_type(curr_map) ...
                            = 2;
        end
        do_map
    end

    function choose_map_type2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        switch get(map_type_group(curr_map), 'selectedobject')
            case map_type_check(1, curr_map)
                curr_map_type(curr_map) ...
                            = 1;
            case map_type_check(2, curr_map)
                curr_map_type(curr_map) ...
                            = 2;
        end
        do_map
    end

    function choose_map_type3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        switch get(map_type_group(curr_map), 'selectedobject')
            case map_type_check(1, curr_map)
                curr_map_type(curr_map) ...
                            = 1;
            case map_type_check(2, curr_map)
                curr_map_type(curr_map) ...
                            = 2;
        end
        do_map
    end

    function choose_map_subtype1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        switch get(map_subtype_group(curr_map), 'selectedobject')
            case map_subtype_check(1, curr_map)
                curr_map_subtype(curr_map) ...
                            = 1;
            case map_subtype_check(2, curr_map)
                curr_map_subtype(curr_map) ...
                            = 2;
            case map_subtype_check(3, curr_map)
                curr_map_subtype(curr_map) ...
                            = 3;
        end
        do_map
    end

    function choose_map_subtype2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        switch get(map_subtype_group(curr_map), 'selectedobject')
            case map_subtype_check(1, curr_map)
                curr_map_subtype(curr_map) ...
                            = 1;
            case map_subtype_check(2, curr_map)
                curr_map_subtype(curr_map) ...
                            = 2;
            case map_subtype_check(3, curr_map)
                curr_map_subtype(curr_map) ...
                            = 3;
        end
        do_map
    end

    function choose_map_subtype3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        switch get(map_subtype_group(curr_map), 'selectedobject')
            case map_subtype_check(1, curr_map)
                curr_map_subtype(curr_map) ...
                            = 1;
            case map_subtype_check(2, curr_map)
                curr_map_subtype(curr_map) ...
                            = 2;
            case map_subtype_check(3, curr_map)
                curr_map_subtype(curr_map) ...
                            = 3;
        end
        do_map
    end

    function do_map(source, eventdata)
        
        if (curr_map_type(curr_map) == 2)
            switch curr_map_subtype(curr_map)
                case 1
                    if ~exist([plots{curr_plot(curr_map), 1} '_grd'], 'var')
                        curr_map_type(curr_map) = 1;
                        set(map_type_group(curr_map), 'selectedobject', map_type_check(1, curr_map))
                        set(map_subtype_group(curr_map), 'selectedobject', map_subtype_check(1, curr_map))
                        return
                    end
                case 2
                    if ~exist([plots{curr_plot(curr_map), 1} '_bound_lo_grd'], 'var')
                        [curr_map_type(curr_map), curr_map_subtype(curr_map)] ...
                            = deal(1);
                        set(map_type_group(curr_map), 'selectedobject', map_type_check(1, curr_map))
                        set(map_subtype_group(curr_map), 'selectedobject', map_subtype_check(1, curr_map))
                        return
                    end
                case 3
                    if ~exist([plots{curr_plot(curr_map), 1} '_bound_hi_grd'], 'var')
                        [curr_map_type(curr_map), curr_map_subtype(curr_map)] ...
                            = deal(1);
                        set(map_type_group(curr_map), 'selectedobject', map_type_check(1, curr_map))
                        set(map_subtype_group(curr_map), 'selectedobject', map_subtype_check(1, curr_map))
                        return
                    end
            end
        elseif (curr_map_subtype(curr_map) > 1)
            if ~exist([plots{curr_plot(curr_map), 1} '_bound'], 'var')
                curr_map_subtype(curr_map) = 1;
                set(map_subtype_group(curr_map), 'selectedobject', map_subtype_check(1, curr_map))
                return
            end
        end
        
        set(plot_list(curr_map), 'foregroundcolor', 'r')
        pause(0.1)
        axes(ax(curr_map))
        
        for ii = 1:num_year %#ok<*FXUP>
            for jj = find(num_fence{ii})
                for kk = 1:2
                    if (any(p_map{curr_map, kk}{ii}{jj}) && any(ishandle(p_map{curr_map, kk}{ii}{jj}))) % delete any previous map handles if they exist
                        delete(p_map{curr_map, kk}{ii}{jj}(logical(p_map{curr_map, kk}{ii}{jj}) & ishandle(p_map{curr_map, kk}{ii}{jj})))
                    end
                end
            end
        end
        if (logical(p_grd(curr_map)) && ishandle(p_grd(curr_map)))
            delete(p_grd(curr_map))
        end
        
        curr_plot(curr_map) = get(plot_list(curr_map), 'value'); % current map number to plot based on sequence in plots
        [curr_min(curr_map), curr_max(curr_map)] ...
                            = deal(curr_clim(curr_plot(curr_map), 1), curr_clim(curr_plot(curr_map), 2));
        caxis([curr_min(curr_map) curr_max(curr_map)]); % use default limits for a new map
        
        switch curr_map_type(curr_map)
            case 1
                for ii = 1:num_year %#ok<*FXUP>
%                     if ~get(yr_check(ii), 'value')
%                         continue
%                     end
                    for jj = find(num_fence{ii})
                        for kk = find(ind_fence{ii}{jj})
                            if (curr_plot(curr_map) < num_normal)
                                if isempty(eval([plots{curr_plot(curr_map), 1} '{ii}{jj}{kk}']))
                                    continue
                                end
                                if (length(eval([plots{curr_plot(curr_map), 1} '{ii}{jj}{kk}'])) == (length(ind_decim{ii}{jj}{kk}) - 1))
                                    ind_decim_mid ...
                                        = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));
                                else
                                    ind_decim_mid ...
                                        = round(ind_decim5{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim5{ii}{jj}{kk}) ./ 2));
                                end
                            elseif (curr_plot(curr_map) == num_normal)
                                    ind_decim_mid ...
                                        = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));                                
                            elseif (curr_plot(curr_map) <= (num_normal + num_depth_norm))
                                if isempty(age_norm1{ii}{jj}{kk})
                                    continue
                                end
                                ind_decim_mid ...
                                        = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));
                            else
                                if isempty(depth_iso1{ii}{jj}{kk})
                                    continue
                                end
                                ind_decim_mid ...
                                        = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));
                            end
%                             p_map{curr_map, 1}{ii}{jj}(kk) ...
%                                 = plot(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 'k', 'linewidth', 1, 'tag', num2str([ii jj kk]));
                            switch curr_map_subtype(curr_map)
                                case 1
                                    if (curr_plot(curr_map) < num_normal)
                                        p_map{curr_map, 2}{ii}{jj}(kk) ...
                                            = scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, eval([plots{curr_plot(curr_map), 1} '{ii}{jj}{kk}']), 'filled', 'tag', num2str([ii jj kk]));
                                    elseif (curr_plot(curr_map) == num_normal)
                                        p_map{curr_map, 2}{ii}{jj}(kk) ...
                                            = scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, (res_shallow{ii}{jj}{kk} ./ num_layer_strain{ii}{jj}{kk}), 'filled', 'tag', num2str([ii jj kk]));                                        
                                    elseif (curr_plot(curr_map) <= (num_normal + num_depth_norm))
                                        p_map{curr_map, 2}{ii}{jj}(kk) ...
                                            = scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, (1e-3 .* age_norm1{ii}{jj}{kk}((curr_plot(curr_map) - num_normal), :)), 'filled', 'tag', num2str([ii jj kk]));
                                    else
                                        p_map{curr_map, 2}{ii}{jj}(kk) ...
                                            = scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, depth_iso1{ii}{jj}{kk}((curr_plot(curr_map) - (num_normal + num_depth_norm)), :), 'filled', 'tag', num2str([ii jj kk]));
                                    end
                                case 2
                                    p_map{curr_map, 2}{ii}{jj}(kk) ...
                                        = scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, eval([plots{curr_plot(curr_map), 1} '_bound{ii}{jj}{kk}(1, :)']), 'filled', 'tag', num2str([ii jj kk]));
                                case 3
                                    p_map{curr_map, 2}{ii}{jj}(kk) ...
                                        = scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, eval([plots{curr_plot(curr_map), 1} '_bound{ii}{jj}{kk}(2, :)']), 'filled', 'tag', num2str([ii jj kk]));
                            end
                        end
                    end
                end
            case 2
                switch curr_map_subtype(curr_map)
                    case 1
                        p_grd(curr_map) ...
                            = pcolor(x_grd, y_grd, eval([plots{curr_plot(curr_map), 1} '_grd']));
                    case 2
                        p_grd(curr_map) ...
                            = pcolor(x_grd, y_grd, eval([plots{curr_plot(curr_map), 1} '_bound_lo_grd']));
                    case 3
                        p_grd(curr_map) ...
                            = pcolor(x_grd, y_grd, eval([plots{curr_plot(curr_map), 1} '_bound_hi_grd']));
                end
                shading interp
        end
        
        if (logical(p_unit(curr_map)) && ishandle(p_unit(curr_map)))
            if ~isempty(plots{curr_plot(curr_map), 4})
                set(p_unit(curr_map), 'string', ['(' plots{curr_plot(curr_map), 4} ')'])
            else
                set(p_unit(curr_map), 'string', '')
            end
        end
        set(min_box(curr_map), 'string', num2str(curr_clim(curr_plot(curr_map), 1)))
        set(max_box(curr_map), 'string', num2str(curr_clim(curr_plot(curr_map), 2)))
        show_basin
        pause(0.1)
        set(plot_list(curr_map), 'foregroundcolor', 'k')
    end

    function do_yr_check(source, ~)
        ii                  = str2double(get(source, 'tag'));
        if get(yr_check(ii), 'value')
            for jj = find(num_fence{ii})
                for kk = find(ind_fence{ii}{jj})
                    if (logical(p_map{curr_map, 2}{ii}{jj}(kk)) && ishandle(p_map{curr_map, 2}{ii}{jj}(kk)))
                        set(p_map{curr_map, 2}{ii}{jj}(kk), 'visible', 'on')
                    end
                end
            end
        else
            for jj = find(num_fence{ii})
                for kk = find(ind_fence{ii}{jj})
                    if (logical(p_map{curr_map, 2}{ii}{jj}(kk)) && ishandle(p_map{curr_map, 2}{ii}{jj}(kk)))
                        set(p_map{curr_map, 2}{ii}{jj}(kk), 'visible', 'off')
                    end
                end
            end
        end
    end

    function plot_trans1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        plot_trans
    end

    function plot_trans2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        plot_trans
    end

    function plot_trans3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        plot_trans
    end

    function plot_trans(source, eventdata) % make transect plot
        
        set(trans_button(curr_map), 'foregroundcolor', 'r')
        pause(0.1)
        
        axes(ax(curr_map))
        [x_pk, y_pk]        = ginput(1);
        dist_min            = Inf;
        for ii = 1:num_year
            for jj = find(num_fence{ii})
                for kk = find(ind_fence{ii}{jj})
                    if (logical(p_map{curr_map, 1}{ii}{jj}(kk)) && ishandle(p_map{curr_map, 1}{ii}{jj}(kk)))
                        if isempty(find((sqrt(((x_pk{ii}{jj}{kk} - x_pk) .^ 2) + ((y_pk{ii}{jj}{kk} - y_pk) .^ 2)) < dist_min), 1))
                            [curr_year, curr_trans, curr_subtrans] ...
                                = deal(ii, jj, kk);
                            dist_min ...
                                = min(sqrt(((x_pk{ii}{jj}{kk} - x_pk) .^ 2) + ((y_pk{ii}{jj}{kk} - y_pk) .^ 2)));
                        end
                    end
                    for ll = 1:3
                        if (logical(p_map{ll, 1}{ii}{jj}(kk)) && ishandle(p_map{ll, 1}{ii}{jj}(kk)))
                            set(p_map{ll, 1}{ii}{jj}(kk), 'linewidth', 1)
                        end
                        if (logical(p_map{ll, 2}{ii}{jj}(kk)) && ishandle(p_map{ll, 2}{ii}{jj}(kk)))
                            set(p_map{ll, 2}{ii}{jj}(kk), 'sizedata', 36)
                        end
                    end
                end
            end
        end
        
        for ii = 1:3
            if (logical(p_map{ii, 1}{curr_year}{curr_trans}(curr_subtrans)) && ishandle(p_map{ii, 1}{curr_year}{curr_trans}(curr_subtrans)))
                set(p_map{ii, 1}{curr_year}{curr_trans}(curr_subtrans), 'linewidth', 3)
            end
            if (logical(p_map{ii, 2}{curr_year}{curr_trans}(curr_subtrans)) && ishandle(p_map{ii, 2}{curr_year}{curr_trans}(curr_subtrans)))
                set(p_map{ii, 2}{curr_year}{curr_trans}(curr_subtrans), 'sizedata', 144)
            end
        end
        pause(0.1)
        
        ind_decim_mid           = round(ind_decim{curr_year}{curr_trans}{curr_subtrans}(1:(end - 1)) + (diff(ind_decim{curr_year}{curr_trans}{curr_subtrans}) ./ 2));
        dist_mid                = dist{curr_year}{curr_trans}{curr_subtrans}(ind_decim_mid);
        
        figure
        sax(1)                  = subplot(311);
        hold on
%         tmp                     = accum_dj_bound{curr_year}{curr_trans}{curr_subtrans};
%         fill([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], (1e2 .* [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')]), [0.75 0.75 1], 'edgecolor', 'none')
%         tmp                     = accum_dj_melt_bound{curr_year}{curr_trans}{curr_subtrans};
%         fill([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], (1e2 .* [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')]), [0.9 0.9 1], 'edgecolor', 'none')
%         tmp                     = accum_nye_bound{curr_year}{curr_trans}{curr_subtrans};
%         fill([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], (1e2 .* [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')]), [1 0.75 0.75], 'edgecolor', 'none')
%         tmp                     = accum_nye_melt_bound{curr_year}{curr_trans}{curr_subtrans};
%         fill([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], (1e2 .* [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')]), [0.75 1 0.75], 'edgecolor', 'none')
%         tmp                     = accum_shallow_bound{curr_year}{curr_trans}{curr_subtrans};
%         fill([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], (1e2 .* [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')]), [1 0.75 1], 'edgecolor', 'none')
        pdj                     = plot(dist_mid, (1e2 .* accum_dj{curr_year}{curr_trans}{curr_subtrans}), 'b', 'linewidth', 2);
        pdjm                    = plot(dist_mid, (1e2 .* accum_dj_melt{curr_year}{curr_trans}{curr_subtrans}), 'c', 'linewidth', 2);
        pn                      = plot(dist_mid, (1e2 .* accum_nye{curr_year}{curr_trans}{curr_subtrans}), 'r', 'linewidth', 2);
        pnm                     = plot(dist_mid, (1e2 .* accum_nye_melt{curr_year}{curr_trans}{curr_subtrans}), 'm', 'linewidth', 2);
        ps                      = plot(dist_mid, (1e2 .* accum_shallow{curr_year}{curr_trans}{curr_subtrans}), 'k', 'linewidth', 2);
        xlim(dist_mid([1 end]))
        set(gca, 'fontsize', 20, 'xgrid', 'on')
        ylabel({'Accumulation rate'; '(cm a^{-1})'})
        legend([pdj pdjm pn pnm ps], 'Dansgaard-Johnsen', 'Dansgaard-Johnsen + melt', 'Nye', 'Nye + basal melt', 'shallow')
        box on
        title(name_trans{curr_year}{curr_trans}, 'interpreter', 'none')
        sax(2)                  = subplot(312);
%         tmp                     = melt_bed_bound{curr_year}{curr_trans}{curr_subtrans};
%         fill([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], (1e3 .* [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')]), [0.75 0.75 1], 'edgecolor', 'none')
%         tmp                     = melt_bed_dj_bound{curr_year}{curr_trans}{curr_subtrans};
%         fill([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], (1e3 .* [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')]), [0.9 0.9 1], 'edgecolor', 'none')
        plot(dist_mid, (1e3 .* melt_bed{curr_year}{curr_trans}{curr_subtrans}), 'b', 'linewidth', 2);
        plot(dist_mid, (1e3 .* melt_bed_dj{curr_year}{curr_trans}{curr_subtrans}), 'c', 'linewidth', 2);        
        axis([dist_mid([1 end]) 0 300])
        set(gca, 'fontsize', 20, 'xgrid', 'on')
        xlabel('Distance (km)')
        ylabel({'Basal melt rate'; '(mm a^{-1})'}, 'color', 'r')
        box off
        sax(3)                  = axes('position', get(gca, 'position'), 'color', 'none', 'xaxislocation', 'top', 'yaxislocation', 'right', 'fontsize', 20, 'xticklabel', {});
%         tmp                     = thick_shear_bound{curr_year}{curr_trans}{curr_subtrans};
%         patch([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')], [1 0.75 0.75], 'edgecolor', 'none')
%         tmp                     = thick_shear_melt_bound{curr_year}{curr_trans}{curr_subtrans};
%         patch([dist_mid(~isnan(tmp(1, :)))'; flipud(dist_mid(~isnan(tmp(2, :)))')], [tmp(1, ~isnan(tmp(1, :)))'; flipud(tmp(2, ~isnan(tmp(2, :)))')], [1 0.75 1], 'edgecolor', 'none')
        line(dist_mid, thick_shear{curr_year}{curr_trans}{curr_subtrans}, 'color', 'r', 'linewidth', 2);
        line(dist_mid, thick_shear_melt{curr_year}{curr_trans}{curr_subtrans}, 'color', 'm', 'linewidth', 2);        
        axis([dist_mid([1 end]) 0 3500])
        ylabel({'Thickness'; '(m)'}, 'color', 'b');
        sax(4)                  = subplot(313);
        hold on
        plot(dist_mid, res_dj{curr_year}{curr_trans}{curr_subtrans}, 'b', 'linewidth', 2)
        plot(dist_mid, res_dj{curr_year}{curr_trans}{curr_subtrans}, 'c', 'linewidth', 2)
        plot(dist_mid, res_nye{curr_year}{curr_trans}{curr_subtrans}, 'r', 'linewidth', 2)
        plot(dist_mid, res_nye_melt{curr_year}{curr_trans}{curr_subtrans}, 'm', 'linewidth', 2)
        plot(dist_mid, res_shallow{curr_year}{curr_trans}{curr_subtrans}, 'color', 'k', 'linewidth', 2)
        xlim(dist_mid([1 end]))
        set(gca, 'fontsize', 20, 'xgrid', 'on')
        ylabel({'\chi^2'; '(non-dimensional)'})
        box on
        linkaxes(sax, 'x')
        
        pause(0.1)
        set(trans_button(curr_map), 'foregroundcolor', 'b')
    end

%% Show drainage basins

    function show_basin1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        show_basin
    end

    function show_basin2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        show_basin
    end

    function show_basin3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        show_basin
    end

    function show_basin(source, eventdata)
        if get(basin_check(curr_map), 'value')
            if (any(p_basin(:, curr_map)) && any(ishandle(p_basin(:, curr_map))))
                set(p_basin((logical(p_basin(:, curr_map)) & ishandle(p_basin(:, curr_map))), curr_map), 'visible', 'on')
                uistack(p_basin((logical(p_basin(:, curr_map)) & ishandle(p_basin(:, curr_map))), curr_map), 'top')
            end
            if (any(p_D1(curr_map)) && any(ishandle(p_D1(curr_map))))
                set(p_D1((logical(p_D1(curr_map)) & ishandle(p_D1(curr_map))), curr_map), 'visible', 'on')
                uistack(p_D1((logical(p_D1(curr_map)) & ishandle(p_D1(curr_map))), curr_map), 'top')
            end
        else
            if (any(p_basin(:, curr_map)) && any(ishandle(p_basin(:, curr_map))))
                set(p_basin((logical(p_basin(:, curr_map)) & ishandle(p_basin(:, curr_map))), curr_map), 'visible', 'off')
            end
            if (any(p_D1(curr_map)) && any(ishandle(p_D1(curr_map))))
                set(p_D1((logical(p_D1(curr_map)) & ishandle(p_D1(curr_map))), curr_map), 'visible', 'off')
            end
        end
        if (any(p_core(:, curr_map)) && any(ishandle(p_core(:, curr_map))))
            uistack(p_core((logical(p_core(:, curr_map)) & ishandle(p_core(:, curr_map))), curr_map), 'top')
        end
    end

%% Pop current map

    function pop1(source, eventdata) %#ok<DEFNU>
        curr_map            = 1;
        pop_map
    end

    function pop2(source, eventdata) %#ok<DEFNU>
        curr_map            = 2;
        pop_map
    end

    function pop3(source, eventdata) %#ok<DEFNU>
        curr_map            = 3;
        pop_map
    end

    function pop_map(source, eventdata)
        set(0, 'DefaultFigureWindowStyle', 'default')
        figure('position', [50 50 600 800])
        hold on
        colormap(jet)
        caxis([curr_min(curr_map) curr_max(curr_map)])
        switch curr_map_type(curr_map)
            case 1
                for ii = 1:num_year %#ok<*FXUP>
                    for jj = find(num_fence{ii})
                        for kk = find(ind_fence{ii}{jj})
                            if isempty(eval([plots{curr_plot(curr_map), 1} '{ii}{jj}{kk}']))
                                continue
                            end
                            ind_decim_mid ...
                                = round(ind_decim{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim{ii}{jj}{kk}) ./ 2));
                            plot(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 'k', 'linewidth', 1, 'tag', num2str([ii jj kk]));
                            switch curr_map_subtype(curr_map)
                                case 1
                                    scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, eval([plots{curr_plot(curr_map), 1} '{ii}{jj}{kk}']), 'filled')
                                case 2
                                    scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, eval([plots{curr_plot(curr_map), 1} '_bound{ii}{jj}{kk}(1, :)']), 'filled')
                                case 3
                                    scatter(x_pk{ii}{jj}{kk}(ind_decim_mid), y_pk{ii}{jj}{kk}(ind_decim_mid), 36, eval([plots{curr_plot(curr_map), 1} '_bound{ii}{jj}{kk}(2, :)']), 'filled')
                            end
                        end
                    end
                end
            case 2
                switch curr_map_subtype(curr_map)
                    case 1
                        pcolor(x_grd, y_grd, eval([plots{curr_plot(curr_map), 1} '_grd']))
                    case 2
                        pcolor(x_grd, y_grd, eval([plots{curr_plot(curr_map), 1} '_bound_lo_grd']))
                    case 3
                        pcolor(x_grd, y_grd, eval([plots{curr_plot(curr_map), 1} '_bound_hi_grd']))
                end
                shading interp
        end
        for ii = 1:num_coast
            plot(x_coast{ii}, y_coast{ii}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for ii = 1:num_basin
            plot(x_basin{ii}, y_basin{ii}, 'k', 'linewidth', 2)
        end
        for ii = 1:num_core
            plot(x_core(ii), y_core(ii), 'k^', 'markerfacecolor', 'm', 'markersize', 12)
        end
        axis equal
        axis([x_min x_max y_min y_max])
        set(gca, 'fontsize', 20, 'layer', 'top')
        xlabel('Polar stereographic X (km)')
        ylabel('Polar stereographic Y (km)')
        title(plots{curr_plot(curr_map), 5}, 'fontweight', 'bold', 'interpreter', 'latex')
        text(850, -550, ['(' plots{curr_plot(curr_map), 4} ')'], 'color', 'k', 'fontsize', 20)
        colorbar('fontsize', 20)
        grid on
        box on
        set(0, 'DefaultFigureWindowStyle', 'docked')        
    end
end