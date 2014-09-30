% GRIS_STRAT_FIG Figures for gris_strat manuscript.
% 
% Joe MacGregor
% Last updated: 09/30/14

clear

plotting                    = false;

load mat/speed_greenland_2007-2010 x_grd y_grd speed_x speed_y
[x_speed_grd, y_speed_grd]  = deal(x_grd, y_grd);
clear x_grd y_grd
load mat/xy_all num_trans num_year x y
load mat/pkstat ind_decim num_pk num_trans_pk x_pk y_pk
ind_decim_pk                = ind_decim;
load mat/merge_all depth_smooth ind_fence num_fence x_pk y_pk
load mat/grl_coast num_coast x_coast y_coast
load mat/grl_basin num_basin x_basin y_basin
load mat/greenland_bed_v3 elev_bed_uncert elev_surf_uncert mask_gris mask_land x y
load mat/greenland_mc_bed_1km thick
thick_grd                   = thick;
clear thick
load mat/trans_outside trans_outside
load mat/core_int name_core num_core x_core y_core
load mat/block_all num_block num_subtrans x_block y_block
load mat/slope_all ind_decim5 res_slope_bed slope_bed slope_bed_predict
load mat/date_all age
% load mat/age_grd2 age_iso age_norm2 age_uncert_norm2 depth_iso2 depth_norm depth_uncert_iso2 x_grd y_grd
load mat/age_grd2_krige age_iso age_uncert_norm2 age_norm2_uncert depth_iso2 depth_iso2_uncert depth_norm depth_uncert_iso2 num_age_iso x_grd y_grd
load mat/age_grd2_krige_alt age_norm2_alt
age_norm2                   = age_norm2_alt;
% clear age_norm2_alt

[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);

% trim speed grid
[speed_x, speed_y]          = deal(speed_x(((y_speed_grd(:, 1) >= y_min) & (y_speed_grd(:, 1) <= y_max)), ((x_speed_grd(1, :) >= x_min) & (x_speed_grd(1, :) <= x_max))), ...
                                   speed_y(((y_speed_grd(:, 1) >= y_min) & (y_speed_grd(:, 1) <= y_max)), ((x_speed_grd(1, :) >= x_min) & (x_speed_grd(1, :) <= x_max))));
[x_speed_grd, y_speed_grd]  = deal(x_speed_grd(((y_speed_grd(:, 1) >= y_min) & (y_speed_grd(:, 1) <= y_max)), ((x_speed_grd(1, :) >= x_min) & (x_speed_grd(1, :) <= x_max))), ...
                                   y_speed_grd(((y_speed_grd(:, 1) >= y_min) & (y_speed_grd(:, 1) <= y_max)), ((x_speed_grd(1, :) >= x_min) & (x_speed_grd(1, :) <= x_max))));
[speed_x, speed_y]          = deal(speed_x(1:2:end, 1:2:end), speed_y(1:2:end, 1:2:end));
speed                       = sqrt((speed_x .^ 2) + (speed_y .^ 2));

[elev_bed_uncert, elev_surf_uncert, mask_land, mask_gris] ...
                            = deal(elev_bed_uncert(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))), elev_surf_uncert(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))), ...
                                   mask_land(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))), mask_gris(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))));
[elev_bed_uncert(elev_bed_uncert == 0), elev_surf_uncert(elev_surf_uncert == 0)] ...
                            = deal(NaN);
                               
mask_combo                  = double(mask_land);
mask_combo(mask_gris)       = 2;

[x_core_label, y_core_label]= deal([-620 -14 -5 264 -140 -48], [-1233 -2630 -1794 -1794 -1252 -1521]);

% order in which to analyze years/campaigns based on their overall data quality
year_ord                    = [17 19 20 6 7 8 4 9 5 2 3 12 1 13 15 16 18 11 10 14];

[depth_iso2, depth_iso2_uncert, depth_uncert_iso2] ...
                            = deal((depth_iso2 .* thick_grd(:, :, ones(1, 1, num_age_iso))), (depth_iso2_uncert .* thick_grd(:, :, ones(1, 1, num_age_iso))), (depth_uncert_iso2 .* thick_grd(:, :, ones(1, 1, num_age_iso))));

[age_norm2_uncert((age_norm2_uncert == 0) | (age_norm2 == 0) | isnan(age_norm2)), age_uncert_norm2((age_uncert_norm2 == 0) | (age_norm2 == 0) | isnan(age_norm2)), ...
 depth_iso2_uncert((depth_iso2_uncert == 0) | (depth_iso2 == 0) | isnan(depth_iso2)), depth_uncert_iso2((depth_uncert_iso2 == 0) | (depth_iso2 == 0) | isnan(depth_iso2))] ...
                            = deal(NaN);
[age_norm2(age_norm2 == 0), depth_iso2(depth_iso2 == 0)] ...
                            = deal(NaN);

thick_eemian_uncert         = sqrt((elev_bed_uncert .^ 2) + (elev_surf_uncert .^ 2) + (squeeze(depth_uncert_iso2(:, :, end)) .^ 2));

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>

%%
    set(0, 'DefaultFigureWindowStyle', 'docked')

%% Figure 1
    figure('position', [100 100 1100 800])
    subplot('position', [0 0.03 0.48 0.91])
    name_year_nice          = {'1993 P3' '1995 P3' '1996 P3' '1997 P3' '1998 P3' '1999 P3' '2001 P3' '2002 P3' '2003 P3' '2005 TO' '2006 TO' '2007 P3' '2008 TO' '2009 TO' '2010 DC8' '2010 P3' '2011 P3' '2011 TO' '2012 P3' '2013 P3' ''};
    hold on
    colors                  = colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(num_year)]);
    axis equal
    axis([-632 846 -3344 -670])
    peb                     = pcolor(x_grd(1, :), y_grd(:, 1), (mask_combo + 1));
    shading interp
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    for ii = 1:length(x_paral)
        plot((1e-3 .* x_paral{ii}), (1e-3 .* y_paral{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:length(x_merid)
        plot((1e-3 .* x_merid{ii}), (1e-3 .* y_merid{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:num_coast
        plot(x_coast{ii}, y_coast{ii}, 'color', 'k', 'linewidth', 1)
    end
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            if ~isempty(intersect(trans_outside, [ii jj], 'rows'))
                continue
            end
            plot(x{ii}{jj}, y{ii}{jj}, '.', 'markersize', 6, 'color', colors((ii + 3), :), 'tag', num2str([ii jj]))
        end
    end
    for ii = 1:num_core
        plot(x_core(ii), y_core(ii), 'k^', 'markersize', 10, 'markerfacecolor', 'm')
        text(x_core_label(ii), y_core_label(ii), name_core{ii}, 'color', 'm', 'fontsize', 18, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
    end
    set(gca, 'fontsize', 20, 'xtick', [], 'ytick', [])
    caxis([1 (num_year + 3)])
    cb1                     = colorbar('fontsize', 20, 'ylim', [3 (num_year + 3)], 'ytick', 3:(num_year + 2), 'yticklabel', name_year_nice);
    set(cb1, 'visible', 'off')
    cb2                     = cbarf([4 24], 4:24, 'vertical', 'linear');
    set(cb2, 'fontsize', 20, 'ytick', 4:24, 'yticklabel', name_year_nice, 'position', get(cb1, 'position'));
    fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
    fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
    text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
    text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
    text(-595, -775, '(a)', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
    text(-542, -3238, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
    text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
    text(552, -2743, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
    text(759, -2975, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
    grid off
    box on
    freezeColors
    subplot('position', [0.50 0.03 0.48 0.91])
    hold on
    colors                  = colormap([([135 206 235] ./ 255); 
                                       ([159 89  39] ./ 255); 
                                        0.9  0.9  0.9;
                                        0    0    1;
                                        0    1    1;
                                        0    1    0;
                                        1    1    0;
                                        1    0    0;
                                        0.5  0    0]);
    axis equal
    axis([-632 846 -3344 -670])
    peb                     = pcolor(x_grd(1, :), y_grd(:, 1), (mask_combo + 1));
    shading interp    
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    for ii = 1:length(x_paral)
        plot((1e-3 .* x_paral{ii}), (1e-3 .* y_paral{ii}), 'k--', 'linewidth', 1)
    end
    for ii = 1:length(x_merid)
        plot((1e-3 .* x_merid{ii}), (1e-3 .* y_merid{ii}), 'k--', 'linewidth', 1)
    end                       
    for ii = 1:num_coast
        plot(x_coast{ii}, y_coast{ii}, 'color', 'k', 'linewidth', 1)
    end
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            if ~isempty(intersect(trans_outside, [ii jj], 'rows'))
                continue
            end
            plot(x{ii}{jj}, y{ii}{jj}, 'k.', 'markersize', 2)
        end
    end
    [x_cat, y_cat, plot_cat]= deal([]);
    range_tmp               = 0:5:30;
    for ii = 1:num_year
        for jj = 1:num_trans(ii)
            if iscell(x_pk{ii}{jj})
                for kk = 1:length(x_pk{ii}{jj})
                    if (isempty(ind_decim_pk{ii}{jj}{kk}) || all(ind_decim_pk{ii}{jj}{kk} == 0))
                        continue
                    end
                    ind_decim_mid ...
                            = round(ind_decim_pk{ii}{jj}{kk}(1:(end - 1)) + (diff(ind_decim_pk{ii}{jj}{kk}) ./ 2));
%                     if ~isempty(num_pk{ii}{jj}{kk}(ind_decim_pk{ii}{jj}{kk}))
%                         scatter(x_pk{ii}{jj}{kk}(ind_decim_pk{ii}{jj}{kk}), y_pk{ii}{jj}{kk}(ind_decim_pk{ii}{jj}{kk}), 12, (3 + interp1(2.5:5:22.5, 1:5, num_pk{ii}{jj}{kk}(ind_decim_pk{ii}{jj}{kk}), 'nearest', 'extrap')), 'filled')
                    num_pk_curr ...
                        = zeros(1, length(ind_decim_mid));
                    for ll = 1:length(ind_decim_mid)
                        if (ind_decim_pk{ii}{jj}{kk}(ll + 1) > ind_decim_pk{ii}{jj}{kk}(ll))
                            num_pk_curr(ll) ...
                                = nanmax(num_pk{ii}{jj}{kk}(ind_decim_pk{ii}{jj}{kk}(ll):ind_decim_pk{ii}{jj}{kk}(ll + 1)));
                        end
                    end
                    [x_cat, y_cat, plot_cat] ...
                            = deal([x_cat x_pk{ii}{jj}{kk}(ind_decim_mid)'], [y_cat y_pk{ii}{jj}{kk}(ind_decim_mid)'], [plot_cat num_pk_curr]);
%                     end
                end
            else
                if (isempty(ind_decim_pk{ii}{jj}) || all(ind_decim_pk{ii}{jj} == 0)) 
                    continue
                end
                ind_decim_mid ...
                            = round(ind_decim_pk{ii}{jj}(1:(end - 1)) + (diff(ind_decim_pk{ii}{jj}) ./ 2));
%                 if ~isempty(num_pk{ii}{jj}(ind_decim_pk{ii}{jj}))
%                     scatter(x_pk{ii}{jj}(ind_decim_pk{ii}{jj}), y_pk{ii}{jj}(ind_decim_pk{ii}{jj}), 12, (3 + interp1(2.5:5:22.5, 1:5, num_pk{ii}{jj}(ind_decim_pk{ii}{jj}), 'nearest', 'extrap')), 'filled')
                    num_pk_curr ...
                        = zeros(1, length(ind_decim_mid));
                    for kk = 1:length(ind_decim_mid)
                        if (ind_decim_pk{ii}{jj}(kk + 1) > ind_decim_pk{ii}{jj}(kk))
                            num_pk_curr(kk) ...
                                = nanmax(num_pk{ii}{jj}(ind_decim_pk{ii}{jj}(kk):ind_decim_pk{ii}{jj}(kk + 1)));
                        end
                    end
                    [x_cat, y_cat, plot_cat] ...
                            = deal([x_cat x_pk{ii}{jj}(ind_decim_mid)'], [y_cat y_pk{ii}{jj}(ind_decim_mid)'], [plot_cat num_pk_curr]);
%                 end
            end
        end
    end
    plot_cat            = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:6, plot_cat, 'nearest', 'extrap');
    xy_plot_cat         = sortrows([x_cat' y_cat' plot_cat'], 3);
    xy_plot_cat         = xy_plot_cat(~isnan(xy_plot_cat(:, 3)), :);
    for ii = 1:size(xy_plot_cat, 1)
        line(xy_plot_cat(ii, 1), xy_plot_cat(ii, 2), 'marker', '.', 'markersize', 12, 'color', colors(xy_plot_cat(ii, 3), :))
    end
    for ii = 1:num_core
        plot(x_core(ii), y_core(ii), 'k^', 'markersize', 10, 'markerfacecolor', 'm')
    end
    axis equal
    axis([-632 846 -3344 -670])
    caxis([1 10])
    set(gca, 'fontsize', 20, 'xtick', [], 'ytick', [])
    colorbar('fontsize', 20, 'ylim', [4 10], 'ytick', 4:10, 'yticklabel', {'0' '5' '10' '15' '20' '25' '>30'})
    text(850, -545, '# reflections traced', 'color', 'k', 'fontsize', 20, 'horizontalalignment', 'center')
    fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
    fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
    text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
    text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
    text(-595, -775, '(b)', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
    text(-542, -3238, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
    text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
    text(552, -2743, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
    text(759, -2975, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
    grid off
    box on

%% Figure 3a/b/c (required independent loading through phaseinterp call)
    
    figure('position', [100 100 1000 1000])
    subplot('position', [0.08 0.69 0.88 0.27])
    hold on
    imagesc(dist_all, (1e6 .* twtt_old), (10 .* log10(abs(data))), [-90 -40])
    plot([0.75 1 1 0.75 0.75], [15 15 25 25 15], 'w--', 'linewidth', 2)
    axis([dist_all([1 end])' 0 40])
    colormap(bone)
    set(gca, 'fontsize', 20)
    ylabel('Traveltime ({\mu}s)')
    colorbar('fontsize', 20)
    text(double(1.01 * dist_all(end)), -3, '(dB)', 'fontsize', 20)
    text(double(0.02 * dist_all(end)), 3, '(a)', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'edgecolor', 'k', 'backgroundcolor', 'w')
    text(1.25, 5, 'surface', 'fontsize', 20, 'color', 'w', 'fontweight', 'bold')
    text(1.25, 32, 'bed', 'fontsize', 20, 'color', 'w', 'fontweight', 'bold')
    box on
    axis ij
    freezeColors
    cbfreeze
    axes('position', [0.6 0.73 0.25 0.2])
    imagesc(dist_all, (1e6 .* twtt_old), (10 .* log10(abs(data))), [-90 -40])
    axis([0.75 1 15 25])
    colormap(bone)
    set(gca, 'fontsize', 20, 'xcolor', 'w', 'ycolor', 'w', 'xtick', [0.75 1])
    box on
    axis ij
    freezeColors
    subplot('position', [0.08 0.38 0.88 0.27])
    hold on
    imagesc(dist_all, (1e6 .* twtt_old), angle(data), [-pi pi])
    plot([0.75 1 1 0.75 0.75], [15 15 25 25 15], 'w--', 'linewidth', 2)
    axis([dist_all([1 end])' 0 40])
    colormap(redblue)
    set(gca, 'fontsize', 20)
    ylabel('Traveltime ({\mu}s)')
    colorbar('fontsize', 20)
    text(double(1.01 * dist_all(end)), -3, '(rad)', 'fontsize', 20)
    text(double(0.02 * dist_all(end)), 3, '(b)', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'edgecolor', 'k', 'backgroundcolor', 'w')
    box on
    axis ij
    axes('position', [0.6 0.42 0.25 0.2])
    imagesc(dist_all, (1e6 .* twtt_old), angle(data), [-pi pi])
    axis([0.75 1 15 25])
    colormap(redblue)
    set(gca, 'fontsize', 20, 'xcolor', 'w', 'ycolor', 'w', 'xtick', [0.75 1])
    box on
    axis ij
    subplot('position', [0.08 0.055 0.88 0.27])
    hold on
    imagesc(dist_all, (1e6 .* twtt_old), phase_diff_filt, [-0.3 0.3])
    plot([0.75 1 1 0.75 0.75], [15 15 25 25 15], 'w--', 'linewidth', 2)
    axis([dist_all([1 end])' 0 40])
    colormap(redblue)
    set(gca, 'fontsize', 20)
    xlabel('Distance (km)')
    ylabel('Traveltime ({\mu}s)')
    cb                      = colorbar('fontsize', 20);
    set(cb, 'ylim', [-0.3 0.3], 'ytick', -0.3:0.1:0.3)
    text(double(0.96 * dist_all(end)), -3, '(rad trace^{-1})', 'fontsize', 20)
    text(double(0.02 * dist_all(end)), 3, '(c)', 'fontsize', 20, 'fontweight', 'bold', 'color', 'k', 'edgecolor', 'k', 'backgroundcolor', 'w')
    box on
    axis ij
    axes('position', [0.6 0.095 0.25 0.2])
    imagesc(dist_all, (1e6 .* twtt_old), phase_diff_filt, [-0.3 0.3])
    axis([0.75 1 15 25])
    colormap(redblue)
    set(gca, 'fontsize', 20, 'xcolor', 'k', 'ycolor', 'k', 'xtick', [0.75 1])
    box on
    axis ij
    
%% Figure 4a (Data_20110506_01_011_wf_2_adc_13_block_12.mat, jj=num_filt)

    figure('position', [100 100 800 800])
    imagesc(wavenum_rad, (1e6 .* twtt_old), (10 .* log10(abs(data_fft))), [-70 -30])
    colormap(bone)
    hold on
    axis([wavenum_rad([1 end]) 0 40])
    vl                      = vline(0, 'w--');
    set(vl, 'linewidth', 2)
    psr                     = vline(wavenum_rad(ind_filt(1)), 'r--');
    set(psr, 'linewidth', 2)
    vl                      = vline(wavenum_rad(ind_filt(end)), 'r--');
    set(vl, 'linewidth', 2)
    pr                      = plot(wavenum_rad(tmp1(~isnan(tmp1)) + ind_filt(1) - 1), (1e6 .* twtt_old(~isnan(tmp1))), 'b.', 'markersize', 16, 'visible', 'off');
    ps                      = plot(phase_diff_filt_curr(:, jj), (1e6 .* twtt_old), 'b', 'linewidth', 2, 'visible', 'off');
    set(gca, 'fontsize', 20)
    xlabel('Wavenumber (rad)')
    ylabel('Traveltime ({\mu}s)')
    text((1.02 * wavenum_rad(end)), -2, '(dB)', 'fontsize', 20, 'color', 'k')
    text(-1.6, 2, 'surface', 'fontsize', 20, 'color', 'w', 'fontweight', 'bold')
    text(-1.9, 15, {'sloped'; 'internal'; 'reflections'}, 'fontsize', 20, 'color', 'w', 'fontweight', 'bold')
    text(-1.3, 36, 'bed', 'fontsize', 20, 'color', 'w', 'fontweight', 'bold')
    text(-3, 1.5, '(a)', 'fontsize', 20, 'color', 'k', 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
    colorbar('fontsize', 20)
    legend([psr pr], 'k_D search limits', 'k_D maxima', 'location', 'southeast')
    box on
    axes('position', [0.6 0.4 0.23 0.5])
    imagesc(wavenum_rad, (1e6 .* twtt_old), (10 .* log10(abs(data_fft))), [-70 -30])
    colormap(bone)
    hold on
    axis([wavenum_rad(ind_filt([1 end])) 5 25])
    vl                      = vline(0, 'w--');
    set(vl, 'linewidth', 2)
    pr                      = plot(wavenum_rad(tmp1(~isnan(tmp1)) + ind_filt(1) - 1), (1e6 .* twtt_old(~isnan(tmp1))), 'b.', 'markersize', 16);
%     ps                      = plot(phase_diff_filt_curr(:, jj), (1e6 .* twtt_old), 'b', 'linewidth', 2);
    set(gca, 'fontsize', 20, 'xcolor', 'w', 'ycolor', 'w')
    box on    
    
%% Figure 8: Age at selected normalized depths

    figure('position', [200 200 1600 1200], 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)]);
    depth_norms             = depth_norm([5 10 15 20]);
    letters                 = ['a':'d'; 'e':'h'];
    ranges                  = zeros(4, 2, 2);
    ranges(:, :, 1)         = [0 5; 0 15; 0 40; 10 90];
    ranges(:, :, 2)         = [0 1; 0 1; 0 2; 1 5];
    incs                    = [0.25 0.75 2 4; 0.05 0.05 0.1 0.2];
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    [cb, ax]                = deal(zeros(2, 4));
    for ii = 1:4
        for kk = 1:2
            ax(kk, ii)      = subplot('position', [(0.005 + (0.25 * (ii - 1))) (0.515 - ((kk - 1) * 0.49)) 0.245 0.44]);
            hold on
            axis equal
            axis([-632 846 -3344 -670])
            range_tmp       = squeeze(ranges(ii, 1, kk)):incs(kk, ii):squeeze(ranges(ii, 2, kk));
            switch kk
                case 1
                    plot_tmp= 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, (1e-3 .* squeeze(age_norm2(:, :, (depth_norm == depth_norms(ii))))), 'nearest', 'extrap');
                case 2
                    plot_tmp= 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, (1e-3 .* squeeze(age_uncert_norm2(:, :, (depth_norm == depth_norms(ii))))), 'nearest', 'extrap');
            end
            plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo(isnan(plot_tmp)) + 1;
            imagesc(x_grd(1, :), y_grd(:, 1), plot_tmp)
            shading interp
            for jj = 1:length(x_paral)
                plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
            end
            for jj = 1:length(x_merid)
                plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
            end
            for jj = 1:num_coast
                plot(x_coast{jj}, y_coast{jj}, 'color', 'k', 'linewidth', 1)
            end
            for jj = 1:num_basin
                plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
            end
            for jj = 1:num_core
                plot(x_core(jj), y_core(jj), 'k^', 'markersize', 10, 'markerfacecolor', 'm')
            end
            set(gca, 'fontsize', 20, 'xtick', [], 'ytick', [])
            caxis([1 24])
            fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
            fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
            text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
            text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
%             text(-595, -775, ['(' letters(kk, ii) ') ' num2str(1e2 * depth_norms(ii)) '%'], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
            text(-608, -3208, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
            text(-339, -3275, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
            text(506, -2763, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
            text(713, -2848, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
            switch kk
                case 1
                    title(['(' letters(kk, ii) ') $a(\hat{z} = ' num2str(1e2 * depth_norms(ii)) '\%)$'], 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold')
                case 2
                    title(['(' letters(kk, ii) ') $\tilde{a}(\hat{z} = ' num2str(1e2 * depth_norms(ii)) '\%)$'], 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold')
            end
            grid off
            box on
            tick_str        = cell(1, 21);
            for jj = 1:21
                tick_str{jj}= '';
            end
            for jj = 1:4:21
                tick_str{jj}= num2str(squeeze(ranges(ii, 1, kk)) + (incs(kk, ii) * (jj - 1)));
            end
            tick_str{end}   = ['>' tick_str{end}];
            if (squeeze(ranges(ii, 1, kk)) > 0)
                tick_str{1} = ['<' tick_str{1}];
            end
            cb1(kk, ii)     = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', squeeze(ranges(ii, :, kk)), 'ytick', squeeze(ranges(ii, 1, kk)):incs(kk, ii):squeeze(ranges(ii, 2, kk)), 'yticklabel', tick_str);
            set(cb1(kk, ii), 'visible', 'off')
            cb2             = cbarf([4 24], 4:24, 'vertical', 'linear');
            set(cb2, 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str, 'position', get(cb1(kk, ii), 'position'))
            text(975, -525, '(ka)', 'color', 'k', 'fontsize', 20, 'horizontalalignment', 'center')
        end
    end
    linkaxes(ax)

%% Figure 9: Depth at selected ages

    figure('position', [200 200 1600 1200], 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)]);
    age_isos                = 1e3 .* [11.7 29 57 115];
    letters                 = ['a':'d'; 'e':'h'];
    ranges                  = zeros(4, 2, 2);
    ranges(:, :, 1)         = [500 2500; 500 2500; 1000 3000; 1000 3000];
    ranges(:, :, 2)         = [0 100; 0 100; 0 100; 0 200];
    incs                    = [100 100 100 100; 5 5 5 10];
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    [cb, ax]                = deal(zeros(2, 4));
    for ii = 1:4
        for kk = 1:2
            ax(kk, ii)      = subplot('position', [(0.005 + (0.25 * (ii - 1))) (0.515 - ((kk - 1) * 0.49)) 0.245 0.44]);
            hold on
            axis equal
            axis([-632 846 -3344 -670])
            range_tmp       = squeeze(ranges(ii, 1, kk)):incs(kk, ii):squeeze(ranges(ii, 2, kk));
            switch kk
                case 1
                    plot_tmp= 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, squeeze(depth_iso2(:, :, (age_iso == age_isos(ii)))), 'nearest', 'extrap');
                case 2
                    plot_tmp= 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, squeeze(depth_uncert_iso2(:, :, (age_iso == age_isos(ii)))), 'nearest', 'extrap');
            end
            plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo(isnan(plot_tmp)) + 1;
            imagesc(x_grd(1, :), y_grd(:, 1), plot_tmp)
            shading interp
            for jj = 1:length(x_paral)
                plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
            end
            for jj = 1:length(x_merid)
                plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
            end
            for jj = 1:num_coast
                plot(x_coast{jj}, y_coast{jj}, 'color', 'k', 'linewidth', 1)
            end
            for jj = 1:num_basin
                plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
            end
            for jj = 1:num_core
                plot(x_core(jj), y_core(jj), 'k^', 'markersize', 10, 'markerfacecolor', 'm')
            end
            set(gca, 'fontsize', 20, 'xtick', [], 'ytick', [])
            caxis([1 24])
            fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
            fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
            text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
            text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
%             text(-595, -775, ['(' letters(kk, ii) ') ' num2str(1e2 * depth_norms(ii)) '%'], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
            text(-608, -3208, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
            text(-339, -3275, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
            text(506, -2763, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
            text(713, -2848, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
            switch kk
                case 1
                    title(['(' letters(kk, ii) ') $z(a = ' num2str(1e-3 * age_isos(ii)) '\,\mathrm{ka})$'], 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold')
                case 2
                    title(['(' letters(kk, ii) ') $\tilde{z}(a = ' num2str(1e-3 * age_isos(ii)) '\,\mathrm{ka})$'], 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold')
            end
            grid off
            box on
            tick_str        = cell(1, 21);
            for jj = 1:21
                tick_str{jj}= '';
            end
            for jj = 1:4:21
                tick_str{jj}= num2str(squeeze(ranges(ii, 1, kk)) + (incs(kk, ii) * (jj - 1)));
            end
            tick_str{end}   = ['>' tick_str{end}];
            if (squeeze(ranges(ii, 1, kk)) > 0)
                tick_str{1} = ['<' tick_str{1}];
            end
            cb1(kk, ii)     = colorbar('fontsize', 20, 'location', 'eastoutside', 'ylim', squeeze(ranges(ii, :, kk)), 'ytick', squeeze(ranges(ii, 1, kk)):incs(kk, ii):squeeze(ranges(ii, 2, kk)), 'yticklabel', tick_str);
            set(cb1(kk, ii), 'visible', 'off')
            cb2             = cbarf([4 24], 4:24, 'vertical', 'linear');
            set(cb2, 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str, 'position', get(cb1(kk, ii), 'position'))
            text(975, -525, '(m)', 'color', 'k', 'fontsize', 20, 'horizontalalignment', 'center')
        end
    end
    linkaxes(ax)
        
%% Figure 10: Eemian thickness

    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    figure('position', [200 200 1200 800], 'renderer', 'zbuffer')
    colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)])
    ranges                  = [0 600; 0 300];
    incs                    = [30 15];
    titles                  = {'$z(a = 115\,\mathrm{ka}) - z_{bed}$' '$\sqrt{\tilde{z}(a = 115\,\mathrm{ka})^2 + \tilde{z}_{bed}^2 + \tilde{z}_{surface}^2}$'};
    letters                 = 'a':'b';
    for ii = 1:2
        subplot('position', [(0.005 + ((ii - 1) * 0.47)) 0.03 0.45 0.90]);
        hold on
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        switch ii
            case 1
                plot_tmp    = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, (thick_grd - squeeze(depth_iso2(:, :, (age_iso == 115e3)))), 'nearest', 'extrap');
            case 2
                plot_tmp    = 3 + interp1((range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2)), 1:20, thick_eemian_uncert, 'nearest', 'extrap');
        end
        plot_tmp(isnan(plot_tmp)) ...
                            = mask_combo(isnan(plot_tmp)) + 1;
        imagesc(x_grd(1, :), y_grd(:, 1), plot_tmp)
        if (ii == 1)
            for jj = 1:num_year %#ok<*FXUP>
                for kk = 1:num_trans(jj)
                    for ll = 1:length(ind_fence{jj}{kk})
                        if isempty(age{jj}{kk}{ll})
                            continue
                        end
                        if isempty(find((age{jj}{kk}{ll} >= 115e3), 1))
                            continue
                        end
                        for mm = find(age{jj}{kk}{ll} >= 115e3)'
                            plot(x_pk{jj}{kk}{ll}(~isnan(depth_smooth{jj}{kk}{ll}(mm, :))), y_pk{jj}{kk}{ll}(~isnan(depth_smooth{jj}{kk}{ll}(mm, :))), '.', 'color', [0.5 0.5 0.5], 'markersize', 10, 'tag', num2str([jj kk ll]));
                        end
                    end
                end
            end
        end
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        for jj = 1:num_core
            plot(x_core(jj), y_core(jj), 'k^', 'markersize', 10, 'markerfacecolor', 'm')
        end
        caxis([1 23])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xtick', [], 'ytick', [])
        text(900, -550, '(m)', 'color', 'k', 'fontsize', 20)
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-542, -3238, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(552, -2743, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(759, -2975, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        title(['(' letters(ii) ') ' titles{ii}], 'interpreter', 'latex', 'fontsize', 20, 'fontweight', 'bold')
        tick_str                = cell(1, 21);
        for jj = 2:2:20
            tick_str{jj}        = '';
        end
        for jj = 1:2:21
            tick_str{jj}        = num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
        end
        tick_str{end}           = ['>' tick_str{end}];
        cb1                     = colorbar('fontsize', 20, 'location', 'eastoutside');
        set(cb1, 'visible', 'off')
        cb2                     = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2, 'position', get(cb1, 'position'), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str)
    end
    
%% Figure 11: Bed slopes a/b/c
    
    figure('position', [150 150 1600 1200])
    colors                  = colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(21); gray(19)]);
    plots                   = {'slope_bed' 'slope_bed_predict' 'res_slope_bed'};
    titles                  = {'$\alpha(H)$' '$\alpha''(H)$' '$\left|\alpha(H)-\alpha''(H)\right|$'};
    letters                 = 'a':'f';
    letters_bonus           = {'observed' 'predicted' 'abs. difference' 'Petermann Gl.' ['Jakobshavn Isbr' char(230)] 'NEGIS'};
    range_tmp1              = 0:5:100;
    x_lim                   = [-300 0; -150 150; 100 400];
    y_lim                   = [-1500 -1000; -2400 -2000; -1500 -1100];
    for ii = 1:3
        subplot('position', [(0 + ((ii - 1) * 0.28)) 0.49 0.28 0.46])
        hold on
        axis equal
        axis([-632 846 -3344 -670])
        pcolor(x_grd(1, :), y_grd(:, 1), (mask_combo + 1))
        shading interp
        [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
        for jj = 1:length(x_paral)
            plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:length(x_merid)
            plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1)
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', 'k', 'linewidth', 1)
        end
        [x_cat, y_cat, plot_cat] ...
                            = deal([]);
        for jj = fliplr(year_ord) %#ok<*FXUP>
            for kk = 1:num_trans(jj)
                for ll = 1:length(ind_fence{jj}{kk})
                    if isempty(slope_bed{jj}{kk}{ll})
                        continue
                    end
                    ind_decim_mid ...
                            = round(ind_decim5{jj}{kk}{ll}(1:(end - 1)) + (diff(ind_decim5{jj}{kk}{ll}) ./ 2));
%                     plot(x_pk{jj}{kk}{ll}(ind_decim_mid), y_pk{jj}{kk}{ll}(ind_decim_mid), 'k', 'linewidth', 1);
                    [x_cat, y_cat, plot_cat] ...
                            = deal([x_cat x_pk{jj}{kk}{ll}(ind_decim_mid)], [y_cat y_pk{jj}{kk}{ll}(ind_decim_mid)], [plot_cat eval([plots{ii} '{jj}{kk}{ll}'])]);
                end
            end
        end
        plot_cat            = 3 + interp1((range_tmp1(1:(end - 1)) + (diff(range_tmp1) ./ 2)), 1:20, abs(plot_cat), 'nearest', 'extrap');
        xy_plot_cat         = sortrows([x_cat' y_cat' plot_cat'], 3);
        xy_plot_cat         = xy_plot_cat(~isnan(xy_plot_cat(:, 3)), :);
        for jj = 1:size(xy_plot_cat, 1)
            line(xy_plot_cat(jj, 1), xy_plot_cat(jj, 2), 'marker', '.', 'markersize', 12, 'color', colors(xy_plot_cat(jj, 3), :))
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        set(gca, 'fontsize', 20, 'xtick', [], 'ytick', [])
        title(titles{ii}, 'interpreter', 'latex', 'fontweight', 'bold')
        caxis([1 44])
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-595, -775, ['(' letters(ii) ') ' letters_bonus{ii}], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        text(-590, -3217, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(513, -2772, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(720, -2893, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        grid off
        box on
        if (ii == 3)
            text(1050, -510, '(m km^{-1})', 'color', 'k', 'fontsize', 20, 'horizontalalignment', 'center')
            cb              = colorbar('fontsize', 20, 'position', [0.80 0.49 0.015 0.46]);
            set(cb, 'visible', 'off')
            cb2             = cbarf([4 24], 4:24, 'vertical', 'linear');
            set(cb2, 'fontsize', 20, 'ytick', 4:24, 'yticklabel', {'0' '' '10' '' '20' '' '30' '' '40' '' '50' '' '60' '' '70' '' '80' '' '90' '' '100'}, 'position', get(cb, 'position'))
            for jj = 1:3
                plot(x_lim(jj, [1 2 2 1 1]),  y_lim(jj, [1 1 2 2 1]), 'm', 'linewidth', 2)
            end
            text(-450, -1100, '(d)', 'color', 'm', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
            text(-250, -2500, '(e)', 'color', 'm', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
            text(425, -1200, '(f)', 'color', 'm', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        end
    end
%% Figure 11 d/e/f
    [pf, p0, p50]           = deal(zeros(1, 3));
    for ii = 1:3
        subplot('position', [(0 + ((ii - 1) * 0.28)) 0.01 0.28 0.46])
        hold on
        axis equal
        axis([x_lim(ii, :) y_lim(ii, :)])
        range_tmp2          = log10([1:10 20:10:100]);
        plot_tmp            = mask_combo + 1;
        plot_tmp((plot_tmp == 3) & ~isnan(speed)) ...
                            = 24 + interp1((range_tmp2(1:(end - 1)) + (diff(range_tmp2) ./ 2)), 1:18, log10(speed((plot_tmp == 3) & ~isnan(speed))), 'nearest', 'extrap');
        pcolor(x_grd(1, :), y_grd(:, 1), plot_tmp)
        shading interp
%         for jj = 1:num_coast
%             plot(x_coast{jj}, y_coast{jj}, 'color', 'k', 'linewidth', 1)
%         end
        for jj = fliplr(year_ord) %#ok<*FXUP>
            for kk = 1:(num_trans(jj))
                for ll = 1:length(ind_fence{jj}{kk})
                    if isempty(slope_bed{jj}{kk}{ll})
                        continue
                    end
                    ind_decim_mid ...
                            = round(ind_decim5{jj}{kk}{ll}(1:(end - 1)) + (diff(ind_decim5{jj}{kk}{ll}) ./ 2));
%                     plot(x_pk{jj}{kk}{ll}(ind_decim_mid), y_pk{jj}{kk}{ll}(ind_decim_mid), 'k', 'linewidth', 1);
                    [x_cat, y_cat, plot_cat] ...
                            = deal([x_cat x_pk{jj}{kk}{ll}(ind_decim_mid)], [y_cat y_pk{jj}{kk}{ll}(ind_decim_mid)], [plot_cat res_slope_bed{jj}{kk}{ll}]);
                end
            end
        end
        plot_cat            = 3 + interp1((range_tmp1(1:(end - 1)) + (diff(range_tmp1) ./ 2)), 1:20, abs(plot_cat), 'nearest', 'extrap');
        xy_plot_cat         = sortrows([x_cat' y_cat' plot_cat'], 3);
        xy_plot_cat         = xy_plot_cat(~isnan(xy_plot_cat(:, 3)), :);
        for jj = 1:size(xy_plot_cat, 1)
            line(xy_plot_cat(jj, 1), xy_plot_cat(jj, 2), 'marker', '.', 'markersize', 20, 'color', colors(xy_plot_cat(jj, 3), :))
        end
        for jj = 1:num_basin
            plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 2)
        end
        set(gca, 'fontsize', 20, 'xtick', [], 'ytick', [])
        caxis([1 44])
        text((x_lim(ii, 1) + (0.03 * diff(x_lim(ii, :)))), (y_lim(ii, 2) - (0.05 * diff(y_lim(ii, :)))), ['(' letters(ii + 3) ') ' letters_bonus{ii + 3}], 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        if (ii == 1)
            pf(ii)          = fill([(x_lim(ii, 1) + 10) (x_lim(ii, 1) + 60) (x_lim(ii, 1) + 60) (x_lim(ii, 1) + 10)], [(y_lim(ii, 1) + 10) (y_lim(ii, 1) + 10) (y_lim(ii, 1) + 15) (y_lim(ii, 1) + 15)], 'k');
            p0(ii)          = text((x_lim(ii, 1) + 8), (y_lim(ii, 1) + 22), '0', 'color', 'k', 'fontsize', 18, 'fontweight', 'bold');
            p50(ii)         = text((x_lim(ii, 1) + 50), (y_lim(ii, 1) + 22), '50 km', 'color', 'k', 'fontsize', 18, 'fontweight', 'bold');
        else
            pf(ii)          = fill([(x_lim(ii, 2) - 100) (x_lim(ii, 2) - 50) (x_lim(ii, 2) - 50) (x_lim(ii, 2) - 100)], [(y_lim(ii, 1) + 10) (y_lim(ii, 1) + 10) (y_lim(ii, 1) + 15) (y_lim(ii, 1) + 15)], 'k');
            p0(ii)          = text((x_lim(ii, 2) - 105), (y_lim(ii, 1) + 22), '0', 'color', 'k', 'fontsize', 18, 'fontweight', 'bold');
            p50(ii)         = text((x_lim(ii, 2) - 60), (y_lim(ii, 1) + 22), '50 km', 'color', 'k', 'fontsize', 18, 'fontweight', 'bold');
        end
        grid off
        box on
        if (ii == 3)
            text(445, -1081, '(m a^{-1})', 'color', 'k', 'fontsize', 20, 'horizontalalignment', 'center')
            cb1             = colorbar('fontsize', 20, 'position', [0.85 0.01 0.015 0.46], 'ylim', [26 46], 'ytick', 26:5:46, 'yticklabel', {'1' '5' '10' '50' '100'});
            set(cb1, 'visible', 'off')
            cb2             = cbarf([26 46], 26:46, 'vertical', 'linear');
            set(cb2, 'fontsize', 20, 'position', [0.83 0.015 0.015 0.445], 'ylim', [26 46], 'ytick', 26:5:46, 'yticklabel', {'1' '5' '10' '50' '100'})
        end
        uistack([pf(ii) p0(ii) p50(ii)], 'top')
    end
%%
end