% LINE_PLOTS
% 
% Joe MacGregor
% Last updated: 08/07/14

clear

plotting                    = false;
age_max                     = 9e3;
[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);

load mat/xy_all num_year name_trans
load mat/merge_all dist ind_decim ind_fence num_fence x_pk y_pk
load mat/core_int x_core y_core num_core
load mat/age_grd2 x_grd y_grd
load(['mat/strain_all_' num2str(1e-3 * age_max) 'ka'], 'melt_bed', 'strain_rate', 'thick_shear')
load mat/grl_coast num_coast x_coast y_coast
load mat/grl_basin num_basin x_basin y_basin
load mat/xy_D1 x_D1 y_D1
load mat/mask_D1 mask_D1 x_D1_ord y_D1_ord

load mat/greenland_bed_v3 mask_gris mask_land x y
[mask_land, mask_gris]      = deal(mask_land(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))), ...
                                   mask_gris(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))));
clear x y
mask_combo                  = double(mask_land);
mask_combo(mask_gris)       = 2;

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
%% line plot
    ranges                  = [0 40; 0 2000; -10 10];
    incs                    = [2 100 1];
    plots                   = {'1e5 .* strain_rate' 'thick_shear' '1e2 .* melt_bed'};
    units                   = {'10^{-5} a^{-1}' 'm' 'cm a^{-1}'};
    [cb1, cb2]              = deal(zeros(1, 3));
    [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland([-632 846], [-3344 -670], 5, 10);
    for ii = 1:3
        figure('position', [200 200 600 800])
        if (ii < 3)
            colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; jet(20)])
        else
            colormap([([135 206 235] ./ 255); ([159 89 39] ./ 255); 0.9 0.9 0.9; redblue(20)])
        end
        subplot('position', [0.005 0.03 0.9 0.90]);
        hold on
        imagesc(x_grd(1, :), y_grd(:, 1), (mask_combo + 1))
        range_tmp           = ranges(ii, 1):incs(ii):ranges(ii, 2);
        range_tmp           = range_tmp(1:(end - 1)) + (diff(range_tmp) ./ 2);
        for jj = 1:num_year %#ok<*FXUP>
            for kk = find(num_fence{jj})
                for ll = find(ind_fence{jj}{kk})
                    if isempty(eval([plots{ii} '{jj}{kk}{ll}']))
                        continue
                    end
                    if (length(eval([plots{ii} '{jj}{kk}{ll}'])) == (length(ind_decim{jj}{kk}{ll}) - 1))
                        ind_decim_mid ...
                            = round(ind_decim{jj}{kk}{ll}(1:(end - 1)) + (diff(ind_decim{jj}{kk}{ll}) ./ 2));
                    else
                        ind_decim_mid ...
                            = round(ind_decim5{jj}{kk}{ll}(1:(end - 1)) + (diff(ind_decim5{jj}{kk}{ll}) ./ 2));
                    end
                    plot(x_pk{jj}{kk}{ll}(ind_decim_mid), y_pk{jj}{kk}{ll}(ind_decim_mid), 'k', 'linewidth', 1, 'tag', num2str([jj kk ll]));
                    ind_tmp = find(~isnan(eval([plots{ii} '{jj}{kk}{ll}'])));
                    scatter(x_pk{jj}{kk}{ll}(ind_decim_mid(ind_tmp)), y_pk{jj}{kk}{ll}(ind_decim_mid(ind_tmp)), 36, (3 + interp1(range_tmp, 1:20, eval([plots{ii} '{jj}{kk}{ll}(ind_tmp)']), 'nearest', 'extrap')), 'filled');
                end
            end
        end
        [pp, pm, pb]        = deal(zeros(1, length(x_paral)), zeros(1, length(x_merid)), zeros(1, length(num_basin)));
        for jj = 1:length(x_paral)
            pp(jj)          = plot((1e-3 .* x_paral{jj}), (1e-3 .* y_paral{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:length(x_merid)
            pm(jj)          = plot((1e-3 .* x_merid{jj}), (1e-3 .* y_merid{jj}), 'k--', 'linewidth', 1);
        end
        for jj = 1:num_coast
            plot(x_coast{jj}, y_coast{jj}, 'color', [0.5 0.5 0.5], 'linewidth', 1)
        end
        for jj = 1:num_basin
            pb(jj)          = plot(x_basin{jj}, y_basin{jj}, 'color', 'k', 'linewidth', 3);
        end
        plot(x_D1_ord, y_D1_ord, 'm', 'linewidth', 2)
        caxis([1 23])
        axis equal
        axis([x_min x_max y_min y_max])
        box on
        set(gca, 'fontsize', 20, 'layer', 'top', 'xticklabel', {}, 'yticklabel', {})
        text(900, -550, ['(' units{ii} ')'], 'color', 'k', 'fontsize', 20)
        fill([325 450 450 325], [-3230 -3230 -3260 -3260], 'k')
        fill([450 575 575 450], [-3230 -3230 -3260 -3260], 'w', 'edgecolor', 'k')
        text(300, -3180, '0', 'color', 'k', 'fontsize', 18)
        text(520, -3180, '250 km', 'color', 'k', 'fontsize', 18)
        text(-542, -3238, '60\circN', 'color', 'k', 'fontsize', 18, 'rotation', -10)
        text(-334, -3250, '50\circW', 'color', 'k', 'fontsize', 18, 'rotation', 82)
        text(552, -2743, '65\circN', 'color', 'k', 'fontsize', 18, 'rotation', 15)
        text(759, -2975, '30\circW', 'color', 'k', 'fontsize', 18, 'rotation', -75)
        tick_str            = cell(1, 21);
        for jj = 2:2:20
            tick_str{jj}= '';
        end
        for jj = 1:2:21
            tick_str{jj}= num2str(ranges(ii, 1) + (incs(ii) * (jj - 1)));
        end
        tick_str{1}         = ['<' tick_str{1}];
        tick_str{end}       = ['>' tick_str{end}];
        cb1(ii)              = colorbar('fontsize', 20, 'location', 'eastoutside');
        set(cb1(ii), 'visible', 'off')
        cb2(ii)             = cbarf([4 24], 4:24, 'vertical', 'linear');
        set(cb2(ii), 'position', get(cb1(ii), 'position'), 'fontsize', 20, 'ytick', 4:24, 'yticklabel', tick_str)
    end
%%
end