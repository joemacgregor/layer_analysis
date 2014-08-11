% GRIS_CAKE Layer cake animations.
%
% Joe MacGregor (UTIG)
% Last updated: 08/07/14

clear all

plotting                    = false;

load mat/grl_coast num_coast x_coast y_coast
load mat/greenland_bed_v3 mask_greenland mask_gris x y
load mat/greenland_mc_bed_1km elev_bed elev_surf thick
load mat/age_grd2_krige depth_norm num_depth_norm x_grd y_grd
load mat/age_grd2_krige_alt age_norm2_alt
age_norm2                   = age_norm2_alt;
clear age_norm2_alt
load mat/core_int name_core num_core x_core y_core
mog                         = geotiffread('../../data/greenland/mog100_hp1_v10_1km');

int_year                    = 1e3 .* [linspace(0, 11.7, 5) linspace(11.7, 115, 5) linspace(115, 130, 3)];
int_year                    = int_year([1:5 7:10 12:13]);
int_year_mid                = int_year(1:(end - 1)) + (diff(int_year) ./ 2);
num_int_year                = length(int_year);

int_diff                    = -100:20:100;
int_diff_mid                = int_diff(1:(end - 1)) + (diff(int_diff) ./ 2);
num_int_diff_mid            = length(int_diff_mid);

[elev_bed_grd, elev_surf_grd, thick_grd] ...
                            = deal(elev_bed, elev_surf, thick);
clear elev_bed elev_surf thick

[x_min, x_max, y_min, y_max]= deal(-632, 846, -3344, -670);
[mask_greenland, mask_gris] = deal(mask_greenland(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))), mask_gris(((y(:, 1) >= y_min) & (y(:, 1) <= y_max)), ((x(1, :) >= x_min) & (x(1, :) <= x_max))));
clear x y

elev_bed_grd(~mask_greenland) ...
                            = NaN;
elev_surf_grd(~mask_greenland | ~mask_gris) ...
                            = NaN;
thick_grd(~mask_greenland | ~mask_gris) ...
                            = NaN;
clear mask_greenland mask_gris

[num_y, num_x]              = size(x_grd);

load mat/greenland_cism accum x y
[accum, x_cism, y_cism]     = deal(double(accum), x, y);
clear x y
[accum, x_cism, y_cism]     = deal(accum(((y_cism(:, 1) >= y_min) & (y_cism(:, 1) <= y_max)), ((x_cism(1, :) >= x_min) & (x_cism(1, :) <= x_max))), x_cism(((y_cism(:, 1) >= y_min) & (y_cism(:, 1) <= y_max)), ((x_cism(1, :) >= x_min) & (x_cism(1, :) <= x_max))), ...
                                   y_cism(((y_cism(:, 1) >= y_min) & (y_cism(:, 1) <= y_max)), ((x_cism(1, :) >= x_min) & (x_cism(1, :) <= x_max))));

age_nye                     = -(thick_grd(:, :, ones(1, 1, num_depth_norm)) ./ accum(:, :, ones(1, 1, num_depth_norm))) .* log(1 - permute(depth_norm(:, ones(1, size(x_grd, 1)), ones(1, size(x_grd, 2))), [2 3 1]));
age_diff                    = 1e2 .* ((age_norm2 - age_nye) ./ age_nye);

[x_core_label, y_core_label]= deal([-751 -14 -5 264 -191 -48], [-1233 -2619 -1794 -1794 -1252 -1521]);

% x/y cake ranges to display
x_cake                      = x_grd(1, find(sum(sum(~isnan(age_norm2), 3)), 1):find(sum(sum(~isnan(age_norm2), 3)), 1, 'last'));
y_cake                      = y_grd(find(sum(sum(~isnan(age_norm2), 3), 2), 1):find(sum(sum(~isnan(age_norm2), 3), 2), 1, 'last'), 1);
num_cake_x                  = length(x_cake);
num_cake_y                  = length(y_cake);

% clean up MOG
[x_mog, y_mog]              = deal(-1199:900, -3399:-600);
mog                         = flipud(double(mog(((y_mog >= y_min) & (y_mog <= y_max)), ((x_mog >= x_min) & (x_mog <= x_max)))));

% ranges for MOG and bed display
[mog_min, mog_max]          = deal(15500, 16500);
[elev_bed_min, elev_bed_max]= deal(-500, 2000);

% on/off indices
[ind_core_x_on, ind_core_y_on] ...
                            = deal(NaN(1, num_core));
for ii = 1:num_core
    ind_core_x_on(ii)       = find((x_cake > x_core(ii)), 1);
    ind_core_y_on(ii)       = find((y_cake > y_core(ii)), 1);
end
core_delay                  = 100; % km
[ind_core_x_off, ind_core_y_off] ...
                            = deal((ind_core_x_on + core_delay), (ind_core_y_on + core_delay));

% indexify bed colors
ind_elev_bed                = elev_bed_grd;
ind_elev_bed                = ind_elev_bed - elev_bed_min;
ind_elev_bed                = 1 + round(1023 .* (ind_elev_bed ./ (elev_bed_max - elev_bed_min)));
ind_elev_bed                = 1025 - ind_elev_bed;
ind_elev_bed(ind_elev_bed < 1) ...
                            = 1;
ind_elev_bed(ind_elev_bed > 1024) ...
                            = 1024;

% indexify MOG
ind_mog                     = mog;
ind_mog(ind_mog < mog_min)  = mog_min;
ind_mog(ind_mog > mog_max)  = mog_max;
ind_mog                     = ind_mog - mog_min;
ind_mog                     = 1 + round(1023 .* (ind_mog ./ (mog_max - mog_min)));
ind_mog                     = 1024 + ind_mog;

% colormap setup
age_plots                   = {'age_norm2' 'age_nye' 'age_diff'};
copper_tmp                  = copper(round(1024 * (4 / 3)));
redblue_tmp                 = redblue(11);
color_year                  = [1 0 0; 1 0.5 0.5; 0 0 0.9; 0.25 0.25 1; 0.5 0.5 1; 0.75 0.75 1; 0 0.9 0; 0.25 1 0.25; 0.5 1 0.5; 0.75 1 0.75];
colors_single               = [flipud(copper_tmp((end - 1023):end, :)); bone(1024); color_year];
colors_all                  = [flipud(copper_tmp((end - 1023):end, :)); bone(1024); color_year; redblue_tmp(1:5, :); redblue_tmp(7:11, :)];

% pre-assign a bunch of age/position slices for faster rendering
age_slice_plots_x           = {'age_slice_x' 'age_slice_x_nye' 'age_slice_x_diff'};
age_slice_plots_y           = {'age_slice_y' 'age_slice_y_nye' 'age_slice_y_diff'};
[age_slice_x, age_slice_x_nye, age_slice_x_diff, depth_rep_x] ...
                            = deal(NaN((num_depth_norm + 1), num_y, num_cake_x));
[age_slice_y, age_slice_y_nye, age_slice_y_diff, depth_rep_y] ...
                            = deal(NaN((num_depth_norm + 1), num_x, num_cake_y));

depth_norm_rep_y            = [zeros(1, num_y); depth_norm(:, ones(1, num_y))];
depth_norm_rep_x            = [zeros(1, num_x); depth_norm(:, ones(1, num_x))];

for ii = 1:3
    for jj = 1:num_cake_x
        if ~mod(jj, 50)
            disp([num2str(ii) ':' num2str(jj) 'x'])
        end
        ind_cake_curr       = find(x_grd(1, :) == x_cake(jj));
        if (ii < 3)
            age_slice       = 2048 + [(10 .* ones(1, num_y)); interp1(int_year_mid, fliplr(1:(num_int_year - 1)), squeeze(eval([age_plots{ii} '(:, ind_cake_curr, :)']))', 'nearest', 'extrap')];
        else
            age_slice       = 2058 + [NaN(1, num_y); interp1(int_diff_mid, 1:num_int_diff_mid, squeeze(eval([age_plots{ii} '(:, ind_cake_curr, :)']))', 'nearest', 'extrap')];
        end
        for kk = find(~isnan(elev_surf_grd(:, ind_cake_curr)))'
            if (all(isnan(age_slice(2:end, kk))) || isnan(age_slice(2, kk)))
                age_slice(1, kk) ...
                            = 2048;
            end
            age_slice(isnan(age_slice(:, kk)), kk) ...
                            = 2048;
        end
        age_slice(1, ~sum(~isnan(age_slice(2:end, :)))) ...
                            = NaN;
        eval([age_slice_plots_x{ii} '(:, :, jj) = age_slice;']);
        if (ii == 1)
            depth_rep_x(:, :, jj) ...
                            = (elev_surf_grd(:, ind_cake_curr(ones((num_depth_norm + 1), 1)))' - (depth_norm_rep_y .* thick_grd(:, ind_cake_curr(ones((num_depth_norm + 1), 1)))'));
        end
    end
    for jj = 1:num_cake_y
        if ~mod(jj, 50)
            disp([num2str(ii) ':' num2str(jj) 'y'])
        end
        ind_cake_curr       = find(y_grd(:, 1) == y_cake(jj));
        if (ii < 3)
            age_slice       = 2048 + [(10 .* ones(1, num_x)); interp1(int_year_mid, fliplr(1:(num_int_year - 1)), squeeze(eval([age_plots{ii} '(ind_cake_curr, :, :)']))', 'nearest', 'extrap')];
        else
            age_slice       = 2058 + [NaN(1, num_x); interp1(int_diff_mid, 1:num_int_diff_mid, squeeze(eval([age_plots{ii} '(ind_cake_curr, :, :)']))', 'nearest', 'extrap')];
        end
        for kk = find(~isnan(elev_surf_grd(ind_cake_curr, :)))
            if (all(isnan(age_slice(2:end, kk))) || isnan(age_slice(2, kk)))
                age_slice(1, kk) ...
                            = 2048;
            end
            age_slice(isnan(age_slice(:, kk)), kk) ...
                            = 2048;
        end
        age_slice(1, ~sum(~isnan(age_slice(2:end, :)))) ...
                            = NaN;
        eval([age_slice_plots_y{ii} '(:, :, jj) = age_slice;']);
        if (ii == 1)
            depth_rep_y(:, :, jj) ...
                            = (elev_surf_grd(ind_cake_curr(ones((num_depth_norm + 1), 1)), :) - (depth_norm_rep_x .* thick_grd(ind_cake_curr(ones((num_depth_norm + 1), 1)), :)));
        end
    end
end

% pre-calculate display position matrices for surf
% x-relevant
x_cake_set                  = NaN((num_depth_norm + 1), num_y, num_cake_x);
for ii = 1:num_cake_x
    x_cake_set(:, :, ii)    = x_cake(ii) .* ones((num_depth_norm + 1), num_y);
end
y_grd_rep                   = y_grd(:, ones((num_depth_norm + 1), 1))';

% y-relevant
x_grd_rep                   = x_grd(ones((num_depth_norm + 1), 1), :);
y_cake_set                  = NaN((num_depth_norm + 1), num_x, num_cake_x);
for ii = 1:num_cake_y
    y_cake_set(:, :, ii)    = y_cake(ii) .* ones((num_depth_norm + 1), num_x);
end

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
    
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
    
%% Layer cake age X (observed only)

    fxo                     = figure('position', [200 200 1920 1080], 'color', 'w', 'renderer', 'zbuffer');
    axes('position', [0.02 0.10 0.9 0.85])
    hold on
    colormap(colors_single)
    surf(x_grd, y_grd, elev_bed_grd, 'cdatamapping', 'direct', 'cdata', ind_elev_bed);
    shading flat
    pes                     = surf(x_grd, y_grd, elev_surf_grd, 'cdatamapping', 'direct', 'cdata', ind_mog);
    shading flat
    annotation('textbox', 'position', [0.86 0.92 0.11 0.04], 'string', 'Ice age (ka)', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.87 0.11 0.04], 'string', '0', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.55 0.11 0.04], 'string', '11.7', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.23 0.11 0.04], 'string', '115', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.07 0.11 0.04], 'string', '>130', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.71 0.11 0.04], 'string', 'Holocene', 'color', [0 0.9 0], 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.40 0.11 0.04], 'string', {'Last Glacial'; 'Period'}, 'color', [0 0 0.9], 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.15 0.11 0.04], 'string', 'Eemian', 'color', 'r', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    [pcl, pct]              = deal(zeros(1, num_core), cell(1, num_core));
    for ii = 1:num_core
        pcl(ii)             = plot3(repmat(x_core(ii), 1, 2), repmat(y_core(ii), 1, 2), [interp2(x_grd, y_grd, elev_bed_grd, x_core(ii), y_core(ii)) interp2(x_grd, y_grd, elev_surf_grd, x_core(ii), y_core(ii))], 'k', 'linewidth', 3, 'visible', 'off');
        if (ii == 3)
            pct{ii}         = text(repmat(x_core(ii), 1, 2), repmat((y_core(ii) - 10), 1, 2), repmat(1500, 1, 2), name_core{ii}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 22, 'fontweight', 'bold', 'visible', 'off');
        else
            pct{ii}         = text(repmat(x_core(ii), 1, 2), repmat((y_core(ii) - 10), 1, 2), repmat(2000, 1, 2), name_core{ii}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 22, 'fontweight', 'bold', 'visible', 'off');
        end
    end
    for ii = 1:num_coast
        plot3(x_coast{ii}, y_coast{ii}, zeros(1, length(x_coast{ii})), 'color', 'k', 'linewidth', 1.5)
    end
    axis tight
    view([-90 30])
    set(gca, 'fontsize', 20, 'dataaspectratio', [1 1 10], 'visible', 'off')
    caxis([1 size(colors_single, 1)])
    cb                      = cbarf([2049 (size(colors_single, 1) + 1)], 2049:(size(colors_single, 1) + 1), 'vertical', 'linear');
    set(cb, 'position', [0.88 0.1 0.015 0.8], 'yticklabel', {})
    annotation('line', [0.894791666666667 0.8796875], [0.9 0.9]);
    grid off
    box off

%% animation X (observed only)

    tic
    ind_mog_curr            = ind_mog;
    rect                    = get(fxo, 'position');
    mp4                     = VideoWriter('fig/gris_cake_x.mp4', 'MPEG-4');
    mp4.FrameRate           = 30;
    mp4.Quality             = 100;
    open(mp4)
    writeVideo(mp4, getframe(fxo, [0 0 rect(3:4)]));
    for ii = 1:10%num_cake_x
        if (ii > 1)
            delete(pc)
            ind_mog_curr(:, find((x_grd(1, :) == x_cake(ii)), 1)) ...
                            = NaN;
        else
            ind_mog_curr(:, 1:find((x_grd(1, :) == x_cake(ii)), 1)) ...
                            = NaN;
        end
        set(pes, 'cdata', ind_mog_curr)
        pc                  = surf(squeeze(x_cake_set(:, :, ii)), y_grd_rep, squeeze(depth_rep_x(:, :, ii)), squeeze(age_slice_x(:, :, ii)), 'facecolor', 'flat', 'edgecolor', 'none', 'facelighting', 'none', 'cdatamapping', 'direct');
        for jj = 1:num_core
            if (ii == ind_core_x_on(jj))
                set(pcl(jj), 'visible', 'on')
                set(pct{jj}, 'visible', 'on')
            end
            if (ii == ind_core_x_off(jj))
                set(pcl(jj), 'visible', 'off')
                set(pct{jj}, 'visible', 'off')
            end
        end
        writeVideo(mp4, getframe(fxo, [0 0 rect(3:4)]));
    end
    close(mp4)
    toc
    
%% Layer cake age X (all 3)

    [ax, pes]               = deal(zeros(1, 3));
    [pcl, pct]              = deal(zeros(3, num_core), cell(3, num_core));
    fxa                     = figure('position', [200 200 1920 1080], 'color', 'w', 'renderer', 'zbuffer');
    colormap(colors_all)
    for ii = 1:3
        ax(ii)              = axes('position', [0.02 (0.60 - ((ii - 1) * 0.32)) 0.85 0.44]);
        hold on
        surf(x_grd, y_grd, elev_bed_grd, 'cdatamapping', 'direct', 'cdata', ind_elev_bed);
        shading flat
        pes(ii)             = surf(x_grd, y_grd, elev_surf_grd, 'cdatamapping', 'direct', 'cdata', ind_mog);
        shading flat
        if (ii == 1)
            annotation('textbox', 'position', [0.85 0.94 0.11 0.04], 'string', 'Ice age (ka)', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.90 0.11 0.04], 'string', '0', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.68 0.11 0.04], 'string', '11.7', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.46 0.11 0.04], 'string', '115', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.35 0.11 0.04], 'string', '>130', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.79 0.11 0.04], 'string', {'Holocene'; ''}, 'color', [0 0.9 0], 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.57 0.11 0.04], 'string', {'Last Glacial'; 'Period'}, 'color', [0 0 0.9], 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.40 0.11 0.04], 'string', 'Eemian', 'color', 'r', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.10 0.94 0.11 0.04], 'string', 'Nye model', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.10 0.62 0.11 0.04], 'string', 'Observed', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.10 0.28 0.15 0.04], 'string', 'Relative difference', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.83 0.31 0.15 0.04], 'string', 'Relative difference (%)', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.27 0.11 0.04], 'string', '>100', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.14 0.11 0.04], 'string', '0', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.90 0.02 0.11 0.04], 'string', '<-100', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
        end
        for jj = 1:num_coast
            plot3(x_coast{jj}, y_coast{jj}, zeros(1, length(x_coast{jj})), 'color', 'k', 'linewidth', 1.5)
        end
        for jj = 1:num_core
            pcl(ii, jj)     = plot3(repmat(x_core(jj), 1, 2), repmat(y_core(jj), 1, 2), [interp2(x_grd, y_grd, elev_bed_grd, x_core(jj), y_core(jj)) interp2(x_grd, y_grd, elev_surf_grd, x_core(jj), y_core(jj))], 'k', 'linewidth', 3, 'visible', 'off');
            if (jj == 3)
                pct{ii, jj} = text(repmat(x_core(jj), 1, 2), repmat((y_core(jj) - 10), 1, 2), repmat(1500, 1, 2), name_core{jj}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 18, 'fontweight', 'bold', 'visible', 'off');
            else
                pct{ii, jj} = text(repmat(x_core(jj), 1, 2), repmat((y_core(jj) - 10), 1, 2), repmat(2000, 1, 2), name_core{jj}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 18, 'fontweight', 'bold', 'visible', 'off');
            end
        end
        axis tight
        view([-90 30])
        set(gca, 'fontsize', 20, 'dataaspectratio', [1 1 10], 'visible', 'off')
        caxis([0 size(colors_all, 1)])
        switch ii
            case 1
                colorbar('ylim', [2048 2058], 'ytick', 2048:2058, 'yticklabel', {}, 'fontsize', 20, 'position', [0.88 0.38 0.015 0.55])
            case 3
                colorbar('ylim', [2058 size(colors_all, 1)], 'ytick', 2058:size(colors_all, 1), 'yticklabel', {}, 'fontsize', 20, 'position', [0.88 0.05 0.015 0.25])
        end
        grid off
        box off
    end
    
%% animation X (all)

    tic
    ind_mog_curr            = ind_mog;
    rect                    = get(fxa, 'position');
    mp4                     = VideoWriter('fig/gris_cake_x_all.mp4', 'MPEG-4');
    mp4.FrameRate           = 30;
    mp4.Quality             = 100;
    open(mp4)
    writeVideo(mp4, getframe(fxa, [0 0 rect(3:4)]));
    for ii = 1:num_cake_x
        if (ii > 1)
            delete(pc)
            ind_mog_curr(:, find((x_grd(1, :) == x_cake(ii)), 1)) ...
                            = NaN;            
        else
            ind_mog_curr(:, 1:find((x_grd(1, :) == x_cake(ii)), 1)) ...
                            = NaN;
        end
        set(pes, 'cdata', ind_mog_curr)
        pc                  = zeros(1, 3);
        for jj = 1:3
            axes(ax(jj))
            pc(jj)          = surf(squeeze(x_cake_set(:, :, ii)), y_grd_rep, squeeze(depth_rep_x(:, :, ii)), squeeze(eval([age_slice_plots_x{jj} '(:, :, ii)'])), 'facecolor', 'flat', 'edgecolor', 'none', 'facelighting', 'none', 'cdatamapping', 'direct');
        end
        for jj = 1:num_core
            if (ii == ind_core_x_on(jj))
                set(pcl(:, jj), 'visible', 'on')
                set(pct{:, jj}, 'visible', 'on')
            end
            if (ii == ind_core_x_off(jj))
                set(pcl(:, jj), 'visible', 'off')
                set(pct{:, jj}, 'visible', 'off')
            end
        end
        writeVideo(mp4, getframe(fxa, [0 0 rect(3:4)]));
    end
    close(mp4)
    
%% Layer cake age Y (observed only)

    fyo                     = figure('position', [200 200 1920 1080], 'color', 'w', 'renderer', 'zbuffer');
    axes('position', [0.02 -0.05 0.85 1.1])
    hold on
    colormap(colors_single)
    surf(x_grd, y_grd, elev_bed_grd, 'cdatamapping', 'direct', 'cdata', ind_elev_bed);
    shading flat
    pes                     = surf(x_grd, y_grd, elev_surf_grd, 'cdatamapping', 'direct', 'cdata', ind_mog);
    shading flat
    annotation('textbox', 'position', [0.86 0.92 0.11 0.04], 'string', 'Ice age (ka)', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.87 0.11 0.04], 'string', '0', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.55 0.11 0.04], 'string', '11.7', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.23 0.11 0.04], 'string', '115', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.07 0.11 0.04], 'string', '>130', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.71 0.11 0.04], 'string', 'Holocene', 'color', [0 0.9 0], 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.40 0.11 0.04], 'string', {'Last Glacial'; 'Period'}, 'color', [0 0 0.9], 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    annotation('textbox', 'position', [0.90 0.15 0.11 0.04], 'string', 'Eemian', 'color', 'r', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none')
    pve                     = annotation('textbox', 'position', [0.70 0.07 0.2 0.04], 'string', '100:1 vertical exaggeration', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold', 'linestyle', 'none');
    for ii = 1:num_coast
        plot3(x_coast{ii}, y_coast{ii}, zeros(1, length(x_coast{ii})), 'color', 'k', 'linewidth', 1.5)
    end
    [pcl, pct]              = deal(zeros(1, num_core), cell(1, num_core));
    for ii = 1:num_core
        pcl(ii)             = plot3(repmat(x_core(ii), 1, 2), repmat(y_core(ii), 1, 2), [interp2(x_grd, y_grd, elev_bed_grd, x_core(ii), y_core(ii)) interp2(x_grd, y_grd, elev_surf_grd, x_core(ii), y_core(ii))], 'k', 'linewidth', 3, 'visible', 'off');
        if (ii == 3)
            pct{ii}         = text(repmat((x_core(ii) + 10), 1, 2), repmat(y_core(ii), 1, 2), repmat(1500, 1, 2), name_core{ii}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 22, 'fontweight', 'bold', 'visible', 'off');
        else
            pct{ii}         = text(repmat((x_core(ii) + 10), 1, 2), repmat(y_core(ii), 1, 2), repmat(2000, 1, 2), name_core{ii}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 22, 'fontweight', 'bold', 'visible', 'off');
        end
    end
    axis tight
    view([0 30])
    set(gca, 'fontsize', 22, 'dataaspectratio', [1 1 10], 'visible', 'off')
    caxis([1 size(colors_single, 1)])
    cb                      = cbarf([2049 (size(colors_single, 1) + 1)], 2049:(size(colors_single, 1) + 1), 'vertical', 'linear');
    set(cb, 'position', [0.88 0.1 0.015 0.8], 'yticklabel', {})
    annotation('line', [0.894791666666667 0.8796875], [0.9 0.9]);
    grid off
    box off
    
%% animation Y (observed only)
    tic
    figure(fyo)
    ind_mog_curr            = ind_mog;
    rect                    = get(fyo, 'position');
    mp4                     = VideoWriter('fig/gris_cake_y.mp4', 'MPEG-4');
    mp4.FrameRate           = 30;
    mp4.Quality             = 100;
    open(mp4)
    writeVideo(mp4, getframe(fyo, [0 0 rect(3:4)]));
    for ii = 1:num_cake_y
        if (ii > 1)
            delete(pc)
            ind_mog_curr(find((y_grd(:, 1) == y_cake(ii)), 1), :) ...
                            = NaN;
        else
            delete(pve)
            ind_mog_curr(1:find((y_grd(:, 1) == y_cake(ii)), 1), :) ...
                            = NaN;
        end
        set(pes, 'cdata', ind_mog_curr)
        pc                  = surf(x_grd_rep, squeeze(y_cake_set(:, :, ii)), squeeze(depth_rep_y(:, :, ii)), age_slice_y(:, :, ii), 'facecolor', 'flat', 'edgecolor', 'none', 'facelighting', 'none', 'cdatamapping', 'direct');
        for jj = 1:num_core
            if (ii == ind_core_y_on(jj))
                set(pcl(jj), 'visible', 'on')
                set(pct{jj}, 'visible', 'on')
            end
            if (ii == ind_core_y_off(jj))
                delete(pcl(jj))
                delete(pct{jj})
            end
        end
        writeVideo(mp4, getframe(fyo, [0 0 rect(3:4)]));
    end
    close(mp4)
    toc
        
%% Layer cake age Y (all 3)

    [ax, pes]               = deal(zeros(1, 3));
    [pcl, pct]              = deal(zeros(3, num_core), cell(3, num_core));
    fya                     = figure('position', [200 200 1920 1080], 'color', 'w', 'renderer', 'zbuffer');
    colormap(colors_all)
    for ii = 1:3
        ax(ii)              = axes('position', [(0.01 + ((ii - 1) * 0.32)) -0.05 0.34 1.1]);
        hold on
        surf(x_grd, y_grd, elev_bed_grd, 'cdatamapping', 'direct', 'cdata', ind_elev_bed);
        shading flat
        pes(ii)             = surf(x_grd, y_grd, elev_surf_grd, 'cdatamapping', 'direct', 'cdata', ind_mog);
        shading flat
        if (ii == 1)
            annotation('textbox', 'position', [0.32 0.12 0.11 0.04], 'string', 'Ice age (ka)', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.575 0.05 0.11 0.04], 'string', '0', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.39 0.05 0.11 0.04], 'string', '11.7', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.21 0.05 0.11 0.04], 'string', '115', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.11 0.05 0.11 0.04], 'string', '>130', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.47 0.05 0.11 0.04], 'string', 'Holocene', 'color', [0 0.9 0], 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.26 0.05 0.11 0.04], 'string', 'Last Glacial Period', 'color', [0 0 0.9], 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.15 0.05 0.11 0.04], 'string', 'Eemian', 'color', 'r', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.13 0.83 0.11 0.04], 'string', '(a) Nye model', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.46 0.83 0.11 0.04], 'string', '(b) Observed', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.75 0.83 0.18 0.04], 'string', '(c) Observed - Nye model', 'color', 'k', 'fontsize', 24, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.78 0.12 0.15 0.04], 'string', 'Relative difference (%)', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.93 0.05 0.11 0.04], 'string', '>100', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.835 0.05 0.11 0.04], 'string', '0', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            annotation('textbox', 'position', [0.72 0.05 0.11 0.04], 'string', '<-100', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none')
            pve             = annotation('textbox', 'position', [0.55 0.15 0.2 0.04], 'string', '100:1 vertical exaggeration', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold', 'linestyle', 'none');
        end
        for jj = 1:num_coast
            plot3(x_coast{jj}, y_coast{jj}, zeros(1, length(x_coast{jj})), 'color', 'k', 'linewidth', 1.5)
        end
        for jj = 1:num_core
            pcl(ii, jj)     = plot3(repmat(x_core(jj), 1, 2), repmat(y_core(jj), 1, 2), [interp2(x_grd, y_grd, elev_bed_grd, x_core(jj), y_core(jj)) interp2(x_grd, y_grd, elev_surf_grd, x_core(jj), y_core(jj))], 'k', 'linewidth', 3, 'visible', 'off');
            if (jj == 3)
                pct{ii, jj} = text(repmat((x_core(jj) + 10), 1, 2), repmat(y_core(jj), 1, 2), repmat(1500, 1, 2), name_core{jj}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 18, 'fontweight', 'bold', 'visible', 'off');
            else
                pct{ii, jj} = text(repmat((x_core(jj) + 10), 1, 2), repmat(y_core(jj), 1, 2), repmat(2000, 1, 2), name_core{jj}, 'color', 'k', 'linewidth', 3, 'edgecolor', 'k', 'backgroundcolor', 'w', 'fontsize', 18, 'fontweight', 'bold', 'visible', 'off');
            end
        end
        axis tight
        view([0 30])
        set(gca, 'fontsize', 20, 'dataaspectratio', [1 1 10], 'visible', 'off')
        caxis([1 size(colors_all, 1)])
        switch ii
            case 1
                cb1         = cbarf([2049 2059], 2049:2059, 'horiz', 'linear');
                set(cb1, 'position', [0.13 0.10 0.45 0.02], 'xticklabel', {}, 'fontsize', 20)
%                 colorbar('location', 'southoutside', 'xlim', [2048 2058], 'xtick', 2048:2058, 'xticklabel', {}, 'fontsize', 20, 'position', [0.13 0.10 0.45 0.02])
            case 2
                cbf         = cbarf([2049 2059], 2049:2059, 'horiz', 'linear');
                set(cbf, 'position', [0.5 0.51 0.05 0.02], 'xticklabel', {}, 'fontsize', 20, 'visible', 'off')
            case 3
                cb2         = cbarf([2059 (size(colors_all, 1) + 1)], 2059:(size(colors_all, 1) + 1), 'horiz', 'linear');
                set(cb2, 'position', [0.74 0.10 0.2 0.02], 'xticklabel', {}, 'fontsize', 20)
%                 colorbar('location', 'southoutside', 'xlim', [2058 size(colors_all, 1)], 'xtick', 2058:size(colors_all, 1), 'xticklabel', {}, 'fontsize', 20, 'position', [0.74 0.10 0.2 0.02])
        end
        annotation('line', [0.58 0.58], [0.10 0.12]);
        grid off
        box off
    end
    
%% animation Y (all)

    tic
    figure(fya)
    ind_mog_curr            = ind_mog;
    rect                    = get(fya, 'position');
    mp4                     = VideoWriter('fig/gris_cake_y_all.mp4', 'MPEG-4');
    mp4.FrameRate           = 30;
    mp4.Quality             = 100;
    open(mp4)
    writeVideo(mp4, getframe(fya, [0 0 rect(3:4)]));
    for ii = 1:num_cake_y
        if (ii > 1)
            delete(pc)
            ind_mog_curr(find((y_grd(:, 1) == y_cake(ii)), 1), :) ...
                            = NaN;
        else
            delete(pve)
            ind_mog_curr(1:find((y_grd(:, 1) == y_cake(ii)), 1), :) ...
                            = NaN;
        end
        set(pes, 'cdata', ind_mog_curr)
        pc                  = zeros(1, 3);
        for jj = 1:3
            axes(ax(jj))
            pc(jj)          = surf(x_grd_rep, squeeze(y_cake_set(:, :, ii)), squeeze(depth_rep_y(:, :, ii)), squeeze(eval([age_slice_plots_y{jj} '(:, :, ii)'])), 'facecolor', 'flat', 'edgecolor', 'none', 'facelighting', 'none', 'cdatamapping', 'direct');
        end
        for jj = 1:num_core
            if (ii == ind_core_y_on(jj))
                set(pcl(:, jj), 'visible', 'on')
                set(pct{1, jj}, 'visible', 'on')
                set(pct{2, jj}, 'visible', 'on')
                set(pct{3, jj}, 'visible', 'on')
            end
            if (ii == ind_core_y_off(jj))
                delete(pcl(:, jj))
                delete(pct{1, jj})
                delete(pct{2, jj})
                delete(pct{3, jj})
%                 set(pcl(:, jj), 'visible', 'off')
%                 set(pct{1, jj}, 'visible', 'off')
%                 set(pct{2, jj}, 'visible', 'off')
%                 set(pct{3, jj}, 'visible', 'off')
            end
        end
        writeVideo(mp4, getframe(fya, [0 0 rect(3:4)]));
    end
    close(mp4)
    toc
    
%%
end