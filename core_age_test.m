% CORE_AGE_TEST Evaluate quasi-Nye method using Greenland ice-core depth-age scales.
%
% Joe MacGregor (UTIG)
% Last updated: 09/30/14

clear all

plotting                    = false;

load mat/core_int name_core_short num_core

% sort out cores
core                        = cell(1, num_core);
for ii = 1:num_core
    core{ii}                = load(['mat/depth_age/' name_core_short{ii} '_depth_age']);
end

colors                      = {'r' 'g' 'b' 'c' 'm' 'k'};
colors_uncert               = [1    0.75 0.75;
                               0.75 1    0.75;
                               0.75 0.75 1;
                               0.75 1    1;
                               1    0.75 1;
                               0.75 0.75 0.75];

% setup parameters
num_test                    = 1e3; % number of tests to perform
num_layer                   = 10; % number of core-intersecting layers
depth_rel_min               = 0.05; % minimum relative depth of core-intersecting layers
depth_rel_max               = 0.85; % maximum relative depth of core-intersecting layers
depth_rel_all               = (0:0.01:1)'; % full relative depth vector
num_depth_rel_all           = length(depth_rel_all);
tol                         = 0.01; % quasi-Nye tolerance, m
iter_max                    = 100; % maximum number of iterations for quasi-Nye dating
depth_rel_shallow_min       = 0.03; % minimum shallow layer relative thickness to date
depth_rel_shallow_max       = 0.2; % maximum shallow layer relative thickness to date
layer_diff_max              = 1; % maximum difference in layer thickness
thick_diff_max              = 0.2; % fraction of ice thickness within which overlapping layer must lie
depth_uncert                = 10; % depth uncertainty, m

% interpolate core depth-age scales to full relative depth vector
age_core_all                = NaN((num_depth_rel_all - 1), num_core);
for ii = 1:num_core
    age_core_all(:, ii)     = interp1((core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), core{ii}.age(~isnan(core{ii}.age)), depth_rel_all(2:end));
end

% initialize quasi age variables
[age_diff, age_mean, age_min, age_qdj, age_qdj_diff_rel, age_qnye, age_qnye_diff_rel, depth_diff] ...
                            = deal(NaN(num_test, num_core, length(depth_rel_all)));
age_layer                   = NaN(num_test, num_core, (num_layer + 1));
age_layer(:, :, 1)          = zeros(num_test, num_core, 1);

% uniformly distributed random layer depths within prescribed depth interval
depth_rel_layer             = sort(cat(3, zeros(num_test, num_core, 1), (depth_rel_min + ((depth_rel_max - depth_rel_min) * rand(num_test, num_core, num_layer)))), 3);

% layer ages
for ii = 1:num_core
    age_layer(:, ii, 2:end) = interp1((core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), core{ii}.age(~isnan(core{ii}.age)), depth_rel_layer(:, ii, 2:end)); % age at synthetic "layers"    
end


%%
% loop through each test
for ii = 1:num_test
    
    if ~mod(ii, 50)
        disp(ii)
    end
       
    % loop through each core
    for jj = 1:num_core
        
        % adjust depths by randomized depth uncertainty
        depth_rel_layer(ii, jj, 2:end) ...
                            = depth_rel_layer(ii, jj, 2:end) + ((depth_uncert * randn(1, 1, num_layer)) ./ core{jj}.depth(end));
        
        % loop through each layer pair (interval)
        for kk = 2:num_layer
            
            if any(isnan(age_layer(ii, jj, kk:(kk + 1))))
                continue
            end
            
            ind_all_curr    = find((depth_rel_all >= depth_rel_layer(ii, jj, kk)) & (depth_rel_all <= depth_rel_layer(ii, jj, (kk + 1))));
            depth_rel_layer_curr ...
                            = squeeze(depth_rel_layer(ii, jj, kk:(kk + 1)));
            depth_layer_curr= depth_rel_layer_curr .* core{jj}.depth(end);
            age_layer_curr  = squeeze(age_layer(ii, jj, kk:(kk + 1)));
            strain_rate_curr= -diff(log(1 - depth_rel_layer_curr)) / diff(age_layer_curr);
            
            % loop through all normalized depths between current layer pair
            for ll = ind_all_curr';
                try
                    age_qnye(ii, jj, ll) ...
                            = date_quasi_nye(depth_layer_curr, age_layer_curr, strain_rate_curr, (depth_rel_all(ll) * core{jj}.depth(end)), tol, iter_max, 'age');
%                     age_qdj(ii, jj, ll) ...
%                             = date_quasi_dj(depth_layer_curr, age_layer_curr, strain_rate_curr, (depth_rel_all(ll) * core{jj}.depth(end)), tol, iter_max, 'age');
                catch ME
                    disp([num2str([ii jj kk]) ' failed...'])
                    disp(ME.identifier)
                    disp(ME.message)
                    disp(ME.stack(:).line)
                end
            end
            
            % age and depth differences of core-dated layers used to date the normalized depth
            age_diff(ii, jj, ind_all_curr) ...
                            = diff(age_layer_curr);
            age_mean(ii, jj, ind_all_curr) ...
                            = mean(age_layer(ii, jj, kk:(kk + 1)));
            age_min(ii, jj, ind_all_curr) ...
                            = age_layer(ii, jj, kk);
            depth_diff(ii, jj, ind_all_curr) ...
                            = diff(depth_layer_curr);
        end
        
        % date shallow layers within restrictions
        ind_shallow         = find((depth_rel_all < depth_rel_layer(ii, jj, 2)) & (depth_rel_all >= depth_rel_shallow_min) & (depth_rel_all <= depth_rel_shallow_max));
        age_qnye(ii, jj, ind_shallow) ...
                            = age_layer(ii, jj, 2) * (log(1 - depth_rel_all(ind_shallow)) / log(1 - depth_rel_layer(ii, jj, 2)));
        age_diff(ii, jj, ind_shallow) ...
                            = age_layer(ii, jj, 2);
        age_mean(ii, jj, ind_shallow) ...
                            = 0.5 * age_layer(ii, jj, 2);
        age_min(ii, jj, ind_shallow) ...
                            = 0;
        depth_diff(ii, jj, ind_shallow) ...
                            = diff(depth_rel_layer(ii, jj, 1:2)) * core{jj}.depth(end);
        
        % date deep layers within restrictions
        ind_layer_last      = find(~isnan(squeeze(age_layer(ii, jj, :))), 1, 'last');
        depth_rel_deep_curr = squeeze(depth_rel_layer(ii, jj, (ind_layer_last - 1):ind_layer_last));
        ind_deep            = find((depth_rel_all > depth_rel_layer(ii, jj, ind_layer_last)) & ((depth_rel_all - depth_rel_layer(ii, jj, ind_layer_last)) <= (layer_diff_max * diff(depth_rel_deep_curr))) & ((depth_rel_all - depth_rel_layer(ii, jj, ind_layer_last)) <= thick_diff_max));
        age_deep_curr       = squeeze(age_layer(ii, jj, (ind_layer_last - 1):ind_layer_last));
        strain_rate_curr    = -diff(log(1 - depth_rel_deep_curr)) / diff(age_deep_curr);
        depth_deep_curr     = depth_rel_deep_curr * core{jj}.depth(end);
        for kk = ind_deep'
            age_qnye(ii, jj, kk) ...
                            = date_quasi_nye(depth_deep_curr, age_deep_curr, strain_rate_curr, (depth_rel_all(kk) * core{jj}.depth(end)), tol, iter_max, 'age');
%             age_qdj(ii, jj, kk) ...
%                             = date_quasi_dj(depth_deep_curr, age_deep_curr, strain_rate_curr, (depth_rel_all(kk) * core{jj}.depth(end)), tol, iter_max, 'age');
        end
        age_diff(ii, jj, ind_deep) ...
                            = diff(age_deep_curr);
        age_mean(ii, jj, ind_deep) ...
                            = mean(age_deep_curr);
        age_min(ii, jj, ind_deep) ...
                            = age_deep_curr(1);
        depth_diff(ii, jj, ind_deep) ...
                            = diff(depth_deep_curr);
        
        % relative age difference from core
        age_qnye_diff_rel(ii, jj, 2:end) ...
                            = (age_qnye(ii, jj, 2:end) - shiftdim(age_core_all(:, jj), -2)) ./ shiftdim(age_core_all(:, jj), -2);
        age_qdj_diff_rel(ii, jj, 2:end) ...
                            = (age_qdj(ii, jj, 2:end) - shiftdim(age_core_all(:, jj), -2)) ./ shiftdim(age_core_all(:, jj), -2);
    end
end
%%
%%

% vertical age gradient and relative difference of upper bounding age
age_depth_grad              = age_diff ./ depth_diff;
age_diff_rel                = age_diff ./ age_mean;
age_qnye_var                = squeeze(nanvar(age_qnye));
age_qnye_diff_rel_mean      = squeeze(nanmean(age_qnye_diff_rel));
age_qdj_var                 = squeeze(nanvar(age_qdj));
age_qdj_diff_rel_mean       = squeeze(nanmean(age_qdj_diff_rel));

% correlation coefficients
r_depth_diff{1}             = corrcoef(depth_diff(~isnan(age_qnye_diff_rel) & ~isnan(depth_diff)), abs(age_qnye_diff_rel(~isnan(age_qnye_diff_rel) & ~isnan(depth_diff))));
r_age_diff{1}               = corrcoef(age_diff(~isnan(age_qnye_diff_rel) & ~isnan(age_diff)), abs(age_qnye_diff_rel(~isnan(age_qnye_diff_rel) & ~isnan(age_diff))));
r_age_depth_grad{1}         = corrcoef(age_depth_grad(~isnan(age_qnye_diff_rel) & ~isnan(age_depth_grad)), abs(age_qnye_diff_rel(~isnan(age_qnye_diff_rel) & ~isnan(age_depth_grad))));
r_age_diff_rel{1}           = corrcoef(age_diff_rel(~isnan(age_qnye_diff_rel) & ~isnan(age_diff_rel)), abs(age_qnye_diff_rel(~isnan(age_qnye_diff_rel) & ~isnan(age_diff_rel))));
r_age_mean{1}               = corrcoef(age_mean(~isnan(age_qnye_diff_rel) & ~isnan(age_mean)), abs(age_qnye_diff_rel(~isnan(age_qnye_diff_rel) & ~isnan(age_mean))));
r_var_age{1}                = corrcoef(age_qnye_var(~isnan(age_qnye_diff_rel_mean) & ~isnan(age_qnye_var)), abs(age_qnye_diff_rel_mean(~isnan(age_qnye_diff_rel_mean) & ~isnan(age_qnye_var))));
% r_depth_diff{2}             = corrcoef(depth_diff(~isnan(age_qdj_diff_rel) & ~isnan(depth_diff)), abs(age_qdj_diff_rel(~isnan(age_qdj_diff_rel) & ~isnan(depth_diff))));
% r_age_diff{2}               = corrcoef(age_diff(~isnan(age_qdj_diff_rel) & ~isnan(age_diff)), abs(age_qdj_diff_rel(~isnan(age_qdj_diff_rel) & ~isnan(age_diff))));
% r_age_depth_grad{2}         = corrcoef(age_depth_grad(~isnan(age_qdj_diff_rel) & ~isnan(age_depth_grad)), abs(age_qdj_diff_rel(~isnan(age_qdj_diff_rel) & ~isnan(age_depth_grad))));
% r_age_diff_rel{2}           = corrcoef(age_diff_rel(~isnan(age_qdj_diff_rel) & ~isnan(age_diff_rel)), abs(age_qdj_diff_rel(~isnan(age_qdj_diff_rel) & ~isnan(age_diff_rel))));
% r_age_mean{2}               = corrcoef(age_mean(~isnan(age_qdj_diff_rel) & ~isnan(age_mean)), abs(age_qdj_diff_rel(~isnan(age_qdj_diff_rel) & ~isnan(age_mean))));
% r_var_age{2}                = corrcoef(age_qdj_var(~isnan(age_qdj_diff_rel_mean) & ~isnan(age_qdj_var)), abs(age_qdj_diff_rel_mean(~isnan(age_qdj_diff_rel_mean) & ~isnan(age_qdj_var))));

% samples for each correlation
N_full                      = num_test * num_core * num_depth_rel_all;
N_mean                      = num_core * num_depth_rel_all;

% standard deviation of Z-statistic
sigma_full                  = 1 / sqrt(N_full - 3);
sigma_mean                  = 1 / sqrt(N_mean - 3);

% Z statistics of correlation coefficients
Z_depth_diff(1)             = 0.5 * log((1 + r_depth_diff{1}(1, 2)) / (1 - r_depth_diff{1}(1, 2)));
Z_age_diff(1)               = 0.5 * log((1 + r_age_diff{1}(1, 2)) / (1 - r_age_diff{1}(1, 2)));
Z_age_depth_grad(1)         = 0.5 * log((1 + r_age_depth_grad{1}(1, 2)) / (1 - r_age_depth_grad{1}(1, 2)));
Z_age_diff_rel(1)           = 0.5 * log((1 + r_age_diff_rel{1}(1, 2)) / (1 - r_age_diff_rel{1}(1, 2)));
Z_age_mean(1)               = 0.5 * log((1 + r_age_mean{1}(1, 2)) / (1 - r_age_mean{1}(1, 2)));
Z_var_age(1)                = 0.5 * log((1 + r_var_age{1}(1, 2)) / (1 - r_var_age{1}(1, 2)));
% Z_depth_diff(2)             = 0.5 * log((1 + r_depth_diff{2}(1, 2)) / (1 - r_depth_diff{2}(1, 2)));
% Z_age_diff(2)               = 0.5 * log((1 + r_age_diff{2}(1, 2)) / (1 - r_age_diff{2}(1, 2)));
% Z_age_depth_grad(2)         = 0.5 * log((1 + r_age_depth_grad{2}(1, 2)) / (1 - r_age_depth_grad{2}(1, 2)));
% Z_age_diff_rel(2)           = 0.5 * log((1 + r_age_diff_rel{2}(1, 2)) / (1 - r_age_diff_rel{2}(1, 2)));
% Z_age_mean(2)               = 0.5 * log((1 + r_age_mean{2}(1, 2)) / (1 - r_age_mean{2}(1, 2)));
% Z_var_age(2)                = 0.5 * log((1 + r_var_age{2}(1, 2)) / (1 - r_var_age{2}(1, 2)));

% correlation confidence bounds (1 sigma)
lim_depth_diff{1}           = tanh(Z_depth_diff(1) + [-sigma_full sigma_full]);
lim_age_diff{1}             = tanh(Z_age_diff(1) + [-sigma_full sigma_full]);
lim_age_depth_grad{1}       = tanh(Z_age_depth_grad(1) + [-sigma_full sigma_full]);
lim_age_diff_rel{1}         = tanh(Z_age_diff_rel(1) + [-sigma_full sigma_full]);
lim_age_mean{1}             = tanh(Z_age_mean(1) + [-sigma_full sigma_full]);
lim_var_age{1}              = tanh(Z_var_age(1) + [-sigma_mean sigma_mean]);
% lim_depth_diff{2}           = tanh(Z_depth_diff(2) + [-sigma_full sigma_full]);
% lim_age_diff{2}             = tanh(Z_age_diff(2) + [-sigma_full sigma_full]);
% lim_age_depth_grad{2}       = tanh(Z_age_depth_grad(2) + [-sigma_full sigma_full]);
% lim_age_diff_rel{2}         = tanh(Z_age_diff_rel(2) + [-sigma_full sigma_full]);
% lim_age_mean{2}             = tanh(Z_age_mean(2) + [-sigma_full sigma_full]);
% lim_var_age{2}              = tanh(Z_var_age(2) + [-sigma_mean sigma_mean]);

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
    
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
    
%% Depth-age scale comparison

    plots                   = {'depth_diff' '1e-3 .* age_diff' '1e2 .* age_diff_rel' '1e-3 .* age_depth_grad'};
    xlabels                 = {'Depth difference (m)' 'Age difference (ka)' 'Relative age difference (%)' 'Vertical age gradient (ka m^{-1})'};
    ranges                  = [0 1e3;
                               0 5e1;
                               0 2e2;
                               0 2e-1];
    incs                    = [1e1 5e-1 2e0 2e-3];
    cutoff                  = [6e2 2e1  1e2 1e-1];
    letters                 = 'c':'f';
    figure('position', [300 50 1200 1200])
    subplot('position', [0.08 0.71 0.41 0.27])
    hold on
    axis ij
    axis([0 120 0 1])
%     fill([get(gca, 'xlim') fliplr(get(gca, 'xlim'))], [repmat(depth_rel_min, 1, 2) repmat(depth_rel_max, 1, 2)], [0.95 0.95 0.95], 'edgecolor', 'none')
    p_core                  = zeros(1, num_core);
    for ii = 1:num_core
        plot((1e-3 .* squeeze(age_qnye(:, ii, :)))', depth_rel_all(:, ones(1, num_test)), 'linewidth', 1, 'color', colors_uncert(ii, :));
    end
    for ii = 1:num_core
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) + core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) - core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
    end
    for ii = 1:num_core
        plot((1e-3 .* core{ii}.age), (core{ii}.depth ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii})
    end
    p_ref                   = zeros(1, 3);
    p_ref(1)                = plot(0, 0, 'color', colors{end}, 'linewidth', 3);
    p_ref(2)                = plot(0, 0, '--', 'color', colors{end}, 'linewidth', 2);
    p_ref(3)                = plot(0, 0, 'color', colors_uncert(end, :), 'linewidth', 1);
    curr_ax                 = axis;
    text((curr_ax(1) + (0.03 * diff(curr_ax(1:2)))), (curr_ax(3) + (0.06 * diff(curr_ax(3:4)))), '(a)', 'fontsize', 20, 'fontweight', 'bold')
    set(gca, 'fontsize', 20)
    axis([0 120 0 1])
    xlabel('Age (ka)')
    ylabel('Depth (normalized)')
    lg                      = legend(p_ref, 'core', 'core uncertainty', 'quasi Nye');
    box on
    axes('position', [0.20 0.85 0.12 0.12])
    hold on
    axis ij
    axis([0 10 0 0.5])
%     fill([get(gca, 'xlim') fliplr(get(gca, 'xlim'))], [repmat(depth_rel_min, 1, 2) repmat(depth_rel_max, 1, 2)], [0.95 0.95 0.95], 'edgecolor', 'none')
    for ii = 1:num_core
        plot((1e-3 .* squeeze(age_qnye(:, ii, :)))', depth_rel_all(:, ones(1, num_test)), 'linewidth', 1, 'color', colors_uncert(ii, :));
    end
    for ii = 1:num_core
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) + core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) - core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
    end
    for ii = 1:num_core
        plot((1e-3 .* core{ii}.age), (core{ii}.depth ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii})
    end
    set(gca, 'fontsize', 20, 'xtick', 0:2:10, 'xticklabel', {'0' '' '' '' '' '10'}, 'ytick', 0:0.1:0.5, 'yticklabel', {'0' '' '' '' '' '0.5'})
    box on
    uistack(lg, 'top')
    subplot('position', [0.57 0.71 0.41 0.27])
    hold on
    for ii = 1:num_core
        p_core(ii)          = plot(conv_ind_stairs(-50:50), conv_dep_stairs((1e2 / (num_test * num_depth_rel_all)) .* hist((1e2 .* reshape(squeeze(age_qnye_diff_rel(:, ii, :)), (num_test * num_depth_rel_all), 1)), -49.5:49.5)), 'linewidth', 3, 'color', colors{ii});
%         fit((-49.5:49.5)', hist((1e2 .* reshape(squeeze(age_qnye_diff_rel(:, ii, :)), (num_test * num_depth_rel_all), 1)), -49.5:49.5)', 'gauss1')
%         disp([ii nanmean(age_qnye_diff_rel_curr) nanstd(age_qnye_diff_rel_curr)])
    end
    axis([-50 50 0 15])
    curr_ax                 = axis;
    text((curr_ax(1) + (0.03 * diff(curr_ax(1:2)))), (curr_ax(3) + (0.94 * diff(curr_ax(3:4)))), '(b)', 'fontsize', 20, 'fontweight', 'bold')
    set(gca, 'fontsize', 20, 'xtick', -50:10:50)
    xlabel('Difference from core age (%)')
    ylabel({'Fraction of'; '1000 simulations (%)'})
    legend(p_core, 'Camp Century', 'DYE-3', 'GISP2', 'GRIP', 'NEEM', 'NorthGRIP', 'location', 'northeast');
    box on
%     axes('position', [0.80 0.25 0.15 0.3])
%     hold on
%     for ii = 1:num_core
%         p_core(ii)          = plot(conv_ind_stairs(-50:50), conv_dep_stairs((1e2 / (num_test * num_depth_rel_all)) .* hist((1e2 .* reshape(squeeze(age_qnye_diff_rel(:, ii, :)), (num_test * num_depth_rel_all), 1)), -49.5:49.5)), 'linewidth', 3, 'color', colors{ii});
%     end
%     axis([-10 10 0 5])
%     set(gca, 'fontsize', 20)
%     box on
    colormap(jet(24))
    for ii = 1:4
        if (ii < 3)
            subplot('position', [(0.08 + ((ii - 1) * 0.49)) 0.38 0.41 0.27])
        else
            subplot('position', [(0.08 + ((ii - 3) * 0.49)) 0.055 0.41 0.27])
        end
        hold on
        var_curr            = eval(plots{ii});
        tmp                 = (1e2 ./ (num_test * num_depth_rel_all * num_core)) .* hist3([var_curr(:) (1e2 .* age_qnye_diff_rel(:))], {ranges(ii, 1):incs(ii):ranges(ii, 2) -50:50});
        tmp(~tmp)           = NaN;
        pcolor(ranges(ii, 1):incs(ii):ranges(ii, 2), -50:50, tmp')
        shading flat
        axis([ranges(ii, :) -50 50])
        caxis([0 0.25])
        plot(get(gca, 'xlim'), zeros(1, 2), '--', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        plot(repmat(cutoff(ii), 1, 2), get(gca, 'ylim'), 'm--', 'linewidth', 2)        
        curr_ax              = axis;
        text((curr_ax(1) + (0.03 * diff(curr_ax(1:2)))), (curr_ax(3) + (0.94 * diff(curr_ax(3:4)))), ['(' letters(ii) ')'], 'fontsize', 20, 'fontweight', 'bold', 'edgecolor', 'k', 'backgroundcolor', 'w')
        set(gca, 'fontsize', 20, 'ytick', -50:10:50)
        xlabel(xlabels{ii})
        ylabel({'Difference from'; 'core age (%)'})
        box on
        if (ii == 2)
            text(13, -31, {'Fraction of'; '1000 simulations (%)'}, 'color', 'k', 'fontsize', 20, 'horizontalalignment', 'center')
            cb              = colorbar('fontsize', 20, 'location', 'south', 'xtick', 0:0.05:0.25, 'xticklabel', {'0' '' '' '' '' '0.25'});
            set(cb, 'position', [0.60 0.385 0.15 0.01])
        end
    end
    
%% Depth-age scale comparison

    figure('position', [300 300 1800 600])
    subplot('position', [0.04 0.11 0.2833 0.85])
    hold on
    axis ij
    axis([0 120 0 1])
    p_core                  = zeros(1, num_core);
    for ii = 1:num_core
        plot((1e-3 .* squeeze(age_qnye(:, ii, :)))', depth_rel_all(:, ones(1, num_test)), 'linewidth', 1, 'color', colors_uncert(ii, :));
%         plot((1e-3 .* squeeze(age_qdj(:, ii, :)))', depth_rel_all(:, ones(1, num_test)), 'linewidth', 1, 'color', colors_uncert(ii, :));
    end
    for ii = 1:num_core
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) + core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) - core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
    end
    for ii = 1:num_core
        plot((1e-3 .* core{ii}.age), (core{ii}.depth ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii})
    end
    p_ref                   = zeros(1, 3);
    p_ref(1)                = plot(0, 0, 'color', colors{end}, 'linewidth', 3);
    p_ref(2)                = plot(0, 0, '--', 'color', colors{end}, 'linewidth', 2);
    p_ref(3)                = plot(0, 0, 'color', colors_uncert(end, :), 'linewidth', 1);
    curr_ax                 = axis;
    text((curr_ax(1) + (0.03 * diff(curr_ax(1:2)))), (curr_ax(3) + (0.05 * diff(curr_ax(3:4)))), '(a)', 'fontsize', 20, 'fontweight', 'bold')
    set(gca, 'fontsize', 20)
    axis([0 120 0 1])
    xlabel('Age (ka)')
    ylabel('Depth (normalized)')
    lg                      = legend(p_ref, 'core', 'core uncertainty', 'quasi Nye');
    box on
    axes('position', [0.10 0.60 0.10 0.32])
    hold on
    axis ij
    axis([0 10 0 0.5])
    for ii = 1:num_core
        plot((1e-3 .* squeeze(age_qnye(:, ii, :)))', depth_rel_all(:, ones(1, num_test)), 'linewidth', 1, 'color', colors_uncert(ii, :));
%         plot((1e-3 .* squeeze(age_qdj(:, ii, :)))', depth_rel_all(:, ones(1, num_test)), 'linewidth', 1, 'color', colors_uncert(ii, :));
    end
    for ii = 1:num_core
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) + core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
        plot((1e-3 .* (core{ii}.age(~isnan(core{ii}.age)) - core{ii}.age_uncert(~isnan(core{ii}.age)))), (core{ii}.depth(~isnan(core{ii}.age)) ./ core{ii}.depth(end)), '--', 'linewidth', 2, 'color', colors{ii});
    end
    for ii = 1:num_core
        plot((1e-3 .* core{ii}.age), (core{ii}.depth ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii})
    end
    set(gca, 'fontsize', 20, 'xtick', 0:2:10, 'xticklabel', {'0' '' '' '' '' '10'}, 'ytick', 0:0.1:0.5, 'yticklabel', {'0' '' '' '' '' '0.5'})
    box on
    uistack(lg, 'top')
    subplot('position', [0.37 0.11 0.2833 0.85])
    hold on
    for ii = 1:num_core
        p_core(ii)          = plot(conv_ind_stairs(-50:50), conv_dep_stairs((1e2 / (num_test * num_depth_rel_all)) .* hist((1e2 .* reshape(squeeze(age_qnye_diff_rel(:, ii, :)), (num_test * num_depth_rel_all), 1)), -49.5:49.5)), 'linewidth', 3, 'color', colors{ii});
%         p_core(ii)          = plot(conv_ind_stairs(-50:50), conv_dep_stairs((1e2 / (num_test * num_depth_rel_all)) .* hist((1e2 .* reshape(squeeze(age_qdj_diff_rel(:, ii, :)), (num_test * num_depth_rel_all), 1)), -49.5:49.5)), 'linewidth', 3, 'color', colors{ii});
    end
    axis([-50 50 0 15])
    curr_ax                 = axis;
    text((curr_ax(1) + (0.03 * diff(curr_ax(1:2)))), (curr_ax(3) + (0.95 * diff(curr_ax(3:4)))), '(b)', 'fontsize', 20, 'fontweight', 'bold')
    set(gca, 'fontsize', 20, 'xtick', -50:10:50)
    xlabel('Model difference from core age (%)')
    ylabel('Fraction of 1000 simulations (%)')
    legend(p_core, 'Camp Century', 'DYE-3', 'GISP2', 'GRIP', 'NEEM', 'NorthGRIP', 'location', 'northeast');
    box on
    subplot('position', [0.70 0.11 0.2833 0.85])
    hold on
    tmp                 = (1e2 ./ (num_test * num_depth_rel_all * num_core)) .* hist3([(1e2 .* age_diff_rel(:)) (1e2 .* age_qnye_diff_rel(:))], {0:2e0:2e2 -50:50});
    tmp(~tmp)           = NaN;
    pcolor(0:2e0:2e2, -50:50, tmp')
    shading flat
    axis([0 2e2 -50 50])
    caxis([0 0.25])
    plot(get(gca, 'xlim'), zeros(1, 2), '--', 'color', [0.5 0.5 0.5], 'linewidth', 2)
    curr_ax              = axis;
    text((curr_ax(1) + (0.03 * diff(curr_ax(1:2)))), (curr_ax(3) + (0.94 * diff(curr_ax(3:4)))), '(c)', 'fontsize', 20, 'fontweight', 'bold')
    set(gca, 'fontsize', 20, 'ytick', -50:10:50)
    xlabel('Relative age difference of bounding reflection pair (%)')
    ylabel('Model difference from core age (%)')
    box on
    text(45, -40, '# simulations', 'color', 'k', 'fontsize', 20, 'horizontalalignment', 'center')
    cb                  = colorbar('fontsize', 20, 'location', 'south', 'xtick', 0:0.05:0.25, 'xticklabel', {'1' '' '' '' '' '>25'});
    set(cb, 'position', [0.71 0.13 0.1 0.02])
    
%%
end