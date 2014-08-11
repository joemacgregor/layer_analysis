% % CORE_AGE_COMP Compare core-located depth-age scales.
%
% Joe MacGregor (UTIG)
% Last updated: 08/07/14

clear all

plotting                    = false;

load mat/core_int name_core_short num_core
load mat/merge_all depth_smooth dist
load mat/date_all_ngrip_lin age age_type age_uncert dist_int_max
% load mat/date_all_ngrip_strain age age_type age_uncert dist_int_max
load mat/date_all int_core_merge
core                        = cell(1, num_core);
for ii = 1:num_core
    core{ii}                = load(['mat/depth_age/' name_core_short{ii} '_depth_age']);
end

colors                      = {'r' 'g' 'b' 'c' 'm', 'k'};
colors_fill                 = {'r' 'g' 'b' 'c' 'm', 'w'};

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
    
%%
    set(0, 'DefaultFigureWindowStyle', 'docked')
    
%% Depth-age scale comparison

    figure('position', [100 100 800 800])
    hold on
    axis ij
    p_core                  = zeros(1, num_core);
    for ii = 1:num_core
        switch ii
            case 1
                ind_break   = interp1(core{ii}.depth, 1:length(core{ii}.depth), 1136, 'nearest', 'extrap');
                p_core(ii)  = plot((1e-3 .* core{ii}.age(1:ind_break)), (core{ii}.depth(1:ind_break) ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii});
                plot((1e-3 .* core{ii}.age(ind_break:end)), (core{ii}.depth(ind_break:end) ./ core{ii}.depth(end)), '--', 'linewidth', 3, 'color', colors{ii})
            case 2
                ind_break   = interp1(core{ii}.depth, 1:length(core{ii}.depth), 1786, 'nearest', 'extrap');
                p_core(ii)  = plot((1e-3 .* core{ii}.age(1:ind_break)), (core{ii}.depth(1:ind_break) ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii});
                plot((1e-3 .* core{ii}.age(ind_break:end)), (core{ii}.depth(ind_break:end) ./ core{ii}.depth(end)), '--', 'linewidth', 3, 'color', colors{ii})
            otherwise
                p_core(ii)  = plot((1e-3 .* core{ii}.age), (core{ii}.depth ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii});
        end
    end
    for ii = 1:size(int_core_merge, 1)
        [jj, kk, ll]        = deal(int_core_merge(ii, 1), int_core_merge(ii, 2), int_core_merge(ii, 3));
        if (length(age{jj}{kk}) < ll)
            continue
        elseif isempty(age{jj}{kk}{ll})
            continue
        end
        curr_layer          = find(~isnan(age{jj}{kk}{ll}))';
        if isempty(curr_layer)
            continue
        end
        [dist_curr, ind_sort] ...
                            = unique(dist{jj}{kk}{ll});
        if ~isempty(find((ind_sort == int_core_merge(ii, 5)), 1))
            depth_smooth_curr ...
                            = depth_smooth{jj}{kk}{ll}(:, ind_sort);
            depth_smooth_curr ...
                            = nanmean(depth_smooth_curr(:, interp1(dist_curr, 1:length(dist_curr), (dist_curr(find((ind_sort == int_core_merge(ii, 5)), 1)) - dist_int_max), 'nearest', 'extrap'):...
                                      interp1(dist_curr, 1:length(dist_curr), (dist_curr(find((ind_sort == int_core_merge(ii, 5)), 1)) + dist_int_max), 'nearest', 'extrap')), 2);
        else
            depth_smooth_curr ...
                            = depth_smooth{jj}{kk}{ll}(:, int_core_merge(ii, 5));
        end
        curr_layer          = intersect(find(~isnan(depth_smooth_curr)), curr_layer);
        if isempty(curr_layer)
            continue
        end
        for mm = curr_layer'
            plot((1e-3 .* [(age{jj}{kk}{ll}(mm) - age_uncert{jj}{kk}{ll}(mm)) (age{jj}{kk}{ll}(mm) + age_uncert{jj}{kk}{ll}(mm))]), ...
                 repmat((depth_smooth_curr(mm) / core{int_core_merge(ii, 4)}.depth(end)), 1, 2), 'color', colors{int_core_merge(ii, 4)}, 'linewidth', 1)
            plot((1e-3 .* repmat((age{jj}{kk}{ll}(mm) - age_uncert{jj}{kk}{ll}(mm)), 1, 2)), ([(depth_smooth_curr(mm) - 20) (depth_smooth_curr(mm) + 20)] ./ core{int_core_merge(ii, 4)}.depth(end)), ...
                 'color', colors{int_core_merge(ii, 4)}, 'linewidth', 1)
            plot((1e-3 .* repmat((age{jj}{kk}{ll}(mm) + age_uncert{jj}{kk}{ll}(mm)), 1, 2)), ([(depth_smooth_curr(mm) - 20) (depth_smooth_curr(mm) + 20)] ./ core{int_core_merge(ii, 4)}.depth(end)), ...
                 'color', colors{int_core_merge(ii, 4)}, 'linewidth', 1)
            p_sym           = plot((1e-3 .* age{jj}{kk}{ll}(mm)), (depth_smooth_curr(mm) / core{int_core_merge(ii, 4)}.depth(end)), 'k^', 'markersize', 12, 'markerfacecolor', colors_fill{int_core_merge(ii, 4)}, 'linestyle', 'none');
            if (age_type{jj}{kk}{ll}(mm) <= 1)
                set(p_sym, 'marker', 'o')
            end
        end
    end
    set(gca, 'fontsize', 20)
    lg                      = legend(p_core, 'Camp Century', 'DYE-3', 'GISP2', 'GRIP', 'NEEM', 'NorthGRIP');
    set(lg, 'position', [0.21 0.68 0.21 0.22])
    axis([0 120 0 1])
    xlabel('Age (ka)')
    ylabel('Depth (normalized)')
    title('(a) linear inter/extrapolation', 'fontweight', 'bold')
%     title('(b) quasi-Nye dating', 'fontweight', 'bold')
    grid on
    box on
    axes('position', [0.48 0.41 0.4 0.49])
    hold on
    axis ij
    for ii = 1:num_core
        plot((1e-3 .* core{ii}.age), (core{ii}.depth ./ core{ii}.depth(end)), 'linewidth', 3, 'color', colors{ii})
    end
    for ii = 1:size(int_core_merge, 1)
        [jj, kk, ll]        = deal(int_core_merge(ii, 1), int_core_merge(ii, 2), int_core_merge(ii, 3));
        if (length(age{jj}{kk}) < ll)
            continue
        elseif isempty(age{jj}{kk}{ll})
            continue
        end
        curr_layer          = find(~isnan(age{jj}{kk}{ll}))';
        if isempty(curr_layer)
            continue
        end
        [dist_curr, ind_sort] ...
                            = unique(dist{jj}{kk}{ll});
        if ~isempty(find((ind_sort == int_core_merge(ii, 5)), 1))
            depth_smooth_curr ...
                            = depth_smooth{jj}{kk}{ll}(:, ind_sort);
            depth_smooth_curr ...
                            = nanmean(depth_smooth_curr(:, interp1(dist_curr, 1:length(dist_curr), (dist_curr(find((ind_sort == int_core_merge(ii, 5)), 1)) - dist_int_max), 'nearest', 'extrap'):...
                                      interp1(dist_curr, 1:length(dist_curr), (dist_curr(find((ind_sort == int_core_merge(ii, 5)), 1)) + dist_int_max), 'nearest', 'extrap')), 2);
        else
            depth_smooth_curr ...
                            = depth_smooth{jj}{kk}{ll}(:, int_core_merge(ii, 5));
        end
        curr_layer          = intersect(find(~isnan(depth_smooth_curr)), curr_layer);
        if isempty(curr_layer)
            continue
        end
        for mm = curr_layer'
            plot((1e-3 .* [(age{jj}{kk}{ll}(mm) - age_uncert{jj}{kk}{ll}(mm)) (age{jj}{kk}{ll}(mm) + age_uncert{jj}{kk}{ll}(mm))]), ...
                 repmat((depth_smooth_curr(mm) / core{int_core_merge(ii, 4)}.depth(end)), 1, 2), 'color', colors{int_core_merge(ii, 4)}, 'linewidth', 1)
            plot((1e-3 .* repmat((age{jj}{kk}{ll}(mm) - age_uncert{jj}{kk}{ll}(mm)), 1, 2)), ...
                 ([(depth_smooth_curr(mm) - 20) (depth_smooth_curr(mm) + 20)] ./ core{int_core_merge(ii, 4)}.depth(end)), 'color', colors{int_core_merge(ii, 4)}, 'linewidth', 1)
            plot((1e-3 .* repmat((age{jj}{kk}{ll}(mm) + age_uncert{jj}{kk}{ll}(mm)), 1, 2)), ([(depth_smooth_curr(mm) - 20) (depth_smooth_curr(mm) + 20)] ./ core{int_core_merge(ii, 4)}.depth(end)), ...
                 'color', colors{int_core_merge(ii, 4)}, 'linewidth', 1)
            p_sym           = plot((1e-3 .* age{jj}{kk}{ll}(mm)), (depth_smooth_curr(mm) / core{int_core_merge(ii, 4)}.depth(end)), 'k^', 'markersize', 12, 'markerfacecolor', colors_fill{int_core_merge(ii, 4)}, 'linestyle', 'none');
            if (age_type{jj}{kk}{ll}(mm) <= 1)
                set(p_sym, 'marker', 'o')
            end
        end
    end
    p_sym(1)                = plot(3.9, 0.025, 'ko', 'markersize', 12);
    p_sym(2)                = plot(3.9, 0.07, 'k^', 'markersize', 12);
    text(4.35, 0.0225, 'NorthGRIP-dated', 'fontsize', 20, 'color', 'k')
    text(4.35, 0.0675, 'reflector-inferred', 'fontsize', 20, 'color', 'k')
    axis([0 10 0 0.5])
    set(gca, 'fontsize', 20)
    grid on
    box on
%%
end