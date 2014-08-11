% STRAIN_MODEL_SENSITIVITY Evaluate model-parameter sensitivity of various strain-rate models.
% 
% Joe MacGregor (UTIG)
% Last updated: 12/08/13

clear

plotting                    = false;

% ice column properties
thick                       = 2500;
depth                       = (0:10:thick)';
num_depth                   = length(depth);

% parameter reference values
accum_ref                   = 0.25;
strain_rate_ref             = 1e-5;
thick_shear_ref             = 200;
melt_bed_ref                = 0.01;
frac_slide_ref              = 0.1;

% parameter ranges
accum                       = 0.05:0.05:0.55;
strain_rate                 = 0:1e-5:1e-4;
thick_shear                 = 0:200:2000;
melt_bed                    = 0:0.01:0.1;
frac_slide                  = 0:0.1:1;

models                      = {'nye' 'shallow_strain' 'nye_melt' 'dj' 'dj_melt'};
models_nice                 = {'Nye' 'Shallow strain' 'Nye + basal melt' 'Dansgaard-Johnsen' 'Dansgaard-Johnsen + basal melt'};
params                      = {'accum' 'thick_shear' 'melt_bed' 'frac_slide' 'strain_rate'};
params_nice                 = {'$\dot{b}$' '$h$' '$\dot{m}$' '$s$' '$\dot{\epsilon}$'};
for ii = 1:length(params)
    params_nice{ii}         = [params_nice{ii} '; ' sprintf('%3.2g', eval([params{ii} '(1)'])) ':' sprintf('%3.2g', mean(diff(eval(params{ii})))) ':' sprintf('%3.2g', eval([params{ii} '(end)'])) ...
                               ' (' sprintf('%3.2g', eval([params{ii} '_ref'])) ')\,\,\,\,\,\,'];
end
param_set                   = [true false false false false;
                               true false false false true;
                               true false true  false false;
                               true true  false false false
                               true true  true  true  false];

% range colors
colors                      = [0 0 0;
                               1 0 0;
                               0 0 1;
                               0 1 0;
                               1 0 1];
colors_all                  = cell(1, 5);
for ii = 1:5
    colors_all{ii}          = zeros(11, 3);
    for jj = 1:11
        colors_all{ii}(jj, :) ...
                            = colors(ii, :) + ((jj - 1) .* [0.085 0.085 0.085]);
    end
    colors_all{ii}(colors_all{ii} > 1) ...
                            = 1;
    colors_all{ii}          = flipud(colors_all{ii});
end

letters                     = 'A':'J';

% initialize modeled ages
for ii = 1:length(models)
    for jj = 1:length(params)
        if param_set(ii, jj)
            eval([models{ii} '_' params{jj} ' = NaN(num_depth, 10);'])
        end
    end
end

% calculate model ages
for ii = 1:11
    dj_accum(:, ii)                     = dj([accum(ii) thick_shear_ref], thick, depth); %#ok<*SAGROW>
    dj_thick_shear(:, ii)               = dj([accum_ref thick_shear(ii)], thick, depth);
    dj_melt_accum(:, ii)                = dj_melt([accum(ii) thick_shear_ref melt_bed_ref frac_slide_ref], thick, depth);
    dj_melt_frac_slide(:, ii)           = dj_melt([accum_ref thick_shear_ref melt_bed_ref frac_slide(ii)], thick, depth);
    dj_melt_melt_bed(:, ii)             = dj_melt([accum_ref thick_shear_ref melt_bed(ii) frac_slide_ref], thick, depth);
    dj_melt_thick_shear(:, ii)          = dj_melt([accum_ref thick_shear(ii) melt_bed_ref frac_slide_ref], thick, depth);
    nye_accum(:, ii)                    = nye(accum(ii), thick, depth);
    nye_melt_accum(:, ii)               = nye_melt([accum(ii) melt_bed_ref], thick, depth);
    nye_melt_melt_bed(:, ii)            = nye_melt([accum_ref melt_bed(ii)], thick, depth);
    shallow_strain_accum(:, ii)         = shallow_strain([accum(ii) strain_rate_ref], depth);
    shallow_strain_strain_rate(:, ii)   = shallow_strain([accum_ref strain_rate(ii)], depth);
end

% NaN out Infs
for ii = 1:length(models)
    for jj = 1:length(params)
        if param_set(ii, jj)
            eval([models{ii} '_' params{jj} '(isinf(' models{ii} '_' params{jj} ')) = NaN;'])
        end
    end
end

% plot ranges
age_range                   = [0 120; 0 9];
depth_max                   = [thick 1500];
incs                        = [{0:40:120} {0:3:9}];

%%
if plotting
%%
    set(0, 'DefaultFigureWindowStyle', 'default') %#ok<UNRCH>
%%
%     set(0, 'DefaultFigureWindowStyle', 'docked')
%%
    figure('position', [10 10 1800 800])
    for ii = 1:5
        for jj = 1:2
            subplot('position', [(0.07 + ((ii - 1) * 0.18)) (0.53 - ((jj - 1) * 0.46)) 0.18 0.44])
            hold on
            axis([age_range(jj, :) 0 depth_max(jj)])
            axis ij square
            p               = zeros(1, length(params));
            for kk = find(param_set(ii, :))
                for ll = 1:11
                    p(kk)   = plot((1e-3 .* eval([models{ii} '_' params{kk} '(:, ll)'])), depth, 'linewidth', 3, 'color', colors_all{kk}(ll, :));
                end
            end
            axis([age_range(jj, :) 0 depth_max(jj)])
            set(gca, 'fontsize', 20)
            if (ii == 1)
                ylabel('Depth (m)')
                set(gca, 'xtick', incs{jj})
                if (jj == 2)
                    text(21.5, 1700, 'Age (ka)', 'fontsize', 20, 'color', 'k')
                end
            else
                set(gca, 'xtick', incs{jj}(2:end), 'yticklabel', {})
            end
            if (jj == 1)
                if (ii == 1)
                    legend(p(logical(p)), {[strtok(params_nice{1}, '(') '\,\,\,\,\,\,']}, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', 18)
                else
                    legend(p(logical(p)), params_nice(param_set(ii, :)), 'location', 'northeast', 'interpreter', 'latex', 'fontsize', 18)
                end
                title(models_nice{ii}, 'fontweight', 'bold')
            end
            t               = text((0.025 * age_range(jj, 2)), (0.95 * depth_max(jj)), letters(ii + ((jj - 1) * 5)), 'fontsize', 20, 'color', 'k', 'fontweight', 'bold');
            if ((ii == 2) && (jj == 1))
                set(t, 'color', 'w')
            end
            grid on
            box on
        end
    end
%%    
figure('position', [10 10 640 640])
hold on
plot((1e-3 .* nye_accum(:, (accum == accum_ref))), depth, 'k', 'linewidth', 3)
plot((1e-3 .* dj_thick_shear(:, end)), depth, 'r', 'linewidth', 3)
plot((1e-3 .* nye_melt_melt_bed(:, end)), depth, 'b', 'linewidth', 3)
axis ij square
axis([age_range(1, :) 0 depth_max(1)])
text(43, 2375, 'accumulation rate = 25 cm a^{-1}', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')
text(3, 1200, 'basal melt rate = 10 cm a^{-1}', 'color', 'b', 'fontsize', 20, 'fontweight', 'bold', 'rotation', -85)
text(45, 1950, 'basal shear layer = 2000 m', 'color', 'r', 'fontsize', 20, 'fontweight', 'bold', 'rotation', -10)
set(gca, 'fontsize', 20)
xlabel('Age (ka)')
ylabel('Depth (m)')
box on
axes('position', [0.42 0.42 0.45 0.45], 'fontsize', 20)
hold on
plot((1e-3 .* nye_accum(:, (accum == accum_ref))), depth, 'k', 'linewidth', 3)
plot((1e-3 .* dj_thick_shear(:, end)), depth, 'r', 'linewidth', 3)
plot((1e-3 .* nye_melt_melt_bed(:, end)), depth, 'b', 'linewidth', 3)
axis ij square
axis([age_range(2, :) 0 depth_max(2)])
set(gca, 'fontsize', 20, 'xtick', incs{2})
legend('Nye', 'Dansgaard-Johnsen', 'Nye + melt')
box on
%%
figure('position', [2600 100 640 640])
hold on
plot((1e-3 .* nye_accum(:, (accum == accum_ref))), depth, 'k', 'linewidth', 3)
plot((1e-3 .* dj_thick_shear(:, end)), depth, 'r', 'linewidth', 3)
plot((1e-3 .* nye_melt_melt_bed(:, end)), depth, 'b', 'linewidth', 3)
plot((1e-3 .* shallow_strain_strain_rate(:, 10)), depth, 'g', 'linewidth', 3)
axis ij square
axis([age_range(1, :) 0 depth_max(1)])
text(43, 2375, 'accumulation rate = 25 cm a^{-1}', 'color', 'k', 'fontsize', 20, 'fontweight', 'bold')
text(3, 1200, 'basal melt rate = 10 cm a^{-1}', 'color', 'b', 'fontsize', 20, 'fontweight', 'bold', 'rotation', -85)
text(45, 1950, 'basal shear layer = 2000 m', 'color', 'r', 'fontsize', 20, 'fontweight', 'bold', 'rotation', -10)
text(14, 1650, '9 \times 10^{-5} a^{-1}', 'color', 'g', 'fontsize', 20, 'fontweight', 'bold', 'rotation', -65)
set(gca, 'fontsize', 20)
xlabel('Age (ka)')
ylabel('Depth (m)')
box on
axes('position', [0.42 0.42 0.45 0.45], 'fontsize', 20)
hold on
plot((1e-3 .* nye_accum(:, (accum == accum_ref))), depth, 'k', 'linewidth', 3)
plot((1e-3 .* dj_thick_shear(:, end)), depth, 'r', 'linewidth', 3)
plot((1e-3 .* nye_melt_melt_bed(:, end)), depth, 'b', 'linewidth', 3)
plot((1e-3 .* shallow_strain_strain_rate(:, 10)), depth, 'g', 'linewidth', 3)
axis ij square
axis([age_range(2, :) 0 depth_max(2)])
set(gca, 'fontsize', 20, 'xtick', incs{2})
legend('Nye', 'Dansgaard-Johnsen', 'Nye + melt', 'Shallow strain')
box on
%%    
end