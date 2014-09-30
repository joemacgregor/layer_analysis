function [speed_bal, speed_bal_surf, accum, strain_rate, thick_shear, shape_fact, depth, vel_vert] ...
                            = bal_vel(age, x, y, accum, depth, strain_rate, thick_shear, mask_D1, mask_gris, thick, az_sin_abs_rel, az_cos_abs_rel, az_sin_sign, az_cos_sign, ind_sort, length_grd, decim, do_filt, scale_type, filt_decim, filt_decim_ref)
% BAL_VEL Age-constrained balance velocity.
% 
% Joe MacGregor (UTIG)
% Last updated: 06/19/14

% decimation indices

ind_x                       = find(~mod(x(1, :), decim), 1):decim:find(~mod(x(1, :), decim), 1, 'last');
ind_y                       = find(~mod(y(:, 1), decim), 1):decim:find(~mod(y(:, 1), decim), 1, 'last');

% grid sizes
[num_y, num_x]              = size(x);
[num_decim_y, num_decim_x]  = deal(length(ind_y), length(ind_x));

% % decimate grids
% if ~isnan(age)
%     [accum, depth, strain_rate] ...
%                             = deal(accum(ind_y, ind_x), depth(ind_y, ind_x), strain_rate(ind_y, ind_x));
% end
% thick_shear                 = thick_shear(ind_y, ind_x);

% filter grids
if do_filt
    vars                    = {'thick_shear'};
    if ~isnan(age)
        vars                = [vars 'accum' 'depth' 'strain_rate'];
    end
    num_var                 = length(vars);
    
    switch scale_type
        case 'ind'
            for ii = 1:num_var
                eval([vars{ii} ' = filter2(filt_decim, ' vars{ii} ');'])
            end
        case 'thick'
            for ii = 1:num_var
                eval([vars{ii} '_tmp = NaN(num_decim_y, num_decim_x);'])
            end
            for jj = 1:num_decim_x
                if ~mod(jj, 50)
                    disp(jj)
                end
                for ii = 1:num_decim_y
                    filt_decim_curr ...
                            = filt_decim{filt_decim_ref(ii, jj)};
                    num_filt_curr ...
                            = (size(filt_decim_curr, 1) - 1) / 2;
                    
                    ind_tmp = -num_filt_curr:num_filt_curr;
                    ii_tmp  = repmat((ind_y(ii) + ind_tmp)', ((2 * num_filt_curr) + 1), 1);
                    jj_tmp  = ind_x(jj) + ind_tmp(ones(((2 * num_filt_curr) + 1), 1), :);
                    jj_tmp  = jj_tmp(:);
                    ind_good= find((ii_tmp > 0) & (jj_tmp > 0) & (ii_tmp <= num_y) & (jj_tmp <= num_x));
                    ind_curr= sub2ind([num_y num_x], ii_tmp(ind_good), jj_tmp(ind_good)); %#ok<NASGU>
                    
                    for kk = 1:num_var
                        eval([vars{kk} '_tmp(ii, jj) = nansum(filt_decim_curr(ind_good) .* ' vars{kk} '(ind_curr)) * ((numel(filt_decim_curr) ^ 2) / (length(ind_good) * length(find(~isnan(' vars{kk} '(ind_curr))))));'])
                    end
                end
            end
            for ii = 1:num_var
                eval([vars{ii} ' = ' vars{ii} '_tmp;'])
            end
            clear *_tmp
    end
end

% mask D<=1 and GrIS
[accum(~mask_D1 | ~mask_gris), depth(~mask_gris), strain_rate(~mask_D1 | ~mask_gris), thick_shear(~mask_D1 | ~mask_gris)] ...
                            = deal(NaN);

% shape factor for ice column down to isochrone
shape_fact                  = NaN(size(accum));
shape_fact(depth <= (thick - thick_shear)) ...
                            = 1;
if ~isempty(find((depth > (thick - thick_shear)), 1))
    ind_deep                = find(depth > (thick - thick_shear));
    shape_fact(ind_deep)    = 1 - (((depth(ind_deep) - (thick(ind_deep) - thick_shear(ind_deep))) .^ 2) ./ (2 .* thick_shear(ind_deep) .* depth(ind_deep)));
end
shape_fact(shape_fact < 0)  = NaN;

% vertical velocity at age_ref or zero if full thickness
if ~isnan(age)
    vel_vert                = accum .* exp(-strain_rate .* age);
else
    vel_vert                = zeros(size(accum));
end

% vertical input minus output rate locally and sanity check
vert_rate                   = accum - vel_vert;
vert_rate(vert_rate < 0)    = 0;
vert_rate(vert_rate > accum)= accum(vert_rate > accum);

% initial balance flux (vertical only), flux from surface accumulation subtracted by vertical velocity through the isochrone, m^3 ice/yr
flux_bal                    = (length_grd ^ 2) .* vert_rate;

% matrix indices of good input fluxes
ind_sort                    = ind_sort(~isnan(flux_bal(ind_sort)));
[jj, kk]                    = ind2sub([num_decim_y num_decim_x], ind_sort);

% sum balance flux by descending elevation
for ii = 1:length(ind_sort)
    % add flux to adjacent cell (column-wise)
    if ~isnan(az_sin_sign(jj(ii), kk(ii)))
        flux_bal((jj(ii) + az_sin_sign(jj(ii), kk(ii))), kk(ii)) ...
                            = flux_bal((jj(ii) + az_sin_sign(jj(ii), kk(ii))), kk(ii)) + (az_sin_abs_rel(jj(ii), kk(ii)) * flux_bal(jj(ii), kk(ii)));
    end
    % add flux to adjacent cell (row-wise)
    if ~isnan(az_cos_sign(jj(ii), kk(ii)))
        flux_bal(jj(ii), (kk(ii) + az_cos_sign(jj(ii), kk(ii)))) ...
                            = flux_bal(jj(ii), (kk(ii) + az_cos_sign(jj(ii), kk(ii)))) + (az_cos_abs_rel(jj(ii), kk(ii)) * flux_bal(jj(ii), kk(ii)));
    end
end

% balance velocity (also surface-corrected)
speed_bal                   = flux_bal ./ (length_grd .* depth);
speed_bal_surf              = speed_bal ./ shape_fact;