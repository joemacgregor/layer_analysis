% HF_DENSITY Herron & Langway (1980) density profile.
% 
% This function models a 2-stage density profile based on the empirical work  
% of Herron & Langway (1980) and known surface glaciological conditions.
% Review that paper for full details.
% 
% densities in Mg/m^3
% temperatures in K
% accumulation rate in m/a ice eq.
% 
% Joe MacGregor
% Last updated: 06/24/09

function density            = hl_density(depth, temp_surf, accum_rate, density_surf)

ideal_gas_const             = 8.315; % in J/K/mol
density_ice                 = 0.917; % full ice density
density_crit                = 0.550; % critical density according to Herron and Langway (1980)

log_density_diff_pre_crit   = log(density_surf ./ (density_ice - density_surf));
log_density_diff_post_crit  = log(density_crit ./ (density_ice - density_crit));
rate_const_pre_crit         = 11 * exp(-10160 ./ (ideal_gas_const .* temp_surf));
rate_const_post_crit        = 575 * exp(-21400 ./ (ideal_gas_const .* temp_surf)); % can you say "empirical"?

% depth of critical density
depth_crit                  = (1 ./ (density_ice .* rate_const_pre_crit)) .* (log_density_diff_post_crit - log_density_diff_pre_crit);

% depth indices before (pre) and after (post) the critical density
inds_precrit                = depth <= depth_crit;
inds_postcrit               = depth > depth_crit;

exp_density                 = zeros(length(depth), 1);
exp_density(inds_precrit)   = exp(((density_ice * rate_const_pre_crit) .* depth(inds_precrit)) + log_density_diff_pre_crit);
exp_density(inds_postcrit)  = exp((((density_ice * rate_const_post_crit) .* (depth(inds_postcrit) - depth_crit)) ./ sqrt(accum_rate)) + log_density_diff_post_crit);

density                     = (density_ice .* exp_density) ./ (1 + exp_density);