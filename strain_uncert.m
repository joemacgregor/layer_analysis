function varargout          = strain_uncert(frac_test, conf_bound, model, thick, depth, age, age_uncert, model_ref, res_ref, lon_strain_rate, lat_strain_rate)
% STRAIN_UNCERT Calculate model-parameter confidence bounds for a strain-rate model.
% 
% [PARAM1,PARAM2,...] =
% STRAIN_UNCERT(FRAC_TEST,CONF_BOUND,MODEL,THICK,DEPTH,AGE,AGE_UNCERT,MODEL_REF,RES_REF,LON_STRAIN_RATE,LAT_STRAIN_RATE)
% calculates the confidence bounds for 1D ice-flow model parameters.
% FRAC_TEST is the vector of fractions of the best-fit (i.e., reference)
% model parameters about which to test, CONF_BOUND is the scalar desired
% confidence-bound fraction (e.g., 95% or 0.95), MODEL is the string
% representing the strain-rate model, THICK is the scalar ice thickness
% (m), DEPTH is the column vector of dated layer depths (m), AGE is the
% column vector of layer ages (a), AGE_UNCERT is the column vector of layer
% age uncertainties (a), MODEL_REF is the vector of best-fit model
% parameters, RES_REF is the scalar chi-squared residual of the best-fit
% model. If MODEL='dj_lat' or 'dj_surf', the longitudinal and lateral
% strain rates (LON_STRAIN_RATE and LAT_STRAIN_RATE, respectively) must
% also be provided.
% 
% STRAIN_UNCERT is called by STRAIN_LAYERS when DO_UNCERT=TRUE.
% 
% See also DJ_FIT, DJ_LAT_FIT, DJ_MELT_FIT, DJ_SURF_FIT,NYE_FIT,
% NYE_MELT_FIT, SHALLOW_STRAIN_FIT and STRAIN_LAYERS.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/19/15

if ~any(nargin == [9 11])
    error('strain_uncert:nargin', 'Incorrect number of input arguments (should be 9 or 11).')
end
if ~exist('delta_chisq', 'file')
    error('strain_uncert:delta_chisq', 'Function DELTA_CHISQ is not available within this user''s path.')
end
if (~isnumeric(frac_test) || ~isvector(frac_test))
    error('strain_uncert:frac_test', 'FRAC_TEST is not a numeric vector.')
end
if (~isnumeric(conf_bound) || ~isscalar(conf_bound))
    error('strain_uncert:tol', 'CONF_BOUND is not a numeric scalar.')
end
if ((conf_bound <= 0) || (conf_bound >= 1))
    error('strain_uncert:tol', 'CONF_BOUND is not between 0 and 1.')
end
if ~ischar(model)
    error('strain_uncert:modelchar', 'MODEL is not a string.')
end
if ~any(strcmp(model, {'dj' 'dj_lat' 'dj_melt' 'dj_surf', 'nye' 'nye_melt' 'shallow_strain'}))
    error('strain_uncert:model', 'MODEL is not one of the recognized strain-rate models: ''dj'', ''dj_lat'', ''dj_melt'', ''dj_surf'', ''nye'', ''nye_melt'' or ''shallow_strain'').')
end
if (~isnumeric(thick) || ~isscalar(thick))
    error('strain_uncert:thick', 'THICK is not a numeric scalar.')
end
if (thick <= 0)
    error('strain_uncert:thickpos', 'THICK is not positive.')
end
if (~isnumeric(depth) || ~iscolumn(depth))
    error('strain_uncert:depth', 'DEPTH is not a numeric column vector.')
end
if (~isnumeric(age) || ~iscolumn(age))
    error('strain_uncert:age', 'AGE is not a numeric column vector.')
end
if (~isnumeric(age_uncert) || ~iscolumn(age_uncert))
    error('strain_uncert:age_uncert', 'AGE_UNCERT is not a numeric column vector.')
end
if (~isnumeric(model_ref) || ~isvector(model_ref))
    error('strain_uncert:model_ref', 'MODEL_REF is not a numeric vector.')
end
if (length(model_ref) ~= nargout)
    error('strain_uncert:model_refnargout', 'Number of elements in MODEL_REF does not match number of parameter confidence ranges to be calculated.')
end
if (~isnumeric(res_ref) || ~isscalar(res_ref))
    error('strain_uncert:res_ref', 'RES_REF is not a numeric scalar.')
end
if (res_ref <= 0)
    error('strain_uncert:res_refpos', 'RES_REF is not positive.')
end
if (any(strcmp(model, {'dj_lat' 'dj_surf'})) && (nargin ~= 11))
    error('strain_uncert:djlatnargin', 'MODEL=''dj_lat/surf'' but LON_SPREAD_RATE and LAT_SPREAD_RATE are not arguments.')
end
if (any(strcmp(model, {'dj_lat' 'dj_surf'})) && (~isscalar(lon_strain_rate) || ~isnumeric(lon_strain_rate)))
    error('strain_uncert:djlonnargin', 'LON_STRAIN_RATE is not a numeric scalar.')
end
if (any(strcmp(model, {'dj_lat' 'dj_surf'})) && (~isscalar(lat_strain_rate) || ~isnumeric(lat_strain_rate)))
    error('strain_uncert:djlatnargin', 'LAT_STRAIN_RATE is not a numeric scalar.')
end
num_param                   = length(model_ref);
model_bound                 = NaN(2, num_param);
num_test                    = length(frac_test);

chisq_diff                  = delta_chisq(conf_bound, ((2 * length(age)) - num_param)); % chi-squared difference given number of degrees of freedom

for ii = 1:num_param
    
    % restart parameter fraction range to examine
    jj                      = 1;
    
    % keep searching for min/max confidence bounds while fraction ranges left to test
    while (any(isnan(model_bound(:, ii))) && (jj <= num_test))
        
        % parameter space to test based on fraction range
        param_test          = linspace((model_ref(ii) - (frac_test(jj) * model_ref(ii))), (model_ref(ii) + (frac_test(jj) * model_ref(ii))));
        param_test(sign(param_test) ~= sign(model_ref(ii))) ...
                            = NaN; % NaN out values with wrong sign
        model_curr          = model_ref; % start out with reference model (in case of other parameters)
        res                 = NaN(1, 100);
        for kk = find(~isnan(param_test)) % loop through current parameter space
            model_curr(ii)  = param_test(kk); % adjust current strain-rate model
            switch model
                case {'dj_lat' 'dj_surf'}
                    res(kk) = eval([model '_fit(model_curr, thick, depth, age, age_uncert, lon_strain_rate, lat_strain_rate)']);
                case 'shallow_strain'
                    res(kk) = eval([model '_fit(model_curr, depth, age, age_uncert)']); % chi-squared of current model
                otherwise
                    res(kk) = eval([model '_fit(model_curr, thick, depth, age, age_uncert)']);
            end
        end
        
        % keep only the residuals that make sense
        res_good1           = find(isreal(res(1:50)) & ~isnan(res(1:50)) & ~isinf(res(1:50)));
        res_good2           = 50 + find(isreal(res(51:end)) & ~isnan(res(51:end)) & ~isinf(res(51:end)));
        
        % determine if residual range crosses the confidence bound target
        if (isnan(model_bound(1, ii)) && (length(res_good1) > 1))
            model_bound(1, ii) ...
                            = model_ref(ii) + interp1((res(res_good1) - res_ref), (param_test(res_good1) - model_ref(ii)), chisq_diff);
        end
        if (isnan(model_bound(2, ii)) && (length(res_good2) > 1))
            model_bound(2, ii) ...
                            = model_ref(ii) + interp1((res(res_good2) - res_ref), (param_test(res_good2) - model_ref(ii)), chisq_diff);
        end
        
        % increment fraction range
        jj                  = jj + 1;
    end
end

varargout                   = cell(1, nargout);
for ii = 1:nargout
    varargout{ii}           = model_bound(:, ii);
end