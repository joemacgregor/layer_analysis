function chisq              = nye_melt_fit(nye_melt_model, thick, depth, age, age_uncert)
% NYE_MELT_FIT Residual between radar-observed and Nye melt-corrected depth-age scales.
%   
%   CHISQ = NYE_MELT_FIT(NYE_MELT_MODEL,THICK,DEPTH,AGE,AGE_UNCERT)
%   calculates the chi-squared residual CHISQ between the ages AGE (and age
%   uncertainties AGE_UNCERT) of radar-observed layers at depths DEPTH with
%   those ages calculated by a one-dimensional Nye vertical strain-rate
%   model that allows for basal melting. This model requires a scalar ice
%   thickness THICK and a two-element vector NYE_MELT_MODEL. The first
%   element of NYE_MELT_MODEL is the steady-state accumulation rate and the
%   second element is the basal melt rate.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for surface-accumulation and basal melt rates, meters (m)
%   for THICK and DEPTH and years (a) for AGE and AGE_UNCERT.
%   
%   References:
%   
%   Fahnestock, M.A. et al, 2001, Internal layer-tracing and
%   age-depth-accumulation relationships for the northern Greenland ice
%   sheet, J. Geophys. Res., 106(D24), 33789-33797.
%   
%   Fahnestock, M.A. et al., 2001, High geothermal heat flow, basal melt, 
%   and the origin of rapid ice flow in central Greenland, Science, 294,
%   2338-2342.
%
% Joe MacGregor (UTIG), Mark Fahnestock (UAF)
% Last updated: 10/30/13

if (nargin ~= 5)
    error('nye_melt_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 5.'])
end
if (~isvector(nye_melt_model) || ~isnumeric(nye_melt_model) || (numel(nye_melt_model) ~= 2))
    error('nye_melt_fit:nyemeltmodel', 'NYE_MELT_MODEL is not a two-element vector.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('nye_melt_fit:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('nye_melt_fit:depth', 'DEPTH is not numeric.')
end
if ~isnumeric(age)
    error('nye_melt_fit:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('nye_melt_fit:ageuncert', 'AGE_UNCERT is not numeric.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('nye_melt_fit:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (nargout > 1)
    error('nye_melt_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

age_nye_melt                = -(thick / (nye_melt_model(1) - nye_melt_model(2))) .* log((((1 - (nye_melt_model(2) / nye_melt_model(1))) ./ thick) .* (thick - depth)) + (nye_melt_model(2) / nye_melt_model(1)));
chisq                       = sum(((age - age_nye_melt) ./ age_uncert) .^ 2);