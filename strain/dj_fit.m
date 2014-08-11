function chisq              = dj_fit(dj_model, thick, depth, age, age_uncert)
% DJ_FIT Residual between radar-observed and Dansgaard-Johnsen depth-age scales.
%   
%   CHISQ = DJ_FIT(DJ_MODEL,THICK,DEPTH,AGE,AGE_UNCERT) calculates the
%   Chi-squared residual CHISQ between the ages AGE (and age uncertainties
%   AGE_UNCERT) of radar-observed layers at depths DEPTH with those ages
%   calculated by a one-dimensional Dansgaard-Johnsen vertical strain-rate
%   model. This model requires a scalar ice thickness THICK and a
%   two-element vector DJ_MODEL. The first element of DJ_MODEL is the
%   steady-state accumulation rate and the second element is the thickness
%   of the basal shear layer, aka "little h" or the height of the shear
%   layer above the bed.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK, DEPTH
%   shear-layer thickness, and years (a) for AGE and AGE_UNCERT.
%   
%   References:
%   
%   Fahnestock, M.A. et al., 2001, Internal layer-tracing and
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
    error('dj_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 5.'])
end
if (~isvector(dj_model) || ~isnumeric(dj_model) || (numel(dj_model) ~= 2))
    error('dj_fit:dj_model', 'DJ_MODEL is not a two-element vector.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj_fit:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('dj_fit:depth', 'DEPTH is not numeric.')
end
if ~isnumeric(age)
    error('dj_fit:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('dj_fit:ageuncert', 'AGE_UNCERT is not numeric.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('dj_fit:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (nargout > 1)
    error('dj_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

age_dj                      = zeros(length(depth), 1);
age_dj(depth <= (thick - dj_model(2))) ...
                            = (((2 * thick) - dj_model(2)) / (2 * dj_model(1))) .* log(((2 * thick) - dj_model(2)) ./ ((2 .* (thick - depth(depth <= (thick - dj_model(2))))) - dj_model(2)));
age_h                       = (((2 * thick) - dj_model(2)) / (2 * dj_model(1))) .* log(((2 * thick) - dj_model(2)) ./ dj_model(2));
age_dj(depth > (thick - dj_model(2))) ...
                            = ((((2 * thick) - dj_model(2)) / dj_model(1)) .* ((dj_model(2) ./ (thick - depth(depth > (thick - dj_model(2))))) - 1)) + age_h;

chisq                       = sum(((age - age_dj) ./ age_uncert) .^ 2);