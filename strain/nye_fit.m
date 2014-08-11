function chisq              = nye_fit(accum, thick, depth, age, age_uncert)
% NYE_FIT Residual between radar-observed and standard Nye depth-age scales.
%   
%   CHISQ = NYE_FIT(ACCUM,THICK,DEPTH,AGE,AGE_UNCERT) calculates the
%   Chi-squared residual CHISQ between the ages AGE (and age uncertainties
%   AGE_UNCERT) of radar-observed layers at depths DEPTH with those ages
%   calculated by a one-dimensional Nye vertical strain-rate model. This
%   model requires a scalar steady-state accumulation rate ACCUM and ice
%   thickness THICK.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK and DEPTH and
%   years (a) for AGE and AGE_UNCERT.
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
    error('nye_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 5.'])
end
if (~isscalar(accum) || ~isnumeric(accum))
    error('nye_fit:accum', 'ACCUM is not a numeric scalar.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('nye_fit:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('nye_fit:depth', 'DEPTH is not numeric.')
end
if ~isnumeric(age)
    error('nye_fit:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('nye_fit:ageuncert', 'AGE_UNCERT is not numeric.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('nye_fit:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (nargout > 1)
    error('nye_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

age_nye                     = -(thick / accum) .* log((thick - depth) ./ thick);
chisq                       = sum(((age - age_nye) ./ age_uncert) .^ 2);