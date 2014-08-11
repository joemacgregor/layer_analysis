function chisq              = shallow_strain_fix(accum, depth, age, age_uncert, strain_rate)
% SHALLOW_STRAIN_FIX Residual between radar-observed and shallow-strain depth-age scales.
%   
%   CHISQ = SHALLOW_STRAIN_FIX(ACCUM,DEPTH,AGE,AGE_UNCERT,STRAIN_RATE)
%   calculates the Chi-squared residual CHISQ between the ages AGE (and age
%   uncertainties AGE_UNCERT) of radar-observed layers at depths DEPTH with
%   those ages calculated by a one-dimensional shallow-strain vertical
%   strain-rate model defined by the two-element vector
%   SHALLOW_STRAIN_MODEL. This model consists of a steady-state
%   accumulation rate (first element) and uniform vertical strain rate
%   (second element).
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, inverse years (1/a) for strain rate,
%   meters (m) for DEPTH and years (a) for AGE and AGE_UNCERT.
%
% Joe MacGregor (UTIG), Mark Fahnestock (UAF)
% Last updated: 05/13/14

if (nargin ~= 5)
    error('shallow_strain_fix:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 5.'])
end
if (~isnumeric(accum) || ~isscalar(accum))
    error('shallow_strain_fix:model', 'ACCUM is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('shallow_strain_fix:depth', 'DEPTH is not numeric.')
end
if ~isnumeric(age)
    error('shallow_strain_fix:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('shallow_strain_fix:ageuncert', 'AGE_UNCERT is not numeric.')
end
if (~isnumeric(strain_rate) || ~isscalar(strain_rate))
    error('shallow_strain_fix:strainrate', 'STRAIN_RATE is not a numeric scalar.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('shallow_strain_fix:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (nargout > 1)
    error('shallow_strain_fix:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

age_shallow_strain          = (1 / strain_rate) .* log(accum ./ (accum - (strain_rate .* depth)));
chisq                       = sum(((age - age_shallow_strain) ./ age_uncert) .^ 2);