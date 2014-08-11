function chisq              = shallow_strain_fit(shallow_strain_model, depth, age, age_uncert)
% SHALLOW_STRAIN_FIT Residual between radar-observed and shallow-strain depth-age scales.
%   
%   CHISQ = SHALLOW_STRAIN_FIT(SHALLOW_STRAIN_MODEL,DEPTH,AGE,AGE_UNCERT)
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
% Last updated: 10/30/13

if (nargin ~= 4)
    error('shallow_strain_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 4.'])
end
if (~isnumeric(shallow_strain_model) || ~isvector(shallow_strain_model) || (numel(shallow_strain_model) ~= 2))
    error('shallow_strain_fit:model', 'SHALLOW_STRAIN_MODEL is not a two-element numeric vector.')
end
if ~isnumeric(depth)
    error('shallow_strain_fit:depth', 'DEPTH is not numeric.')
end
if ~isnumeric(age)
    error('shallow_strain_fit:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('shallow_strain_fit:ageuncert', 'AGE_UNCERT is not numeric.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('shallow_strain_fit:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (nargout > 1)
    error('shallow_strain_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

age_shallow_strain          = (1 / shallow_strain_model(2)) .* log(shallow_strain_model(1) ./ (shallow_strain_model(1) - (shallow_strain_model(2) .* depth)));
chisq                       = sum(((age - age_shallow_strain) ./ age_uncert) .^ 2);