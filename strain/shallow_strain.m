function age              = shallow_strain(shallow_model, depth)
% SHALLOW_STRAIN Shallow-strain depth-age scales.
%   
%   AGE = SHALLOW_STRAIN_FIT(SHALLOW_MODEL,DEPTH) calculates
%   one-dimensional shallow-strain vertical strain-rate model defined by
%   the two-element SHALLOW_MODEL. This model requires a scalar
%   steady-state accumulation rate (first element) and uniform vertical
%   strain rate (second element).
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, inverse years (1/a) for strain rate
%   and meters (m) for DEPTH.
%
% Joe MacGregor (UTIG), Mark Fahnestock (UAF)
% Last updated: 10/17/13

if (nargin ~= 2)
    error('shallow_strain:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 2.'])
end
if (~isnumeric(shallow_model) || (numel(shallow_model) ~= 2))
    error('shallow_strain:model', 'SHALLOW_MODEL is not a two-element numeric vector.')
end
if ~isnumeric(depth)
    error('shallow_strain:depth', 'DEPTH is not numeric.')
end

age                         = (1 / shallow_model(2)) .* log(shallow_model(1) ./ (shallow_model(1) - (shallow_model(2) .* depth)));
