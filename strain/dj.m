function age                = dj(dj_model, thick, depth)
% DJ Dansgaard-Johnsen depth-age scale.
%   
%   AGE = DJ(DJ_MODEL,THICK,DEPTH) calculates a one-dimensional
%   Dansgaard-Johnsen vertical strain-rate model. This model requires a
%   scalar ice thickness THICK and a two-element vector DJ_MODEL. The first
%   element of DJ_MODEL is the steady-state accumulation rate and the
%   second element is the thickness of the basal shear layer, aka "little
%   h" or the height of the shear layer above the bed.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK, DEPTH
%   and shear-layer thickness.
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
% Last updated: 08/11/15

if (nargin ~= 3)
    error('dj:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 3.'])
end
if (~isvector(dj_model) || ~isnumeric(dj_model) || (length(dj_model) ~= 2))
    error('dj:dj_model', 'DJ_MODEL is not a two-element numeric vector.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj:thick', 'THICK is not a scalar.')
end
if ~isnumeric(depth)
    error('dj:depth', 'DEPTH is not numeric.')
end

age                         = zeros(length(depth), 1);
age(depth <= (thick - dj_model(2))) ...
                            = (((2 * thick) - dj_model(2)) / (2 * dj_model(1))) .* log(((2 * thick) - dj_model(2)) ./ ((2 .* (thick - depth(depth <= (thick - dj_model(2))))) - dj_model(2)));
age_h                       = (((2 * thick) - dj_model(2)) / (2 * dj_model(1))) .* log(((2 * thick) - dj_model(2)) ./ dj_model(2));
age(depth > (thick - dj_model(2))) ...
                            = ((((2 * thick) - dj_model(2)) / dj_model(1)) .* ((dj_model(2) ./ (thick - depth(depth > (thick - dj_model(2))))) - 1)) + age_h;