function age                = nye(accum, thick, depth)
% NYE Nye depth-age scale.
%   
%   AGE = NYE(ACCUM,THICK,DEPTH) calculates a one-dimensional Nye vertical
%   strain-rate model. This model requires a scalar steady-state
%   accumulation rate ACCUM and ice thickness THICK.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK and DEPTH.
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
    error('nye:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 3.'])
end
if (~isscalar(accum) || ~isnumeric(accum))
    error('nye:accum', 'ACCUM is not a scalar.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('nye:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('nye:depth', 'DEPTH is not numeric.')
end

age                         = -(thick / accum) .* log((thick - depth) ./ thick);