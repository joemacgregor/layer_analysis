function age                = nye_melt(nye_melt_model, thick, depth)
% NYE_MELT Nye melt-corrected depth-age scale.
%   
%   AGE = NYE_MELT(NYE_MELT_MODEL,THICK,DEPTH) calculates a one-dimensional
%   Nye vertical strain-rate model that allows for basal melting. This
%   model requires a scalar ice thickness THICK and a two-element vector
%   NYE_MELT_MODEL. The first element of NYE_MELT_MODEL is the steady-state
%   accumulation rate and the second element is the basal melt rate.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for surface-accumulation and basal melt rates, meters (m)
%   for THICK and DEPTH.
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
% Last updated: 10/17/13

if (nargin ~= 3)
    error('nye_melt:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 3.'])
end
if ((length(nye_melt_model) ~= 2) || ~isnumeric(nye_melt_model))
    error('nye_melt:nyemeltmodel', 'NYE_MELT_MODEL is not a two-element vector.')
end
if ((length(thick) ~= 1) || ~isnumeric(thick))
    error('nye_melt:thick', 'THICK is not a scalar.')
end
if ~isnumeric(depth)
    error('nye_melt:depth', 'DEPTH is not numeric.')
end

age                         = -(thick / (nye_melt_model(1) - nye_melt_model(2))) .* log((((1 - (nye_melt_model(2) / nye_melt_model(1))) ./ thick) .* (thick - depth)) + (nye_melt_model(2) / nye_melt_model(1)));