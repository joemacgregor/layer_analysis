function age                = dj_melt(dj_melt_model, thick, depth)
% DJ_MELT Dansgaard-Johnsen + basal melting depth-age scale.
%   
%   AGE = DJ_MELT_FIT(DJ_MELT_MODEL,THICK,DEPTH)
%   calculates a one-dimensional Dansgaard-Johnsen vertical
%   strain-rate model with basal melting. This model requires a scalar ice
%   thickness THICK and a four-element vector DJ_MELT_MODEL. The first
%   element of DJ_MELT_MODEL is the steady-state accumulation rate, the
%   second element is the thickness of the basal shear layer, aka "little
%   h" or the height of the shear layer above the bed, the third element is
%   the basal melt rate and the fourth element is the sliding fraction,
%   i.e., the ratio of the basal sliding velocity to the surface velocity.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation and basal melt rates, meters (m) for THICK,
%   DEPTH and shear-layer thickness.
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
%   Dahl-Jensen, D. et al., 2003, Basal melt at NorthGRIP modeled from
%   borehole, ice-core and radio-echo sounder observations, Ann. Glaciol.,
%   37, 207-212.
%
% Joe MacGregor (UTIG), Mark Fahnestock (UAF)
% Last updated: 08/11/15

if (nargin ~= 3)
    error('dj_melt:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 3.'])
end
if (~isvector(dj_melt_model) || ~isnumeric(dj_melt_model) || (length(dj_melt_model) ~= 4))
    error('dj_melt:dj_melt_model', 'DJ_MELT_MODEL is not a four-element numeric vector.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj_melt:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('dj_melt:depth', 'DEPTH is not numeric.')
end

k                           = (dj_melt_model(1) - dj_melt_model(3)) / (((dj_melt_model(2) * (dj_melt_model(4) - 1)) / 2) + thick); % factor rr in Dahl-Jensen et al., 2003, Equation 4, (bdot-mdot)/(((h(s-1))/2)+H)
depth_0                     = -((dj_melt_model(3) / k) + ((dj_melt_model(2) * (dj_melt_model(4) - 1)) / 2)); % intercept depth 
age                         = zeros(length(depth), 1);
age(depth <= (thick - dj_melt_model(2))) ...
                            = ((depth_0 - thick) / dj_melt_model(1)) .* log(((thick - depth(depth <= (thick - dj_melt_model(2)))) - depth_0) ./ (thick - depth_0));
age_h                       = ((depth_0 - thick) / dj_melt_model(1)) * log((dj_melt_model(2) - depth_0) / (thick - depth_0));
a                           = k * ((dj_melt_model(4) - 1) / (2 * dj_melt_model(2)));
b                           = -k * dj_melt_model(4);
c                           = -dj_melt_model(3);
quad                        = (4 * a * c) - (b ^ 2);
if (quad > 0)
    quad_root               = sqrt(quad);
    age(depth > (thick - dj_melt_model(2))) ...
                            = age_h + ((2 / quad_root) .* (atan((((2 * a) .* (thick - depth(depth > (thick - dj_melt_model(2))))) + b) ./ quad_root) - atan((((2 * a) * dj_melt_model(2)) + b) ./ quad_root)));
elseif (quad < 0)
    quad_root               = sqrt((b ^ 2) - (4 * a * c));
    age(depth > (thick - dj_melt_model(2))) ...
                            = age_h + ((2 / quad_root) .* (atanh((((2 * a) * dj_melt_model(2)) + b) ./ quad_root) - atanh((((2 * a) .* (thick - depth(depth > (thick - dj_melt_model(2))))) + b) ./ quad_root)));    
elseif (quad == 0)
    age(depth > (thick - dj_melt_model(2))) ...
                            = age_h + ((2 / (((2 * a) * dj_melt_model(2)) + b)) - (2 / (((2 * a) * (thick - depth(depth > (thick - dj_melt_model(2))))) + b)));
end