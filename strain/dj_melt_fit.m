function chisq              = dj_melt_fit(dj_melt_model, thick, depth, age, age_uncert)
% DJ_MELT_FIT Residual between radar-observed and Dansgaard-Johnsen + basal melting depth-age scales.
%   
%   CHISQ = DJ_MELT_FIT(DJ_MELT_MODEL,THICK,DEPTH,AGE,AGE_UNCERT)
%   calculates the chi-squared residual CHISQ between the ages AGE (and age
%   uncertainties AGE_UNCERT) of radar-observed layers at depths DEPTH with
%   those ages calculated by a one-dimensional Dansgaard-Johnsen vertical
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
%   DEPTH and shear-layer thickness, and years (a) for AGE and AGE_UNCERT.
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
% Last updated: 08/10/15

if (nargin ~= 5)
    error('dj_melt_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 5.'])
end
if (~isnumeric(dj_melt_model) || ~isvector(dj_melt_model) || (numel(dj_melt_model) ~= 4))
    error('dj_melt_fit:dj_model', 'DJ_MELT_MODEL is not a four-element vector.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj_melt_fit:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('dj_melt_fit:depth', 'DEPTH is not numeric.') 
end
if ~isnumeric(age)
    error('dj_melt_fit:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('dj_melt_fit:ageuncert', 'AGE_UNCERT is not numeric.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('dj_melt_fit:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (nargout > 1)
    error('dj_melt_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

k                           = (dj_melt_model(1) - dj_melt_model(3)) / (((dj_melt_model(2) * (dj_melt_model(4) - 1)) / 2) + thick); % factor rr in Dahl-Jensen et al., 2003, Equation 4, (bdot-mdot)/(((h(s-1))/2)+H)
depth_0                     = -((dj_melt_model(3) / k) + ((dj_melt_model(2) * (dj_melt_model(4) - 1)) / 2)); % intercept depth 
age_dj_melt                 = zeros(length(depth), 1);
age_dj_melt(depth <= (thick - dj_melt_model(2))) ...
                            = ((depth_0 - thick) / dj_melt_model(1)) .* log(((thick - depth(depth <= (thick - dj_melt_model(2)))) - depth_0) ./ (thick - depth_0));
age_h                       = ((depth_0 - thick) / dj_melt_model(1)) * log((dj_melt_model(2) - depth_0) / (thick - depth_0));
a                           = k * ((dj_melt_model(4) - 1) / (2 * dj_melt_model(2)));
b                           = -k * dj_melt_model(4);
c                           = -dj_melt_model(3);
quad                        = (4 * a * c) - (b ^ 2);
if (quad > 0)
    quad_root               = sqrt(quad);
    age_dj_melt(depth > (thick - dj_melt_model(2))) ...
                            = age_h + ((2 / quad_root) .* (atan((((2 * a) .* (thick - depth(depth > (thick - dj_melt_model(2))))) + b) ./ quad_root) - atan((((2 * a) * dj_melt_model(2)) + b) ./ quad_root)));
elseif (quad < 0)
    quad_root               = sqrt((b ^ 2) - (4 * a * c));
    age_dj_melt(depth > (thick - dj_melt_model(2))) ...
                            = age_h + ((2 / quad_root) .* (atanh((((2 * a) * dj_melt_model(2)) + b) ./ quad_root) - atanh((((2 * a) .* (thick - depth(depth > (thick - dj_melt_model(2))))) + b) ./ quad_root)));    
elseif (quad == 0)
    age_dj_melt(depth > (thick - dj_melt_model(2))) ...
                            = age_h + ((2 / (((2 * a) * dj_melt_model(2)) + b)) - (2 / (((2 * a) * (thick - depth(depth > (thick - dj_melt_model(2))))) + b)));
end
chisq                       = sum(((age - age_dj_melt) ./ age_uncert) .^ 2);