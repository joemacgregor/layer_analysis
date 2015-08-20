function age                = dj_surf(thick_shear, thick, depth, lon_strain_rate, lat_strain_rate)
% DJ_SURF Dansgaard-Johnsen surface-strain-corrected depth-age scale.
%   
%   AGE = DJ_SURF(THICK_SHEAR,THICK,DEPTH,LON_STRAIN_RATE,LAT_STRAIN_RATE)
%   calculates a one-dimensional Dansgaard-Johnsen ice-flow model that is
%   corrected for the surface longitudinal and lateral strain rates
%   (LON_STRAIN_RATE and LAT_STRAIN_RATE, respectively). This model
%   requires a scalar ice thickness THICK and THICK_SHEAR, which is the
%   thickness of the basal shear layer, aka "little h" or the height of the
%   shear layer above the bed.
%   
%   Values must have self-consistent units, e.g., meters (m) for THICK,
%   DEPTH and THICK_SHEAR, and inverse years (1/a) for horizontal strain
%   rates.
%
% Joe MacGregor (UTIG)
% Last updated: 08/20/15

if (nargin ~= 5)
    error('dj_surf:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 5.'])
end
if (~isscalar(thick_shear) || ~isnumeric(thick_shear))
    error('dj_surf:thickshear', 'THICK_SHEAR is not a numeric scalar.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj_surf:thick', 'THICK is not a scalar.')
end
if ~isnumeric(depth)
    error('dj_surf:depth', 'DEPTH is not numeric.')
end
if (~isscalar(lon_strain_rate) || ~isnumeric(lon_strain_rate))
    error('dj_surf:lonstrainrate', 'LON_STRAIN_RATE is not a numeric scalar.')
end
if (~isscalar(lat_strain_rate) || ~isnumeric(lat_strain_rate))
    error('dj_surf:latstrainrate', 'LAT_STRAIN_RATE is not a numeric scalar.')
end
if (nargout > 1)
    error('dj_surf:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

accum                       = (lon_strain_rate + lat_strain_rate) * (((2 * thick) - thick_shear) / 2);
depth_0                     = thick_shear / 2;
age_h                       = ((depth_0 - thick) / accum) * log((thick_shear - depth_0) / (thick - depth_0));
age                         = zeros(length(depth), 1);
age(depth <= (thick - thick_shear)) ...
                            = ((depth_0 - thick) / accum) .* log(((thick - depth(depth <= (thick - thick_shear))) - depth_0) ./ (thick - depth_0));
age(depth > (thick - thick_shear)) ...
                            = age_h + ((2 / (lon_strain_rate + lat_strain_rate)) .* ((thick_shear ./ (thick - depth(depth > (thick - thick_shear)))) - 1));