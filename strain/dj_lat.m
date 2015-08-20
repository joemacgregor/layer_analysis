function age                = dj_lat(dj_lat_model, thick, depth, lon_strain_rate, lat_strain_rate)
% DJ_LAT Dansgaard-Johnsen lateral-strain-corrected depth-age scale.
%   
%   AGE = DJ_LAT(DJ_LAT_MODEL,THICK,DEPTH,LON_STRAIN_RATE,LAT_STRAIN_RATE)
%   calculates a one-dimensional Dansgaard-Johnsen ice-flow model that is
%   corrected for the longitudinal and lateral strain rates
%   (LON_STRAIN_RATE and LAT_STRAIN_RATE, respectively). This model
%   requires a scalar ice thickness THICK and a two-element vector
%   DJ_LAT_MODEL. The first element of DJ_LAT_MODEL is the steady-state
%   accumulation rate and the second element is the thickness of the basal
%   shear layer, aka "little h" or the height of the shear layer above the
%   bed.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK, DEPTH and
%   shear-layer thickness, and inverse years (1/a) for horizontal strain
%   rates.
%
% Joe MacGregor (UTIG)
% Last updated: 08/20/15

if (nargin ~= 5)
    error('dj_lat:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 5.'])
end
if (~isvector(dj_lat_model) || ~isnumeric(dj_lat_model) || (length(dj_lat_model) ~= 2))
    error('dj_lat:dj_model', 'DJ_LAT_MODEL is not a numeric two-element vector.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj_lat:thick', 'THICK is not a scalar.')
end
if ~isnumeric(depth)
    error('dj_lat:depth', 'DEPTH is not numeric.')
end
if (~isscalar(lon_strain_rate) || ~isnumeric(lon_strain_rate))
    error('dj_lat:lonstrainrate', 'LON_STRAIN_RATE is not a numeric scalar.')
end
if (~isscalar(lat_strain_rate) || ~isnumeric(lat_strain_rate))
    error('dj_lat:latstrainrate', 'LAT_STRAIN_RATE is not a numeric scalar.')
end
if (nargout > 1)
    error('dj_lat:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

alpha                       = lat_strain_rate / lon_strain_rate;
k                           = (1 + alpha) * ((2 * (dj_lat_model(1))) / ((2 * thick) - dj_lat_model(2)));
depth_0                     = (dj_lat_model(2) / 2) / (1 + alpha);
age_h                       = ((depth_0 - thick) / dj_lat_model(1)) * log((dj_lat_model(2) - depth_0) / (thick - depth_0));
age                         = zeros(length(depth), 1);
age(depth <= (thick - dj_lat_model(2))) ...
                            = ((depth_0 - thick) / dj_lat_model(1)) .* log(((thick - depth(depth <= (thick - dj_lat_model(2)))) - depth_0) ./ (thick - depth_0));
age(depth > (thick - dj_lat_model(2))) ...
                            = age_h + ((2 / k) .* ((dj_lat_model(2) ./ (thick - depth(depth > (thick - dj_lat_model(2))))) - 1));