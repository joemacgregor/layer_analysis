function age                = dj_lat(dj_lat_model, thick, depth, lat_spread_rate)
% DJ_LAT Dansgaard-Johnsen lateral-spreading-corrected depth-age scale.
%   
%   AGE = DJ_LAT(DJ_LAT_MODEL,THICK,DEPTH,LAT_SPREAD_RATE) calculates a
%   one-dimensional Dansgaard-Johnsen vertical strain-rate model that is
%   corrected for the lateral spreading rate LAT_SPREAD_RATE. This model
%   requires a scalar ice thickness THICK and a two-element vector
%   DJ_LAT_MODEL. The first element of DJ_LAT_MODEL is the steady-state
%   accumulation rate and the second element is the thickness of the basal
%   shear layer, aka "little h" or the height of the shear layer above the
%   bed.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK, DEPTH and
%   shear-layer thickness, and inverse years (1/a) for lateral spreading
%   rate.
%
% Joe MacGregor (UTIG)
% Last updated: 08/11/15

if (nargin ~= 4)
    error('dj_lat:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 4.'])
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
if (~isscalar(lat_spread_rate) || ~isnumeric(lat_spread_rate))
    error('dj_lat:latspreadrate', 'LAT_SPREAD_RATE is not a numeric scalar.')
end
if (nargout > 1)
    error('dj_lat:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

k                           = (2 * (dj_lat_model(1) - (lat_spread_rate * thick))) / ((2 * thick) - dj_lat_model(2));
depth_0                     = (k * dj_lat_model(2)) / (2 * (k + lat_spread_rate));
age_h                       = ((depth_0 - thick) / dj_lat_model(1)) * log((dj_lat_model(2) - depth_0) / (thick - depth_0));
a                           = -k / (2 * dj_lat_model(2));
b                           = -lat_spread_rate;
age                         = zeros(length(depth), 1);
age(depth <= (thick - dj_lat_model(2))) ...
                            = ((depth_0 - thick) / dj_lat_model(1)) .* log(((thick - depth(depth <= (thick - dj_lat_model(2)))) - depth_0) ./ (thick - depth_0));
if (b == 0)
    age(depth > (thick - dj_lat_model(2))) ...
                            = age_h + ((1 / a) .* ((1 / dj_lat_model(2)) - (1 ./ (thick - depth(depth > (thick - dj_lat_model(2)))))));   
else
    age(depth > (thick - dj_lat_model(2))) ...
                            = age_h + ((2 / b) .* (atanh((((2 * a) .* dj_lat_model(2)) + b) / b) - atanh((((2 * a) .* (thick - depth(depth > (thick - dj_lat_model(2))))) + b) ./ b)));
end