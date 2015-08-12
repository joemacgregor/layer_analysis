function chisq              = dj_lat_fit(dj_lat_model, thick, depth, age, age_uncert, lat_spread_rate)
% DJ_LAT_FIT Residual between radar-observed and Dansgaard-Johnsen lateral-spreading-corrected depth-age scales.
%   
%   CHISQ = DJ_LAT_FIT(DJ_LAT_MODEL,THICK,DEPTH,AGE,AGE_UNCERT,LAT_SPREAD_RATE)
%   calculates the Chi-squared residual CHISQ between the ages AGE (and age
%   uncertainties AGE_UNCERT) of radar-observed layers at depths DEPTH with
%   those ages calculated by a one-dimensional Dansgaard-Johnsen vertical
%   strain-rate model that is corrected for the lateral spreading rate
%   LAT_SPREAD_RATE. This model requires a scalar ice thickness THICK and a
%   two-element vector DJ_MODEL. The first element of DJ_LAT_MODEL is the
%   steady-state accumulation rate and the second element is the thickness
%   of the basal shear layer, aka "little h" or the height of the shear
%   layer above the bed.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK, DEPTH,
%   shear-layer thickness, years (a) for AGE and AGE_UNCERT and inverse
%   years (1/a) for lateral spreading rate.
%   
% Joe MacGregor (UTIG)
% Last updated: 08/11/15

if (nargin ~= 6)
    error('dj_lat_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 6.'])
end
if (~isvector(dj_lat_model) || ~isnumeric(dj_lat_model) || (numel(dj_lat_model) ~= 2))
    error('dj_lat_fit:dj_lat_model', 'DJ_LAT_MODEL is not a numeric two-element vector.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj_lat_fit:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('dj_lat_fit:depth', 'DEPTH is not numeric.')
end
if ~isnumeric(age)
    error('dj_lat_fit:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('dj_lat_fit:ageuncert', 'AGE_UNCERT is not numeric.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('dj_lat_fit:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (~isscalar(lat_spread_rate) || ~isnumeric(lat_spread_rate))
    error('dj_lat_fit:latspreadrate', 'LAT_SPREAD_RATE is not a numeric scalar.')
end
if (nargout > 1)
    error('dj_lat_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

k                           = (2 * (dj_lat_model(1) - (lat_spread_rate * thick))) / ((2 * thick) - dj_lat_model(2));
depth_0                     = (k * dj_lat_model(2)) / (2 * (k + lat_spread_rate));
age_h                       = ((depth_0 - thick) / dj_lat_model(1)) * log((dj_lat_model(2) - depth_0) / (thick - depth_0));
a                           = -k / (2 * dj_lat_model(2));
b                           = -lat_spread_rate;
age_dj_lat                  = zeros(length(depth), 1);
age_dj_lat(depth <= (thick - dj_lat_model(2))) ...
                            = ((depth_0 - thick) / dj_lat_model(1)) .* log(((thick - depth(depth <= (thick - dj_lat_model(2)))) - depth_0) ./ (thick - depth_0));
if (b == 0)
    age_dj_lat(depth > (thick - dj_lat_model(2))) ...
                            = age_h + ((1 / a) .* ((1 / dj_lat_model(2)) - (1 ./ (thick - depth(depth > (thick - dj_lat_model(2)))))));   
else
    age_dj_lat(depth > (thick - dj_lat_model(2))) ...
                            = age_h + ((2 / b) .* (atanh((((2 * a) .* dj_lat_model(2)) + b) / b) - atanh((((2 * a) .* (thick - depth(depth > (thick - dj_lat_model(2))))) + b) ./ b)));
end
chisq                       = sum(((age - age_dj_lat) ./ age_uncert) .^ 2);