function chisq              = dj_lat_fit(dj_lat_model, thick, depth, age, age_uncert, lon_strain_rate, lat_strain_rate)
% DJ_LAT_FIT Residual between radar-observed and Dansgaard-Johnsen lateral-strain-corrected depth-age scales.
%   
%   CHISQ = DJ_LAT_FIT(DJ_LAT_MODEL,THICK,DEPTH,AGE,AGE_UNCERT,LON_STRAIN_RATE,LAT_STRAIN_RATE)
%   calculates the Chi-squared residual CHISQ between the ages AGE (and age
%   uncertainties AGE_UNCERT) of radar-observed layers at depths DEPTH with
%   those ages calculated by a one-dimensional Dansgaard-Johnsen vertical
%   strain-rate model that is corrected for longitudinal and lateral
%   strain rates (LON_STRAIN_RATE and LAT_STRAIN_RATE, respectively).
%   This model requires a scalar ice thickness THICK and a two-element
%   vector DJ_MODEL. The first element of DJ_LAT_MODEL is the steady-state
%   accumulation rate and the second element is the thickness of the basal
%   shear layer, aka "little h" or the height of the shear layer above the
%   bed.
%   
%   Values must have self-consistent units, e.g., ice-equivalent meters per
%   year (m/a) for accumulation rate, meters (m) for THICK, DEPTH,
%   shear-layer thickness, years (a) for AGE and AGE_UNCERT and inverse
%   years (1/a) for horizontal strain rates.
%   
% Joe MacGregor (UTIG)
% Last updated: 08/20/15

if (nargin ~= 7)
    error('dj_lat_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 7.'])
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
if (~isscalar(lon_strain_rate) || ~isnumeric(lon_strain_rate))
    error('dj_lat_fit:lonstrainrate', 'LON_STRAIN_RATE is not a numeric scalar.')
end
if (~isscalar(lat_strain_rate) || ~isnumeric(lat_strain_rate))
    error('dj_lat_fit:latstrainrate', 'LAT_STRAIN_RATE is not a numeric scalar.')
end
if (nargout > 1)
    error('dj_lat_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

alpha                       = lat_strain_rate / lon_strain_rate;
k                           = (1 + alpha) * ((2 * (dj_lat_model(1))) / ((2 * thick) - dj_lat_model(2)));
depth_0                     = (dj_lat_model(2) / 2) / (1 + alpha);
age_h                       = ((depth_0 - thick) / dj_lat_model(1)) * log((dj_lat_model(2) - depth_0) / (thick - depth_0));
age_dj_lat                  = zeros(length(depth), 1);
age_dj_lat(depth <= (thick - dj_lat_model(2))) ...
                            = ((depth_0 - thick) / dj_lat_model(1)) .* log(((thick - depth(depth <= (thick - dj_lat_model(2)))) - depth_0) ./ (thick - depth_0));
age_dj_lat(depth > (thick - dj_lat_model(2))) ...
                            = age_h + ((2 / k) .* ((dj_lat_model(2) ./ (thick - depth(depth > (thick - dj_lat_model(2))))) - 1));
chisq                       = sum(((age - age_dj_lat) ./ age_uncert) .^ 2);