function chisq              = dj_surf_fit(thick_shear, thick, depth, age, age_uncert, lon_strain_rate, lat_strain_rate)
% DJ_SURF_FIT Residual between radar-observed and Dansgaard-Johnsen surface-strain-corrected depth-age scales.
%   
%   CHISQ = DJ_SURF_FIT(THICK_SHEAR,THICK,DEPTH,AGE,AGE_UNCERT,LON_STRAIN_RATE,LAT_STRAIN_RATE)
%   calculates the Chi-squared residual CHISQ between the ages AGE (and age
%   uncertainties AGE_UNCERT) of radar-observed layers at depths DEPTH with
%   those ages calculated by a one-dimensional Dansgaard-Johnsen ice-flow
%   model that is corrected for longitudinal and lateral strain rates
%   (LON_STRAIN_RATE and LAT_STRAIN_RATE, respectively). This model
%   requires a scalar ice thickness THICK and THICK_SHEAR, the thickness of
%   the basal shear layer, aka "little h" or the height of the shear layer
%   above the bed.
%   
%   Values must have self-consistent units, e.g.,meters (m) for
%   THICK_SHEAR, THICK and DEPTH, years (a) for AGE and AGE_UNCERT and
%   inverse years (1/a) for horizontal strain rates.
%   
% Joe MacGregor (UTIG)
% Last updated: 08/20/15

if (nargin ~= 7)
    error('dj_surf_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 7.'])
end
if (~isscalar(thick_shear) || ~isnumeric(thick_shear))
    error('dj_surf_fit:thickshear', 'THICK_SHEAR is not a numeric scalar.')
end
if (~isscalar(thick) || ~isnumeric(thick))
    error('dj_surf_fit:thick', 'THICK is not a numeric scalar.')
end
if ~isnumeric(depth)
    error('dj_surf_fit:depth', 'DEPTH is not numeric.')
end
if ~isnumeric(age)
    error('dj_surf_fit:age', 'AGE is not numeric.')
end
if ~isnumeric(age_uncert)
    error('dj_surf_fit:ageuncert', 'AGE_UNCERT is not numeric.')
end
if ~isequal(size(depth), size(age), size(age_uncert))
    error('dj_surf_fit:depthagematch', 'DEPTH/AGE/AGE_UNCERT sizes are not equal.')
end
if (~isscalar(lon_strain_rate) || ~isnumeric(lon_strain_rate))
    error('dj_surf_fit:lonstrainrate', 'LON_STRAIN_RATE is not a numeric scalar.')
end
if (~isscalar(lat_strain_rate) || ~isnumeric(lat_strain_rate))
    error('dj_surf_fit:latstrainrate', 'LAT_STRAIN_RATE is not a numeric scalar.')
end
if (nargout > 1)
    error('dj_surf_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 1.'])
end

accum                       = (lon_strain_rate + lat_strain_rate) * (((2 * thick) - thick_shear) / 2);
depth_0                     = thick_shear / 2;
age_h                       = ((depth_0 - thick) / accum) * log((thick_shear - depth_0) / (thick - depth_0));
age_dj_surf                 = zeros(length(depth), 1);
age_dj_surf(depth <= (thick - thick_shear)) ...
                            = ((depth_0 - thick) / accum) .* log(((thick - depth(depth <= (thick - thick_shear))) - depth_0) ./ (thick - depth_0));
age_dj_surf(depth > (thick - thick_shear)) ...
                            = age_h + ((2 / (lon_strain_rate + lat_strain_rate)) .* ((thick_shear ./ (thick - depth(depth > (thick - thick_shear)))) - 1));
chisq                       = sum(((age - age_dj_surf) ./ age_uncert) .^ 2);