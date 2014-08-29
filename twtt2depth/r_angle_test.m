function residual           = r_angle_test(r_angle_start, ind_refract, half_antenna_sepn, thick_density)
% R_ANGLE_TEST Residual between modeled horizontal distance traveled by refracted radio wave and antenna half-separation. 
% 
% RESIDUAL = R_ANGLE_TEST(R_ANGLE_START,IND_REFRACT,HALF_ANTENNA_SEPN,THICK_DENSITY)
% calculates the horizontal distance traveled by the refracted raypath
% through the firn and ice, and how close that horizontal distance is to
% half of the antenna separation (HALF_ANTENNA_SEPN), as measured by the
% difference between the modeled and measured antenna half-separations
% (RESIDUAL, m). R_ANGLE_START is the starting incidence angle (rad) and
% IND_REFRACT is the depth profile of the index of refraction.
% 
% Joe MacGregor (UTIG, joemac@ig.utexas.edu)
% Last updated: 08/05/14

if (nargin ~= 4)
    error('r_angle_test:nargin', 'Number of arguments must be equal to 4.')
end
if (~isnumeric(r_angle_start) || ~isscalar(r_angle_start) || (r_angle_start < 0))
    error('r_angle_test:r_angle_start', 'R_ANGLE_START is not a numeric positive scalar.')
end
if (~isnumeric(ind_refract) || ~iscolumn(ind_refract) || any(ind_refract < 1))
    error('r_angle_test:ind_refract', 'IND_REFRACT is not a numeric positive column vector.')
end
if (~isnumeric(half_antenna_sepn) || ~isscalar(half_antenna_sepn) || (half_antenna_sepn < 0))
    error('r_angle_test:half_antenna_sepn', 'ANTENNA_SEPN is not a numeric positive scalar.')
end
if (~isnumeric(thick_density) || ~iscolumn(thick_density) || any(thick_density < 0))
    error('r_angle_test:thick_density', 'THICK_DENSITY is not a numeric positive column vector.')
end
if ~isequal(length(ind_refract), length(thick_density))
    error('r_angle_test:density_length', 'IND_REFRACT and THICK_DENSITY are not equal lengths.')
end
if (nargout > 1)
    error('r_angle_test:nargout', 'Number of outputs cannot be more than 1.')
end

r_angle                     = zeros(length(thick_density), 1);
r_angle(1)                  = r_angle_start; % angle intially and gets iterated
for ii = 2:length(thick_density)
    r_angle(ii)             = asin((ind_refract(ii - 1) ./ ind_refract(ii)) .* sin(r_angle(ii - 1)));
end
residual                    = sum(thick_density .* tan(r_angle)) - half_antenna_sepn; % residual between what was traveled and what should have been traveled