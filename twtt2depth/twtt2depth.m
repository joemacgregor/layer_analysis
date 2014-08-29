function [depth, firn_corr, depth_lin, refract_gain, raypath] ...
                            = twtt2depth(twtt, depth_density, thick_density, density, density_ice, permitt_ice, antenna_sepn, antenna_height, do_trigger)
% TWTT2DEPTH Conversion from two-way traveltime to depth for ice-penetrating radar data.
% 
% [DEPTH,FIRN_CORR,DEPTH_LIN,REFRACT_GAIN,RAYPATH] = TWTT2DEPTH(TWTT,DEPTH_DENSITY,THICK_LAYER,DENSITY,DENSITY_ICE,PERMITT_ICE,ANTENNA_SEPN,ANTENNA_HEIGHT,DO_TRIGGER)
% converts two-way traveltime (aka "fast time") (TWTT; s) to depth (DEPTH;
% m) given a depth profile of firn/ice density. DEPTH_DENSITY (m) is the
% depth of the intervals (not midpoints) between density measurements or or
% model, THICK_DENSITY (m) is the thickness of those intervals and DENSITY
% is their density (g/cm^3). DENSITY_ICE is the pure ice density (e.g.,
% 0.917 g/cm^3) and PERMITT_ICE is the real part of the complex relative
% permittivity of ice (e.g., 3.2). ANTENNA_SEPN (m) is the scalar horizontal
% separation between the transmitting (Tx) and receiving (Rx) antennas, and
% ANTENNA_HEIGHT (m) is their scalar height above the subaerial ice-sheet
% surface. DO_TRIGGER is a logical scalar that should be set to true if
% TWTT=0 when the receiver is triggered, rather than when the radar signal
% is transmitted.
% 
% TWTT should be a column vector, and TWTTDEPTH will return DEPTH as a
% column vector of equal length to TWTT. THICK_DENSITY and DENSITY should
% be column vectors of the same length, whereas DEPTH_DENSITY should be a
% column vector whose length is one element greater than THICK_DENSITY or
% DENSITY.
% 
% When ANTENNA_SEPN~=0, TWTT2DEPTH accounts for the curved
% raypath geometry due to refraction through the firn.
% 
% When ANTENNA_HEIGHT~=0, negative depths are heights in air above the
% surface.
% 
% FIRN_CORR (m) is the firn correction that represents the depth difference
% between a purely linear relationship between traveltime and depth at
% TWTT(end), as compared to that calculated by TWTT2DEPTH when accounting
% for faster radio-wave velocities in the lower-density firn. DEPTH_LIN (m)
% is the firn-corrected linear relationship between TWTT and depth.
% REFRACT_GAIN is the refraction gain (dB) as a function of depth due to
% refraction within the firn. RAYPATH (m) is the length of the raypath
% within each density layer.
% 
% See also: LOOYENGA, R_ANGLE_TEST
% 
% Joe MacGregor (UTIG, joemac@ig.utexas.edu)
% Last updated: 08/05/14

if ~exist('looyenga.m', 'file')
    error('twtt2depth:looyenga', 'Necessary function LOOYENGA is not available within this user''s path.')
end
if ~exist('r_angle_test.m', 'file')
    error('twtt2depth:r_angle_test', 'Necessary function R_ANGLE_TEST is not available within this user''s path.')
end
if (nargin ~= 9)
    error('twtt2depth:nargin', 'Number of arguments must be equal to 9.')
end
if (~isnumeric(twtt) || ~iscolumn(twtt))
    error('twtt2depth:twtt', 'TWTT is not a numeric column vector.')
end
if (~isnumeric(depth_density) || ~iscolumn(depth_density))
    error('twtt2depth:depth_density', 'DEPTH_DENSITY is not a numeric column vector.')
end
if (~isnumeric(thick_density) || ~iscolumn(thick_density) || any(thick_density < 0))
    error('twtt2depth:thick_density', 'THICK_DENSITY is not a numeric positive column vector.')
end
if (~isnumeric(density) || ~iscolumn(density) || any(density < 0))
    error('twtt2depth:density', 'DENSITY is not a numeric positive column vector.')
end
if ~isequal((length(depth_density) - 1), length(thick_density), length(density))
    error('twtt2depth:density_length', 'DEPTH_DENSITY/THICK_DENSITY/DENSITY not the correct lengths.')
end
if (~isnumeric(density_ice) || ~isscalar(density_ice) || (density_ice < 0))
    error('twtt2depth:density_ice', 'DENSITY_ICE is not a numeric positive scalar.')
end
if (~isnumeric(permitt_ice) || ~isscalar(permitt_ice) || (permitt_ice < 0))
    error('twtt2depth:permitt_ice', 'PERMITT_ICE is not a numeric positive scalar.')
end
if (~isnumeric(antenna_sepn) || ~isscalar(antenna_sepn) || (antenna_sepn < 0))
    error('twtt2depth:antenna_sepn', 'ANTENNA_SEPN is not a numeric positive scalar.')
end
if (~isnumeric(antenna_height) || ~isscalar(antenna_height))
    error('twtt2depth:antenna_height', 'ANTENNA_HEIGHT is not a numeric scalar.')
end
if (~isscalar(do_trigger) || ~islogical(do_trigger))
    error('twtt2depth:do_trigger', 'DO_TRIGGER is not a logical scalar.')
end
if (nargout > 5)
    error('twtt2depth:nargout', 'Number of outputs cannot be more than 5.')
end

speed_vacuum                = 299792458; % speed of light in the vacuum in m/s, c
ind_refract                 = sqrt(looyenga((density ./ density_ice), permitt_ice(ones(length(density))), ones(length(density), 1))); % index of refraction (sqrt of permittivity) in each layer

% account for "air time" if antenna height is non-zero
if antenna_height 
    thick_density           = [antenna_height; thick_density];
    depth_density           = [-antenna_height; depth_density];
    ind_refract             = [1; ind_refract];
end

% account for increased traveltime due to recording triggered by arrival of airwave
if do_trigger
    twtt                    = twtt + (antenna_sepn / speed_vacuum); % twtt=0 is when Rx is triggered, so L/c is the twtt from Tx to Rx
end

speed                       = speed_vacuum ./ ind_refract; % radio-wave velocity within each layer
num_layer                   = length(thick_density); % number of layers with a uniform density / index of refraction
r_angle                     = zeros(num_layer); % square matrix because raypath incidence angles are different for each layer

% only proceed if refraction angle is non-zero, i.e., antenna are separated
if antenna_sepn
    % starting refraction angles based on geometry
    if antenna_height
        r_angle(1, 1)       = atan(antenna_sepn / (2 * antenna_height));
    else
        r_angle(1, 1)       = atan(antenna_sepn / (2 * depth_density(2)));
    end;
    for ii = 2:num_layer
        r_angle(1, ii)      = fminsearch(@r_angle_test, r_angle(1, (ii - 1)), [], ind_refract(1:ii), (antenna_sepn / 2), thick_density(1:ii)); % start angle for each layer based on previous layer's angle
    end
    % refraction-angle profiles
    for ii = 2:num_layer
        for jj = 2:ii
            r_angle(jj, ii) = asin((ind_refract(jj - 1) / ind_refract(jj)) * sin(r_angle((jj - 1), ii))); % r_angle-depth profiles (each column is for each depth)
        end
    end
end

% length of raypath through each layer and its equivalent twtt
length_ray                  = zeros(num_layer);
for ii = 1:num_layer
    length_ray(1:ii, ii)    = thick_density(1:ii) ./ cos(r_angle(1:ii, ii)); % diagonal raypath length through each depth increment
end
twtt_layer                  = 2 .* (length_ray ./ speed(:, ones(1, num_layer))); % twtt through each layer

% modeled twtt profile based on density
twtt_model                  = zeros((num_layer + 1), 1);
for ii = 2:(num_layer + 1)
    twtt_model(ii)          = sum(twtt_layer(:, (ii - 1))); % sum the twtt through each increment (total in column)
end

% best depth calculation based on interpolation between modeled twtt, modeled depth, and measured twtt
depth                       = interp1(twtt_model, depth_density, twtt, 'linear', 'extrap');

% bonus stuff
if antenna_height
    depth_lin               = interp1(twtt_model([1:2 end]), depth_density([1:2 end]), twtt, 'linear', 'extrap');
else
    depth_lin               = (speed(end) / 2) .* twtt; % simple linear twtt-depth profile using maximum density in profile
end
firn_corr                   = depth(end) - depth_lin(end); % firn correction based on difference between linear depth profile and the more complex one
depth_lin                   = depth_lin + firn_corr; % firn-corrected linear depth vector
refract_gain                = 20 .* log10((2 .* depth_density(2:end) .* tan(r_angle(1, :))') ./ antenna_sepn); % refraction gain vs depth
if (length(depth_density) > 2)
    refract_gain            = interp1(depth_density(2:end), refract_gain, depth, 'linear', 'extrap'); % interpolate to modeled depths
    raypath                     = interp1(depth_density(2:end), sum(length_ray)', depth, 'linear', 'extrap'); % length of raypath through each layer
end