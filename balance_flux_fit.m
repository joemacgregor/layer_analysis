function [res, flux_in, flux_out] ...
                            = balance_flux_fit(speed_bal_fg, width_fg, depth_iso_fg, shape_fg, speed_mod_fg, age_iso_fg, x_fb, y_fb, dist_fb, speed_fb, x_grd, y_grd, accum)
% BALANCE_FLUX_FIT Residual between input and output balance ice fluxes.
%
%   [RES,FLUX_IN,FLUX_OUT] = BALANCE_FLUX_FIT(SPEED_BAL_FG,WIDTH_FG,DEPTH_ISO_FG,SHAPE_FG,SPEED_MOD_FG,AGE_FG,X_FB,Y_FB,DIST_FB,SPEED_FB,X_GRD,Y_GRD,ACCUM)
%   calculates the residual RES between the input and output balance fluxes
%   (FLUX_IN and FLUX_OUT, respectively) for a given balance velocity
%   SPEED_BAL_FG through a fluxgate of width WIDTH_FB and depth
%   DEPTH_ISO_FG whose modern surface speed is SPEED_MOD_FG. The age of the
%   ice column at DEPTH_ISO_FG is AGE_ISO_FG and the shape factor for the
%   ice column to DEPTH_FG is SHAPE_FG. The x/y positions of the edges of
%   the upstream flowband and distances along it are given by X_FB, Y_FB
%   and DIST_FB respectively, which are all 1x2 cells. The modern surface
%   speed along these flowband edges is SPEED_FB (also a 1x2 cell). The
%   accumulation-rate grid ACCUM is necessary to calculate the input flux,
%   and its x/y positions are X_GRD and Y_GRD, respectively.
%
%   Assumed units are m/a for SPEED_BAL_FG, SPEED_MOD_FG, SPEED_FB and
%   ACCUM, km for X_GRD and Y_GRD, m for WIDTH_FG and DEPTH_FG, yr for
%   AGE_ISO_FG and dimensionless for SHAPE_FG. Outputs will be in Gt/yr.
%   
% Joe MacGregor (UTIG)
% Last updated: 11/14/13

if (nargin ~= 13)
    error('balance_flux_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 13.'])
end
if (~isscalar(speed_bal_fg) || ~isnumeric(speed_bal_fg))
    error('balance_flux_fit:speedbalfg', 'SPEED_BAL_FG is not a numeric scalar.')
end
if (~isscalar(width_fg) || ~isnumeric(width_fg))
    error('balance_flux_fit:widthfg', 'WIDTH_FG is not a numeric scalar.')
end
if (~isscalar(depth_iso_fg) || ~isnumeric(depth_iso_fg))
    error('balance_flux_fit:depthfg', 'DEPTH_ISO_FG is not a numeric scalar.')
end
if (~isscalar(shape_fg) || ~isnumeric(shape_fg))
    error('balance_flux_fit:shapefg', 'SHAPE_FG is not a numeric scalar.')
end
if (~isscalar(speed_mod_fg) || ~isnumeric(speed_mod_fg))
    error('balance_flux_fit:speedmodfg', 'SPEED_MOD_FG is not a numeric scalar.')
end
if (~isscalar(age_iso_fg) || ~isnumeric(age_iso_fg))
    error('balance_flux_fit:ageisofg', 'AGE_ISO_FG is not a numeric scalar.')
end
if (~iscell(x_fb) || (length(x_fb) ~= 2))
    error('balance_flux_fit:xfb', 'X_FB is not a two-element cell.')
end
if (~iscell(y_fb) || (length(y_fb) ~= 2))
    error('balance_flux_fit:yfb', 'Y_FB is not a two-element cell.')
end
if (~iscell(dist_fb) || (length(dist_fb) ~= 2))
    error('balance_flux_fit:distfb', 'DIST_FB is not a two-element cell.')
end
if (~iscell(speed_fb) || (length(speed_fb) ~= 2))
    error('balance_flux_fit:speedfb', 'SPEED_FB is not a two-element cell.')
end
if ~isnumeric(x_grd)
    error('balance_flux_fit:xgrd', 'X_GRD is not numeric.')
end
if ~isnumeric(y_grd)
    error('balance_flux_fit:ygrd', 'Y_GRD is not numeric.')
end
if ~isnumeric(accum)
    error('balance_flux_fit:accum', 'ACCUM is not numeric.')
end
if ~isequal(size(x_grd), size(y_grd), size(accum))
    error('balance_flux_fit:grd', 'X_GRD, Y_GRD and ACCUM are not the same size.')
end
if (nargout > 3)
    error('balance_flux_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 3.'])
end

% flowband limits based on fluxgate-adjusted velocity
for ii = 1:2
    speed_fb_adj            = (speed_bal_fg / (speed_mod_fg * shape_fg)) .* speed_fb{ii}; % adjusted speed along flowband edge (fluxgate balance speed corrected to surface by shape factor and compared to modern speed)
    age_fb                  = cumsum([0 (diff(1e3 .* dist_fb{ii}) ./ speed_fb_adj)]); % age along flowband edge, yr
    ind_age_max             = find((age_fb <= age_iso_fg), 1, 'last'); % index of last point on flowband edge with age less than isochrone age
    [x_fb{ii}, y_fb{ii}]    = deal([x_fb{ii}(1:ind_age_max) interp1(age_fb, x_fb{ii}, age_iso_fg, 'linear', 'extrap')], ...
                                   [y_fb{ii}(1:ind_age_max) interp1(age_fb, y_fb{ii}, age_iso_fg, 'linear', 'extrap')]); % trim flowband
end

% flowband outline
[x_fb_ord, y_fb_ord]        = deal([x_fb{1} fliplr(x_fb{2})], [y_fb{1} fliplr(y_fb{2})]);

% fluxes, Gt/yr
flux_in                     = 1e-12 * 917 * 1e6 * polyarea(x_fb_ord, y_fb_ord) * nanmean(accum(inpolygon(x_grd, y_grd, x_fb_ord, y_fb_ord))); % input flux
flux_out                    = 1e-12 * 917 * width_fg * depth_iso_fg * speed_bal_fg; % output flux
res                         = sqrt((flux_out - flux_in) ^ 2); % residual