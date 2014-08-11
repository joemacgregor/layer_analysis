function [paral, merid, x_paral, y_paral, x_merid, y_merid] ...
                            = graticule_greenland(x_lim, y_lim, inc_paral, inc_merid)
% GRATICULE_GREENLAND Construct graticule vectors from bounds of a polar stereographic image using the GIMP projection.
% 
%   [PARAL,MERID,X_PARAL,Y_PARAL,X_MERID,Y_MERID] = GRATICULE_GREENLAND(X_LIM,Y_LIM,INC_PARAL,INC_MERID,LAT_STD)
%   calculates a set of vectors in polar stereographic coordinates
%   (X_PARAL,Y_PARAL,X_MERID,Y_MERID) that are equivalent to a graticule
%   for a polar stereographic image along parallels PARAL and meridians
%   MERID. The projection is best suited to Greenland and is the same as
%   for the GIMP DEM. GRATICULE makes it easier to plot a graticule over an
%   existing image that is plotted using a regular polar stereographic grid
%   without requiring the Mapping Toolbox. X_LIM and Y_LIM are the x-axis
%   and y-axis limits (units: kilometers) of the polar stereographic image
%   (e.g., get(gca,'xlim') and get(gca,'ylim')). INC_PARAL and INC_MERID
%   are the desired intervals for the parallels (lines of equal latitude)
%   and meridians (lines of equal longitude), respectively (units: decimal
%   degrees).
% 
%   GRATICULE_GREENLAND requires the Mapping Toolbox.
%
% Joe MacGregor
% Last updated: 08/07/14

if (nargin ~= 4)
    error('graticule_greenland:nargin', 'Incorrect number of inputs (must be 5).')
end
if (~isnumeric(x_lim) || (length(x_lim) ~= 2))
    error('graticule_greenland:xlim', 'X_LIM must be a two-element vector.')
end
if (~isnumeric(y_lim) || (length(y_lim) ~= 2))
    error('graticule_greenland:ylim', 'Y_LIM must be a two-element vector.')
end
if (~isnumeric(inc_paral) || (length(inc_paral) ~= 1))
    error('graticule_greenland:incparal', 'INC_PARAL must be a scalar.')
end
if (~isnumeric(inc_merid) || (length(inc_merid) ~= 1))
    error('graticule_greenland:incmerid', 'INC_MERID must be a scalar.')
end
if (nargout ~= 6)
    error('graticule_greenland:nargout', 'Incorrect number of outputs (must be 6).')
end

load mat/gimp_proj gimp_proj

% polar stereographic bounding box
[x_box, y_box]              = deal([x_lim(1) x_lim(1) x_lim(2) x_lim(2) x_lim(1)], [y_lim fliplr(y_lim) y_lim(1)]);
% lat/lon bounding box
[lat_box, lon_box]          = projinv(gimp_proj, (1e3 .* x_box), (1e3 .* y_box));
% lat/lon box edges
[lat_edge, lon_edge]        = deal(cell(1, 4));
for ii = 1:4
    [lat_edge{ii}, lon_edge{ii}] ...
                            = projinv(gimp_proj, (1e3 .* linspace(x_box(ii), x_box(ii + 1))), (1e3 .* linspace(y_box(ii), y_box(ii + 1))));
end

% intersections with box edges at desired inc_lat/inc_lon spacing
[lat, lon]                  = deal(cell(1, 4));
for ii = 1:4
    if (lat_box(ii) > lat_box(ii + 1))
        lat{ii}             = (lat_box(ii) - mod(lat_box(ii), inc_paral)):-inc_paral:lat_box(ii + 1);
    else
        lat{ii}             = (lat_box(ii) - mod(lat_box(ii), -inc_paral)):inc_paral:lat_box(ii + 1);
    end
    if (lon_box(ii) > lon_box(ii + 1))
        lon{ii}             = (lon_box(ii) - mod(lon_box(ii), inc_merid)):-inc_merid:lon_box(ii + 1);
    else
        lon{ii}             = (lon_box(ii) - mod(lon_box(ii), -inc_merid)):inc_merid:lon_box(ii + 1);
    end
end

% parallels at desired spacing
lat_cat                     = catcell(lat);                 % all parallel intersections
paral                       = unique(lat_cat);              % unique parallels
num_paral                   = length(paral);                % number of parallels
inds_sum                    = cumsum(sizec(lat));           % index tracker for all parallels
[x_paral, y_paral]          = deal(cell(1, num_paral));

for ii = 1:num_paral
    inds                    = find(lat_cat == paral(ii));   % pair of indices where current parallel intersects
    if (inds(1) <= inds_sum(1))                             % assign edge indices (1-4) to parallel pair
        ind1                = 1;
    elseif (inds(1) <= inds_sum(2))
        ind1                = 2;
    else
        ind1                = 3;
    end
    if (inds(2) <= inds_sum(2))
        ind2                = 2;
    elseif (inds(2) <= inds_sum(3))
        ind2                = 3;
    else
        ind2                = 4;
    end
    lon1                    = interp1(lat_edge{ind1}, lon_edge{ind1}, paral(ii));               % longitude at parallel intersection for each edge
    lon2                    = interp1(lat_edge{ind2}, lon_edge{ind2}, paral(ii));
    [x_paral{ii}, y_paral{ii}] ...
                            = projfwd(gimp_proj, repmat(paral(ii), 1, 100), linspace(lon1, lon2));  % ps equivalent of current parallel
end

% meridians at desired spacing
lon_cat                     = catcell(lon);
merid                       = unique(lon_cat);
num_merid                   = length(merid);
inds_sum                    = cumsum(sizec(lon));
[x_merid, y_merid]          = deal(cell(1, num_merid));
for ii = 1:num_merid
    inds                    = find(lon_cat == merid(ii));
    if (inds(1) <= inds_sum(1))
        ind1                = 1;
    elseif (inds(1) <= inds_sum(2))
        ind1                = 2;
    else
        ind1                = 3;
    end
    if (inds(2) <= inds_sum(2))
        ind2                = 2;
    elseif (inds(2) <= inds_sum(3))
        ind2                = 3;
    else
        ind2                = 4;
    end
    lat1                    = interp1(lon_edge{ind1}, lat_edge{ind1}, merid(ii));
    lat2                    = interp1(lon_edge{ind2}, lat_edge{ind2}, merid(ii));
    [x_merid{ii}, y_merid{ii}] ...
                            = projfwd(gimp_proj, linspace(lat1, lat2), repmat(merid(ii), 1, 100));
end