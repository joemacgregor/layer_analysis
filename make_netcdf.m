% MAKE_NETCDF Make NetCDF layer database for distribution
% 
% Joe MacGregor
% Last updated: 03/17/15

clear

nc_ver                      = 1.2;

load mat/greenland_mc_bed_1km thick
thick_grd                   = thick;
clear thick
load mat/age_grd2_krige age_iso age_uncert_norm2 age_norm2_uncert depth_iso2 depth_iso2_uncert depth_norm depth_uncert_iso2 num_age_iso num_depth_norm x_grd y_grd
load mat/age_grd2_krige_alt age_norm2_alt
age_norm2                   = age_norm2_alt;
clear age_norm2_alt

[depth_iso2, depth_iso2_uncert, depth_uncert_iso2] ...
                            = deal((depth_iso2 .* thick_grd(:, :, ones(1, 1, num_age_iso))), (depth_iso2_uncert .* thick_grd(:, :, ones(1, 1, num_age_iso))), (depth_uncert_iso2 .* thick_grd(:, :, ones(1, 1, num_age_iso))));

[age_norm2_uncert((age_norm2_uncert == 0) | (age_norm2 == 0) | isnan(age_norm2)), age_uncert_norm2((age_uncert_norm2 == 0) | (age_norm2 == 0) | isnan(age_norm2)), ...
 depth_iso2_uncert((depth_iso2_uncert == 0) | (depth_iso2 == 0) | isnan(depth_iso2)), depth_uncert_iso2((depth_uncert_iso2 == 0) | (depth_iso2 == 0) | isnan(depth_iso2))] ...
                            = deal(NaN);
[age_norm2(age_norm2 == 0), depth_iso2(depth_iso2 == 0)] ...
                            = deal(NaN);

path_save                   = '/Users/joemac/Desktop/';
file_nc                     = 'Greenland_age_grid.nc';

%%

% variable writing to NetCDF file
nccreate([path_save file_nc], 'depth_iso', 'Dimensions', {'number of grid points in y-direction' size(x_grd, 1) 'number of grid points in x-direction' size(x_grd, 2) 'number of isochrones' (num_age_iso - 1)}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'depth_iso', depth_iso2(:, :, 2:end))
nccreate([path_save file_nc], 'depth_iso_uncert', 'Dimensions', {'number of grid points in y-direction' size(x_grd, 1) 'number of grid points in x-direction' size(x_grd, 2) 'number of isochrones' (num_age_iso - 1)}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'depth_iso_uncert', depth_uncert_iso2(:, :, 2:end))
nccreate([path_save file_nc], 'age_norm', 'Dimensions', {'number of grid points in y-direction' size(x_grd, 1) 'number of grid points in x-direction' size(x_grd, 2) 'number of vertical layers' num_depth_norm}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'age_norm', age_norm2)
nccreate([path_save file_nc], 'age_norm_uncert', 'Dimensions', {'number of grid points in y-direction' size(x_grd, 1) 'number of grid points in x-direction' size(x_grd, 2) 'number of vertical layers' num_depth_norm}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'age_norm_uncert', age_uncert_norm2)
nccreate([path_save file_nc], 'x', 'Dimensions', {'number of grid points in y-direction' size(x_grd, 1) 'number of grid points in x-direction' size(x_grd, 2)}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'x', x_grd)
nccreate([path_save file_nc], 'y', 'Dimensions', {'number of grid points in y-direction' size(x_grd, 1) 'number of grid points in x-direction' size(x_grd, 2)}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'y', y_grd)
nccreate([path_save file_nc], 'age_iso', 'Dimensions', {'number of isochrones' (num_age_iso - 1)}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'age_iso', age_iso(2:end))
nccreate([path_save file_nc], 'num_age_iso', 'Datatype', 'double')
ncwrite([path_save file_nc], 'num_age_iso', (num_age_iso - 1))
nccreate([path_save file_nc], 'depth_norm', 'Dimensions', {'number of vertical layers' num_depth_norm}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'depth_norm', depth_norm)
nccreate([path_save file_nc], 'num_depth_norm', 'Datatype', 'double')
ncwrite([path_save file_nc], 'num_depth_norm', num_depth_norm)
nccreate([path_save file_nc], 'thick', 'Dimensions', {'number of grid points in y-direction' size(x_grd, 1) 'number of grid points in x-direction' size(x_grd, 2)}, 'Datatype', 'double')
ncwrite([path_save file_nc], 'thick', thick_grd)
% nccreate([path_save file_nc], 'Name', 'Dimensions', {'length_filename' length(file_nc)}, 'Datatype', 'char')
% ncwrite([path_save file_nc], 'Name', file_nc)

% global attributes
ncwriteatt([path_save file_nc], '/', 'title', 'Gridded age structure of the Greenland Ice Sheet')
ncwriteatt([path_save file_nc], '/', 'original file name', file_nc)
ncwriteatt([path_save file_nc], '/', 'version', num2str(nc_ver))
ncwriteatt([path_save file_nc], '/', 'citation', ['MacGregor, J.A., M.A. Fahnestock, G.A. Catania, J.D. Paden, S.P. Gogineni, S.K. Young, S.C. Rybarski, A.N. Mabrey, B.M. Wagman and M. Morlighem, 2015, Radiostratigraphy and age structure of the Greenland Ice Sheet, ' ...
                                                  'Journal of Geophysical Research Earth Surface, 120'])
ncwriteatt([path_save file_nc], '/', 'date generated', datestr(now))
ncwriteatt([path_save file_nc], '/', 'point of contact', 'Joseph MacGregor, joemac@ig.utexas.edu')
ncwriteatt([path_save file_nc], '/', 'grid projection', 'EPSG:3413')

% variable-specific attributes
ncwriteatt([path_save file_nc], 'depth_iso', 'description', 'depths of selected isochrones')
ncwriteatt([path_save file_nc], 'depth_iso', 'units', 'm')
ncwriteatt([path_save file_nc], 'depth_iso_uncert', 'description', 'depth uncertainty of selected isochrones')
ncwriteatt([path_save file_nc], 'depth_iso_uncert', 'units', 'm')
ncwriteatt([path_save file_nc], 'age_norm', 'description', 'age at ice-thickness-normalized depths, evenly spaced vertically')
ncwriteatt([path_save file_nc], 'age_norm', 'units', 'yr')
ncwriteatt([path_save file_nc], 'age_norm_uncert', 'description', 'age uncertainty at ice-thickness-normalized depths, evenly spaced vertically')
ncwriteatt([path_save file_nc], 'age_norm_uncert', 'units', 'yr')
ncwriteatt([path_save file_nc], 'x', 'description', 'x-dimension grid centered on Greenland')
ncwriteatt([path_save file_nc], 'x', 'units', 'km')
ncwriteatt([path_save file_nc], 'y', 'description', 'y-dimension grid centered on Greenland')
ncwriteatt([path_save file_nc], 'y', 'units', 'km')
ncwriteatt([path_save file_nc], 'age_iso', 'description', 'age of selected isochrones')
ncwriteatt([path_save file_nc], 'age_iso', 'units', 'yr')
ncwriteatt([path_save file_nc], 'num_age_iso', 'description', 'number of isochrones')
ncwriteatt([path_save file_nc], 'depth_norm', 'description', 'ice-thickness-normalized depth of vertical layers')
ncwriteatt([path_save file_nc], 'depth_norm', 'units', 'depth as fraction of ice thickness')
ncwriteatt([path_save file_nc], 'num_depth_norm', 'description', 'number of vertical layers')
ncwriteatt([path_save file_nc], 'thick', 'description', 'ice thickness')
ncwriteatt([path_save file_nc], 'thick', 'units', 'm')
ncwriteatt([path_save file_nc], 'thick', 'citation', '1-km natural-neighbor-interpolated version of the ice-thickness grid described by Morlighem, M. et al., 2014, Deeply incised submarine glacier valleys beneath the Greenland Ice Sheet, Nature Geosci., 7, 418-422')