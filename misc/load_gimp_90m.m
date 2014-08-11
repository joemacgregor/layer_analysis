% LOAD_GIMP_90M Load GIMP 90-m DEM from full-resolution GeoTIFF files.
% 
% Joe MacGregor (UTIG, joemac@ig.utexas.edu)
% Last updated: 06/11/14

clear

file_tiff                   = '/Volumes/polar/Arctic/data/GIMP_DEM/gimpdem_90m.tif';
elev_surf_gimp              = geotiffread(file_tiff);
elev_surf_gimp(elev_surf_gimp == -9999) ...
                            = NaN;
gimp_info                   = geotiffinfo(file_tiff);

[x_min, y_min]              = deal(gimp_info.BoundingBox(1, 1), gimp_info.BoundingBox(1, 2));
[num_x, num_y]              = deal(gimp_info.Width, gimp_info.Height);
gimp_inc                    = gimp_info.GeoTIFFTags.ModelPixelScaleTag(1);

x_gimp                      = (x_min + (gimp_inc / 2)) + (0:gimp_inc:((num_x - 1) * gimp_inc));
y_gimp                      = (y_min + (gimp_inc / 2)) + (0:gimp_inc:((num_y - 1) * gimp_inc))';

elev_surf_gimp              = flipud(elev_surf_gimp);

% OLD / WRONG?
% x_gimp                      = (-639985 + 45) + (0:90:((16620 - 1) * 90));
% y_gimp                      = (-3355565 + 45) + (0:90:((30000 - 1) * 90))';

save('mat/gimp_90m', '-v7.3', 'elev_surf_gimp', 'gimp_info', 'x_gimp', 'y_gimp')