% MAKE_HDF5 Make HDF5 database of Greenland radiostratigraphy.
% 
% Joe MacGregor
% Last updated: 03/17/15

clear

mat_ver                      = 1.2;

load mat/xy_all name_trans num_trans num_year
load mat/merge_all depth_smooth dist elev_bed_gimp elev_smooth_gimp elev_surf_gimp int_smooth lat lon num_layer num_subtrans num_trace_tot thick time twtt_smooth twtt_surf x_pk y_pk
load mat/date_all age age_uncert
load mat/gimp_proj gimp_proj

name_year_nice              = {'1993 P3' '1995 P3' '1996 P3' '1997 P3' '1998 P3' '1999 P3' '2001 P3' '2002 P3' '2003 P3' '2005 TO' '2006 TO' '2007 P3' '2008 TO' '2009 TO' '2010 DC8' '2010 P3' '2011 P3' '2011 TO' '2012 P3' '2013 P3'};

path_save                   = '/Users/joemac/Desktop/';
file_mat                    = 'Greenland_radiostratigraphy.mat';

%%

% variable writing to MATLAB structure

gris_strat                  = struct;
gris_strat.num_campaign     = num_year;
gris_strat.campaign         = struct;
for ii = 1:num_year
    gris_strat.campaign(ii).name ...
                            = name_year_nice{ii};
    gris_strat.campaign(ii).num_segment ...
                            = num_trans(ii);
    gris_strat.campaign(ii).segment ...
                            = struct;
    for jj = 1:num_trans(ii)
        gris_strat.campaign(ii).segment(jj).division ...
                            = struct;
        gris_strat.campaign(ii).segment(jj).name ...
                            = name_trans{ii}{jj};
        gris_strat.campaign(ii).segment(jj).num_division ...
                            = num_subtrans{ii}(jj);
        for kk = 1:num_subtrans{ii}(jj)
            [gris_strat.campaign(ii).segment(jj).division(kk).distance, gris_strat.campaign(ii).segment(jj).division(kk).elevation_bed, gris_strat.campaign(ii).segment(jj).division(kk).elevation_surface, gris_strat.campaign(ii).segment(jj).division(kk).latitude, ...
             gris_strat.campaign(ii).segment(jj).division(kk).longitude, gris_strat.campaign(ii).segment(jj).division(kk).num_layer, gris_strat.campaign(ii).segment(jj).division(kk).num_trace, gris_strat.campaign(ii).segment(jj).division(kk).thickness, ...
             gris_strat.campaign(ii).segment(jj).division(kk).time_gps, gris_strat.campaign(ii).segment(jj).division(kk).traveltime_surface, gris_strat.campaign(ii).segment(jj).division(kk).x, gris_strat.campaign(ii).segment(jj).division(kk).y] ...
                            = deal(dist{ii}{jj}{kk}, elev_bed_gimp{ii}{jj}{kk}, elev_surf_gimp{ii}{jj}{kk}, lat{ii}{jj}{kk}, lon{ii}{jj}{kk}, num_layer{ii}{jj}(kk), num_trace_tot{ii}{jj}(kk), thick{ii}{jj}{kk}, time{ii}{jj}{kk}, twtt_surf{ii}{jj}{kk}, x_pk{ii}{jj}{kk}, y_pk{ii}{jj}{kk});
            gris_strat.campaign(ii).segment(jj).division(kk).layer ...
                            = struct;
            for ll = 1:num_layer{ii}{jj}(kk)
                [gris_strat.campaign(ii).segment(jj).division(kk).layer(ll).age, gris_strat.campaign(ii).segment(jj).division(kk).layer(ll).age_uncert, gris_strat.campaign(ii).segment(jj).division(kk).layer(ll).depth, ...
                 gris_strat.campaign(ii).segment(jj).division(kk).layer(ll).echo_intensity, gris_strat.campaign(ii).segment(jj).division(kk).layer(ll).elevation, gris_strat.campaign(ii).segment(jj).division(kk).layer(ll).traveltime] ...
                            = deal(age{ii}{jj}{kk}(ll), age_uncert{ii}{jj}{kk}(ll), depth_smooth{ii}{jj}{kk}(ll, :), int_smooth{ii}{jj}{kk}(ll, :), elev_smooth_gimp{ii}{jj}{kk}(ll, :), twtt_smooth{ii}{jj}{kk}(ll, :));
            end
%             gris_strat.campaign(ii).segment(jj).division(kk) ...
%                             = orderfields(gris_strat.campaign(ii).segment(jj).division(kk));
        end
    end
end

% global attributes
gris_strat.title            = 'Traced radiostratigraphy for the Greenland Ice Sheet from ice-penetrating radar data collected by The University of Kansas between 1993 and 2013';
gris_strat.filename         = file_mat;
gris_strat.version          = mat_ver;
gris_strat.citation         = 'MacGregor, J.A., M.A. Fahnestock, G.A. Catania, J.D. Paden, S.P. Gogineni, S.K. Young, S.C. Rybarski, A.N. Mabrey, B.M. Wagman and M. Morlighem, 2015, Radiostratigraphy and age structure of the Greenland Ice Sheet, Journal of Geophysical Research Earth Surface, 120';
gris_strat.date_generated   = datestr(now);
gris_strat.point_of_contact = 'Joseph MacGregor, joemac@ig.utexas.edu';
gris_strat.projection_name  = 'EPSG:3413';
gris_strat.projection       = gimp_proj;
gris_strat.variables        = {'VARIABLE NAME'      'UNITS' 'DESCRIPTION';
                               'age'                'yr'    'reflection age';
                               'age_uncert'         'yr'    'reflection age uncertainty';
                               'distance'           'km'    'distance along segment';
                               'depth'              'm'     'ice depth to reflection';
                               'echo_intensity'     'dB'    'layer relative echo intensity';
                               'elevation'          'm'     'layer elevation referenced to GIMP';
                               'elevation_bed'      'm'     'bed elevation referenced to GIMP';
                               'elevation_surface'  'm'     'surface elevation referenced to GIMP';
                               'latitude'           'deg'   'latitude';
                               'longitude'          'deg'   'longitude';
                               'time_gps'           's'     'measurement time since 1 January 0000';
                               'thickness'          'm'     'ice thickness';
                               'traveltime'         's'     'traveltime to traced reflection';
                               'traveltime_surface' 's'     'traveltime to surface';
                               'x'                  'km'    'x value in EPSG:3413';
                               'y'                  'km'    'y value in EPSG:3413'};

gris_strat                  = orderfields(gris_strat);

save([path_save file_mat], '-v7.3', 'gris_strat');