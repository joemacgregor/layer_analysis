% MAKE_HDF5_FRAME Make HDF5 database of Greenland deep radiostratigraphy aligned with CReSIS frames.
% 
% Joe MacGregor
% Last updated: 06/22/15

clear

mat_ver                      = 1.2; % manually updated version number

load mat/xy_all name_trans name_year num_frame num_trans num_year % name and number of transects
load mat/date_all age age_uncert % inferred layer ages

% nicer campaign names
name_year_nice              = {'1993 P3' '1995 P3' '1996 P3' '1997 P3' '1998 P3' '1999 P3' '2001 P3' '2002 P3' '2003 P3' '2005 TO' '2006 TO' '2007 P3' '2008 TO' '2009 TO' '2010 DC8' '2010 P3' '2011 P3' '2011 TO' '2012 P3' '2013 P3'};

% dir_ref                     = '/Volumes/icebridge/data/';
dir_ref                     = '~/icebridge/';
path_save                   = 'mat/';
% path_save                   = '/Users/joemac/Desktop/';
file_mat                    = 'Greenland_radiostratigraphy_frame.mat';

% letters for sub-transects
letter_str                  = 'a':'z';
letter_cell                 = cell(1, 26);
for ii = 1:26
    letter_cell{ii}         = letter_str(ii);
end

%% variable writing to MATLAB structure

% initialize structure
gris_strat                  = struct;
gris_strat.num_campaign     = num_year;
gris_strat.campaign         = struct;

% loop through each campaign
for ii = 1:num_year
    
    disp(name_year{ii})
    
    % initialize campaign structure
    gris_strat.campaign(ii).name ...
                            = name_year_nice{ii};
    gris_strat.campaign(ii).num_segment ...
                            = num_trans(ii);
    gris_strat.campaign(ii).segment ...
                            = struct;
    
    % loop through each transect/segment
    for jj = 1:num_trans(ii)
        
        disp(name_trans{ii}{jj})
        
        % set up structure for each frame
        gris_strat.campaign(ii).segment(jj).frame ...
                            = struct;
        gris_strat.campaign(ii).segment(jj).name ...
                            = name_trans{ii}{jj};
        gris_strat.campaign(ii).segment(jj).num_frame ...
                            = num_frame{ii}(jj); % number of frames in each segment
        
        % extract frame filenames and merged picks filenames
        name_frame_curr     = dir([dir_ref name_year{ii} '/data/' name_trans{ii}{jj}  '/*.mat']);
        name_frame_curr     = {name_frame_curr.name};
        name_merge_curr     = dir([dir_ref name_year{ii} '/merge/' name_trans{ii}{jj} '*.mat']);
        name_merge_curr     = {name_merge_curr.name};
        
        % do not proceed if no merged picks files
        if isempty(name_merge_curr)
            continue
        end
        
        % loop through each frame
        for kk = 1:num_frame{ii}(jj)
            
            mm              = 1; % dummy counter
            pk              = []; % empty pk file to start
            
            % load each merge file and extract its data onto relevant frame positions
            while true
                pk          = pickinterp([dir_ref name_year{ii} '/merge/'], name_merge_curr{mm}, [dir_ref name_year{ii} '/data/' name_trans{ii}{jj} '/'], name_frame_curr{kk}(1:(end - 4)), 'frame', false); % pickinterp is where the magic happens
                if (isempty(pk) && (mm < length(name_merge_curr))) % move on to next merge file if pk empty and some merge files left for this campaign
                    mm      = mm + 1;
                else
                    break
                end
            end
            
            % always move on if pk empty
            if isempty(pk)
                continue
            end
            
            % extract position variables
            [gris_strat.campaign(ii).segment(jj).frame(kk).elevation_bed, gris_strat.campaign(ii).segment(jj).frame(kk).elevation_surface, gris_strat.campaign(ii).segment(jj).frame(kk).num_layer] ...
                            = deal(pk.elev_bed_gimp, pk.elev_surf_gimp, length(pk.layer));
            
            % re-assign layer structure
            gris_strat.campaign(ii).segment(jj).frame(kk).layer ...
                            = pk.layer;
            
            % determine which sub-transect we're on
            [~, tmp]        = strtok(name_merge_curr{mm}, '_');
            tmp             = strtok(tmp, '_');
            if ~isempty(find(strcmp(letter_cell, tmp(end)), 1))
                mm          = find(strcmp(letter_cell, tmp(end)));
            end
            
            % loop through each layer and assign age
            for ll = 1:pk.num_layer
                if isempty(age{ii}{jj}) % oops no ages
                    [gris_strat.campaign(ii).segment(jj).frame(kk).layer(ll).age, gris_strat.campaign(ii).segment(jj).frame(kk).layer(ll).age_uncert] ...
                            = deal(NaN);
                else % assign ages based on pk.ind_match (provides master layer ID of layers within the frame, determined pickinterp)
                    [gris_strat.campaign(ii).segment(jj).frame(kk).layer(ll).age, gris_strat.campaign(ii).segment(jj).frame(kk).layer(ll).age_uncert] ...
                            = deal(age{ii}{jj}{mm}(pk.ind_match(ll)), age_uncert{ii}{jj}{mm}(pk.ind_match(ll)));
                end
                % assign layer echo intensity, elevation and traveltime
                [gris_strat.campaign(ii).segment(jj).frame(kk).layer(ll).echo_intensity, gris_strat.campaign(ii).segment(jj).frame(kk).layer(ll).elevation, gris_strat.campaign(ii).segment(jj).frame(kk).layer(ll).traveltime] ...
                            = deal(pk.layer(ll).int_smooth, pk.layer(ll).elev_smooth_gimp, pk.layer(ll).twtt_smooth);
            end
            % remove unnecessary variables
            gris_strat.campaign(ii).segment(jj).frame(kk).layer ...
                            = rmfield(gris_strat.campaign(ii).segment(jj).frame(kk).layer, {'depth_smooth' 'elev' 'elev_gimp' 'elev_smooth' 'elev_smooth_gimp' 'ind_y' 'ind_y_smooth' 'int' 'int_smooth' 'twtt' 'twtt_ice' 'twtt_ice_smooth' 'twtt_smooth'});
        end
    end
end

% miscellaneous variables
gris_strat.title            = 'Traced radiostratigraphy for the Greenland Ice Sheet from ice-penetrating radar data collected by The University of Kansas between 1993 and 2013';
gris_strat.filename         = file_mat;
gris_strat.version          = mat_ver;
gris_strat.citation         = 'MacGregor, J.A., M.A. Fahnestock, G.A. Catania, J.D. Paden, S.P. Gogineni, S.K. Young, S.C. Rybarski, A.N. Mabrey, B.M. Wagman and M. Morlighem, 2015, Radiostratigraphy and age structure of the Greenland Ice Sheet, Journal of Geophysical Research Earth Surface, 120';
gris_strat.date_generated   = datestr(now);
gris_strat.point_of_contact = 'Joseph MacGregor, joemac@ig.utexas.edu';
gris_strat.variables        = {'VARIABLE NAME'      'UNITS' 'DESCRIPTION';
                               'age'                'yr'    'reflection age';
                               'age_uncert'         'yr'    'reflection age uncertainty';
                               'depth'              'm'     'ice depth to reflection';
                               'echo_intensity'     'dB'    'layer relative echo intensity';
                               'elevation'          'm'     'layer elevation referenced to GIMP';
                               'elevation_bed'      'm'     'bed elevation referenced to GIMP';
                               'elevation_surface'  'm'     'surface elevation referenced to GIMP';
                               'traveltime'         's'     'traveltime to traced reflection'};

% order fields alphabetically
gris_strat                  = orderfields(gris_strat);

save([path_save file_mat], '-v7.3', 'gris_strat');