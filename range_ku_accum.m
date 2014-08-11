% RANGE_KU_ACCUM Calculate range resolution for each accumulation-radar depth sounder used by KU/CReSIS.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/08/14

clear

load mat/xy_all_accum num_year name_year

speed_vacuum                = 299792458;
permitt_ice                 = 3.15;
fact_window                 = 1.53; % no windowing

bandwidth                   = 1e6 .* [320 320 320];

range_resolution            = (fact_window * speed_vacuum) ./ ((2 * sqrt(permitt_ice)) .* bandwidth);

save mat/range_resolution_accum -v7.3 range_resolution