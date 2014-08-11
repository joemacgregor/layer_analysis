% RANGE_KU Calculate range resolution for each radar depth sounder used by KU/CReSIS.
% 
% Joe MacGregor (UTIG)
% Last updated: 08/08/14

clear

load mat/xy_all num_year name_year

speed_vacuum                = 299792458;
permitt_ice                 = 3.15;
fact_window                 = 0.88; % no windowing

bandwidth                   = 1e6 .* [17 17 17 17 17 17 17 17 20 20 20 20 20 20 20 30 30 30 30 30];

range_resolution            = (fact_window * speed_vacuum) ./ ((2 * sqrt(permitt_ice)) .* bandwidth);

save mat/range_resolution -v7.3 range_resolution