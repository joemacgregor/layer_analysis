function [atten_rate, atten_rate_uncert] ...
                            = atten_gris(temp, H, Cl, NH4, hfcp, temp_uncert, H_uncert, Cl_uncert, NH4_uncert)
% ATTEN_GRIS Radar attenuation-rate model for the Greenland Ice Sheet.
% 
%   [ATTEN_RATE,ATTEN_RATE] = ATTEN_GRIS(TEMP,H,CL,NH4,HFCP,TEMP_UNCERT,H_UNCERT,CL_UNCERT,NH4_UNCERT)
%   calculates the temperature- and impurity-dependent radar attenuation
%   rate ATTEN_RATE (dB/m) for Greenlandic polar ice, based on the ice
%   temperature TEMP (K) and concentrations of the the soluble impurities
%   H+ (H), Cl- (CL) and NH4+ (NH4), whose units are all mol/L. *_UNCERT
%   are their uncertainties. HFCP is a structure defined by HFCONDPROP that
%   provides the necessary dielectric properties.
% 
%   The temperature dependence of this model is based on Stillman et al.,
%   2013, The role of acids in electrical conduction through ice, J.
%   Geophys. Res., 118. The impurity-concentration dependence of this model
%   is the GRIP DEP-chemistry relationship summarized by Wolff et al.,
%   1997, Factors controlling the electrical conductivity of ice from the
%   polar regions - a summary, J. Phys. Chem. B, 101, 6090-6094.
%   
% Joe MacGregor (UTIG)
% Last updated: 07/17/14

if ~exist('hfcond', 'file')
    error('atten_gris:hfcond', 'Function HFCOND is not available within this user''s path.')
end
if ~any(nargin == [5 9])
    error('atten_gris:nargin', ['Number of arguments (' num2str(nargin) ') not equal to four or eight.'])
end
if (~isnumeric(temp) || any(temp < 0))
    error('atten_gris:temp', 'TEMP is not a numeric positive matrix.')
end
if (~isnumeric(H) || any(H < 0))
    error('atten_gris:H', 'H is not a numeric positive matrix.')
end
if (~isnumeric(Cl) || any(Cl < 0))
    error('atten_gris:Cl', 'Cl is not a numeric positive matrix.')
end
if (~isnumeric(NH4) || any(NH4 < 0))
    error('atten_gris:NH4', 'NH4 is not a numeric positive matrix.')
end
if ~isstruct(hfcp)
    error('atten_gris:hfcp', 'HFCP is not a structure.')
end
if (nargin == 9)
    if (~isnumeric(temp_uncert) || any(temp_uncert < 0))
        error('atten_gris:temp_uncert', 'TEMP_UNCERT is not a numeric positive matrix.')
    end
    if (~isnumeric(H_uncert) || any(H_uncert < 0))
        error('atten_gris:H_uncert', 'H_UNCERT is not a numeric positive matrix.')
    end
    if (~isnumeric(Cl_uncert) || any(Cl_uncert < 0))
        error('atten_gris:Cl_uncert', 'Cl_UNCERT is not a numeric positive matrix.')
    end
    if (~isnumeric(NH4_uncert) || any(NH4_uncert < 0))
        error('atten_gris:NH4_uncert', 'NH4_UNCERT is not a numeric positive matrix.')
    end
end
if ~any(nargout == [1 2])
    error('atten_gris:nargout', ['Number of outputs (' num2str(nargout) ') is not 1 or 2.'])
end

permitt_vacuum              = 8.8541878176e-12; % permittivity of the vacuum, F/m
speed_light                 = 299792458; % speed of the light in the vacuum, m/s

% radar attenuation rate (dB/m) calculated from conversion factor * conductivity at reference temperature * Arrhenius-form temperature correction
if (nargin == 5)
    atten_rate              = ((10 * log10(exp(1))) / (sqrt(hfcp.permitt_ice) * permitt_vacuum * speed_light)) .* hfcond(temp, H, Cl, NH4, ones(size(H)), hfcp);
else
    [conduct, conduct_uncert] ...
                            = hfcond(temp, H, Cl, NH4, ones(size(H)), hfcp, temp_uncert, H_uncert, Cl_uncert, NH4_uncert);
    [atten_rate, atten_rate_uncert] ...
                            = deal((((10 * log10(exp(1))) / (sqrt(hfcp.permitt_ice) * permitt_vacuum * speed_light)) .* conduct), (((10 * log10(exp(1))) / (sqrt(hfcp.permitt_ice) * permitt_vacuum * speed_light)) .* conduct_uncert));
end