function [chisq, atten_model] ...
                            = atten_fit(temp, atten_obs, atten_uncert, depth_inc, H, Cl, NH4, hfcp)
% LOSS_FIT Residual between radar-observed and modeled patterns of englacial returned power.
%   
%   [CHISQ,ATTEN_MODEL] = ATTEN_FIT(TEMP,ATTEN_OBS,ATTEN_UNCERT,DEPTH_INC,H,CL,NH4,HFCP)
%   calculates the chi-squared residual CHISQ between the
%   spreading-corrected echo-intensity loss pattern LOSS with uncertainties
%   LOSS_UNCERT (both in dB) of radar reflectors separated by intervals
%   DEPTH_INC (m) at midpoint depths DEPTH_MID (m) and the pattern
%   calculated by a temperature- and chemistry-dependent radar-attenuation
%   model (ATTEN_MODEL, as described by ATTEN_GRIS). The mean annual
%   surface temperature TEMP_SURF (K) and the concentrations of the soluble
%   impurities H+ (H), Cl- (Cl) and NH4+ (NH4) (mol/L) are also required.
%   HFCP is a structure defined by HFCONDPROP that provides the necessary
%   dielectric properties.
% 
% Joe MacGregor (UTIG)
% Last updated: 07/21/14

if ~exist('atten_gris', 'file')
    error('atten_fit:attengris', 'Function ATTEN_GRIS is not available within this user''s path.')
end
if (nargin ~= 8)
    error('atten_fit:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 8.'])
end
if (~isnumeric(temp) || ~isvector(temp))
    error('atten_fit:temp', 'TEMP is not a numeric vector.')
end
if (~isnumeric(atten_obs) || ~isscalar(atten_obs))
    error('atten_fit:atten_obs', 'ATTEN_OBS is not a numeric scalar.')
end
if (~isnumeric(atten_uncert) || ~isscalar(atten_uncert))
    error('atten_fit:atten_uncert', 'ATTEN_UNCERT is not a numeric scalar.')
end
if (~isnumeric(depth_inc) || ~isvector(depth_inc) || any(depth_inc <= 0))
    error('atten_fit:depthinc', 'DEPTH_INC is not a numeric positive vector.')
end
if (~isnumeric(H) || ~isvector(H) || any(H < 0))
    error('atten_fit:H', 'H is not a numeric positive vector.')
end
if (~isnumeric(Cl) || ~isvector(Cl) || any(Cl < 0))
    error('atten_fit:Cl', 'Cl is not a numeric positive vector.')
end
if (~isnumeric(NH4) || ~isvector(NH4)|| any(NH4 < 0))
    error('atten_fit:NH4', 'NH4 is not a numeric positive vector.')
end
if ~isstruct(hfcp)
    error('atten_fit:hfcp', 'HFCP is not a structure.')
end
if ~isequal(length(depth_inc), length(H), length(Cl), length(NH4))
    error('atten_fit:match', 'DEPTH_INC/H/CL/NH4 vector lengths not equal.')
end
if (nargout > 2)
    error('atten_fit:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 2.'])
end

% temperature model depending on model_type
atten_model                 = sum(atten_gris(temp, H, Cl, NH4, hfcp) .* depth_inc) ./ sum(depth_inc); % mean modeled attenuation rate
chisq                       = ((atten_obs - atten_model) / atten_uncert) ^ 2; % chi-squared residual between observed and modeled attenuation rates