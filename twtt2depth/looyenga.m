function permitt_mix        = looyenga(vol_frac, permitt_clast, permitt_matrix)
% LOOYENGA Permittivity of a mixture from Looyenga's equation.
% 
%   PERMITT_MIX = LOOYENGA(VOL_FRAC,PERMITT_CLAST,PERMITT_MATRIX)
%   calculates the permittivity PERMIT_MIX of a mixture of two materials of
%   different permittivities using Looyenga (1965)'s classic relationship,
%   the permittivities of clasts and matrix materials, PERMITT_CLAST and
%   PERMITT_MATRIX, respectively, and the volume fraction VOL_FRAC, which
%   is the fraction of the total volume that the clast material occupies.
%   For example, VOL_FRAC can be calculated as the ratio of the density of
%   an ice sample to that of pure ice.
%   
%   Reference:
%   Looyenga, H., 1965, Dielectric constants of heterogenous mixtures,
%   Physica, 31(3), 401-406
%   
%   permitt_mix = ((vol_frac * ((permitt_clast ^ (1 / 3)) - (permitt_matrix ^ (1 / 3)))) + (permitt_matrix ^ (1 / 3))) ^ 3
%   
% Joe MacGregor (UTIG, joemac@ig.utexas.edu)
% Last updated: 08/05/14

if (nargin ~= 3)
    error('looyenga:nargin', 'Number of inputs must be equal to 3.')
end
if (~isnumeric(vol_frac) || any(vol_frac(:) <= 0))
    error('looyenga:vol_frac', 'VOL_FRAC is not positive numeric.')
end
if (~isnumeric(permitt_clast) || any(permitt_clast(:) <= 0))
    error('looyenga:permitt_clast', 'PERMITT_CLAST is not positive numeric.')
end
if (~isnumeric(permitt_matrix) || any(permitt_matrix(:) <= 0))
    error('looyenga:permitt_matrix', 'PERMITT_MATRIX is not positive numeric.')
end
if ~isequal(size(vol_frac), size(permitt_clast), size(permitt_matrix))
    error('looyenga:sizes', 'Inputs are not the same size.')
end
if (nargout > 1)
    error('looyenga:nargout', 'Number of outputs cannot be more than 1.')
end

permitt_mix                 = ((vol_frac .* ((permitt_clast .^ (1 / 3)) - (permitt_matrix .^ (1 / 3)))) + (permitt_matrix .^ (1 / 3))) .^ 3;