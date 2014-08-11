function id_out             = idtrans(varargin)
% IDTRANS Identify transect by three-number code or name.
% 
% Joe MacGregor (UTIG)
% Last updated: 11/21/13

load mat/xy_all name_trans num_trans num_year
letters                     = 'a':'z';

switch nargin
    case 1
        % determine current year/transect
        tmp1                = varargin{1};
        tmp2                = 'fail';
        if (isnan(str2double(tmp1(end))) || ~isreal(str2double(tmp1(end))))
            tmp3            = tmp1(end);
            tmp1            = tmp1(1:(end - 1));
        else
            tmp3            = 0;
        end
        for ii = 1:num_year
            for jj = 1:num_trans(ii)
                if strcmp(tmp1, name_trans{ii}{jj})
                    break
                end
            end
            if strcmp(tmp1, name_trans{ii}{jj})
                tmp2        = 'success';
                break
            end
        end
        
        if strcmp(tmp2, 'fail')
            disp('Transect not identified. Try again.')
            return
        end
        
        % fix for 2011 P3/TO ambiguity
        if (ii == 17)
            tmp1            = input('2011 TO (press T)?', 's');
            if strcmpi(tmp1, 'T')
                tmp2        = 'fail';
                ii          = 18;
                for jj = 1:num_trans(ii)
                    if strcmp(tmp1, name_trans{ii}{jj})
                        tmp2= 'success';
                        break
                    end
                end
                if strcmp(tmp2, 'fail')
                    disp('Transect incorrectly identified as 2011 TO. Try again.')
                    return
                end
            end
        end
        if ischar(tmp3)
            for kk = 1:length(letters)
                if strcmp(tmp3, letters(kk))
                    break
                end
            end
        else
            kk              = 0;
        end
        if strcmp(tmp2, 'success')
            id_out          = [ii jj kk];
        end
    case 2
        id_out              = name_trans{varargin{1}}{varargin{2}};
    case 3
        if varargin{3}
            id_out          = [name_trans{varargin{1}}{varargin{2}} letters(varargin{3})];
        else
            id_out          = name_trans{varargin{1}}{varargin{2}};
        end
end