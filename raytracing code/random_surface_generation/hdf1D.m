function [hdf,bc,h] = hdf1D(f,b,opt)
%
% [hdf,bc,h] = hdf1D(f,b,opt)
%
% calculates the height distribution function and rms height of a 
% 1-d surface profile f(x) using b bins.
%
% Input:    f   - surface heights
%           b   - number of bins
%           opt - optional parameter (type 'hist' for drawing histogram,
%                 type 'plot' for continuous plot) 
%
% Output:   hdf - height distribution function
%           bc  - bin centers
%           h   - rms height
%
% Last updated: 2009-03-11 (David Bergström)
%

format long

h = std(f); % rms height

[hdf,bc] = hist(f,b); 
hdf = hdf/sum(hdf); % normalization to get height distribution function

% optional plotting
if nargin<3 || isempty(opt)
    return;
end;
if nargin==3
    if ischar(opt)
        if strcmp(opt,'hist')
            bar(bc,hdf);
            xlabel('surface height')
            ylabel('probability')
            title('Histogram of height distribution function for height of surface profile y=f(x)')
        elseif strcmp(opt,'plot');
            plot(bc,hdf);
            xlabel('surface height')
            ylabel('probability')
            title('Plot of height distribution function for height of surface profile y=f(x)')
        else fprintf('%s is not a valid option. Type \''help pmf1D\'' for further details.\n',opt);
        end;
    else fprintf('Option must be a string. Type \''help pmf1D\'' for further details.\n');
    end;
end;