% MINMAXNORM Minimum-Maximum normalization.
%   Y = MINMAXNORM(X) normalizes data in X (N samples-by-D features) with
%   the minimum and maximum values in X. The output data range is [-1,+1].
%
%   [Y,MN,MX] = MINMAXNORM(X) returns the minimum and maximum of each column 
%   in X.
%
%   Y = MINMAXNORM(X,[MN;MX]) normalizes data in X with the given minimum and 
%   maximum values.
%
%   Example:
%   -------
%   X1 = rand(10,2);
%   [Y1,mn,mx] = minmaxnorm(X1);
%   X2 = rand(10,2);
%   Y2 = minmaxnorm(X2,[mn;mx]);
%
%   See also SOFTMAXNORM SOFTMAXNORM
%
%
%   Reference:
%   ---------
%   K. L. Priddy, P. E. Keller, Artificial Neural Networks: An Introduction.
%   Bellingham, WA: SPIE-The Int. Soc. Optical Eng., 2005.

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico) - LUS/PEB/COPPE/UFRJ (Brazil)
%   MINMAXNORM Version 1.0 (Matlab R2014a Unix)
%   November 2016
%   Copyright (c) 2016, Wilfrido Gomez Flores
% ------------------------------------------------------------------------

function [vect,dmn,dmx] = minmaxnorm(vect,stats)
if nargin == 1
    dmx = max(vect,[],1);
    dmn = min(vect,[],1);
elseif nargin == 2
    dmn = stats(1,:);
    dmx = stats(2,:);
end
N = size(vect,1);
mx = repmat(dmx,N,1);
mn = repmat(dmn,N,1);
%vect = (2.*(vect - mn)./(mx-mn))-1;
vect = (vect - mn)./(mx-mn);