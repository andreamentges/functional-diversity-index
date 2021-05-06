%% Computes functional biodiversity based on Rao's index
%
% [rao] = get_functional_diversity(D, int, varargin)
%
% Where:
%   D       column vector of the trait (e.g. mass, H/C ratios) 
%   int     column vector of normalized abundances/intensities (sum==1)
%
%
% Source of Rao index:
% Botta-Dukat, Z. (2005). Rao’s quadratic entropy as a measure of 
% functional diversity on mulitple traits. J. Veg. Sci., 16, 533–540.
% (Note the typo in equation 1b, it should be dij!)
%
% Annotation: 
% The code is not vectorized, so it can only handle a single sample and thus
% only one column of intensities.You'll need to call the function for each 
% of your samples, or tweak it to to handle a whole matrix column-wise.
% Also, the code only handles a single trait. 
% 
% The column vector of a samples' abundances should comprise the same
% molecular formulas in the same order as the column vector of trait values.

function [rao] = get_functional_diversity(D, int)

% normalize the intensities again (just in case)
int = int./ repmat(sum(int), size(int));

% securtiy checks
assert(abs((sum(int))-1)<10^-5 | isnan(sum(int)), 'Intensities are not normalized to sum ==1!')
assert(size(D,1)==size(int,1),...
    'D and int are expected to be column vectors of same length.')
assert(size(int,1)>size(int,2),...
    'Int is expected to be a column vector.')

%% Generate all combinations of molecule pairs

% set up combinations for i = 1 : n-1, j = i+1 : n 
% Note: built-in MATLAB function nchoosek is too slow
n = length(int); 
i     = NaN(1, sum(1:n-1)); 
j     = NaN(1, sum(1:n-1)); 
marker = [0 cumsum(fliplr(1:n-1))];
for m = 1:n-1
    i(marker(m)+1 : marker(m+1))        = m;
    j(marker(m)+1 : marker(m+1))        = m+1:n;
end

%% Calculate distance matrix dij

dij = pdist(D, 'euclidean')';

%% Calculate rao index 

rao = sum(int(i) .* int(j) .*sum(dij,2));    


end
