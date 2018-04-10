function [z] = conditional_entropy(x,y,is_continuous,num_discrete_bins)
% computes H(Y|X).  If X is discrete, this computation is trivial, if it is
% continuous, then this function discretizes according to the number of
% bins specified by the user.
% This function assumes that Y is discrete

% if(~is_continuous)
%     c = unique(x);
% else
%     c = 1;  % TODO: define this 
% end

if(is_continuous)
    if(nargin>3)
        ndb = num_discrete_bins;
        [N,edges] = histcounts(x, ndb, 'Normalization', 'probability');
    else
        % auto-detect
        [N,edges] = histcounts(x, 'Normalization', 'probability');
    end
else
    error('Not yet implemented!');
end

z = 0;
for nn_idx=1:length(N)
    nn = N(nn_idx);
    I = (x>=edges(nn_idx) & x<edges(nn_idx+1));
    y_subset = y(I);
    z = z + discrete_entropy(y_subset)*nn;
end

end