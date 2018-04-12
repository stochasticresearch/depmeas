function [z] = h_mi_interface(X,Y,X_continuous,num_discrete_bins)
% computes MI(X,Y) using the relation H(Y) - H(Y|X).  This may be useful
% measure if X is continuous and Y is discrete.

% Inputs:
% Outputs:

if(nargin>3)
    z = discrete_entropy(Y) - conditional_entropy(X,Y,X_continuous,num_discrete_bins);
else
    z = discrete_entropy(Y) - conditional_entropy(X,Y,X_continuous);
end

end