function [metric] = tau_mi_interface(X,Y)

tau_val = corr(X,Y,'type','kendall');
metric = dep2mi(tau_val);

end