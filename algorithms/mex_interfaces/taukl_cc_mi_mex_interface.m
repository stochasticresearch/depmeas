function [metric] = taukl_cc_mi_mex_interface(X,Y,autoDetectHybrid,isHybrid,continuousRvIndicator)

[u,v] = pobs_sorted_cc_mex(X,Y); 
taukl_val = taukl_cc_mex(u,v,autoDetectHybrid,isHybrid,continuousRvIndicator);
metric = dep2mi(taukl_val);

end