function [taukl_val] = taukl_cc_mex_interface(X,Y,autoDetectHybrid,isHybrid,continuousRvIndicator)

[u,v] = pobs_sorted_cc_mex(X,Y); 
taukl_val = taukl_cc_mex(u,v,autoDetectHybrid,isHybrid,continuousRvIndicator);

end