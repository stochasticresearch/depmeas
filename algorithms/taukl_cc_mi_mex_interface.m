function [metric] = taukl_cc_mi_mex_interface(X,Y,autoDetectHybrid,isHybrid,continuousRvIndicator)

[u,v] = pobs_sorted_cc_mex(X,Y); 
taukl_val = taukl_cc_mex(u,v,int32(autoDetectHybrid),int32(isHybrid),int32(continuousRvIndicator));
metric = dep2mi(taukl_val);

end