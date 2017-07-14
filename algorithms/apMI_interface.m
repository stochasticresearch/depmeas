function [m] = apMI_interface(x,y)

ds = [1;1]; mult = 1;
iSHAP = IShannon_AP_initialization(mult);
m = IShannon_AP_estimation([x y]',ds,iSHAP);
                
end