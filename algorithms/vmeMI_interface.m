function [m] = vmeMI_interface(x,y)

ds = [1;1]; mult = 1;
iSHvME = IShannon_vME_initialization(mult);
m = IShannon_vME_estimation([x y]',ds,iSHvME);
                
end