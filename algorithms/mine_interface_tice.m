function [m] = mine_interface_tice(x,y,alpha,c,cfgStr)

minestats = mine(x',y',alpha,c,cfgStr);
m = minestats.tic;
                
end