function [m] = mine_interface_mic(x,y,alpha,c,cfgStr)

minestats = mine(x',y',alpha,c,cfgStr);
m = minestats.mic;
                
end