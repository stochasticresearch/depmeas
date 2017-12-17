function all_colors = helper_power_colormap()
rows = 100;
up_grad = linspace(0, 1, rows)';
dn_grad = flipud(up_grad);
orange = [1 .5 0];
blue = [.5 .7 1]; %[0 0 1];
white = [1 1 1];


% Create a gradient colormap running from red to blue
all_colors = [...
    dn_grad*blue + up_grad*white;
    dn_grad*white + up_grad*orange];
        
end
