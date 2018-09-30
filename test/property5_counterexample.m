
%% Counter-Example of Property 5 by Reviewer #3 ISC

clear;
clc;

x = rand(500,1);
y = zeros(500,1);
for ii=1:500
    if(x(ii)>=0 && x(ii)<=1/2)
        y(ii) = 0.5+x(ii);
    else
        y(ii) = -0.5+x(ii);
    end
end

[v,rr] = cim(x,y)

scatter(x,y)
hold on;
% plot the region boundary
darkGreenRGB = [0 0.5 0];
plot([0.5 0.5], [0 1], '--','LineWidth',3,'Color',darkGreenRGB);
grid on;
xlabel('x','FontSize',20);
ylabel('g(x)','FontSize',20);