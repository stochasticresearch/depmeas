clear;
clc;

% dynamic range estimator tests
x = linspace(0.01,1,100);
y1 = x;
y2 = x.^3;
y3 = x.^12;

var_level = 0.05;

cmap = downsample(colormap,round(64/5));

%plot(x,y1,'+','LineWidth',2);
boundedline(x,y1,var_level*ones(1,length(x)),'+','cmap', cmap(1,:), 'transparency', 0.5);
hold on;

%plot(x,y2,'o','LineWidth',2);
boundedline(x,y2,var_level*ones(1,length(x)),'o','cmap', cmap(2,:), 'transparency', 0.5);

%plot(y2,x,'-','LineWidth',2);
boundedline(y2,x,var_level*ones(1,length(x)),'-','cmap', cmap(3,:), 'transparency', 0.5);

%plot(x,y3,'o','LineWidth',2);
boundedline(x,y3,var_level*ones(1,length(x)),'o','cmap', cmap(4,:), 'transparency', 0.5);

%plot(y3,x,'-','LineWidth',2);
boundedline(y3,x,var_level*ones(1,length(x)),'-','cmap', cmap(5,:), 'transparency', 0.5);

grid on;
xlabel('Dependence Strength','FontSize',16);
ylabel('Estimator Response','FontSize',20);
axis([0,1,0,1])