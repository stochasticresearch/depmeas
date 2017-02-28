clear;
clc;

rng(123);

fontSize = 20;

M = 500;
x = rand(M,1)*10+2;
y1 = x; y2 = exp(x);
sigma = 2;
noise = sigma*randn(M,1);

u = pobs(x);
v1 = pobs(y1); v1_n = pobs(y1+noise);
v2 = pobs(y2); v2_n = pobs(y2+noise);

subplot(2,2,1); scatter(x,y1); grid on; title('(a)', 'FontSize', fontSize);
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
subplot(2,2,2); scatter(u,v1_n); 
xlabel('u', 'FontSize', 20); ylabel('v', 'FontSize', 20);
title('(b)', 'FontSize', fontSize); grid on;

subplot(2,2,3); scatter(x,y2); grid on; title('(c)', 'FontSize', fontSize);
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
subplot(2,2,4); scatter(u,v2_n); 
xlabel('u', 'FontSize', 20); ylabel('v', 'FontSize', 20);
title('(d)', 'FontSize', fontSize); grid on;