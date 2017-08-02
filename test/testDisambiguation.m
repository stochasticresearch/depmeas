%% script to test the disambiguation
clear;
clc;

r1 = [0 0.5 0 0.5; 0.5 1 0.5 1; 0 0 0.5 0.5; 0.5 0.5 1 1];
r1Out = disambiguateRectangleMat(r1);

% r2 = [ 0 0 0.5 0.5; 0.5 0.5 1 1; 0 0.5 0 0.5; 0.5 1 0.5 1];
% r2Out = disambiguateRectangleMat(r2);

% r3 = [[0 0.45 0 1]' [0.45 1 0 1]'];
% r3Out = disambiguateRectangleMat(r3)