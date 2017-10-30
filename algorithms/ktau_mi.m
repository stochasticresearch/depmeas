function [metric] = ktau_mi(x, y)

metric = corr(x,y,'type','kendall');
metric = dep2mi(metric);