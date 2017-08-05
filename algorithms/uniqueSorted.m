function [values,instances] = uniqueSorted(u)
coder.inline('always');

p = find([numel(u);diff(u);numel(u)]);
values = u(p(1:end-1));
instances = diff(p);

end