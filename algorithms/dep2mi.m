function [y] = dep2mi(x)

MAX_VAL = 0.999999;
x(x>MAX_VAL) = MAX_VAL;
x(x<-MAX_VAL) = -MAX_VAL;

y = -0.5*log(1-x.^2);

end