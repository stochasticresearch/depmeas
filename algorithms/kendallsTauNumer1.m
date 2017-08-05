function [K] = kendallsTauNumer1(X,Y)
K = 0;
len = length(X);
for k = 1:len-1
    K = K + sum( sign(X(k)-X(k+1:len)) .* sign(Y(k)-Y(k+1:len)) );
end

end