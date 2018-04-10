function [y] = entropy2d(c,u,v)
    % Computes the entropy of a 2-D distribution estimated based on the
    % intervals of u & v
    
    % assume uniform sampling of u and v
    delta_u = mean(diff(u(1,:)));
    delta_v = mean(diff(v(:,1)));
    deltaA  = delta_u*delta_v;
    y = 0;
    for u_idx=1:length(u)
        for v_idx=1:length(v)
            c_val = c(u_idx,v_idx);
            if(c_val>0)
                y = y + c_val*log(c_val)*deltaA;
            end
        end
    end
    y = y*-1;  % entropy is negative integral
end