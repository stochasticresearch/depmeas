function [metric] = cosf(x, y)
%COSF - a copula based test of dependence called the copula statistic. See
%       <TODO: insert paper reference here>
% Inputs:
%  x - realizations of the first random variable
%  y - realizations of the second random variable
% Outputs
%  metric - a real-valued number between 0 and 1 which represents the
%           statistical dependence between x and y
%
% TODO:
%   [ ] - Make this multivariate.
%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

u = pobs(x);
v = pobs(y);

n = length(x);
E_copula = zeros(1,n);

for ii=1:n
    E_copula(ii) = sum(u(ii)>=u & v(ii)>=v)/(n+1);
end
W_copula = min(u,v);
M_copula = max(u+v-1,0);
Pi_copula = u.*v;

score = 0; wanc = 1;

if(E_copula(1)==Pi_copula(1))
    wanc = 0;
end
if( (E_copula(1)==W_copula(1)) || (E_copula(1)==M_copula(1)) || (E_copula(1)==E_copula(2)) )
    wanc = 1;
end

ss_minIdx = 1;
for ii=2:n
    ss = E_copula(ss_minIdx:ii);
    if(~issorted(ss) && ~issorted(fliplr(ss)))      % just a way of testing if they are tied?
        
        if(E_copula(ii-1)==E_copula(ii-2))
            jj = ii-2;
        else
            jj = ii-1;
        end
        
        % compute the relative distance function
        if (E_copula(jj)>Pi_copula(jj))
            w = (E_copula(jj)-Pi_copula(jj))/(W_copula(jj)-Pi_copula(jj));
        else
            if (E_copula(jj)<Pi_copula(jj))
                w=(Pi_copula(jj)-E_copula(jj))/(Pi_copula(jj)-M_copula(jj));
            else
                w=0;
            end
        end
      
        if (( (round(E_copula(ii-1),2)==round(W_copula(ii-1),2)) || ...
              (round(E_copula(ii-1),2)==round(M_copula(ii-1),2)) || ...
              (E_copula(ii-1)==(1/(n+1)))) && (length(ss)>4) )
            w=1;
        end
      
        if (( (E_copula(ii)==E_copula(ii-2)) || ...
              (E_copula(ii-1)==E_copula(ii-2))) && (length(ss)>4) ) 
            w=1;
        end
      
        if (( (round(Pi_copula(ii-1),2)==round(W_copula(ii-1),2)) || ...
              (round(Pi_copula(ii-1),2)==round(M_copula(ii-1),2))) && (length(ss)>4))
            w=1;
        end
        
        score = score+(length(ss)-1)*(w+wanc)/2;
        wanc = w;
        ss_minIdx = ii;
    else
        if(ii==n)
            % compute the relative distance function
            if (E_copula(ii)>Pi_copula(ii))
                w=(E_copula(ii)-Pi_copula(ii))/(W_copula(ii)-Pi_copula(ii));
            else
                if (E_copula(ii)<Pi_copula(ii))
                    w=(Pi_copula(ii)-E_copula(ii))/(Pi_copula(ii)-M_copula(ii));
                else
                    w=0;
                end
            end
            
            if (( (E_copula(ii-1)==W_copula(ii-1)) || ...
                  (E_copula(ii-1)==M_copula(ii-1))) && (length(ss)>=4 ))
                w=1;
            end
            
            score=score+(length(ss)-1)*(w+wanc)/2;
        end
    end
end

metric = (score+1)/n;

end