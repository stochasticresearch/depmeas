%% Test the ktau_numer functionality

clear;
clc;

x = rand(500,1); y = rand(500,1);
[u,v] = pobs_sorted_cc(x,y);

numer1 = kendallsTauNumer1(u,v);
numer2 = kendallsTauNumer2(u,v);
numer3 = ktau_numer(u,v);

fprintf('numer1=%d numer2=%d numer3=%d\n',numer1,numer2,numer3);

%% Test old taukl_cc vs. the new taukl_cc w/ the C offload of numerator calculation
clear;
clc;

M = 5000;

numSims = 100;

t1 = 0; t2 = 0;
t1_vec = zeros(1,numSims);
t2_vec = zeros(1,numSims);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

tauVec = 0.1:.1:.9;
%copulasToTest = {'rand', 'gaussian','gumbel','frank','clayton'};
copulasToTest = {'rand'};

for cIdx=1:length(copulasToTest)
    c = copulasToTest{cIdx};
    for tIdx=1:length(tauVec)
        t = tauVec(tIdx);
        for simNum=1:numSims
            dispstat(sprintf('%d/%d',simNum, numSims),'timestamp');
            % generate data for this copula and tau value
            if(strcmpi(c,'rand'))
                x = rand(M,1); y = rand(M,1);
            else
                depVal = copulaparam(c,t);
                U = copularnd(c,depVal,M);
                x = U(:,1); y = U(:,2);
            end

            [u,v] = pobs_sorted_cc(x,y);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tic;
%             t1_vec(simNum) = taukl_cc_mex(u,v,int32(1),int32(0),int32(0));
            t1_vec(simNum) = taukl_cc_deprecated(u,v,1,0,0);
            z = toc;
            t1 = t1 + z;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tic;
            t2_vec(simNum) = double(taukl_cc(u,v,1,0,0));
            z = toc;
            t2 = t2 + z;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        dispstat(sprintf('err[%s(%0.02f)]=%0.02f', c,t,mean((t2_vec-t1_vec).^2)),'keepthis','timestamp');
    end
end

fprintf('t1=%0.02f t2=%0.02f\n',t1,t2);

%% Test it for hybrid cases
clear;
clc;
% dbstop if error;

M = 500;
tau = 0.7;
cop = 'gaussian';

iTau = copulaparam(cop,tau,'type','kendall');

U = copularnd(cop,iTau,M);
distObj1 = makedist('Normal');
X = icdf(distObj1,U(:,1));
distObj2 = makedist('Multinomial','probabilities',[0.5,0.5]);
Y = icdf(distObj2,U(:,2));

[X,Y] = pobs_sorted_cc(X,Y);

if(length(unique(X))~=M || length(unique(Y))~=2)
    warning('data messed uP?');
end
% X = [1 2 3 4 5 6 7 8];
% Y = [0 1 1 1 1 1 1 1];
% 
% X
% Y

taukl_cc(X,Y,0,1,0)
taukl_cc_v2(X,Y,0,1,0)

if(length(unique(X))~=M || length(unique(Y))~=2)
    warning('data messed uP?');
end