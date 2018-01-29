%% Speed Test taukl_cc vs. taukl_cc_mex
clear;
clc;

numMCSim = 1000;
M = 10000;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

taukl_cc_vec = zeros(1,numMCSim);  t1 = 0;
taukl_cc_mex_vec = zeros(1,numMCSim); t2 = 0;
for mcSimNum=1:numMCSim
    dispstat(sprintf('%d/%d',mcSimNum,numMCSim),'timestamp');
    u = rand(M,1); v = rand(M,1);
    
    tic;
    taukl_cc_vec(mcSimNum) = taukl_cc(u,v,1,0,0);
    z = toc;
    t1 = t1 + z;
    
    tic;
    taukl_cc_mex_vec(mcSimNum) = taukl_cc_mex(u,v,1,0,0);
    z = toc;
    t2 = t2 + z;
end

err = sum((taukl_cc_vec-taukl_cc_mex_vec).^2);
fprintf('err=%0.02f t1=%0.02f t2=%0.02f\n', err, t1, t2);

% Conclusion --> mexing the TAUKL_CC function does *not* speed it up!