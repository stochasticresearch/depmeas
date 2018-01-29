%% Speed Test taukl_cc vs. taukl_cc_mex
clear;
clc;

numMCSim = 100;
M = 500;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

msi = 0.015625;
alpha = 0.2;

cim_cc_vec = zeros(1,numMCSim);  t1 = 0;
cim_cc_mex_vec = zeros(1,numMCSim); t2 = 0;
cim_cc_deprecated_vec = zeros(1,numMCSim); t3 = 0;
for mcSimNum=1:numMCSim
    dispstat(sprintf('%d/%d',mcSimNum,numMCSim),'timestamp');
    u = rand(M,1); v = rand(M,1);
    
    tic;
    cim_cc_vec(mcSimNum) = cim_cc(u,v,msi,alpha);
    z = toc;
    t1 = t1 + z;
    
    tic;
    cim_cc_mex_vec(mcSimNum) = cim_cc_mex(u,v,msi,alpha);
    z = toc;
    t2 = t2 + z;
    
    tic;
    cim_cc_deprecated_vec(mcSimNum) = cim_cc_deprecated(u,v,msi);
    z = toc;
    t3 = t3 + z;
end

err1 = sum((cim_cc_vec-cim_cc_mex_vec).^2);
err2 = sum((cim_cc_mex_vec-cim_cc_deprecated_vec).^2);
fprintf('err1=%0.02f err2=%0.02f t1=%0.02f t2=%0.02f t3=%0.02f\n',...
         err1, err2, t1/numMCSim, t2/numMCSim, t3/numMCSim);

% Conclusion --> mexing the cim_CC function does speed it up!