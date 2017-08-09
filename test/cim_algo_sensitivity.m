function [sensitivityResults] = cim_algo_sensitivity(cimfunc, M, scanincrsToTest)

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;
numMCSim = 500;

sensitivityResults = cell(1,length(scanincrsToTest));

for scanincrIdx=1:length(scanincrsToTest)
    minscanincrVal = scanincrsToTest(scanincrIdx);
    % Test Empirical copula distance as a distance measure under various
    % amounts of noise
    linearDep     = zeros(3,length(noiseVec));
    quadraticDep  = zeros(3,length(noiseVec));
    cubicDep      = zeros(3,length(noiseVec));
    sinusoidalDep = zeros(3,length(noiseVec));
    hiFreqSinDep  = zeros(3,length(noiseVec));
    fourthRootDep = zeros(3,length(noiseVec));
    circleDep     = zeros(3,length(noiseVec));
    stepDep       = zeros(3,length(noiseVec));
    indep         = zeros(3,length(noiseVec));
                    % 1 - using tau directly
                    % 2 - segmenting into regions
    for l=noiseVec
        dispstat(sprintf('Computing for noise level=%d >> M=%d',l,M),'keepthis', 'timestamp');
        t1 = 0; c1 = 0; c11 = 0;
        t2 = 0; c2 = 0; c22 = 0;
        t3 = 0; c3 = 0; c33 = 0;
        t4 = 0; c4 = 0; c44 = 0;
        t6 = 0; c6 = 0; c66 = 0;
        t5 = 0; c5 = 0; c55 = 0;
        t7 = 0; c7 = 0; c77 = 0;
        t8 = 0; c8 = 0; c88 = 0;
        t9 = 0; c9 = 0; c99 = 0;
        parfor mcSimNum=1:numMCSim
            x = rand(M,1)*(xMax-xMin)+xMin;
            y1 = x + noise*(l/num_noise)*randn(M,1);
            y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
            y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
            y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
            y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
            y7=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
            y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
            y9 = rand(M,1)*(xMax-xMin)+xMin;

            U1 = pobs([x y1]);
            U2 = pobs([x y2]);
            U3 = pobs([x y3]);
            U4 = pobs([x y4]);
            U5 = pobs([x y5]);
            U6 = pobs([x y6]);
            U7 = pobs([x y7]);
            U8 = pobs([x y8]);
            U9 = pobs([x y9]);

            t1 = t1 + abs(corr(x,y1,'type','kendall'));
            c1 = c1 + abs(corr(x,y1,'type','kendall'));
            c11 = c11 + cimfunc(x,y1,minscanincrVal);

            t2 = t2 + abs(corr(x,y2,'type','kendall'));
            r1 = inBoundedPts(U2(:,1),U2(:,2),0,0.5,0,1); 
            r2 = inBoundedPts(U2(:,1),U2(:,2),0.5,1,0,1);
            c2 = c2 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.5 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.5 );
            c22 = c22 + cimfunc(x,y2,minscanincrVal);

            t3 = t3 + abs(corr(x,y3,'type','kendall'));
            r1 = inBoundedPts(U3(:,1),U3(:,2),0,0.1138,0,1); 
            r2 = inBoundedPts(U3(:,1),U3(:,2),0.1138,0.5768,0,1);
            r3 = inBoundedPts(U3(:,1),U3(:,2),0.5768,1,0,1);
            c3 = c3 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.1138 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.4368 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.4232);
            c33 = c33 + cimfunc(x,y3,minscanincrVal);

            t4 = t4 + abs(corr(x,y4,'type','kendall'));
            r1 = inBoundedPts(U4(:,1),U4(:,2),0,0.1138,0,1); 
            r2 = inBoundedPts(U4(:,1),U4(:,2),0.1138,0.3413,0,1);
            r3 = inBoundedPts(U4(:,1),U4(:,2),0.3413,.6208,0,1);
            r4 = inBoundedPts(U4(:,1),U4(:,2),0.6208,.8483,0,1);
            r5 = inBoundedPts(U4(:,1),U4(:,2),0.8483,1,0,1);
            c4 = c4 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.1138 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.2275 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.2795 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*0.2275 + ...
                        abs(corr(r5(:,1),r5(:,2),'type','kendall'))*0.1517);
            c44 = c44 + cimfunc(x,y4,minscanincrVal);

            t5 = t5 + abs(corr(x,y5,'type','kendall'));
            zz1 = 0; zz2 = 0.02900;   r1 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w1 = zz2-zz1;
            zz1 = zz2; zz2 = 0.09158; r2 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w2 = zz2-zz1;
            zz1 = zz2; zz2 = 0.1528;  r3 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w3 = zz2-zz1;
            zz1 = zz2; zz2 = 0.2108;  r4 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w4 = zz2-zz1;
            zz1 = zz2; zz2 = 0.2715;  r5 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w5 = zz2-zz1;
            zz1 = zz2; zz2 = 0.3339;  r6 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w6 = zz2-zz1;
            zz1 = zz2; zz2 = 0.3999;  r7 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w7 = zz2-zz1;
            zz1 = zz2; zz2 = 0.4613;  r8 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w8 = zz2-zz1;
            zz1 = zz2; zz2 = 0.5227;  r9 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w9 = zz2-zz1;
            zz1 = zz2; zz2 = 0.5893;  r10 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w10 = zz2-zz1;
            zz1 = zz2; zz2 = 0.6547;  r11 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w11 = zz2-zz1;
            zz1 = zz2; zz2 = 0.7135;  r12 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w12 = zz2-zz1;
            zz1 = zz2; zz2 = 0.7734;  r13 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w13 = zz2-zz1;
            zz1 = zz2; zz2 = 0.8438;  r14 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w14 = zz2-zz1;
            zz1 = zz2; zz2 = 0.9074;  r15 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w15 = zz2-zz1;
            zz1 = zz2; zz2 = 0.9694;  r16 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w16 = zz2-zz1;
            zz1 = zz2; zz2 = 1;       r17 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w17 = zz2-zz1;
            c5 = c5 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*w1 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*w2 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*w3 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*w4 + ...
                        abs(corr(r5(:,1),r5(:,2),'type','kendall'))*w5 + ...
                        abs(corr(r6(:,1),r6(:,2),'type','kendall'))*w6 + ...
                        abs(corr(r7(:,1),r7(:,2),'type','kendall'))*w7 + ...
                        abs(corr(r8(:,1),r8(:,2),'type','kendall'))*w8 + ...
                        abs(corr(r9(:,1),r9(:,2),'type','kendall'))*w9 + ...
                        abs(corr(r10(:,1),r10(:,2),'type','kendall'))*w10 + ...
                        abs(corr(r11(:,1),r11(:,2),'type','kendall'))*w11 + ...
                        abs(corr(r12(:,1),r12(:,2),'type','kendall'))*w12 + ...
                        abs(corr(r13(:,1),r13(:,2),'type','kendall'))*w13 + ...
                        abs(corr(r14(:,1),r14(:,2),'type','kendall'))*w14 + ...
                        abs(corr(r15(:,1),r15(:,2),'type','kendall'))*w15 + ...
                        abs(corr(r16(:,1),r16(:,2),'type','kendall'))*w16 + ...
                        abs(corr(r17(:,1),r17(:,2),'type','kendall'))*w17);
            c55 = c55 + cimfunc(x,y5,minscanincrVal);

            t6 = t6 + abs(corr(x,y6,'type','kendall'));
            c6 = c6 + abs(corr(x,y6,'type','kendall'));
            c66 = c66 + cimfunc(x,y6,minscanincrVal);

            t7 = t7 + abs(corr(x,y7,'type','kendall'));
            r1 = inBoundedPts(U7(:,1),U7(:,2),0,0.5,0,0.5); 
            r2 = inBoundedPts(U7(:,1),U7(:,2),0.5,1,0,0.5);
            r3 = inBoundedPts(U7(:,1),U7(:,2),0,0.5,0.5,1);
            r4 = inBoundedPts(U7(:,1),U7(:,2),0.5,1,0.5,1);
            c7 = c7 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.25 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.25 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.25 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*0.25);
            c77 = c77 + cimfunc(x,y7,minscanincrVal);

            t8 = t8 + abs(corr(x,y8,'type','kendall'));
            c8 = c8 + abs(corr(x,y8,'type','kendall'));
            c88 = c88 + cimfunc(x,y8,minscanincrVal);

            t9 = t9 + abs(corr(x,y9,'type','kendall'));
            c9 = c9 + abs(corr(x,y9,'type','kendall'));
            c99 = c99 + cimfunc(x,y9,minscanincrVal);
        end
        linearDep(1,l+1) = t1/numMCSim; linearDep(2,l+1) = c1/numMCSim; linearDep(3,l+1) = c11/numMCSim;
        quadraticDep(1,l+1) = t2/numMCSim; quadraticDep(2,l+1) = c2/numMCSim; quadraticDep(3,l+1) = c22/numMCSim;
        cubicDep(1,l+1) = t3/numMCSim; cubicDep(2,l+1) = c3/numMCSim; cubicDep(3,l+1) = c33/numMCSim;
        sinusoidalDep(1,l+1) = t4/numMCSim; sinusoidalDep(2,l+1) = c4/numMCSim; sinusoidalDep(3,l+1) = c44/numMCSim;
        hiFreqSinDep(1,l+1) = t5/numMCSim; hiFreqSinDep(2,l+1) = c5/numMCSim; hiFreqSinDep(3,l+1) = c55/numMCSim;
        fourthRootDep(1,l+1) = t6/numMCSim; fourthRootDep(2,l+1) = c6/numMCSim; fourthRootDep(3,l+1) = c66/numMCSim;
        circleDep(1,l+1) = t7/numMCSim; circleDep(2,l+1) = c7/numMCSim; circleDep(3,l+1) = c77/numMCSim;
        stepDep(1,l+1) = t8/numMCSim; stepDep(2,l+1) = c8/numMCSim; stepDep(3,l+1) = c88/numMCSim;
        indep(1,l+1) = t9/numMCSim; indep(2,l+1) = c9/numMCSim; indep(3,l+1) = c99/numMCSim;
    end
    x = struct;
    x.linearDep = linearDep; x.quadraticDep = quadraticDep; x.cubicDep = cubicDep;
    x.sinusoidalDep = sinusoidalDep; x.hiFreqSinDep = hiFreqSinDep;
    x.fourthRootDep = fourthRootDep; x.circleDep = circleDep;
    x.stepDep = stepDep; x.indep = indep;
    sensitivityResults{scanincrIdx} = x;
end

end