function [powerCurve] = cim_power_sensitivity(cimfunc, M, scanincrsToTest)

% standard test configurations
nsim_null = 500;
nsim_alt  = 500;
num_noise = 30;
noise = 3;
numDepTests = 8;

cimNull = zeros(1,nsim_null);
cimAlt = zeros(1,nsim_alt);
powerCurve = zeros(length(scanincrsToTest),numDepTests,num_noise);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;
xMin = 0;
xMax = 1;
for scanincrIdx=1:length(scanincrsToTest)
    minscanincrVal = scanincrsToTest(scanincrIdx);
    for lIdx=1:length(noiseVec)
        l = noiseVec(lIdx);
        for typ=1:numDepTests
            dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
            % simulate data under the null w/ correct marginals
            parfor ii=1:nsim_null
                x = rand(M,1)*(xMax-xMin)+xMin;
                switch(typ)
                    case 1
                        % linear
                        y = x + noise*(l/num_noise)*randn(M,1); 
                    case 2
                        % parabolic
                        y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                    case 3
                        % cubic
                        y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                    case 4
                        % low-freq sin
                        y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                    case 5
                        % high-freq sin
                        y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                    case 6
                        % fourth root
                        y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                    case 7
                        % circle
                        y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                    case 8
                        % step function
                        y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                    otherwise
                        error('unknown dep type!');
                end
                % resimulate x so we have null scenario
                x = rand(M,1)*(xMax-xMin)+xMin;

                % calculate the metrics
                cimNull(ii)   = cimfunc(x,y,minscanincrVal);
            end

            % compute the rejection cutoffs
            cim_cut = quantile(cimNull, 0.95);

            % resimulate the data under the alternative hypothesis
            parfor ii=1:nsim_alt
                x = rand(M,1)*(xMax-xMin)+xMin;
                switch(typ)
                    case 1
                        % linear
                        y = x + noise*(l/num_noise)*randn(M,1); 
                    case 2
                        % parabolic
                        y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                    case 3
                        % cubic
                        y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3) + 10*noise*(l/num_noise)*randn(M,1);
                    case 4
                        % low-freq sin
                        y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                    case 5
                        % high-freq sin
                        y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                    case 6
                        % fourth root
                        y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                    case 7
                        % circle
                        y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                    case 8
                        % step function
                        y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                    otherwise
                        error('unknown dep type!');
                end

                % calculate the metrics
                cimAlt(ii)   = cimfunc(x,y,minscanincrVal);
            end

            % compute the power
            powerCurve(scanincrIdx, typ, lIdx)   = sum(cimAlt > cim_cut)/nsim_alt;
        end
    end
end

end
