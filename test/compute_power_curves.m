function [powerCurve] = compute_power_curves(M,functionHandlesCell, functionArgsCell, ...
                        nsim_null,nsim_alt)

% standard test configurations
if(nargin<4)
    nsim_null = 500;
    nsim_alt  = 500;
elseif(nargin<5)
    nsim_alt = 500;
end

num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;
num_noise = 30;
noise = 3;
numDepTests = 8;

numDepMeasures = length(functionHandlesCell);

nullDataVec = zeros(numDepMeasures,nsim_null);
altDataVec = zeros(numDepMeasures,nsim_alt);
powerCurve = zeros(numDepMeasures,numDepTests,length(noiseVec));

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
xMin = 0;
xMax = 1;
for lIdx=1:length(noiseVec)
    l = noiseVec(lIdx);
    for typ=1:numDepTests
        dispstat(sprintf('Computing for M=%d noise_level=%d Dependency=%d',M,l, typ),'keepthis', 'timestamp');
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
            for jj=1:numDepMeasures
                f = functionHandlesCell{jj};
                fargs = functionArgsCell{jj};
                nullDataVec(jj,ii)   = f(x,y,fargs{:});
            end
        end

        % compute the rejection cutoffs
        cut = quantile(nullDataVec, 0.95, 2);  % compute row-wise

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
            for jj=1:numDepMeasures
                f = functionHandlesCell{jj};
                fargs = functionArgsCell{jj};
                altDataVec(jj,ii)   = f(x,y,fargs{:});
            end
        end

        % compute the power
        powerCurve(:, typ, lIdx)   = sum(altDataVec > cut, 2)/nsim_alt;  % compute sum row-wise
    end
end

end
