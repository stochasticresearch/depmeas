function [resid, residAssocIdxs] = ...
    computeResidual_MISC(rectangleCellAggr, ax1pts, ax2pts)
% the cell array that is passed in rectangles for each of the scan-incr's
% we compute the residual for each of these box configurations, and take
% the one that has the lowest residual (i.e. best fit).

seePlots = 1;
if(seePlots)
    f = figure(1);
end

resid = {}; residAssocIdxs = {};
residIdx = 1;
for kk=1:length(rectangleCellAggr)
    rectangleCell = rectangleCellAggr{kk};
    residualCell = cell(1,length(rectangleCell));
    residAssocIdxsCell = cell(1,length(rectangleCell));
    for ii=1:length(rectangleCell)
        rectangles = rectangleCell{ii};
        residualMat = cell(1,size(rectangles,2)); 
        residAssocIdxsMat = cell(1,size(rectangles,2));
        residualMatIdx = 1;
        for jj=1:size(rectangles,2)
            rectangle = rectangles(:,jj);
            ax1min = rectangle(1); ax1max = rectangle(2);
            ax2min = rectangle(3); ax2max = rectangle(4);

            ax1_match = find(ax1pts>ax1min & ax1pts<=ax1max);
            ax2_match = find(ax2pts>ax2min & ax2pts<=ax2max);
            matchIdxs = intersect(ax1_match,ax2_match);
            matchPts = [ax1pts(matchIdxs) ax2pts(matchIdxs)];

            XX = [ones(length(matchIdxs),1) matchPts(:,1)];
            zz = XX\matchPts(:,2);
            offset = zz(1); regressionSlope = zz(2);

            yPredict = regressionSlope*matchPts(:,1)+offset;
            residualVec = matchPts(:,2) - yPredict;
            
            % some debug plotting
            if(seePlots)
                clf(f);
                figure(f);
                scatter(ax1pts, ax2pts); hold on;
                scatter(matchPts(:,1),matchPts(:,2), 'k'); hold on;
                scatter(matchPts(:,1),yPredict,'r'); 
                grid on;
                rectangle
                pause;
            end

            residualMat{residualMatIdx} = residualVec;
            residAssocIdxsMat{residualMatIdx} = matchIdxs;
            residualMatIdx = residualMatIdx + 1;
        end
        residualCell{ii} = residualMat;
        residAssocIdxsCell{ii} = residAssocIdxsMat;
    end

    % compute the sum squared of all the residuals, and take the minimum
    minErr = Inf; minErrIdx = 1;
    for ii=1:length(residualCell)
        residualMat = residualCell{ii};
        err = 0;
        for jj=1:length(residualMat)
            err = err + sum(residualMat{jj}.^2);
        end
        if(err<minErr)
            minErr = err;
            minErrIdx = ii;
        end
    end

    resid{residIdx} = residualCell{minErrIdx};
    residAssocIdxs{residIdx} = residAssocIdxsCell{minErrIdx};
    residIdx = residIdx + 1;
end

end