function [rectMatOut] = disambiguateRectangleMat(rectMatIn)

allIdxs = 1:4;

numRects = size(rectMatIn,2);
if(mod(numRects,2)==0)
    rectDiff = diff(rectMatIn,1,2);
    % find "groupings" of rectangles
    orientationIdxs = find(rectDiff(:,1)==0)';
    checkIdxs = setdiff(allIdxs,orientationIdxs);
    if(diff(orientationIdxs)==1)
        breakPts = find(rectDiff(orientationIdxs(1),:)~=0);
        if(isempty(breakPts))
            rectMatOut = rectMatIn;
        else
            % merge greedily and recursively
            mergeStartIdx = 1;
            g1StartIdx = mergeStartIdx;
            g1StopIdx = breakPts(mergeStartIdx);
            g2StartIdx = breakPts(mergeStartIdx)+1;
            if(length(breakPts)==1)
                g2StopIdx = numRects;
            else
                g2StopIdx = breakPts(mergeStartIdx+1);
            end
            if(isequal(rectMatIn(checkIdxs(1),g1StartIdx:g1StopIdx),...
                       rectMatIn(checkIdxs(1),g2StartIdx:g2StopIdx)) && ...
               isequal(rectMatIn(checkIdxs(2),g1StartIdx:g1StopIdx),...
                       rectMatIn(checkIdxs(2),g2StartIdx:g2StopIdx)) )
                   % merge
                   
            else
                
            end
            
            % recurse to get all merges
        end
    else
        rectMatOut = rectMatIn;
    end
    
else
    % for now, this means that there is no disambiguation that can be done
    % ...
    rectMatOut = rectMatIn;
end

end