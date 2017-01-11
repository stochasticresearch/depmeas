% Code to verify we merge rectangles properly

function [] = testMergeRectangles()

fprintf('********** TEST 1 **********\n');
rectangles = [0 1 0 1]'; tauVec = 1;
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

fprintf('********** TEST 2 **********\n');
rectangles = [ [0 0.5 0 1]' [0.5 1 0 1]']; 
tauVec = [0.5 0.5];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

fprintf('********** TEST 3 **********\n');
rectangles = [ [0 0.5 0 1]' [0.5 1 0 1]']; 
tauVec = [0.5 -0.5];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

fprintf('********** TEST 4 **********\n');
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 -0.5 0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

fprintf('********** TEST 5 **********\n');
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 0.5 0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

fprintf('********** TEST 6 **********\n');
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 0.5 0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = 7;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 0.5 0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 0.5 -0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 0.5 -0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 -0.5 0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 -0.5 0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 -0.5 -0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [0.5 -0.5 -0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = 7;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 0.5 0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = 7;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 0.5 0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 0.5 -0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 0.5 -0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 -0.5 0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 -0.5 0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 -0.5 -0.3 0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

testNum = testNum + 1;
fprintf('********** TEST %d **********\n', testNum);
rectangles = [ [0 0.25 0 1]' [0.25 0.5 0 1]' [0.5 0.75 0 1]' [0.75 1 0 1]']; 
tauVec = [-0.5 -0.5 -0.3 -0.3];
fprintf('*** Input ***\n');
rectangles
tauVec
rectangles = mergeRectangles(rectangles, tauVec);
fprintf('*** Output ***\n');
rectangles
fprintf('****************************\n');

end

function [rectangles] = mergeRectangles(rectangles, tauVec)

if(size(rectangles,2)>1)
    % merge rectangles w/ same sign of tau_kl, this might increase power of
    % test?
    beginIdx = 1;
    endIdx = 2;
    insertIdx = 1;
    rectanglesMerged = zeros(size(rectangles));
    % TODO: the below loop is inefficient ... look to optimize this
    for ii=1:size(rectangles,2)-1
        taukl_begin = tauVec(beginIdx);
        taukl_end = tauVec(endIdx);
        if(sign(taukl_begin)==sign(taukl_end))
            endIdx = endIdx + 1;
        else
            mergedRectangle = [rectangles(1,beginIdx) rectangles(2,endIdx-1) ...
                               rectangles(3,beginIdx) rectangles(4,endIdx-1)];
            rectanglesMerged(:,insertIdx) = mergedRectangle;
            beginIdx = endIdx;
            endIdx = beginIdx + 1;
            insertIdx = insertIdx + 1;
        end
    end
    % insert the remaining rectangles
    mergedRectangle = [rectangles(1,beginIdx) rectangles(2,endIdx-1) ...
                       rectangles(3,beginIdx) rectangles(4,endIdx-1)];
    rectanglesMerged(:,insertIdx) = mergedRectangle;
    % remove unused entries
    rectanglesMerged(:,insertIdx+1:end) = [];
    rectangles = rectanglesMerged;
end

end