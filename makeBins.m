function [xBinCenters,yMeans] = makeBins(dData,yData,numBins)
% Bins dData (distances) and returns means of yData in each bin
% dData and yData can be matrices
% Only includes data from positive distances
%-------------------------------------------------------------------------------
if nargin < 3 || isempty(numBins)
    numBins = 10;
end
numThresholds = numBins + 1;

%-------------------------------------------------------------------------------
% Keep d > 0
isPositiveDistance = (dData > 0);
dData = dData(isPositiveDistance);
yData = yData(isPositiveDistance);

%-------------------------------------------------------------------------------
xThresholds = linspace(min(dData),max(dData)+eps,numThresholds);
yMeans = arrayfun(@(x)mean(yData(dData>=xThresholds(x) & dData < xThresholds(x+1))),1:numThresholds-1);
xBinCenters = mean([xThresholds(1:end-1);xThresholds(2:end)],1);

end
