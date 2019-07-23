function [xBinCenters,yMeans] = makeBins(xData,yData,numBins)
% Bins xData and returns means of yData in each bin
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(numBins)
    numBins = 10;
end
numThresholds = numBins + 1;

%-------------------------------------------------------------------------------
% Filter out NaNs:
goodBoth = (~isnan(xData) & ~isnan(yData));
if ~any(goodBoth)
    error('No good data');
elseif any(~goodBoth)
    xData = xData(goodBoth);
    yData = yData(goodBoth);
    fprintf(1,'Removed %u bad samples from x/y data\n',sum(~goodBoth));
end

xThresholds = linspace(min(xData),max(xData)+eps,numThresholds);
yMeans = arrayfun(@(x)mean(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);
xBinCenters = mean([xThresholds(1:end-1);xThresholds(2:end)],1);

end
