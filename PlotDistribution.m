function PlotDistribution(xData,numBins)
% Plots a histogram with a given number of bins
if nargin < 2
    numBins = 20;
end

f = figure('color','w');
h = histogram(xData,numBins); % Plot the distribution of `kTot` as a histogram
h.FaceColor = 'w'; % White bars
h.EdgeColor = 'k'; % Black borders
h.LineWidth = 1;

ylabel('Frequency (# neurons)')

end
