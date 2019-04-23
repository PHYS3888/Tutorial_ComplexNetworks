function neuronLabels = GiveMeNeuronLabels(makeCategorical)

if nargin < 1
    makeCategorical = false;
end

% Load in data:
dataFile = fullfile('Data','CElegansConnectivityData.mat');
load(dataFile,'C')

% Logical for each type:
isHead = logical(C.RegionM(:,strcmp(C.neuronAnatomyNames,'head neuron')));
isTail = logical(C.RegionM(:,strcmp(C.neuronAnatomyNames,'tail neuron')));
isBody = logical(~isHead & ~isTail);

% Form categorical
numNeurons = C.numNeurons;
neuronLabels = zeros(numNeurons,1);
neuronLabels(isHead) = 1;
neuronLabels(isBody) = 2;
neuronLabels(isTail) = 3;

if makeCategorical
    neuronLabels = categorical(neuronLabels,1:3,{'head','body','tail'});
end

end
