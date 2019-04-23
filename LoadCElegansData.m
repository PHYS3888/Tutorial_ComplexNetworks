function [adjAll,neuronNames,positionXY] = LoadCElegansData()
% Script to load in C. elegans data

dataFile = fullfile('Data','CElegansConnectivityData.mat');
load(dataFile,'C')
adjAll = C.Adj_B{3}; % this is the relevant data
neuronNames = C.NeuronNames;
positionXY = -C.Pos; % flip the sign by convention

end
