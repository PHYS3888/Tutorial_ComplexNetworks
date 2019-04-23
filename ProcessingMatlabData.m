
load(CElegansConnectivityData.mat)

% Write binary connectivity information (chemical + electrical)
csvwrite('CElegansBinaryAll.csv',C.Adj_B{3})

% Write neuron names to a text file:
fid = fopen('NeuronNames.txt','w','n');
numNeurons = length(C.NeuronNames);
for i = 1:numNeurons
    fprintf(fid,'%s\n',C.NeuronNames{i});
end
fclose(fid)
