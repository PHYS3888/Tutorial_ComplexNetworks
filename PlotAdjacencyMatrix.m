function PlotAdjacencyMatrix(adjMatrix)

figure('color','w')
imagesc(adjMatrix)
colormap('gray')
axis('square')
xlabel('Target neuron')
ylabel('Source neuron')
title('Neuronal connectivity in the nematode worm')

end
