## Outputs:

![Adjacency matrix](figs/adjMatBinary.png)

![Directed graph](figs/Directedgraph.png)

![Degree distribution](figs/degreeDistribution.png)

## Q1: Reordering an adjacency matrix

```matlab
[~,ix] = sort(positionXY(:,1),'descend'); % head-to-tail
adjSorted = adjMatrix(ix,ix);
imagesc(adjSorted)
colormap('gray')
axis('square')
xlabel('Source neuron')
ylabel('Target neuron')
```

## Q2: Distinguishing head, body, and tail neurons

```matlab
figure('color','w')
G.Nodes.NodeColors = GiveMeNeuronLabels;
p = plot(G,'XData',positionXY(:,1),'YData',positionXY(:,2));
p.NodeCData = G.Nodes.NodeColors;
p.MarkerSize = 6;
colormap([233, 163, 201;
            10, 10, 10;
            161, 215, 106]/255)
cB = colorbar();
cB.Ticks = [1,2,3];
cB.TickLabels = {'head','body','tail'};
```
![Nodes colored by type](figs/nodesColored.png)

## Computing network properties

How many neurons?:
```matlab
size(adjMatrix)
length(neuronNames)
```

How many edges?:
```matlab
sum(adjMatrix(:))
```

__Q__:
The degree distribution contains a heavy tail with some strongly connected neurons at the upper end of connectivity. These high-degree neurons are known as 'hubs'.


__Q__: Use the `isBodyNeuron` indicator to reduce the full adjacency matrix down to include information about body neurons only, as a new adjacency matrix, `adjMatrixBody`. Upload your code.

```matlab
neuronLabels = GiveMeNeuronLabels();
isBodyNeuron = (neuronLabels==2);
adjMatrixBody = adjMatrix(isBodyNeuron,isBodyNeuron);
```


__Q__ Distance matrix:
```matlab
distMatrix = squareform(pdist(positionXY));
distMatrixBody = distMatrix(isBodyNeuron,isBodyNeuron);
```
