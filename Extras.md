
## Symmetrizing a graph

The _C. elegans_ connectome contains directed information---not all connections go in both directions.
Because of this, some connections exist i -> j, for which the reverse connection, j -> i, does not exist.
Thus, the matrix, `adjMatrix`, is not symmetric.

We can form a symmetrized version of the network (an undirected graph) efficiently by keeping an edge when there is either connection i->j OR j->i using the OR operator `|` with the transpose of the adjacency matrix:

```matlab
adjMatrixSym = (adjMatrix | adjMatrix');
```

We can plot the symmetrized graph:
```matlab
G = graph(adjMatrixSym); % construct a graph object
p = plot(G); % plot the graph
% Adjust some of the plotting properties:
p.Marker = 'o';
p.NodeColor = 'r';
p.MarkerSize = 8;
```


## Shortest paths

Being highly connected, many shortest paths travel through hub neurons (assuming that shortest paths are a relevant communication mechanism for this system).
