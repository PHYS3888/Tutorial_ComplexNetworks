# Complex Networks Analysis of Neuronal Connectivity in _C. Elegans_

This tutorial will walk you through a basic understanding of complex networks.

Networks are a representation of objects (nodes) and the connections between pairs of nodes (edges).

As an example dataset we will use the network of neurons, and their connectivity in the nematode worm, _C. elegans_ :bug:.

![](figs/CElegans.png)

This data has been analyzed by many physicists (and other scientists).
This tutorial follows results from [this recent paper](http://dx.plos.org/10.1371/journal.pcbi.1005989).

---

## Part 1: Representing and visualizing networks

First load in the _C. elegans_ connectivity data using the function `LoadCElegansData`:
```matlab
[adjMatrix,neuronNames,positionXY] = LoadCElegansData();
```

The variable `adjMatrix` represents the adjacency matrix of the network.
This contains information about all connections from every neuron (row) to every other neuron (column).

### Plotting an adjacency matrix

When analyzing any type of data, your habit should be to start by getting a good visualization of it.
Starting an analysis with a comprehensive visualization can help identify any issues with the data, and can motivate the most suitable types of analysis to perform on it.

Let's start with plotting the adjacency matrix as an image:
```matlab
figure('color','w')
imagesc(adjMatrix)
colormap('gray')
axis('square')
xlabel('Target neuron')
ylabel('Source neuron')
title('Neuronal connectivity in the nematode worm')
```

This is the adjacency matrix representation of the _C. elegans_ connectome, where every neuron is a row (and a column), and edges represent the complex connectivity patterns between pairs of neurons.

---

#### Q1: Sorted adjacency matrix

In their paper, Arnatkeviciute et al. plotted the same data in Figure 1A, shown here:

![Arnatkeviciute et al.](figs/Arnatkeviciute.png)

If it is the exact same data, then why does it look so different to what we plotted?

Plots can be deceiving---sometimes random data can be plotted in a way that appears to the human eye to contain non-random structure.
On closer inspection, Arnatkeviciute ordered their matrix by the position from head-to-tail.
The spatial co-ordinates of each neuron is in the variable `positionXY`.
Can you reorder the network so that neurons are ordered according to their position from head-to-tail and verify that the result matches the result from Arnatkeviciute et al. (ignoring coloring)?

_**Q1**: Upload the lines of code you used to plot the sorted network that matches Arnatkeviciute et al._


### Plotting a network in space

What about visualizing the same information as a graph?

```matlab
G = digraph(adjMatrix); % construct a graph object
p = plot(G); % plot the graph
```

Do you know what you're looking at?! Zoom in and have a play.

Note that this visualization is in an abstract space, that is, the coordinates don't correspond to physical space.
But we have two-dimensional coordinates for every neuron, `positionXY`, that we can use to plot the network information in physical space:

```matlab
figure('color','w')
p = plot(G,'XData',positionXY(:,1),'YData',positionXY(:,2));
axis('equal') % scale the x and y axes comparably
p.Marker = 'o'; % make each neuron a circle
p.NodeColor = 'r'; % make circles red
p.MarkerSize = 6; % make circles size 6
```

You can zoom in to explore the cluster of head neurons to the left of the plot, and the cluster of tail neurons towards the right of the plot.
Body neurons are scattered through the length of the worm.

---

#### Q2: Distinguishing head, body, and tail neurons

Adjust the plot above to color head, body, and tail neurons a different color?
Once you have a nice plot, upload the code you used to generate it.

_Note:_ You can retrieve the labeling of neurons as a categorical using the function `GiveMeNeuronLabels` (head = 1, body = 2, tail = 3).
You can set node colors using `p.NodeCData = numericLabels`, for a given set of numeric labels.

---

## Part 2: Computing network properties

We can ask many basic questions by running simple operations on the adjacency matrix.

#### How many neurons?
Use the `size()` function to determine how many neurons there are.

#### Is the network binary or weighted?
Recall the difference between a binary and a weighted network.
Do a simple test to determine whether the network captured in `adjMatrix` is binary or weighted.
(_Hint:_ one way is to use the `unique()` function).

:question: Is the network binary or weighted?

#### Is the network directed or undirected?
Recall the difference between an undirected and a directed network.
Verify that the below tests for this:

```matlab
symmetricMatches = (adjMatrix'==adjMatrix);
```
What should you look for in the `symmetricMatches` matrix?

:question: Is the matrix directed or undirected?

#### Are there any self-connections?

Where do self-connections show up in the adjacency matrix?
Use the `diag()` function to determine whether there `adjMatrix` contains self-connections.

:question: Does the connectome contain self-connections?

#### How many edges?

Use the `sum()` command to count the total number of edges.

:question: How many edges are there?

#### In-degree and Out-degree

How many inward-coming connections does each neuron have, its _in-degree_, `kIn`?
What about the total number of outward-going connections, its _out-degree_, `kOut`?

Since each neuron is a row and a column, we can compute degree using sums.
This adjacency matrix has source neurons as rows, and target neurons as columns.
Thus we can compute the in-degree, `kIn`, as the number of sources (rows) coming into a given target (column) by summing down each column:
```matlab
kIn = sum(adjMatrix,1);
```

Similarly for `kOut`, as the number of targets (columns) for a given source:
```matlab
kOut = sum(adjMatrix,2)';
```

The total number of connections involving a neuron can be computed as the sum of these two quantities, the total degree, `kTot` (the total number of connections, in and out, from a neuron):

```matlab
kTot = kIn + kOut;
```

Do all neurons have a similar number of connections, with some spread around a mean value (a Gaussian distribution), or are there some highly-connected hub neurons?
To test this, we can plot the degree distribution to see how connectivity is distributed across neurons:
```matlab
f = figure('color','w');
h = histogram(kTot);
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth = 1;
xlabel('Total degree, kTot')
ylabel('Frequency (# neurons)')
```

:question: Does this distribution provide evidence for the existence of hubs? Why?

#### Which neurons are the hubs?

The most complex behavior in _C. elegans_ is its locomotion, which is governed by a set of ten _"command interneurons"_, which control both forward (neurons: AVBL, AVBR, PVCL, PVCR) and backward (neurons: AVAL, AVAR, AVDL, AVDR, AVEL, AVER).
(If you're interested you can [read more here](https://www.frontiersin.org/articles/10.3389/fncom.2013.00128/full)).

I wonder if any of these show up in our list of highly-connected neurons.
Let's list the top ten to see:
```matlab
[~,ix] = sort(kTot,'descend');
for i = 1:10
    fprintf(1,'%s, k = %u\n',neuronNames{ix(i)},kTot(ix(i)))
end
```

:question: Is there overlap between the complex locomotion behavior of _C. elegans_ and the neurons that exhibit the strongest connectivity?__

---

## Part 3: Physical embedding

Organisms are physical systems.
As a physicist, you might be interested in characterizing how the connectivity in _C. elegans_ is physically embedded.

Figure 3A in Arnatkeviciute et al. (2018) is reproduced below, and shows a decay in connection probability as a function of distance that is clearest for interconnections between body neurons.

![Connection Probability](figs/ConnectionProbability.png)

Let's try to reproduce the distance-dependence of body-body worm-brain connectivity.

### A body connectome

The first step is to filter down to an adjacency matrix containing only body neurons.

```matlab
neuronLabels = GiveMeNeuronLabels();
isBodyNeuron = (neuronLabels==2);
```

:question: Use the `isBodyNeuron` indicator to reduce the full adjacency matrix down to include information about body neurons only, as a new adjacency matrix, `adjMatrixBody`. Upload your code.

Visualize the body neurons, as we did for the full matrix above (using `imagesc()`).
Visually estimate the probability that if a pair of neurons are connected, that this connection is reciprocal?
Verify your intuition by running the following code:
```matlab
% There is a connection in either direction:
adjMatrixBodyConnEither = (adjMatrixBody==1 | adjMatrixBody'==1);
% There is a connection in both directions:
adjMatrixBodyConnBoth = (adjMatrixBody==1 & adjMatrixBody'==1);
adjMatrixBodyConnEitherUpper = adjMatrixBodyConnEither(triu(true(size(adjMatrixBody))));
adjMatrixBodyConnBothUpper = adjMatrixBodyConnBoth(triu(true(size(adjMatrixBody))));
pRecip = mean(adjMatrixBodyConnBothUpper(adjMatrixBodyConnEitherUpper==1));
fprintf(1,'%.1f%% reciprocal connectivity rate\n',pRecip*100);
```

### Euclidean distances

Ok, so we have connectivity information, now we need distance information.

:question: Convert coordinates in `positionXY` into Euclidean distances, storing the result in a new variable, `distMatrix`.
Filter down to `distMatrixBody` using the `isBodyNeuron` indicator you created above.
_Hint_: `pdist` is a relevant function for computing Euclidean distances.
Upload your code.

### Binning
Ok, so now we have a distance for every connection between pairs of body neurons, and whether that connection exists or not.

We want to include pairwise connections in both directions (i.e., both the upper triangle and the lower triangle of our connection matrix, `adjMatrixBody`).

Calculating probabilities are easier for binary data, because we can compute the probability of being a `1` as the mean of a binary vector.
For example, note that `mean([1,1,0,0]) = 0.5`, or `mean([1,1,1,1,0]) = 0.8`.
We use this property to compute the connection probability of the mean of the binary connection indicators.
Run the code below, which makes 10 equally-spaced bins (`numBins = 10`) through the range of body-body distances, and, for pairs of neurons in that distance range, computes the probability that connections exist between them.
Note that diagonal entries are self-connections and should be excluded:

```matlab
% Important not to include the diagonal entries (self-connections)
notDiag = ~logical(eye(size(distMatrixBody)));
distDataBody = distMatrixBody(notDiag);
connDataBody = adjMatrixBody(notDiag);
numBins = 10;
% Use the makeBins function to hide the dirty work:
[distBinCenters,connProb] = makeBins(distDataBody,connDataBody,numBins);
```

Plot the connection probability, `connProb`, as a function of distance in the body of the worm _C. elegans_.

```matlab
f = figure('color','w');
plot(distBinCenters,connProb,'o-k')
xlabel('Separation distance (mm)')
ylabel('Connection probability')
```

:question: Does connection probability depend on separation distance in _C. elegans_?

---

#### :fire: Challenge :fire:: Connection probability

Challenge: Repeat for the head -> tail, and tail -> head results.

---
