# Complex Networks Analysis of Neuronal Connectivity in _C. Elegans_ :bug:

This tutorial will walk you through a basic understanding of complex networks.

Networks are a representation of objects (nodes) and the connections between pairs of nodes (edges).

As an example dataset we will use the network of neurons, and their connectivity in the nematode worm, _C. elegans_ :bug:, measured through painstaking reconstruction from electron microscopy.
Note that the network of connectivity between elements of the brain is called a _connectome_.

![](figs/CElegans.png)

This data has been analyzed by many physicists (and other scientists).
This tutorial follows results from [this recent paper](http://dx.plos.org/10.1371/journal.pcbi.1005989).

## PRE-WORK

In this section you can brush up on a few basic commands for reordering and subsetting vectors/matrices in Matlab.

Let's start with a 5x5 matrix:
```matlab
M = magic(5);
```

### 0.1 Reordering

Let's reorder the rows according to the ordering defined by: `[1,4,3,5,2]`:

```matlab
ix = [1,4,3,5,2];
M_row = M(ix,:);
```
Verify how the order of the rows of `M` have been switched using the `ix` permutation.
We could do the same to the columns as `M_col = M(:,ix);`.
And we can reorder both together as `M_both = M(ix,ix);`

### 0.2 Subsetting
What if we want to keep only a subset of rows/columns of `M`?
One way is to define a logical indicator for the rows we want to keep:

```matlab
keepMe = [false,false,true,true,false];
```

Then we can keep just these rows, as `M(keepMe,:)`, just the columns, as `M(:,keepMe)`.
We can do both all at once as `M(keepMe,keepMe)`.

Note that the same results would be obtained if we instead defined the indices we want to keep: `keepMe = [3,4]`.

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

#### Plotting a sorted adjacency matrix

In their paper, Arnatkeviciute et al. plotted the same data in Figure 1A, shown here:

![Arnatkeviciute et al.](figs/Arnatkeviciute.png)

If it is the exact same data, then why does it look so different to what we plotted?

Plots can be deceiving---sometimes random data can be plotted in a way that appears to the human eye to contain non-random structure.
On closer inspection, we find that Arnatkeviciute  et al. (2018) ordered their matrix by the position from head-to-tail (called the 'anterior-posterior' axis).
The spatial co-ordinates of each neuron is in the variable `positionXY`.
The first column of `positionXY` is the `x`-coordinates (broadly from head to tail), and the second column contains the `y`-coordinates.

#### :question::question::question: Q1: Sorting nodes by location

Can you reorder the network so that neurons are ordered according to their position from head-to-tail and verify that the result matches the result from Arnatkeviciute et al. (2018) (ignoring coloring)?

Upload these lines of code to Canvas.
_Hint:_ the `sort` function is relevant to get the relevant permutation (you can read about this function using `doc sort`).

Note that the remainder of this tutorial will work with the original (unordered) matrix, `adjMatrix`.

### Plotting a network in physical space

Let's visualize the same information as a graph.

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
p.Marker = 'o'; % plot each neuron as a circle
p.NodeColor = 'r'; % make circles red
p.MarkerSize = 6; % make circles size 6
axis('equal') % make horizontal and vertical scales comparable
```

You can zoom in to explore the cluster of head neurons to the left of the plot, and the cluster of tail neurons towards the right of the plot.
Body neurons are scattered through the length of the worm.


#### :question::question::question: Q2: Distinguishing head, body, and tail neurons

Adjust the plot above to color head, body, and tail neurons a different color?
Once you have a nice plot, upload the code you used to generate it.

_Hint:_ You can retrieve the labeling of neurons using the function `GiveMeNeuronLabels` (head = 1, body = 2, tail = 3).
You can set node colors by setting the `p.NodeCData` property to a given set of numeric labels.

---

## Part 2: Computing network properties

We can ask many basic questions by running simple operations on the adjacency matrix.

#### How many neurons?
Use the `size()` function to determine how many neurons there are.

#### Is the network binary or weighted?

Recall the difference between a binary and a weighted network.
Do a simple test to determine whether the network captured in `adjMatrix` is binary or weighted.
(_Hint:_ one way is to use the `unique()` function).


#### Is the network directed or undirected?

Recall the difference between an undirected and a directed network.
Verify that the code below tests for directedness:

```matlab
symmetricMatches = (adjMatrix'==adjMatrix);
```
What should you look for in the `symmetricMatches` matrix?

#### Does the connectome contain self-connections?

Where do self-connections show up in the adjacency matrix?
Use the `diag()` function to determine whether `adjMatrix` contains self-connections.

#### How many edges are in the connectome?

Use the `sum()` command to count the total number of edges.
Do you need to divide this number by 2?
Why/why not?

#### In-degree, `kIn`, and out-degree, `kOut`

How many inward-coming connections does each neuron have, its _in-degree_, `kIn`?
What about the total number of outward-going connections, its _out-degree_, `kOut`?

Since each neuron is a row and a column, we can compute degree using sums.
This adjacency matrix has source neurons as rows, and target neurons as columns.
Thus we can compute the in-degree, `kIn`, as the number of sources (rows) coming into a given target (column) by summing down each column:
```matlab
kIn = sum(adjMatrix,1); % dimension 1: sum down columns
```

Similarly for `kOut`, as the number of targets (columns) for a given source (row):
```matlab
kOut = sum(adjMatrix,2)'; % dimension 2: sum across rows
```

The total number of connections involving a neuron (both outgoing and incoming) can be computed as the sum of these two quantities, the total degree, `kTot`:

```matlab
kTot = kIn + kOut;
```

Recall that in a random network, there is a tight distribution about the mean degree.
If neurons connect at random, this would mean that most neurons will have a similar number of connections, with a tight spread around a mean value (a Binomial/Poisson distribution).
We can test this by plotting the degree distribution to see how connectivity is distributed across neurons:
```matlab
f = figure('color','w');
numBins = 20;
h = histogram(kTot,numBins); % Plot the distribution of `kTot` as a histogram
h.FaceColor = 'w'; % White bars
h.EdgeColor = 'k'; % Black borders
h.LineWidth = 1;
xlabel('Total degree, kTot')
ylabel('Frequency (# neurons)')
```

__What about this distribution tells us that there are highly connected hubs in the _C. elegans_ connectome?__

#### :question::question::question: What do the hub neurons do?

The most complex behavior in _C. elegans_ is its locomotion, which is governed by a set of ten _"command interneurons"_, which control both forward (neurons: AVBL, AVBR, PVCL, PVCR) and backward (neurons: AVAL, AVAR, AVDL, AVDR, AVEL, AVER).
(If you're interested you can [read more here](https://www.frontiersin.org/articles/10.3389/fncom.2013.00128/full)).

I wonder if any of these show up in our list of highly-connected neurons...
Let's list the top ten to see:
```matlab
% Sort neurons from the most to the least connected
[~,ix] = sort(kTot,'descend');
for i = 1:10
    fprintf(1,'%s, k = %u\n',neuronNames{ix(i)},kTot(ix(i)))
end
```

Is there overlap between the neurons that control the worm's locomotion and the neurons that are most strongly connected in the network?
In the Canvas quiz, select all of the overlapping neurons.

---

## Part 3: Physical embedding

The worm's nervous system is a physical system, across a head, body, and tail.
Representing it as nodes and edges does not capture this physical embedding.
As a physicist, you might be interested in characterizing how the connectivity in _C. elegans_ is physically embedded.

Take a look at Fig. 3A in Arnatkeviciute et al. (2018), shown below, which reveals a decrease in connection probability with the distance between pairs of neurons.
Note that this relationship is clearest for connections between body neurons.

![Connection Probability](figs/ConnectionProbability.png)

Let's try to reproduce this distance-dependence of body-body connectivity in the _C. elegans_ nervous system.

#### :question::question::question: A body connectome

Our first step is to filter down to an adjacency matrix containing only body neurons.

```matlab
neuronLabels = GiveMeNeuronLabels(); % Label neurons in the body as '2'
isBodyNeuron = (neuronLabels==2); % Construct a binary indicator for body neurons
```

How many connections are collectively made between the body neurons of _C. elegans_?

What is the total number of connections made _from_ a body neuron to a head neuron?

Use the `isBodyNeuron` indicator to reduce the full adjacency matrix down to include information about body neurons only, as a new adjacency matrix, `adjMatrixBody`.
Upload your code.

#### Visualize body-body neuron connectivity
Visualize interconnectivity between the worm's body neurons using `imagesc()` (as we did for the full network above).
Visually estimate the probability that if a pair of neurons are connected, that this connection is reciprocal?

Assess your visual estimation by running the function `whatProportionReciprocal`, using `adjMatrixBody` as input.

How close was your estimate?

#### :question::question::question: Converting coordinates to Euclidean distances

Ok, so now we have the connectivity information for body neurons, now we need their separation distances.

Convert coordinates in `positionXY` into Euclidean distances, storing the result in a new variable, `distMatrix`.
Filter down to `distMatrixBody` using the `isBodyNeuron` indicator you created above.
_Hint_: `pdist` is a relevant function for computing Euclidean distances (see also the `squareform` function to convert between vector and matrix representations).
Upload your code.

### Binning
Ok, so now we have a distance for every connection between pairs of body neurons, and whether that connection exists or not.

We want to include pairwise connections in both directions (i.e., both the upper triangle and the lower triangle of our connection matrix, `adjMatrixBody`).

Calculating probabilities are easier for binary data, because we can compute the probability of being a `1` as the mean of a binary vector.
For example, note that `mean([1,1,0,0]) = 0.5`, or `mean([1,1,1,1,0]) = 0.8`.
We use this property to compute the connection probability of the mean of the binary connection indicators.
Run the code below, which makes 10 equally-spaced bins (`numBins = 10`) through the range of body-body distances, and, for pairs of neurons in that distance range, computes the probability that connections exist between them.
Note that diagonal entries are self-connections and should be excluded:

#### Step 1: Exclude diagonal entries
We need to make sure that we don't include self-connections in our computations.

Check out this indicator matrix: `logical(eye(size(distMatrixBody)))`

Use this to turn `distMatrixBody` into a vector, `distDataBody`, that contains only non-diagonal entries.
Do the same for `connDataBody` (as the non-diagonal data from `adjMatrixBody`).

#### Step 2: Compute probability across equally spaced distance bins
Now we want to compute the probability that a connection exists in each of 10 distance bins.
Set `numBins` and run `makeBins` as below to do this.

```matlab
% The makeBins function hides the dirty work:
[distBinCenters,connProb] = makeBins(distDataBody,connDataBody,numBins);
```

#### :question::question::question: Connection probability as a function of distance
Plot the connection probability, `connProb`, as a function of distance in the body of the worm _C. elegans_.

```matlab
f = figure('color','w');
plot(distBinCenters,connProb,'o-k')
xlabel('Separation distance (mm)')
ylabel('Connection probability')
```

How does connection probability depend on separation distance in _C. elegans_?
Upload your plot.

---

#### :fire: Optional challenge :fire:: connection probability

Challenge: Repeat for the head -> tail, and tail -> head results.
