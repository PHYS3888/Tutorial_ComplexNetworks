
#### Task: Reordering an adjacency matrix

```matlab
[~,ix] = sort(C.Pos(:,1),'descend'); % head-to-tail
adjSorted = adjAll(ix,ix);
imagesc(adjSorted)
colormap('gray')
axis('square')
xlabel('Neuron')
ylabel('Neuron')
```
