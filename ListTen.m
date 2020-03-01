function ListTen(neuronNames,kTot,ix)

for i = 1:10
    fprintf(1,'Neuron %s: k = %u\n',neuronNames{ix(i)},kTot(ix(i)))
end

end
