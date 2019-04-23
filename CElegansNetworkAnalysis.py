# Import packages:
import pandas as pd
import os
import networkx as nx

#-------------------------------------------------------------------------------
# Load and process data:
#-------------------------------------------------------------------------------
dataFile = os.path.join('Data','CElegansBinaryAll.csv')
adjMatrix = pd.read_csv(dataFile,header=None)

# Each neuron is a row of this matrix:
dataFile = os.path.join('Data','NeuronNames.txt')
f = open(dataFile,'r')
neuronNames = f.read().split('\n')

# Assign columns


#-------------------------------------------------------------------------------
# We now have a 279 x 279 matrix:
#-------------------------------------------------------------------------------
adjMatrix.shape
