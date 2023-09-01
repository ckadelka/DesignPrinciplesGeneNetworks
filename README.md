# DesignPrinciplesGeneNetworks
A collection of functions and tools to analyze Boolean gene regulatory networks

This repository contains code to perform the analyses described in this paper: [https://arxiv.org/abs/2009.01216](https://arxiv.org/abs/2009.01216). 
Most functionality is also implemented in the interactive project website [https://booleangenenetworks.math.iastate.edu](https://booleangenenetworks.math.iastate.edu).

The three key pieces of code are:
1. load_database13.py,
2. analyse_database13.py, and
3. canalizing_function_toolbox_v13.

# load_database13.py
This file contains all code needed to load the Boolean network models via load_database(folders). The models are stored in a list of folders, as text files in a standardized format:
```text
A = B OR C
B = A OR (C AND D)
C = NOT A
```
This little example represents a model with three genes, A, B and C, and one external parameter D (which only appears on the right side of the equations).

# canalizing_function_toolbox_v13.py
This file contains a variety of functions to analyze Boolean functions and Boolean networks. Each Python function has its own documentation. A Boolean function is considered as a list of 0-1 entries of length 2^n where n is the number of inputs. A Boolean network of N nodes is a list of N Boolean functions. For example,
```python
f_A = [0,1,1,1]
f_B = [0,0,0,1,1,1,1,1]
f_C = [1,0]
F = [f_A,f_B,f_C]
```
describes the Boolean network from above. One can also get this via
```python
import load_database13 as db

with open('example.txt', 'w') as writer:
    writer.write('A = B OR C\nB = A OR (C AND D)\nC = NOT A')
F, I, degree, variables, constants = db.text_to_BN(folder='',textfile='example.txt')
```
which yields in addition the adjacency matrix I, the in-degree of each node, the names of the variables (in order) and the names of potential external parameters.

# analyse_database13.py
This file combines everything. It loads the database and analyzes it. Most Python functions defined in this file are used to generate the plots in [https://arxiv.org/abs/2009.01216](https://arxiv.org/abs/2009.01216), which will soon be published in updated format but with the same title in a journal.
 

