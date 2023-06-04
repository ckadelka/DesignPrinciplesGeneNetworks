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

# canalizing_function_toolbox_v13
This file contains a variety of functions to analyze Boolean functions and Boolean networks. Each function has its own documentation. 

# analyse_database13.py
This file combines everything. It loads the database and analyzes it. Most functions defined in this file are used to generate the plots in [https://arxiv.org/abs/2009.01216](https://arxiv.org/abs/2009.01216).
 

