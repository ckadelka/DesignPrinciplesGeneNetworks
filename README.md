# DesignPrinciplesGeneNetworks
A collection of functions and tools to analyze Boolean gene regulatory networks

This repository contains code to perform the analyses described in this paper: [https://arxiv.org/abs/2009.01216](https://arxiv.org/abs/2009.01216). 
Most functionality is also implemented in the interactive project website [https://booleangenenetworks.math.iastate.edu](https://booleangenenetworks.math.iastate.edu).

The three key pieces of code are:
1. load_database13.py,
2. analyse_database13.py, and
3. canalizing_function_toolbox_v13.

# load_database13.py
This file contains all code needed to load the Boolean network models, which are stored as text files in a standardized format as follows:
```text
A = B OR C
B = A OR (C AND D)
C = NOT A
```
