import canalizing_function_toolbox_v1_9 as can
import networkx as nx

max_loop = 6

def get_G(network):
    F = network[0]
    I = network[1]
    edges = []
    for j,regulators in enumerate(I):
        n_variables = can.n_from_f(F[j])
        if j>=n_variables: #exclude constant self-loops
            break
        for ii in regulators:
            edges.append((ii,j))
    G=nx.DiGraph(edges)
    return G

def get_loops(network):
    G = get_G(network)
    loops = list(can.simple_cycles(G,max_loop))
    return loops
