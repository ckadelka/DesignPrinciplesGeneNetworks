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

def record_loops(network):
    ret_length = []
    ret_type = []
    F = network[0]
    I = network[1]
    D = network[2]
    loops = get_loops(network)
    for loop in loops:
        l_length = len(loop)
        l_type = can.get_type_of_loop(loop,F,I)
        type_num = can.get_loop_type_number(l_type)
        ret_length.append(l_length)
        ret_type.append(type_num)
    return ret_length,ret_type

def record_attractors(F, I, N):
    ret = []
    attractors = can.num_of_attractors_v2(F, I, N)[0]
    for attractor in attractors:
        ret.append(len(attractor))
    return ret

def avg_loop_len(loops, types):
    counter = 0
    lengths = 0
    for i in range(len(loops)):
        if types[i] != -1:
            counter += 1
            lengths += loops[i]
    if counter != 0:
        avg_len = lengths / counter
    else:
        avg_len = 0

    return avg_len

def percent_FBLs_pos(network):
    F = network[0]
    I = network[1]
    D = network[2]
    num_pos = 0
    num_neg = 0
    loops = get_loops(network)
    for loop in loops:
        l_type = can.get_type_of_loop(loop, F, I)
        if can.is_pos_loop(l_type) == True:
            num_pos += 1
        else:
            num_neg += 1
    if num_pos and num_neg != 0:
        ratio = num_pos / (num_pos + num_neg)
    else:
        ratio = -1
    return ratio

