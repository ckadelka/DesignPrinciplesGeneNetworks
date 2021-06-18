0) Basics
tobin(x):
	returns the binary representation (in array form) of a decimal number
 param: double x the decimal number
 return: int[] binary representation array of the number 

dec2bin(x, n=[]):
 returns the binary representation (in array form) of a decimal number.
 Input can be an array itself, in which case each decimal number in the 
 array is separately converted into a binary array.
 The second input n, if specified, describes the number of positions in the
 binary representation array for each number.
 Example: dec2bin(10)=dec2bin(10,4)=[1,0,1,0]
 dec2bin(10,6)=[0,0,1,0,1,0]

 param: double x decimal number OR double[] x array of decimals
 param: (optional) int[] specifies number of positions in binary representation of array for each number 
 return: int[] binary representation array of the number

dim(A): opt
 param: A array 
 return: a tuple of ints representing the size of each dimension of A

find_all_indices(array, el):
 Find all indicies at which el occurs in array
 param: array an array
 param: el a potential element of the array
 return: a list of ints representing the indicies of array at which el occurs

edgelist_to_I(edgelist): 
 Converts edgelist to a tuple containing an adjacency list and a list of which node indexes the adjacency entries correspond to 
 param: edgelist a list where each entry is of length 2, with the first entry the regulator node and the second is the regulated 
 return: A tuple where the first entry is an adjacency matrix linking nodes to a list of things they are regulated by. The first entry of this tuple
 is the adjacency list. The second indicates which node each adjacency list entry corresponds to.
 Example: input: [[0,1],[4,1]] output: ([[], [0, 2], []], [0, 1, 4]))

f_from_expression(expr):
 Converts a string-represented expression to a tuple containing two arrays. If the input is a mathematical expression with numbers and operators (*, /, +, -), 
 the first element of this tuple is an array containing the output of this expression when evaluated, and the second is empty. If the expression is in terms of 
 variables and logical operators (and, or, not, etc.), the first entry of the tuple is an array containing 0 and 1 corresponding to the output of the statement 
 given the presumed variables are true/false. 
 param: expr a string representing an expression 
 return: a tuple 'F' containing the truth table for input variables, and the names of those input variables. The first entry in this tuple represents
	a function.
 Example: input: "x and (y or not y)" output: ([0, 0, 1, 1], ['x', 'y']) 
 How this example is evaluated: [x false y false, x false y true, x true y false, x true y true]

1) Methods to analyze Boolean functions
is_degenerated(F):?
 determines whether tuple containing function 'F' is degenerated. A function is degenerated if it contains any variables which can have
		no impact on the output of the function for all possible input combinations to other variables of the function.
 param: tuple F or list F- a representation of a function. This can be either the entire tuple returned by f_from_expression(expr), or just the 
		list representing the function.
 return: boolean - True if F is degenerated, False otherwise

nr_essential_variables(F):
 returns the number of essential variables in function representation F
 param: tuple F or list F- a representation of a function. This can be either the entire tuple returned by f_from_expression(expr), or just the 
		list representing the function.
 return: int - the number of essential variables in F 
 (a variable is essential if its status as True or False can potentially impact the value of the expression)
 Example: input: f_from_expression("x and (y or not y)") output: 1

is_constant(F):
 determine if the function always evaluates to True or False for any inputs
 param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
 return: boolean - true if the output of F is independent of input variables

is_monotonic(F, GET_DETAILS=False): opt
 determine if function F contains only essential variables. If GET_DETAILS is specified as true, the effect of each variable 
 ("increasing", "decreasing", "not essential", "not monotonic") is returned in a list as the second element of a tuple,
 with the first tuple element being whether F contains only essential variables.
 param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
 return: (boolean,list) - true if F is monotonic, false if F is not monotonic. list contains string for effect of each variable if GET_DETAILS
 was specified. Only the boolean is returned otherwise.

get_symmetry_groups(F, bool_list=[]): 
	Get a list of the symmetry groups in F. Two variables of Fare in the same symmetry group if they could be substituted for eachother in
		the function and not change output.
	param: F - a list representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: bool_list - Pre-computed truth table for the function
	return: a 2d list indicating symmetry groups. Each entry in the outer list is a symmetry group, which lists the integers corresponding to the 
		variables in that symmetry group.
	Example:
		input: can.f_from_expression("x and y or z")[0]
		output: [[0, 1], [2]]

is_canalizing(F,n): opt - why require n passed?
	determine if function F is canalizing. n must be the number of variables in F. A function is canalizing if it has a variable which has 
		a state (true or false) that determines the evaluation of the function regardless of the other variables.
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int n - the number of variables in F

is_collectively_canalizing(F, k, n): 
	determine if function F is collectively canalizing amongst k variables. F is collectively canalizing for k variables if setting k variables
		a certain way (true or false) will always gaurantee a constant output (true or false) of the function, regardless of the other variables. 
		If F has no variables, it will always be considered collectively canalizing for k=0.
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int k - The number of variables to check if F is canalizing for. k must be greater than or equal to 0 and less than or equal to n.
	param: int n - the number of variables in F
	return: boolean - True if F is collectively canalizing amongst k variables, False otherwise
	Examples:
		input: (f_from_expression("(x and y and z) or (not x and y and z) or (x and not y and z) or (x and y and not z)")[0], 3, 3) output: True
		input: (f_from_expression("(x and y and z) or (not x and y and z) or (x and not y and z) or (x and y and not z)")[0], 2, 3) output: False
	
get_proportion_of_collectively_canalizing_input_sets(F, k, n, bool_list=[]):
	determine what proportion of input sets of length n-k collectively canalize f.
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int k - the number of variables for F to be canalizing for
	param: int n - the number of variables in F
	param: list bool_list - Pre-computed truth table for the function
	return: float the k-set canalizing proportion of F

get_canalizing_strength(F,bool_list=[]): - pre-computed truth table : the set corresponding output for each possible input
	determines the canalizing strength of function F. F must have number of variables equal to 2^n for some n>1, and the number of variables must be 
		greater than 1. Canalizing strength is determined from the k-set canalizing proportions of F for all k between 0 and number of variables in F,
		exclusive. 
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: list bool_list - Pre-computed truth table for the function
	return: tuple - The first entry of the tuple is a float representing the canalizing strength of F. The second entry of the tuple is a list
		of floats, each representing the k-set canalizing proportion for 0<k<number of variables in F.
	Example:
		input: can.f_from_expression("(x and y and z) or (not x and y and z) or (x and not y and z) or (x and y and not z)")[0]
		output: (0.3333333333333333, [0.0, 0.5])

is_k_canalizing(F,k,n):
	determines if F is k-canalizing for k. A function is k-canalizing if it has k variables which all act canalizing assuming other more canalizing 
		variables are not in their canalizing states (for example, suppose F is true if variable x is true, and otherwise true if and only if 
		variable y is true. This function would be 2-canalizing, or k-canalizing for k=2). If F is k-canalizing for a k-value greater than k,
		it is also k-canalizing for k
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int k - the canalizing depth to check for in the function
	param: int n - the number of variables in F
	return: boolean - True if F is k-canalizing for k, False otherwise
	Example: 
		input: (can.f_from_expression("x or y")[0], 2, 2) output: True
	
is_k_canalizing_return_inputs_outputs_corefunctions(F,k,n,can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int)):
	determines if F is k-canalizing for k, and also returns the canalizing inputs, canalizing outputs, and the core function. The core function is 
		the function with all canalizing variables removed
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int k - the canalizing depth to check for in the function
	param: int n - the number of variables in F
	param: can_inputs (optional) - the array to store the canalizing inputs
	param: can_outputs (optional) - the array to store the canalizing outputs
	return: (boolean, array, array, array) - boolean represents whether F is k-canalizing for k. The first array contains the canalizing inputs for each
		canalizing variable of F. The second array contains the outputs of F when canalized by each canalizing variable. The third array represents
		the core function of F.
	Example:
		input: (can.f_from_expression("x or y or (a and b)")[0], 2, 4)
		output: (True, array([1, 1]), array([1, 1]), array([0, 0, 0, 1]))
		
is_k_canalizing_return_inputs_outputs_corefunction_order(F,k,n,can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int),can_order=np.array([],dtype=int),variables=[]): ?
	determines if F is k-canalizing for k, and returns the canalizing inputs, canalizing outputs, and the core function. The core function is 
		the function with all canalizing variables removed. The canalizing variables are checked in order specified by the variables parameter.
		(first to last if order is not specified)
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int k - the canalizing depth to check for in the function
	param: int n - the number of variables in F
	param: can_inputs (optional) - the array to store the canalizing inputs
	param: can_outputs (optional) - the array to store the canalizing outputs
	param: can_order (optional) - the array to store the order in which canalizing variables were checked
	param: variables (optional) - the order in which to check the canalizing variables
	return: (boolean, array, array, array, array) - the boolean indicates whether F is k-canalizing. The first array contains the canalizing inputs 
		for each canalizing variable of F, in the order specified by the variables parameter. The second array 
		contains the outputs of F when canalized by each canalizing variable. The third array represents the order in which canalizing 
		variables were checked.
		
get_all_canalizing_variables_of_boolean_function(rule,or_symbol = 'or', and_symbol = 'and', not_symbol='not',can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int),can_variables=np.array([],dtype=int)): ??
	given a string representation of a function in CNF or DNF, ??
	
average_sensitivity_old_wrong(F,nsim=10000): (deprecated?)
	estimate the average sensitivity of a function in the old and wrong manner.
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int nsim (optional) - the number of simulations to run to determine the average sensitivity of F. Higher values lead to
		a more accurate output, but also make the function execute slower. (default: 10000)
	return: float - the average sensitivity of F
	
average_sensitivity(F,nsim=10000,EXACT=False,NORMALIZED=True):
	determine the average sensitivity of F (estimate if EXACT=False).
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int nsim (optional) - the number of simulations to run to determine the average sensitivity of F. Higher values lead to
		a more accurate output, but also make the function execute slower. (default: 10000)
	param: boolean EXACT - whether the exact average sensitivity should be determined. This will guarantee that the estimate is exact,
		but may increase runtime (default: False)
	param: boolean NORMALIZED - Whether the average sensitivity should be normalized to the number of variables in the function.
		The average sensitivity is normalized by dividing by the number of variables in the function. (default: True)
	return: float - the average sensitivity of F
	
absolute_bias(F,n=None):
	return the absolute bias of function F. Higher absolute bias correlates with increased canalization.
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	param: int n - the number of variables in F. If unspecified, it is set to the number of variables in F regardless.
	return: float - the absolute bias of F
	
	
2) Put everything together to obtain canalizing depth, layer structure, canalized outputs, canalizing inputs as well as core function 
	(could also calculate canalizing variables in future versions)
get_canalizing_depth_inputs_outputs_corefunction(F):
	return a tuple containing the number of variables in F, the maximum k value for which F is k-canalizing, the canalizing inputs, 
	canalizing outputs, and the core function. The core function is the function with all canalizing variables removed. 
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	return: (int, int, array, array, array)
		int: the number of variables in F
		int: the k value for which F is k-canalizing
		array: the canalizing inputs of each canalizing variable of F
		array: the canalizing outputs of each canalizing variable of F
		array: The core function of F
	Example:
		input: can.f_from_expression("x and (y or z)")[0]
		output: (3, 3, array([0, 1, 1]), array([0, 1, 1]), array([0]))

get_canalizing_depth_inputs_outputs_corefunction_order(F, variables = []):
	return a tuple containing the number of variables in F, the maximum k value for which F is k-canalizing, the canalizing inputs, 
	canalizing outputs, and the core function. The core function is the function with all canalizing variables removed. 
	The canalization of the variables is evaluated in order specified by the variables parameter. (first to last if order is not specified)
	param: list F - a representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
	return: (int, int, array, array, array, array):
		int: the number of variables in F
		int: the maximum k-value for which F is k-canalizing
		array: the canalizing inputs of each canalizing variable of F
		array: the canalizing outputs of each canalizing variable of F
		array: The core function of F
		array: The order in which the canalization of variables was evaluated
	Example:
		input: (can.f_from_expression("x and (y or z)")[0], [1,0,2])
		output: (3, 3, array([0, 1, 1]), array([0, 1, 1]), array([0]), array([1, 0, 2]))

get_layer_structure_given_outputs_corefunction(can_outputs_F,corefunction_F,n_F):
	returns the layer structure of the nested canalizing portion of function F, represented by can_outputs_F. The layer structure is returned as 
		an array where each entry is the size of the corresponding layer of can_outputs_F. A layer of F is all adjacent numbers of the same value
		(0 or 1) within can_outputs_F - meaning the consequitive variables that all canalize to the same output. The size of a layer is the 
		number of variables in the layer.
	param: list can_outputs_F - The cananalizing outputs of F
	param: list corefunction_F - The corefunction of F
	param: int n_F - The number of variables in F
	return: list - A list where each entry represents the size of the corresponding layer
	Example:
		Python:
			x = "x and (y or z)"
			func = can.get_canalizing_depth_inputs_outputs_corefunction(can.f_from_expression(x)[0])
			print(can.get_layer_structure_given_outputs_corefunction(func[3], func[4], func[0]))
		
		Terminal:
			[1, 2]
			
kindoflayer(k,w):
	For n-canalizing functions, there is a bijection (function that is both one-to-one and onto, so each input maps to exactly one output, and
		every output is mapped to) that maps Hamming weight w (assuming w = 2^n-w where function is n-canalizing) to the layer structure. 
		This method uses this rule to find a list containing the layer sizes of corresponding layers.
	param: int k - the number of essential variables of the function
	param: int w - the odd Hamming weight of the function
	return: (int, list):
		int: the number of layers determined by this rule
		list: a list of integers representing the layer sizes of corresponding layers


3) Methods to randomly generate Boolean functions (uniform distribution) and Boolean networks
random_function(n):
	return a random function with n variables
	param: int n - the number of variables for the random function to have
	return: array - an array representing the function
	
random_non_degenerated_function(n):
	return a random non-degenerated function with n variables. A function is non-degenerated if, for each variable of the function, there exists
		some possible combination of inputs to other variables such that the output of the function depends on this variable.
	param: int n - the number of variables for the function to have
	return: an array representing the non-degenerated function
	
random_non_canalizing_function(n):
	return a random non-canalizing function with n variables. A function is non-canalizing if it has no variable which has 
		a state (true or false) that determines the evaluation of the function regardless of the other variables.
	param: int n - the number of variables for the function to have
	return: an array representing the non-canalizing function
	
random_non_canalizing_non_degenerated_function(n):
	return a random function that is neither canalizing nor degenerated. 
	param: int n - the number of variables for the function to have
	return: an array representing the non-canalizing non-degenerated function

random_k_canalizing(n, k, EXACT_DEPTH_K=False, x=[]):
	returns a random function with n variables that is k-canalizing.
	param: int n - the number of variables for the function to have
	param: int k - the k-value for which the function should be k-canalizing
	param: boolean EXACT_DEPTH_K=False (optional) - if specified as true, the function cannot be k-canalizing for k values greater than k.
	param: list x=[] (optional)? pre-computed something for size n?
	return: array - represents the k-canalizing function
	
random_k_canalizing_with_specific_weight(n, kis, EXACT_DEPTH_K=False, x=[]):
	returns a random k-canalizing function with n variables and layer sizes equivalent to the corresponding entries of kis.?
	param: int n - the number of variables for the function to have
	param: int k - the k-value for which the function should be k-canalizing
	param: list kis - list of integers representing the layer sizes of the function
	param: boolean EXACT_DEPTH_K=False (optional) - if specified as true, the function cannot be k-canalizing for k values greater than k.
	param: list x=[] (optional) - ?
	
random_adj_matrix(N,ns,NO_SELF_REGULATION=True,STRONGLY_CONNECTED=False):
	recursively generates a random adjacency matrix with N nodes and number of nodes regulated by the ith node is ns[i].
	param: int N - the number of nodes in the graph
	param: list ns - the number of nodes for each node to regulate, where ns[i] is the number of nodes regulated by the ith node
	param: boolean NO_SELF_REGULATION=True (optional) - whether nodes are allowed to regulate themselves
	param: boolean STRONGLY_CONNECTED=False (optional) - whether the graph should be strongly connected. A graph is strongly connected if a 
		path of edges can be found from one node to every other node of the graph.
	return: (2d array, 2d array):
		2d array: the adjacency matrix of the graph
		2d array: the adjacency list of the graph

random_edge_list(N,ns,NO_SELF_REGULATION):
	generates a random edge list with N nodes and where number of nodes regulated by ith node is ns[i]
	param: int N - the number of nodes in the graph
	param: list ns - the number of nodes for each node to regulate, where ns[i] is the number of nodes regulated by the ith node
	param: boolean NO_SELF_REGULATION - whether nodes are allowed to regulate themselves (edges going from and to same node)
	return: array - each entry of this array is a tuple representing an edge. The first entry of this tuple is the node being regulated, and
		the second is the node regulating it
	
random_BN(N, n = 2, k = 0, STRONGLY_CONNECTED = True, indegree_distribution = 'constant', list_x=[], kis = None, EXACT_DEPTH=False,NO_SELF_REGULATION=True):
	generates a random boolean network with N nodes
	params:
		int N - the number of nodes in the network
		list or int n=2 (optional) - the number of nodes regulated by each node. If n is an int, the same n value is used for each node.
		list or int k=0 (optional) - the k-canalizing depth for each node. If k is an int, the same k value is used for each node.
		boolean STRONGLY_CONNECTED=True (optional) - whether the network should be strongly connected (a path along edges 
			can be found from any node to any other node
		String indegree_distribution='constant' (optional) - "uniform", "constant", or "poisson". Specifies the distribution of number of nodes
			regulated by each node. If specified as "delta" or "dirac", it will be treated the same as "constant".
		list list_x=[] (optional) - ?
		list kis=None (optional) - the size of each level of the function
		boolean EXACT_DEPTH=False (optional) - whether the functions can be k-canalizing for values greater than their randomly assigned k value
		boolean NO_SELF_REGULATION=True (optional) - Whether nodes are allowed to regulate themselves
	return: (2d array, 2d array, array)
		2d array - Each entry 'F' of this array is the activation function for the corresponding node of the network (ith entry goes to ith node)
		2d array - The ith entry of this array is an array of the nodes regulated by the ith node (adjacency list)
		array - The ith entry of this array is the number of nodes regulated by the ith node
		
		
4) Enumeration methods
nr_non_canalizing_by_weight_exact(n): ??

nr_non_canalizing_by_weight_simulation(n,nsim=10000): ??

stratify_Boolean_fcts_by_canalization_ns(n,nsim=10000): ??


5) Analysis methods
get_constant_nodes(I,degree,N):
	Get an array of all nodes which are regulated by only themselves
	params:
		2d list I - Adjacency list where each entry represents the regulators of that node
		list degree - A list of the number of regulators of each node
		int N - The number of nodes in I
	return:
		array - A numpy array of nodes regulated only by themselves
	Example:
		Input - get_constant_nodes([[0],[0,1]], [1, 2], 2)
		Output - [0]

rnd_edge_deletion(F, I, N, degree, nsim=100, bool_lists=[]):??
	
update(F, I, N, X):
	Given a network and its activity state, return its new activity state after one update step
	params:
		2d list F - A list of lists where each sub-list represents a node's update function
		2d list I - An adjacency list where each entry represents regulators of that node
		int N - The number of nodes in the network
		list X - A list representing the activity of each node - 0 for inactive, 1 for active.
	return:
		list - Representing the activity of each node after an update step. 0 for inactive, 1 for active.

average_sensitivity(F, nsim=10000, EXACT=False, NORMALIZED=True):
	Given a fucntion represented in the form of the right side of a truth table, return its average sensitivity.
	params:
		list F - A representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
		(optional) nsim - The number of simulations to run if EXACT is False. Higher values increase precision and runtime. (default: 10000)
		(optional) EXACT - Whether the exact average sensitivity should be computed rather than running simulations. Setting this to True
			will cause the function to run in exponential time with regard to the number of variables in F. (default: False)
		(optional) NORMALIZED - Whether the value returned should be normilized by dividing by the number of
			varaibles in F. (default: True)
	return:
		float - The average sensitivity of F, normalized to the number of nodes in F if NORMALIZED is True.

absolute_bias(F, n=None):
	Given a fucntion represented in the form of the right side of a truth table, return the absolute bias of these rules.
	params:
		list F - A representation of a function, in the format of the list at index 0 in the tuple returned by f_from_expression(expr)
		(optional) n - The number of variables in F. If None, it is calculated from the length of F. (default: None)
	return:
		float - the absolute bias of F
	
num_of_attractors(F, I, N, nsim=500, EXACT=False, bool_list=[]):
	Return the attractors of a network.
	params:
		2d list F - A list of lists where each sub-list represents a node's update function
		2d list I - An adjacency list for the network
		int N - the number of nodes in the network
		(optional) nsim - The number of simulations to run. Higher values increase percision and decrease speed. (default: 500)
		(optional) EXACT - Whether all of the attractors should be found exactly. If set to True, the function will run in
			exponential time with respect to N. (default: False)
		(optional) bool_list - pre-computed left side of truth table for N
	return (tuple):
		2d list - Each sub-list represents an attractor.
		int - The total number of attractors
		list - List of integers representing the size of each basin
		
num_of_attractors_v2(F, I, N, nb=500):
	Return the attractors of a network. This function has the same output as num_of_attractors, although it is
		faster for networks with long average path lengths.
	params:
		2d list F - A list of lists where each sub-list represents a node's update function
		2d list I - An adjacency list for the network
		(optional) nsim - The number of simulations to run. Higher values increase percision and decrease speed. (default: 500)
	return (tuple):
		2d list - Each sub-list represents an attractor.
		int - The total number of attractors
		list - List of integers representing the size of each basin
		
basin_size_largest(basin_sizes):
	Returns the largest basin size divided by the sum of all basin sizes in the list.
	params:
		list basin_sizes - A list of numbers
	return:
		float - The largest number in the list divided by the sum of all numbers in the list
		
entropy(basin_sizes): <explain what entropy means, not just what it is>??
	For each element in basin_sizes, divide that element by the sum of all elements, then the corresponding output element is
		the negative natural logarithm of that number times that number.
	params:
		list basin_sizes - A list of numbers
	return:
		list - Each element in this list is -ln(S/total)*S/total where S is the corresponding input element and total is the sum of
			all elements in the input list

d(x, y):
	Returns a positive integer if numpy arrays x and y are different, 0 otherwise.
	params:
		numpy array x - first array to compare
		numpy array y - second array to compare
	return:
		int - positive if x and y are different, 0 otherwise.
		
derrida_value(F, I, N, m, nsim=500):
	Caluclate the derrida value of a network. ?? m = size of preturbation
	
adjacency_matrix(I,constants=[],IGNORE_SELFLOOPS=False,IGNORE_CONSTANTS=True):
	Get an adjacency matrix representation of input adjacency list I.
	params:
		2d list I - the adjacency list to convert to an adjacency matrix.
		(optional) list constants - 




































