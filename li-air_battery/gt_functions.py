import numpy as np
from copy import deepcopy


# Functions with no dependencies
def read_pair_list(file_path):
    """

    :param file_path: String containing file path to pair list of a frame
    :return: List of pairs (in lists) formatted from
    'SOL42560-Side \t 0ONB554-Side \t 100.00%\n' to ['id1', 'id2', 1.0]
    and then saved to list [['id1', 'id2', w] ['id4', 'id9', w] ..]
    """
    edge_list = list()

    with open(file=file_path) as f:
        all_lines = f.readlines()

        # Skip first two lines
        for line in all_lines[2:]:

            # Add weight of edge
            if float(line.split()[2].strip('%')) >= 100.0:
                w = 1.0
            else:
                w = float(line.split()[2].strip('%')) / 100

            donor, acceptor = line.split()[:2]

            if donor.endswith('-Side'):
                donor = donor[:-5]

            if acceptor.endswith('-Side'):
                acceptor = acceptor[:-5]

            edge_list.append([donor, acceptor, w])

    return edge_list


def uniqe_components(edge_list):
    """

    :param edge_list: list of pairs of chain molecules
    :return: list of unique molecule ids (resids)
    """
    # Edge list comes in with entries like:
    # [['id1', 'id2', w] ['id4', 'id9' w] ...]

    # Utilize the uniqueness of the set function
    s = set()

    for p in edge_list:
        s.add(p[0])
        s.add(p[1])

    ucl = list(s)
    ucl.sort()

    return ucl


def dir2undir_el(edge_list):
    """

    :param edge_list: nested list - directed edgelist of edges with weights
    :return edge_list: nested list - undirected edgelist of edges with weights
    """

    # Takes directed edgelist and makes it undirected
    temp_list = deepcopy(edge_list)
    for n in range(len(temp_list)):
        temp_list[n][0], temp_list[n][1] = temp_list[n][1], temp_list[n][0]

    # Merge edgelist with its flipped copy to create an undirected edgelist
    for pair in temp_list:
        edge_list.append(pair)

    return edge_list


def is_directed(adj):

    """

    is_directed(adj) checks if the adjacency matrix adj is directed.
    Returns True if it is directed, False otherwise

    Usage:
    Bool = isDirected(adj)

    INPUT: adjacency matrix adj, a NumPy array
    OUTPUT: True/False
    """

    # adj is directed if adj != adj'
    # False if they are equal,
    # True if they are not equal

    return not (np.array_equal(adj, np.transpose(adj)))


def k_neighbors(adj, index, k):

    """

    :param adj: Adjacency matrix
    :param index: index of node to check k-neighbors for
    :param k: integer, number of edges away to check number of neighbors
    :return nl: List of k-neighbor indices

    This function is based of the theorem that states that
    if adj**k is the k-th power of the adjacency matrix, then
    pk_ij = [adj**k]_ij is the number of paths, of length k, linking node i
    to node j.
    Source: BETA Mathematics Handbook for Science and Engineering by
    Lennart Råde and Bertil Westergren
    """

    ak = adj

    for n in range(k-1):
        ak = ak.dot(adj)

    return ak[index, :].nonzero()


# Functions with one dependency
def adj_from_edgelist(edge_list):
    members = uniqe_components(edge_list)
    adj_m = np.zeros(shape=(len(members), len(members)))

    # Row = index of start node
    # Column = index of end node

    for pairs in edge_list:
        r, c, w = members.index(pairs[0]), members.index(pairs[1]), pairs[2]
        adj_m[r, c] = w
    return adj_m


def degree(adj):

    """

    :param adj: Adjacency matrix
    :return: total degree, indegree, outdegree
    """

    if is_directed(adj):
        indeg = np.sum(adj, 0)
        outdeg = np.sum(adj, 1)
        deg = indeg + outdeg

    else:
        deg = np.sum(adj, 0) + np.diag(adj)
        indeg, outdeg = deg, deg

    return deg, indeg, outdeg


# Functions with multiple dependencies
def find_comp_i(adj, ind):

    """

    :param adj: Adjacency matrix (np 2d-array)
    :param ind: index of node i
    :return: (np 1d-array) indices of nodes connected
    to node i, including index of node i
    """

    # Loop through 'neigh', use a "do-while" equivalent since
    # a first iteration is needed.

    # Combine the initial node index with the indices of its 1-neighbors
    neigh = np.unique(np.append(ind, k_neighbors(adj, ind, 1)))
    neigh_temp1 = neigh
    # Repeat until no new neighbors are found
    while True:
        
        # Go through all indices in 'neigh', find neighbors for each and
        # merge list (could be made more efficient by only checking the
        # latest addition of neighbors)
        for n_ind in neigh:
            # Combine current index 'n_ind' with its 1-neighbors in temp2
            neigh_temp2 = np.unique(np.append(n_ind, k_neighbors(adj, n_ind, 1)))
            # Expand temp1 with temp2 
            neigh_temp1 = np.unique(np.append(neigh_temp1, neigh_temp2))

        # Check if neigh == neigh_temp1
        if np.array_equal(neigh, neigh_temp1):
            break

        # If not equal, update neigh with neigh_temp1
        else:
            neigh = neigh_temp1

    return neigh


def find_components(adj):

    """

    :param adj: Adjacency matrix
    :return conn_comps: A list of lists, each one containing the members of each of the grouped nodes

    This function goes over each node (node index) and finds the nodes (node indices)
    connected to it, and returns a set containing tuples
    """

    # Initialize the list
    conn_comps = list()

    # loop the total number of nodes
    for n in range(adj.shape[0]):
        # find connected components of node i
        temp = list(find_comp_i(adj, n))

        # Compare to previous conn comps ('if not in list...')
        if temp not in conn_comps:
            # save into new
            conn_comps.append(temp)

    conn_comps.sort()

    return conn_comps
