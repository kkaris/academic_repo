import numpy as np
from gt_functions import *
from glob import glob

# Translation; resid-string: molecule
# 'SOL': 'H2O',
# 'DMSO': 'dmso',
# '0ONB': 'BF4',
# '_MPO': 'CO2',
# 'EAIF': 'emim'

# List of file names. One file per frame of trajectory.
file_list = glob('test_files/*.dat')

# Open a file to write chain statistics to
outf = open('chain_stats.csv', 'w')
outf.write('size,h2o-h2o_bonds,H2O,BF4,emim,dmso,co2\n')

# Loop over frames/files
for file_name in file_list:

    # Load pair list (starts on 2nd row) and make edgelist:
    # [['id1' 'id2', w1] ['id4' 'id9', w2] ...]
    edge_list = read_pair_list(file_name)
    unique_nodes = uniqe_components(edge_list)
    edgelist_undir = dir2undir_el(edge_list)

    # Create undirected adjacency matrix
    adj_m = adj_from_edgelist(edgelist_undir)

    # Get list of chains in frame
    conn_comp = find_components(adj_m)

    # Loop over chains
    for chain in conn_comp:

        # Initialize count for node types
        nh2ob, nh2o, nbf4, ndmso, nco2, nemim = 0, 0, 0, 0, 0, 0

        # Loop over nodes in chain
        for node_index in chain:

            # Count chain node types on resid
            if 'SOL' in unique_nodes[node_index]:
                nh2o += 1

                # Check all connections with current node using edgelist
                # (that is directed) and add to h2o-h2o bond count
                for edge in edge_list:

                    # matching current id and check if acceptor is
                    if edge[0] == unique_nodes[node_index] and \
                       'SOL' in edge[1]:
                        nh2ob += 1

            elif '0ONB' in unique_nodes[node_index]:
                nbf4 += 1

            elif '_MPO' in unique_nodes[node_index]:
                nco2 += 1

            elif 'EAIF' in unique_nodes[node_index]:
                nemim += 1

            elif 'DMSO' in unique_nodes[node_index]:
                ndmso += 1

        # Save count of chain member types per chain size:
        # size,h2o-h2o_bonds,H2O,BF4,emim,dmso,co2
        outf.write('{},{},{},{},{},{},{}\n'.format(
            nh2o + ndmso + nemim + nco2 + nbf4,
            nh2ob, nh2o, nbf4, nemim, ndmso, nco2))

outf.close()
