# -*- coding: iso-8859-1 -*-
import numpy as np
import hotnet as hn
import multiprocessing as mp
from union_find import UnionFind
from collections import namedtuple, defaultdict

def get_component_sizes(arrs):
    return [len(arr) for arr in arrs]

def delta_too_small(component_sizes, max_size ):
    if max_size > max(component_sizes): 
        return True
    else: 
        return False

Edge = namedtuple("Edge", ["node1", "node2", "weight"])

def find_best_delta(permuted_sim, ks, start=0.05, quiet=True):
    """
    Return a dict mapping each size in ks to the median delta that maximizes the number of
    connected components of size at least k in the graph corresponding the the given similarity
    matrix.
    
    Arguments:
    permuted_sim -- 2D ndarray representing a similarity matrix
    ks -- list of minimum sizes for connected components to be counted. This must be at least 2.
    start -- only deltas in the top start proportion of edge weights will be considered when
             searching for deltas that maximize the number of connected components
    """
    
    if not quiet: 
        print("Finding median delta that maximizes the # of CCs of size >= l")
    edges = get_edges(permuted_sim, start)
    k2delta = {}

    for k in ks:
        _, bestDeltas = find_best_delta_for_given_k(permuted_sim, edges, k)
        k2delta[k] = np.median(bestDeltas)

    return k2delta

def find_best_delta_for_given_k(permuted_sim, edges, k):
    
    if k < 2:
        raise ValueError("k must be at least 2")

    max_num_ccs = 0 #initially, each node is its own CC of size 1, so none is of size >= k for k >= 2
    bestDeltas = [edges[0].weight]
    uf = UnionFind()

    for edge in edges:
        uf.union(edge.node1, edge.node2)
        num_ccs = len([root for root in uf.roots if uf.weights[root] >= k])
        if num_ccs > max_num_ccs:
            max_num_ccs = num_ccs
            bestDeltas = [edge.weight]
        elif num_ccs == max_num_ccs:
            bestDeltas.append(edge.weight)

    return max_num_ccs, bestDeltas

def get_edges(sim, start=.05):
    """Return a list of Edge tuples representing edges in the top start% of edge weights"""
    flattened = np.ndarray.flatten(sim)
    edges = [Edge(i/len(sim), i%len(sim), flattened[i]) for i in range(len(flattened)) if i/len(sim) <= i%len(sim)]
    edges = sorted(edges, key=lambda x: x.weight, reverse=True)
    edges = edges[:int(start*len(edges))]
    return edges

def heat_delta_wrapper(infmat, index2gene, heat_permutation, sizes):  
    M, index2gene = hn.induce_infmat(infmat, index2gene, sorted(heat_permutation.keys()))
    heat = hn.heat_vec(heat_permutation, index2gene)
    sim_mat = hn.similarity_matrix(M, heat)
    return find_best_delta(sim_mat, sizes)

def heat_delta_selection(infmat, index2gene, heat_permutations, sizes, parallel=True):
    """Return a dict mapping each size in sizes to a list of the best deltas for each heat
    permutation for that size. 
    
    Arguments:
    infmat -- 2D ndarray representing an influence matrix
    index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
                  index in the influence matrix
    heat_permutations -- list of heat permutations (dicts mapping gene name to heat score)
    sizes -- list of sizes for largest CC / min size for CCs to be counted (based on selection_fn)
    parallel -- whether finding the best delta for each permuted network should be performed in parallel
    
    """
    
    if parallel:
        pool = mp.Pool()
        map_fn = pool.map
    else:
        map_fn = map
    args = [(infmat, index2gene, heat_permutation, sizes)
            for heat_permutation in heat_permutations]
    deltas = map_fn(heat_delta_wrapper, args)

    if parallel:
        pool.close()
        pool.join()
    
    # Parse the deltas into one dictionary
    sizes2deltas = defaultdict(list)
    for size2delta in deltas:
        for s in sizes: sizes2deltas[s].append(size2delta[s])
         
    return sizes2deltas
