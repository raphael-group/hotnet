#!/usr/bin/env python

import sys
import hnap
import hnio
import scipy.io
import numpy as np
import findThreshold as ft
import heat as hnheat
import hotnet as hn
import permutations as p
import stats
import viz
from constants import *
import os
import json
import shutil

MIN_CC_SIZE = 3
MAX_CC_SIZE = 25
INFMAT_NAME = "Li"

def parse_args(raw_args): 
    description = "Helper script for simple runs of generalized HotNet,.\
                   including automated parameter selection."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mf', '--infmat_file', required=True,
                        help='Path to .mat file containing influence matrix')
    parser.add_argument('-if', '--infmat_index_file', required=True,
                        help='Path to tab-separated file containing an index in the first column\
                              and the name of the gene represented at that index in the second\
                              column of each line.')
    parser.add_argument('-hf', '--heat_file', required=True,
                        help='Path to a tab-separated file containing a gene name in the first\
                              column and the heat score for that gene in the second column of\
                              each line.')
    parser.add_argument('-ms', '--min_heat_score', type=float,
                        help='Minimum heat score for a gene to be eligible for inclusion in a\
                              returned connected component. By default, all genes with positive\
                              heat scores will be included. (To include genes with score zero, set\
                              min_heat_score to 0).')
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=3,
                        help='Minimum size connected components that should be returned.')
    parser.add_argument('-n', '--num_permutations', type=int, default=100,
                        help='Number of permutations that should be used for parameter selection\
                              and statistical significance testing.')
    parser.add_argument('-o', '--output_directory', default='hotnet_output',
                        help='Output directory. Files results.json, components.txt, and\
                              significance.txt will be generated in subdirectories for each delta.')
    parser.add_argument('--parallel', dest='parallel', action='store_true',
                        help='Run permutation tests in parallel. Only recommended for machines\
                              with at least 8 cores.')
    parser.add_argument('--no-parallel', dest='parallel', action='store_false',
                        help='Run permutation tests sequentially. Recommended for machines\
                              with fewer than 8 cores.')
    parser.add_argument('-ef', '--edge_file',
                        help='Path to TSV file listing edges of the interaction network, where\
                              each row contains the indices of two genes that are connected in the\
                              network. This is used to create subnetwork visualizations; if not\
                              provided, visualizations will not be made.')
    parser.add_argument('-nn', '--network_name', default='Network',
                        help='Display name for the interaction network. (Used for subnetwork\
                              visualizations)')
    parser.set_defaults(parallel=False)
    
    return parser.parse_args(raw_args)

JSON_OUTPUT = "results.json"
COMPONENTS_TSV = "components.txt"
SIGNIFICANCE_TSV = "significance.txt"
HEAT_JSON = "heat.json"

def run(args):
    # create output directory if doesn't exist; warn if it exists and is not empty
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    if len(os.listdir(args.output_directory)) > 0:
        print("WARNING: Output directory is not empty. Any conflicting files will be overwritten. ")
        print("(Ctrl-c to cancel).")
    
    infmat = scipy.io.loadmat(args.infmat_file)[INFMAT_NAME]
    infmat_index = hnio.load_index(args.infmat_index_file)
    heat = hnio.load_heat_tsv(args.heat_file)
    
    # filter out genes with heat score less than min_heat_score
    heat, addtl_genes, args.min_heat_score = hnheat.filter_heat(heat, args.min_heat_score)
    
    # find delta that maximizes # CCs of size >= MIN_SIZE for each permuted data set
    deltas = ft.get_deltas_for_heat(infmat, infmat_index, heat, addtl_genes, args.num_permutations,
                                    args.parallel)

    #find the multiple of the median delta s.t. the size of the largest CC in the real data
    #is <= MAX_CC_SIZE
    medianDelta = np.median(deltas[MIN_CC_SIZE])
    M, gene_index = hn.induce_infmat(infmat, infmat_index, sorted(heat.keys()), quiet=False)
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h)
    
    for i in range(1, 11):
        G = hn.weighted_graph(sim, gene_index, i*medianDelta)
        max_cc_size = max([len(cc) for cc in hn.connected_components(G)])
        if max_cc_size <= MAX_CC_SIZE:
            break
    
    # load interaction network edges and determine location of static HTML files for visualization
    edges = hnio.load_ppi_edges(args.edge_file) if args.edge_file else None
    index_file = '%s/viz_files/%s' % (os.path.realpath(__file__).rsplit('/', 1)[0], VIZ_INDEX)
    subnetworks_file = '%s/viz_files/%s' % (os.path.realpath(__file__).rsplit('/', 1)[0], VIZ_SUBNETWORKS)
    gene2index = dict([(gene, index) for index, gene in list(infmat_index.items())])

    #and run HotNet with that multiple and the next 4 multiples
    run_deltas = [i*medianDelta for i in range(i, i+5)]
    for delta in run_deltas: 
        # create output directory
        delta_out_dir = args.output_directory + "/delta_" + str(delta)
        if not os.path.isdir(delta_out_dir):
            os.mkdir(delta_out_dir)
        
        # find connected components
        G = hn.weighted_graph(sim, gene_index, delta)
        ccs = hn.connected_components(G, args.min_cc_size)
        
        # calculate significance (using all genes with heat scores)
        print("* Performing permuted heat statistical significance...")
        heat_permutations = p.permute_heat(heat, args.num_permutations, addtl_genes, args.parallel)
        sizes = list(range(2, 11))
        print("\t- Using no. of components >= k (k \\in")
        print("[%s, %s]) as statistic" % (min(sizes), max(sizes)))
        sizes2counts = stats.calculate_permuted_cc_counts(infmat, infmat_index, heat_permutations,
                                                          delta, sizes, args.parallel)
        real_counts = stats.num_components_min_size(G, sizes)
        size2real_counts = dict(list(zip(sizes, real_counts)))
        sizes2stats = stats.compute_statistics(size2real_counts, sizes2counts, args.num_permutations)
    
        # sort ccs list such that genes within components are sorted alphanumerically, and components
        # are sorted first by length, then alphanumerically by name of the first gene in the component 
        ccs = [sorted(cc) for cc in ccs]
        ccs.sort(key=lambda comp: comp[0])
        ccs.sort(key=len, reverse=True)
    
        # write output
        heat_dict = {"heat": heat, "parameters": {"heat_file": args.heat_file}}
        heat_out = open(os.path.abspath(delta_out_dir) + "/" + HEAT_JSON, 'w')
        json.dump(heat_dict, heat_out, indent=4)
        heat_out.close()
        
        args.heat_file = os.path.abspath(delta_out_dir) + "/" + HEAT_JSON
        args.delta = delta
        output_dict = {"parameters": vars(args), "sizes": hn.component_sizes(ccs),
                       "components": ccs, "statistics": sizes2stats}
        hnio.write_significance_as_tsv(os.path.abspath(delta_out_dir) + "/" + SIGNIFICANCE_TSV,
                                       sizes2stats)
        
        json_out = open(os.path.abspath(delta_out_dir) + "/" + JSON_OUTPUT, 'w')
        json.dump(output_dict, json_out, indent=4)
        json_out.close()
        
        hnio.write_components_as_tsv(os.path.abspath(delta_out_dir) + "/" + COMPONENTS_TSV, ccs)

        # write visualization output if edge file given
        if args.edge_file:
            viz_data = {"delta": delta, 'subnetworks': list()}
            for cc in ccs:
                viz_data['subnetworks'].append(viz.get_component_json(cc, heat, edges, gene2index, args.network_name))
                
            delta_viz_dir = '%s/viz/delta%s' % (args.output_directory, delta)
            if not os.path.isdir(delta_viz_dir):
                os.makedirs(delta_viz_dir)
            viz_out = open('%s/subnetworks.json' % delta_viz_dir, 'w')
            json.dump(viz_data, viz_out, indent=4)
            viz_out.close()
   
            shutil.copy(subnetworks_file, delta_viz_dir)
    
    if args.edge_file:
        viz.write_index_file(index_file, '%s/viz/%s' % (args.output_directory, VIZ_INDEX), run_deltas)
    
if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
