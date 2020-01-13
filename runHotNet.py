#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import hnap
import hnio
import hotnet as hn
import permutations
import stats
import sys
import scipy.io
import json
import os

def parse_args(raw_args): 
    description = "Runs HotNet algorithm.\
                   Note that some or all parameters can be specified via a configuration file by\
                   passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python runHotNet.py --runname TestRun @testConf.txt'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')

    parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parser.add_argument('-mf', '--infmat_file', required=True,
                        help='Path to .mat file containing influence matrix')
    parser.add_argument('-mn', '--infmat_name', default='Li',
                        help='Variable name of the influence matrix in the .mat file')
    parser.add_argument('-if', '--infmat_index_file', required=True,
                        help='Path to tab-separated file containing an index in the first column\
                              and the name of the gene represented at that index in the second\
                              column of each line.')
    parser.add_argument('-hf', '--heat_file', required=True,
                        help='JSON heat score file generated via generateHeat.py')
    parser.add_argument('-d', '--delta', type=float, required=True,
                        help='Weight threshold for edge removal')
    parser.add_argument('-ccs', '--min_cc_size', type=int, default=3,
                        help='Minimum size connected components that should be returned.')
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory. Files results.json, components.txt, and\
                              significance.txt will be generated.')
    
    #parent parser for arguments common to all permutation types
    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-n', '--num_permutations', type=int, required=True,
                               help='Number of permutation tests to run; set to 0 to skip running\
                                     permutation tests.')
    parent_parser.add_argument('-s', '--cc_start_size', type=int, default=2,
                               help='Smallest connected component size to count in permutation tests')
    parent_parser.add_argument('-l', '--cc_stop_size', type=int, default=10,
                               help='Largest connected component size to count in permutation tests')
    parent_parser.add_argument('--parallel', dest='parallel', action='store_true',
                               help='Run permutation tests in parallel. Only recommended for machines\
                                     with at least 8 cores.')
    parent_parser.add_argument('--no-parallel', dest='parallel', action='store_false',
                               help='Run permutation tests sequentially. Recommended for machines\
                                     with fewer than 8 cores.')
    parent_parser.set_defaults(parallel=False)
    
    subparsers = parser.add_subparsers(title='Heat score type', dest='permutation_type')
    
    subparsers.add_parser('none', help='Do not perform statistical significance permutation tests')
    
    heat_parser = subparsers.add_parser('heat', help='Permute heat scores', parents=[parent_parser])
    heat_parser.add_argument('-pgf', '--permutation_genes_file',
                             help='Path to file containing a list of additional genes that can have\
                                   permuted heat values assigned to them in permutation tests')
    
    mutation_parser = subparsers.add_parser('mutations', help='Permute mutation data',
                                             parents=[parent_parser])
    mutation_parser.add_argument('-glf', '--gene_length_file', required=True,
                                 help='Path to tab-separated file containing gene names in the\
                                       first column and the length of the gene in base pairs in\
                                       the second column')
    mutation_parser.add_argument('-gof', '--gene_order_file', required=True,
                                 help='Path to file containing tab-separated lists of genes on\
                                 each chromosme, in order of their position on the chromosome, one\
                                  chromosome per line')
    mutation_parser.add_argument('-b', '--bmr', type=float, required=True,
                                 help='Default background mutation rate')
    mutation_parser.add_argument('-bf', '--bmr_file',
                                 help='File listing gene-specific BMRs. If none, the default BMR\
                                       will be used for all genes.')
    
    return parser.parse_args(raw_args)

JSON_OUTPUT = "results.json"
COMPONENTS_TSV = "components.txt"
SIGNIFICANCE_TSV = "significance.txt"

def run(args):
    # create output directory if doesn't exist; warn if output files already exist
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    dir_contents = os.listdir(args.output_directory)
    if JSON_OUTPUT in dir_contents or COMPONENTS_TSV in dir_contents or SIGNIFICANCE_TSV in dir_contents:
        print("WARNING: Output directory already contains HotNet results file(s), which will be "
              "overwritten. (Ctrl-c to cancel).")
    
    # load data
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]  
    infmat_index = hnio.load_index(args.infmat_index_file)
    heat, heat_params = hnio.load_heat_json(args.heat_file)
  
    # compute similarity matrix and extract connected components
    M, gene_index = hn.induce_infmat(infmat, infmat_index, sorted(heat.keys()), quiet=False)
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h)
    G = hn.weighted_graph(sim, gene_index, args.delta)
    ccs = hn.connected_components(G, args.min_cc_size)
    
    # calculate significance
    if args.permutation_type != "none":
        if args.permutation_type == "heat":
            sizes2stats = heat_permutation_significance(args, heat, infmat, infmat_index, G)
        elif args.permutation_type == "mutations":
            if heat_params["heat_fn"] != "load_mutation_heat":
                raise RuntimeError("Heat scores must be based on mutation data to perform\
                                    significance testing based on mutation data permutation.")
            sizes2stats = mutation_permutation_significance(args, infmat, infmat_index, G, heat_params)
        else:
            raise ValueError("Unrecognized permutation type %s" % (args.permutation_type))
    
    #sort ccs list such that genes within components are sorted alphanumerically, and components
    #are sorted first by length, then alphanumerically by name of the first gene in the component 
    ccs = [sorted(cc) for cc in ccs]
    ccs.sort(key=lambda comp: comp[0])
    ccs.sort(key=len, reverse=True)
    
    # write output
    output_dict = {"parameters": vars(args), "heat_parameters": heat_params,
                   "sizes": hn.component_sizes(ccs), "components": ccs}
    if args.permutation_type != "none":
        output_dict["statistics"] = sizes2stats
        hnio.write_significance_as_tsv(os.path.abspath(args.output_directory) + "/" + SIGNIFICANCE_TSV,
                                       sizes2stats)
    
    json_out = open(os.path.abspath(args.output_directory) + "/" + JSON_OUTPUT, 'w')
    json.dump(output_dict, json_out, indent=4)
    json_out.close()
    
    hnio.write_components_as_tsv(os.path.abspath(args.output_directory) + "/" + COMPONENTS_TSV, ccs)

def heat_permutation_significance(args, heat, infmat, infmat_index, G):
    print("* Performing permuted heat statistical significance...")
    
    addtl_genes = hnio.load_genes(args.permutation_genes_file) if args.permutation_genes_file else None
    heat_permutations = permutations.permute_heat(heat, args.num_permutations, addtl_genes, args.parallel)
    return calculate_significance(args, infmat, infmat_index, G, heat_permutations)

def mutation_permutation_significance(args, infmat, infmat_index, G, heat_params):
    print("* Performing permuted mutation data statistical significance...")
    
    heat_permutations = permutations.generate_mutation_permutation_heat(
                            heat_params["heat_fn"], heat_params["sample_file"],
                            heat_params["gene_file"], list(infmat_index.values()), heat_params["snv_file"],
                            args.gene_length_file, args.bmr, args.bmr_file, heat_params["cna_file"],
                            args.gene_order_file, heat_params["cna_filter_threshold"],
                            heat_params["min_freq"], args.num_permutations, args.parallel)

    return calculate_significance(args, infmat, infmat_index, G, heat_permutations)

def calculate_significance(args, infmat, infmat_index, G, heat_permutations):
    sizes = list(range(args.cc_start_size, args.cc_stop_size+1))
    
    print("\t- Using no. of components >= k (k \\in")
    print("[%s, %s]) as statistic" % (min(sizes), max(sizes)))

    #size2counts is dict(size -> (list of counts, 1 per permutation))
    sizes2counts = stats.calculate_permuted_cc_counts(infmat, infmat_index, heat_permutations,
                                                      args.delta, sizes, args.parallel)
    real_counts = stats.num_components_min_size(G, sizes)
    size2real_counts = dict(list(zip(sizes, real_counts)))
    return stats.compute_statistics(size2real_counts, sizes2counts, args.num_permutations)

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
