#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import hnap
import hnio
import hotnet as hn
import delta
import permutations
import sys
import json
import scipy.io
import numpy as np

MIN_CC_SIZE = 3
MAX_CC_SIZE = 25

def parse_args(raw_args):
    description = "Runs HotNet threshold-finding procedure.\
                   Note that some or all parameters can be specified via a configuration file by\
                   passing '@<ConfigFileName>' as a command-line parameter, e.g.\
                   'python findThreshold.py @testConf.txt --runname TestRun'."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    
    #create parent parser for arguments common to both permutation types
    parent_parser = hnap.HotNetArgParser(add_help=False, fromfile_prefix_chars='@')
    parent_parser.add_argument('-r', '--runname', help='Name of run / disease.')
    parent_parser.add_argument('-mf', '--infmat_file', required=True,
                             help='Path to .mat file containing influence matrix')
    parent_parser.add_argument('-mn', '--infmat_name', default='Li',
                                help='Variable name of the influence matrix in the .mat file')
    parent_parser.add_argument('-if', '--infmat_index_file', required=True, default=None,
                                help='Path to tab-separated file containing an index in the first\
                                      column and the name of the gene represented at that index in\
                                      the second column of each line.')
    parent_parser.add_argument('-hf', '--heat_file', required=True,
                               help='JSON heat score file generated via generateHeat.py')
    parent_parser.add_argument('-n', '--num_permutations', type=int, required=True,
                                help='Number of permuted data sets to generate')
    parent_parser.add_argument('--parallel', dest='parallel', action='store_true',
                               help='Run permutation tests in parallel. Only recommended for machines\
                                     with at least 8 cores.')
    parent_parser.add_argument('--no-parallel', dest='parallel', action='store_false',
                               help='Run permutation tests sequentially. Recommended for machines\
                                     with fewer than 8 cores.')
    parent_parser.add_argument('-o', '--output_file',
                        help='Output file.  If none given, output will be written to stdout.')
    parent_parser.set_defaults(parallel=False)
    
    subparsers = parser.add_subparsers(title='Permutation techniques', dest='perm_type')
    
    #create subparser for options for permuting heat scores
    heat_parser = subparsers.add_parser('heat', help='Permute heat scores', parents=[parent_parser])
    heat_parser.add_argument('-pgf', '--permutation_genes_file', default=None,
                             help='Path to file containing a list of additional genes that can have\
                                   permuted heat values assigned to them in permutation tests')

    #create subparser for options for permuting mutation data
    mutation_parser = subparsers.add_parser('mutations', help='Permute mutation data',
                                            parents=[parent_parser])
    mutation_parser.add_argument('-glf', '--gene_length_file', required=True,
                                 help='Path to tab-separated file containing gene names in the\
                                       first column and the length of the gene in base pairs in\
                                       the second column')
    mutation_parser.add_argument('-gof', '--gene_order_file', required=True,
                                 help='Path to file containing tab-separated lists of genes on\
                                 each chromosome, in order of their position on the chromosome, one\
                                  chromosome per line')
    mutation_parser.add_argument('-b', '--bmr', type=float, required=True,
                                 help='Default background mutation rate')
    mutation_parser.add_argument('-bf', '--bmr_file',
                                 help='File listing gene-specific BMRs. If none, the default BMR\
                                       will be used for all genes.')
                        
    return parser.parse_args(raw_args)

def run(args):
    infmat = scipy.io.loadmat(args.infmat_file)[args.infmat_name]
    infmat_index = hnio.load_index(args.infmat_index_file)
    heat, heat_params = hnio.load_heat_json(args.heat_file)
        
    if args.perm_type == "heat":
        addtl_genes = hnio.load_genes(args.permutation_genes_file) if args.permutation_genes_file else None
        deltas = get_deltas_for_heat(infmat, infmat_index, heat, addtl_genes,
                                         args.num_permutations, args.parallel)
    elif args.perm_type == "mutations":
        deltas = get_deltas_for_mutations(args, infmat, infmat_index, heat_params)
    else:
        raise ValueError("Invalid mutation permutation type: %s" % args.perm_type)
    
    #find the multiple of the median delta s.t. the size of the largest CC in the real data
    #is <= MAX_CC_SIZE
    medianDelta = np.median(deltas[MIN_CC_SIZE])
    M, gene_index = hn.induce_infmat(infmat, infmat_index, sorted(heat.keys()))
    h = hn.heat_vec(heat, gene_index)
    sim = hn.similarity_matrix(M, h)
    
    for i in range(1, 11):
        G = hn.weighted_graph(sim, gene_index, i*medianDelta)
        max_cc_size = max([len(cc) for cc in hn.connected_components(G)])
        if max_cc_size <= MAX_CC_SIZE:
            break
    
    #and recommend running HotNet with that multiple and the next 4 multiples
    recommended_deltas = [i*medianDelta for i in range(i, i+5)]

    output_file = open(args.output_file, 'w') if args.output_file else sys.stdout
    json.dump({"parameters": vars(args), "heat_parameters": heat_params,
               "recommended_deltas": recommended_deltas}, output_file, indent=4)
    if (args.output_file): 
        output_file.close()

def get_deltas_for_heat(infmat, index2gene, gene2heat, addtl_genes, num_permutations, parallel):
    print("* Performing permuted heat delta selection...")
    heat_permutations = permutations.permute_heat(gene2heat, num_permutations, addtl_genes, parallel)
    return get_deltas_from_heat_permutations(infmat, index2gene, heat_permutations, parallel)

def get_deltas_for_mutations(args, infmat, index2gene, heat_params):    
    print("* Performing permuted mutation data delta selection...")
    heat_permutations = permutations.generate_mutation_permutation_heat(
                            heat_params["heat_fn"], heat_params["sample_file"],
                            heat_params["gene_file"], list(index2gene.values()), heat_params["snv_file"],
                            args.gene_length_file, args.bmr, args.bmr_file, heat_params["cna_file"],
                            args.gene_order_file, heat_params["cna_filter_threshold"],
                            heat_params["min_freq"], args.num_permutations, args.parallel)
    return get_deltas_from_heat_permutations(infmat, index2gene, heat_permutations, args.parallel)

def get_deltas_from_heat_permutations(infmat, gene_index, heat_permutations, parallel):   
    return delta.heat_delta_selection(infmat, gene_index, heat_permutations, [MIN_CC_SIZE], parallel)

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))
