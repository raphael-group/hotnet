#!/usr/bin/env python

import hnap
import sys
import json
import hnio
import os
from collections import defaultdict

def parse_args(raw_args):
    description = "For HotNet output from mutation data, generates an annotated output file that\
                   includes the number of samples mutated in each gene and the mutation frequency\
                   of each subnetwork."
    parser = hnap.HotNetArgParser(description=description, fromfile_prefix_chars='@')
    
    parser.add_argument('--hotnet_output_json', required=True,
                        help='Path to JSON output file produced by running "runHotNet.py"')
    
    return parser.parse_args(raw_args)

ANNOTATION_TSV = "annotated_components.txt"

def run(args):
    output_f = open(args.hotnet_output_json)
    output_blob = json.load(output_f)
    output_f.close()
    
    heat_parameters = output_blob["heat_parameters"]
    
    if heat_parameters["heat_fn"] != "load_mutation_heat":
        raise ValueError("Heat scores must have been calculated from mutation data to annotate output.")
    
    components = output_blob["components"]
    genes = hnio.load_genes(heat_parameters["gene_file"])
    samples = hnio.load_samples(heat_parameters["sample_file"])
    snvs = hnio.load_snvs(heat_parameters["snv_file"], genes, samples)
    cnas = hnio.load_cnas(heat_parameters["cna_file"], genes, samples)
    
    if not samples:
        samples = set([snv.sample for snv in snvs] + [cna.sample for cna in cnas])
    
    gene2mutsam = defaultdict(set)
    for mut in snvs + cnas:
        gene2mutsam[mut.gene].add(mut.sample)
    
    annotated_ccs = list()
    for component in components:
        annotated_cc = list()
        annotated_ccs.append(annotated_cc)
        cc_mutated_samples = set()
        for gene in component:
            cc_mutated_samples.update(gene2mutsam[gene])
            annotated_cc.append("%s(%s)" % (gene, len(gene2mutsam[gene])))
        annotated_cc.insert(0, "%s(%s%%)" % (len(cc_mutated_samples),
                                             len(cc_mutated_samples) / float(len(samples)) * 100))
    
    output_directory = output_blob["parameters"]["output_directory"]
    output_file = os.path.abspath(output_directory) + "/" + ANNOTATION_TSV
    hnio.write_components_as_tsv(output_file, annotated_ccs)

if __name__ == "__main__": 
    run(parse_args(sys.argv[1:]))