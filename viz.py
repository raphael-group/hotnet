#!/usr/bin/env python

import hnio

def get_nodes(cc, gene2heat):
    return [{'name': gene, 'heat': gene2heat[gene]} for gene in cc]

def get_edges(cc, edges, gene2index, networkName):
    edgeData = list()
    for i in range(len(cc)):
        for j in range(i+1, len(cc)):
            gene1 = min(cc[i], cc[j])
            gene2 = max(cc[i], cc[j])
            if (gene2index[gene1], gene2index[gene2]) in edges or \
               (gene2index[gene2], gene2index[gene1]) in edges:
                edgeData.append({'source': gene1, 'target': gene2, 'networks': [networkName]})

    return edgeData

def get_component_json(cc, gene2heat, edges, gene2index, networkName):
    nodes = get_nodes(cc, gene2heat)
    cc_edges = get_edges(cc, edges, gene2index, networkName)

    return {'nodes': nodes, 'edges': cc_edges}

def write_index_file(index_file, out_file, deltas):
    index = hnio.load_file(index_file)
    index += '<ul>\n'
    for delta in deltas:
        index += '<li><a href="delta%s/subnetworks.html">&delta; = %s</a></li>\n' % (delta, delta)
    index += '</ul>'
    hnio.write_file(out_file, index)