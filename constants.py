#!/usr/bin/env python

from collections import namedtuple

SNV = "snv"
AMP = "amp"
DEL = "del"

VIZ_INDEX = 'index.html'
VIZ_SUBNETWORKS = 'subnetworks.html'

Mutation = namedtuple("Mutation", ["sample", "gene", "mut_type"])