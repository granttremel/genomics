

#%%

import numpy as np
from matplotlib import pyplot as plt

from collections import Counter
import os
from pprint import pprint

from cyvcf2 import VCF
import pyranges

import ggene
from ggene import msats
from ggene.utils import inspect, write_object, locals
from ggene.draw import draw_gene_structure

from ggene.genomemanager import GenomeManager
geneman = GenomeManager()

import logging
logger = logging.getLogger(__name__)
logging.getLogger().setLevel(logging.INFO)


# locs = locals("test1")

def write_gene_data(name, chrom):
    testgene = geneman.assemble_gene(name, chrom, min_qual = 5)
    locs = locals(f"./analysis/{name}")
    locs.add(name,testgene)
#***************************************************************************
# varmap = geneman.make_variant_map(nchr)
# locals.add('varmap',varmap)

# qstats,quals = geneman.get_quality_stats(nchr)
# locals.add('quals',pe = 'data')

# ldels, lins = geneman.find_long_variants(nchr, 300)

# l,c,n = geneman.genemap.neighbors(nchr,lins[0]['pos'])

# locs = locals(f'chr{nchr}-longvars')
# locs.add('long_dels', ldels)
# locs.add('long_inserts', lins)
# if l:
#     locs.add('l',l)
# if n:
#     locs.add('n',n)
# if c:
#     locs.add('c',c)
# c,res = geneman.gene_map.find_gene('HTR2A')

interesting_genes = [
    ('HTR2A',13),
    ('HTR2B',2),
    ('HTR2C','X'),
    ('SLC6A4',17),
    ('FTO',16),
    ('IRX3',16),
    ('KISS1R',19),
    ('AKNAD1',1)
    ]
# interesting_gene = 'HTR2C'
# for g,c in interesting_genes:
#     gene = geneman.assemble_gene(g,c)

d,ba = ggene.make_splice_re()
logger.info(d)
logger.info(ba)

gene = geneman.assemble_gene(*interesting_genes[3])

first_tx = gene.transcripts[next(iter(gene.transcripts.keys()))]

seq = geneman.get_feature_sequence(first_tx, personal = True)

for sf in first_tx.subfeatures:
    if sf.type == "start_codon":
        sfstart = sf.start - first_tx.start
        sfend = sf.end - first_tx.start
        subseq = seq[sfstart-3:sfend+3]
        logger.info(subseq)
        logger.info(ggene.complement(subseq))

# nbs = len(seq)
# if gene.strand == '-':
#     newseq = ggene.complement(seq, rna = True)
# else:
#     newseq = ggene.to_rna(seq)

# donor_mi = d.finditer(newseq)
# for m in donor_mi:
#     logger.info('donor: '+str(m.span(0)))
    
# ba_mi = ba.finditer(newseq)
# for m in ba_mi:
#     logger.info('branch+acceptor: '+str(m.span(0)))
    
# nbs = first_tx.end - first_tx.start
# for sf in first_tx.subfeatures:
#     if sf.type == "intron":
#      logger.info(f'{nbs - sf.end_relative}, {nbs - sf.start_relative}')
    
# acc_mi = acc_ptrn.finditer(newseq)
# for m in acc_mi:
#     logger.info(str(m))

# geneman.assemble_gene("KISS1", 1)

#*******************************************************************************
# chrom = 1
# name = "KISS1"
# testgene = geneman.get_gene(name, chrom, min_qual = 5)

# locs.add('kiss1',testgene)

# var5 = testgene.get_feature('variant',queries = {'delta':lambda a:a>300})[0]

# varseq = geneman.get_feature_sequence(var5)
# res = geneman.get_feature_sequence(var5, upstream = 200, downstream = 200)
# cres = ggene.complement(res)
# cres_rna = ggene.complement(res, rna = True)
# ccres_rna = ggene.complement(cres_rna, rna = True)

# locs.add('var5seq',varseq)
# locs.add('moreseq',res)
# locs.add('compseq',cres)
# locs.add('compseq_rna',cres_rna)
# locs.add('seq_rna',ccres_rna)
# locs.add('var5',var5)

# tg = locals('KISS1')
# tg.add('KISS1',testgene)

# name = 'KISS1'
# testgene = geneman.assemble_gene(name, nchr, min_qual = 5)
# locs = locals('test5')
# locs.add('testgene',testgene.to_dict())


# ls = geneman.genemap.list_genes(nchr, start=0.9e8, end = 2e8)
# locs = locals('chr1')
# locs.add('chr1',ls)

# name = 'KISS1R'
# testgene = geneman.assemble_gene(name, nchr, min_qual = 5)
# locs = locals('kiss1r')
# locs.add('kiss1r',testgene.to_dict())

# name = 'RNR1'
# testgene = geneman.assemble_gene(name, nchr, min_qual = 5)
# locs = locals('rnr1')
# locs.add('rnr1',testgene.to_dict())

# name = 'NGF'
# testgene = geneman.assemble_gene(name, nchr, min_qual = 5)
# print(testgene)
# locs = locals('ngf')
# locs.add('ngf',testgene.to_dict())


# name = 'BDNF'
# testgene = geneman.assemble_gene(name, nchr, min_qual = 5)
# print(testgene)
# locs = locals(name)
# locs.add(name,testgene.to_dict())

# name = 'PEBP4'
# testgene = geneman.assemble_gene(name, nchr, min_qual = 5)
# print(testgene)
# locs = locals(name)
# locs.add(name,testgene.to_dict())

# name = 'SLC6A4'
# testgene = geneman.assemble_gene(name, nchr, min_qual = 5)
# print(testgene)
# locs = locals(name)
# locs.add(name,testgene.to_dict())

# name = 'PTEN'


# write_gene_data('HTR2A', 13)
# write_gene_data('HTR2B', 2)

# list all genes and count their features
#*******************************************************************************

# genelist = geneman.genemap.list_genes(nchr)
# genestrs = ['\t'.join([str(gitem) for gitem in gl]) for gl in genelist]
# nfeats = [(row[0],row[-1]) for row in genelist]
# nfeats = sorted(nfeats, key = lambda a:-a[1])

# print(nfeats[:10])

# locs = locals(fname = 'genelist')
# locs.add('genelist',genelist)


# summarize features of AKNAD1
#*******************************************************************************


# gene = geneman.genemap.get_gene('AKNAD1',1)['gene'][0]
# print(gene)
# featsumm = geneman.genemap.summarize_features(1, gene['start'], gene['end'])
# print(featsumm)
# locs = locals(fname = 'aknad1-features')
# locs.add('feature_summary',featsumm)


# testgene = 'THAP3'
# _gene = geneman.genemap.make_gene(testgene, nchr)
# # for ts in _gene.transcripts.values():
# #     print(ts.subfeatures)
# locals.add(testgene,_gene.to_dict())

# 4916519 vars total, 386873 in chr1


""" var:
ALT ['T']
FILTER None
FILTERS ['PASS']
FORMAT ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'GQ', 'PL', 'GP', 'PRI', 'SB', 'MB']
ID None
INFO <cyvcf2.cyvcf2.INFO object at 0x76a23cf2b840>
POS 10254
QUAL 4.329999923706055
REF TA
aaf 1.0
call_rate 1.0
end 10255
genotype [[1 1 0]]
genotypes [[1, 1, False]]
gt_alt_depths [10]
gt_alt_freqs [0.83333333]
gt_bases ['T/T']
gt_depths [12]
gt_phases [False]
gt_phred_ll_het [0]
gt_phred_ll_homalt [2]
gt_phred_ll_homref [60]
gt_quals [3.]
gt_ref_depths [2]
gt_types [3]
is_deletion True
is_indel True
is_mnp False
is_snp False
is_sv False
is_transition False
nucl_diversity 0.0
num_called 1
num_het 0
num_hom_alt 1
num_hom_ref 0
num_unknown 0
ploidy 2
start 10253
var_subtype del
var_type indel

"""

""" INFO:
('AC', 2)
('AF', 1.0)
('AN', 2)
('DP', 15)
('FS', 0.0)
('MQ', 25.979999542236328)
('MQRankSum', 0.9670000076293945)
('QD', 2.1700000762939453)
('ReadPosRankSum', 0.7519999742507935)
('SOR', 0.6930000185966492)
('FractionInformativeReads', 0.800000011920929)

"""