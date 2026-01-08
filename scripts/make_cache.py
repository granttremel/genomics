
import numpy as np

from ggene.config import other_paths
from ggene.database.genome_manager import GenomeManager
from ggene.database.cache import GenomeCache
from ggene.database.ufeature import UFeat

from ggene.genome import genes, ncrna

from ggene.seqs import process, bio, vocab, lambdas

from ggene import draw
from ggene.draw.scalar_plot import ScalarPlot


def load_genome():
    return GenomeManager(skip_clinvar = True)


def test_dfam(gm):
    
    dfam = gm.annotations.streams.get("dfam")
    print(dfam, dfam.filepath)
    n = 0
    for rpt in dfam.stream("3", start = int(10e6), end = int(11e6)):
        print(rpt)
        n += 1
        if n >= 128:
            break


def main():
    
    gm = load_genome()
    
    gc = GenomeCache(gm)
    
    # ssps = ['gc', 'at', 'ag', 'ac', 'cpg', 'cpg_to_gc', 'polya', 'polyt', 'polyg', 'polyc', 'polyy', 'polyr', 'polys', 'polyw']
    ssps = ['features', 'genes', 'exons', 'motifs', 'pseudo', 'lncRNA', 'ncRNA', 'nongenes', 'simple_repeats', 'Alu', 'AluJ', 'AluY', 'AluS', 'Line1', 'L1HS', 'L1PA', 'L1P', 'L1M', 'L2', 'LTR', 'MIR', 'SVA']

    gc.cache_all_chromes(ssps)
    # gc.cache_chromes(['Y'], ssps)



if __name__=="__main__":
    main()


