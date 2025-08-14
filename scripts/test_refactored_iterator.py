from cyvcf2 import Variant
from ggene.features import Feature
from ggene.genomemanager import GenomeManager
from ggene.genemap import GeneMap
from ggene.genome_iterator import GenomeIterator

gm = GenomeManager()

gmap = GeneMap()



start = 1000000
end = 2000000
geneiter = GenomeIterator(gm, 1, start, end, iteration_type="sequence",direction=1, window_size = 160, strand='+',as_rna=False)


while True:
    refseq, altseq = geneiter.get_sequence_window()
    print(refseq)
    
    # comp = ''.join(['X' if not a==b else ' ' for a,b in zip(refseq,altseq)])
    # print(comp)
    # if '-' in refseq or '-' in altseq:
    #     print('- found!')
    
    
    