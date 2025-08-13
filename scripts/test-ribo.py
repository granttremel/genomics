
from ggene.genomemanager import GenomeManager
from ggene.genome_browser import browse_genome

geneman = GenomeManager()

k1 = geneman.load_gene("KISS1")
seq1 = geneman.transcribe_gene("MDM4", transcript_index=0, verify_start_stop=False)
seq2 = geneman.transcribe_gene("MDM4", transcript_index=0, verify_start_stop=False, personal=False)
print(seq1['protein_sequence'])
print(seq2['protein_sequence'])

