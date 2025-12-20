
from ggene import get_paths
DEFAULT_VCF_PATH, DEFAULT_GTF_PATH, DEFAULT_FASTA_PATH, DEFAULT_LIBRARY = get_paths()
from ggene.database.unified_stream import UnifiedGenomeAnnotations
from ggene.database.genomemanager import GenomeManager
import os
    
# Create unified annotation system with sequence streaming
annotations = UnifiedGenomeAnnotations(
    fasta_path=DEFAULT_FASTA_PATH,
    vcf_path=DEFAULT_VCF_PATH,
)

annotations.add_gtf(DEFAULT_GTF_PATH, "genes")
print("✓ Added GTF annotations")

seqstream = annotations.sequence_stream

kchr = 1
kstart = 204190241 
kend = 204196591

print(seqstream.get_sequence_with_alignment(kchr,kstart,kend))
# print(seqstream.get_variants_in_region(kchr,kstart,kend))

