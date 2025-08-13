



from . import COMPLEMENT_MAP, CODON_TABLE_DNA, CODON_TABLE_RNA, complement, to_rna
from .features import Feature, Gene
import itertools

import logging
logger = logging.getLogger(__name__)

codon_table_dna = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
start_codon = 'AUG'
stop_codons = ('UAA','UGA','UAG')

codon_table_rna = {to_rna(k):v for k,v in codon_table_dna.items()}

class Ribosome:
    """
    RNA->protein
    """
    
    def __init__(self):
        
        pass
    
    def get_coding_sequence(self, transcript, sequence_generator):
        
        starts = []
        stops = []
        mrna_seqs = []
        
        for feat in transcript.subfeatures:
            if feat.type == "start_codon":
                starts.append(feat.five_prime())
            elif feat.type == "stop_codon":
                stops.append(feat.three_prime())
            elif feat.type == "cds":
                seq = sequence_generator(feat)
                mrna_seqs.append(seq)
        
        return ','.join(mrna_seqs), starts, stops
        
    def from_transcript(self, transcript:Feature, sequence_generator):
        
        mrna_seq, starts, stops = self.get_coding_sequence(transcript, sequence_generator)
        
        if len(starts) > 1:
            logger.debug(f"multiple starts found! ({len(starts)})")
        elif len(starts) < 1:
            #find it manually?
            pass
        if len(stops) > 1:
            logger.debug(f"multiple stops found! ({len(stops)})")
        elif len(stops) < 1:
            #find it manually?
            pass
        
        possibles = [(a,b) for a,b in itertools.product(starts,stops) if (a - b)%3 == 0]
        
        if len(possibles) > 1:
            logger.debug(f'multiple possible start-stop pairs found! ({len(possibles)})')
        elif len(possibles) < 1:
            logger.debug(f'no possible start-stop pairs detected!')
            return ''
        
        if transcript.strand == '-':
            mrna_seq = complement(mrna_seq, rna = True)
        else:
            mrna_seq = to_rna(mrna_seq)
            
        aas = self.translate(mrna_seq, strand=transcript.strand)
            
        return aas
    
    def translate(self, sequence, start = 0, stop = 0):
        
        pseq = []
        i = start
        first = sequence[i:i+3]
        if not first == start_codon:
            logger.debug(f'start codon not found at {start} ({first})')
        pseq.append(codon_table_rna[codon])
        
        go = True
        while go:
            i+=3
            codon = sequence[i:i+3]
            pseq.append(codon_table_rna[codon])
            
        return ','.join(pseq)

    def compare_variant(variant, transcript, sequence_generator, border = 5):
        
        res = False
        variants = transcript.get_feature("variant")
        start = transcript.get_feature("start_codon")
        
        shift = 0
        for v in variants:
            if v.start >= variant.start:
                break
            
            shift += v.delta
        
        shift = shift % 3
        
        varseq = sequence_generator(variant.start-3*border, variant.end+3*border)
        
        
        
        
        pass
    
    def find_reading_frame(self, sequence):
        
        for i in range(len(sequence)-2):
            
            pass
        pass

    def find_neighbors(self, codon):
        
        pass
