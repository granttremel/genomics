
from typing import Any, Optional, List, Dict
from dataclasses import dataclass
import time
import numpy as np

import pyhmmer
from pyhmmer.plan7 import HMM, Pipeline, Domain, TopHits, Hit, Builder, Background, Profile, Alignment
from pyhmmer.easel import DigitalSequence, DigitalSequenceBlock, MSA, DigitalMSA, MSAFile

from ggene.seqs.bio import reverse_complement
from ggene.config import DATA_DIR
from ggene import draw
from .motif import BaseMotif

import h5py

hmm_path = DATA_DIR / "dfam" / "human_repeats.hmm"

AB = pyhmmer.easel.Alphabet.dna()
BG = Background(AB, uniform = True)

@dataclass
class SearchResult:
    target_name:str
    query_name:str
    loge_val: float
    logp_val:float
    score:float
    norm_score:float
    
    chrom:str
    start:int
    end:int
    strand:str
    
    target_from:int
    target_to:int
    target_length:int
    
    query_from:int
    query_to:int
    query_length:int
    
    target_cov: float = None
    query_cov: float = None
    
    _target_seq:str = ""
    _query_seq:str = ""
    _id_seq:str = ""
    
    _best_domain:Any=None
    _alignment:Alignment = None
    
    @classmethod
    def from_hit(cls, hit:Hit, target_seq:str = "", chrom = "", offset = 0, strand = '+'):
        
        algn = hit.best_domain.alignment
        
        idseq = hit.best_domain.alignment.identity_sequence
        idlen = len(idseq.replace("+","").replace(" ",""))
        norm_score = hit.score / hit.best_domain.alignment.target_length
        
        target_cov = idlen/algn.target_length
        query_cov = idlen/algn.hmm_length
        
        return SearchResult(
            hit.name.decode(),
            hit.best_domain.alignment.hmm_name.decode(),
            np.log(hit.evalue),
            np.log(hit.pvalue),
            hit.score,
            norm_score,
            
            chrom,
            algn.target_from + offset,
            algn.target_to + offset,
            strand,
            
            algn.target_from,
            algn.target_to,
            algn.target_length,
            algn.hmm_from,
            algn.hmm_to,
            algn.hmm_length,
            
            target_cov = target_cov,
            query_cov = query_cov,
            
            _target_seq = target_seq,
            _id_seq = idseq,
            _alignment = algn
        )
    
    def __repr__(self):
        
        parts = []
        
        for k, v in self.__dict__.items():
            
            if k.startswith("_"):
                continue
            
            if isinstance(v, float):
                vstr = format(v,"0.3f")
            else:
                vstr = str(v)
            
            parts.append(f"{k}={vstr}")
            
            pass
        
        return f"SearchResult({", ".join(parts)})"

# todo
# class HMMMotif(BaseMotif):
    
#     def __init__(self, name, hmm, scoring_function, allow_rc = True, motif_class = ""):
#         super().__init__(name)
        
#         self.hmm = hmm
#         self.score_func = scoring_function
#         self.allow_rc = allow_rc
#         self.motif_class = motif_class
    
#     @property
#     def length(self):
#         return self.hmm.M
    
    
    
#     def __call__(self, seq):
#         return 0
    
#     def score(self, seq, return_positions=False):
#         return None
    
#     def count_instances(self, seq):
#         return 0
    
#     def find_instances(self, seq, threshold=None):
        
#         insts = []
        
#         res = search_hmm(self.hmm, [seq])
        
#         for hit in res:
            
#             if threshold and hit.score < threshold:
#                 continue
#             algn = hit.best_domain.alignment
#             insts.append((algn.target_from, algn.target_to, hit.score))
        
#         return insts

def get_hmmmotifs(names = [], classes = [], families = []):
    
    motifs = []
    
    hmms = load_hmms()
    
    for hmm in hmms:
        
        name = hmm.name.decode("utf-8")
        cn, fn = get_repeat_class(name)
        if classes and cn not in classes:
            continue
        
        if families and not fn in families:
            continue
        
        if names and name not in names:
            continue
        
        # new_motif = HMMMotif(name, hmm, None)
        # motifs.append(new_motif)
    
    return motifs


def load_hmms(hmm_path = hmm_path, name_filt = None):
    
    hmms = []
    
    name_filt = bytes(name_filt, encoding="utf-8") if name_filt else ""
    
    with pyhmmer.plan7.HMMFile(hmm_path) as hmf:
        for hh in hmf:
            if name_filt and not name_filt in hh.name:
                continue
            
            hmms.append(hh)
    return hmms

hmm_cls = ["Alu","LTR","MER", "L1", "L2", "L3", "L4", "L5", "CR1", "Charlie", "Tigger", "Arthur", "Zaphod", "Ricksha", "Kanga", "MamGyp", "MST","UCON", "HERV", "Eulor", "MLT", "tRNA", "Eutr", "EUTRE", "MARE", "MADE", "ERVL", "U", "HAL", "HSAT", "GSAT", "THE1", "HUERS", "X", "hAT", "SVA", "PRIMA", "MamTip", "MamRep", "MIR", "PABL", "HY", "Eut", "COMP", "BSR"]

def load_hmm_classes(classes = []):
    
    out_hmms = {}
    
    if not classes:
        classes = hmm_cls
    
    hmms = load_hmms()
    
    for hmm in hmms:
        name = hmm.name.decode("utf-8")
        
        cn, fn = get_repeat_class(name)
        
        if not cn in classes:
            continue
        
        if not cn in out_hmms:
            out_hmms[cn] = {}
        
        cnd = out_hmms[cn]
        
        if not fn in cnd:
            cnd[fn] = []
        
        cnd[fn].append(hmm)
    
    return out_hmms

def get_repeat_class(name):
    
    for cls_name in sorted(hmm_cls, key = lambda k:-len(k)):
        
        if name.startswith(cls_name):
            
            fam_name = name.removeprefix(cls_name)
            
            return cls_name, fam_name
    
    return "other", ""

def load_alus():
    return load_hmms(name_filt = "Alu")

def load_line1s():
    return load_hmms(name_filt = "L1")
    
def check_truncation(repeat, hmm):
    
    rpt_len = repeat.end - repeat.start
    rpt_hmm_start = repeat.attributes.get("start_hmm", 0)
    rpt_hmm_end = repeat.attributes.get("end_hmm", 0)
    
    hmm_len = len(hmm.consensus)
    
    if repeat.strand == '-':
        rpt_hmm_start, rpt_hmm_end = rpt_hmm_end, rpt_hmm_start
    
    tr_5p = rpt_hmm_start / hmm_len
    tr_3p = 1 - rpt_hmm_end / hmm_len
    
    # if repeat.strand == '-':
    #     tr_5p, tr_3p = tr_3p, tr_5p
    
    repeat.attributes["motif_len"] = hmm_len
    repeat.attributes["motif_start"] = repeat.start - rpt_hmm_start
    repeat.attributes["motif_end"] = repeat.end + rpt_hmm_end
    repeat.attributes["5p_truncation"] = tr_5p
    repeat.attributes["3p_truncation"] = tr_3p
    
    return repeat

def prepare_seqs(seqs, do_rc = False):
    
    if not isinstance(seqs, dict):
        seqs = {str(i):seqs[i] for i in range(len(seqs))}
    
    dseqs = []
    names = []
    for name, seq in seqs.items():
        dseq = DigitalSequence(AB, sequence = AB.encode(seq), name = name.encode())
        dseqs.append(dseq)
        names.append(name)
    
    db = DigitalSequenceBlock(AB, dseqs)
    
    return db, names

def search_hmm(hmm, seqs):
    
    dseqs = []
    for seq in seqs:
        
        seq_enc = AB.encode(seq)
        
        dseq = pyhmmer.easel.DigitalSequence(AB, sequence = seq_enc)
        dseqs.append(dseq)
    
    
    dsb = pyhmmer.easel.DigitalSequenceBlock(AB, iterable = dseqs)
    ppl = pyhmmer.plan7.Pipeline(hmm.alphabet)
    res = ppl.search_hmm(hmm, dsb)
    
    return res

def print_repeats(rpts, start_offset = 0):
    
    for rpt in rpts:
        
        parts = []
        
        parts.append(f"repeat_len={rpt.end - rpt.start}")
        
        motif_len = rpt.get("motif_len", 0)
        if motif_len:
            parts.append(f"motif_len={motif_len}")
        
        for k in ["bits","e-value","bias","kimura_div", "start_hmm", "end_hmm"]:
            v = rpt.attributes.get(k, None)
            if v:
                parts.append(f"{k}={v}")
        
        for k in ["5p_truncation", "3p_truncation"]:
            
            v = rpt.attributes.get(k, 0)
            if v > 0.1:
                parts.append(f"{k}={v:0.0%}")
        
        start_hmm = rpt.attributes.get("start_hmm",0)
        end_hmm = rpt.attributes.get("end_hmm", 0)
        
        inst_pos_str = f"{rpt.chr}:{rpt.start}-{rpt.end}"
        # inst_pos_off_str = f"{rpt.start-start_offset}-{rpt.end-start_offset}"
        # hmm_pos_str = f"{start_hmm}-{end_hmm}/{motif_len}"
    
        print(f"{rpt.name} at {inst_pos_str}, ({", ".join(parts)})")

def search_hmms(seqs:Dict[str,str], hmms, norm_score_thresh = 0.5):
    
    scores = {}
    
    dseqs = []
    names = []
    for name, seq in seqs.items():
        dseq = DigitalSequence(AB, sequence = AB.encode(seq), name = name.encode())
        dseqs.append(dseq)
        names.append(name)
        
    db = DigitalSequenceBlock(AB, dseqs)
    
    ppl = Pipeline(AB)
    
    for hmm_name, hmm in hmms.items():
        
        res = ppl.search_hmm(hmm, db)
        
        if not res:
            continue
        
        max_res = max(res, key = lambda r:r.score)
        if max_res:
            target_seq = seqs.get(max_res.name.decode(), "")
            sres = SearchResult.from_hit(max_res, target_seq = target_seq)
            
            if sres.norm_score < norm_score_thresh:
                continue
            
            if not hmm.name in scores:
                scores[hmm_name] = []
            scores[hmm_name].append(sres)
        
    return scores
    

def _make_msa(name:str, seqs:Dict[str,str], preserve_gaps:bool = False):
    """Create MSA from sequences.

    Args:
        name: MSA name
        seqs: Dictionary of sequences (can be aligned with gaps or unaligned)
        preserve_gaps: If True, sequences are already aligned and gaps should be preserved

    Returns:
        DigitalMSA object
    """
    if not preserve_gaps:
        # Simple case: unaligned sequences
        dsb, names = prepare_seqs(seqs)
        msa = DigitalMSA(AB, sequences = dsb, name=name.encode())
        return msa

    # For aligned sequences with gaps, write to file and read back
    # pyhmmer's file parsers handle gaps correctly
    from pathlib import Path
    import tempfile

    # Write to Stockholm format (supports gaps)
    temp_file = Path(tempfile.gettempdir()) / f"{name}_temp.sto"

    with open(temp_file, 'w') as f:
        f.write("# STOCKHOLM 1.0\n")

        # Find max name length for formatting
        max_name_len = max(len(n) for n in seqs.keys())

        # Write aligned sequences
        for seq_name, seq in seqs.items():
            f.write(f"{seq_name:<{max_name_len}}  {seq}\n")

        f.write("//\n")

    # Read back with pyhmmer's parser (preserves gaps)
    with MSAFile(temp_file, digital=True, alphabet=AB) as msa_file:
        msa = msa_file.read()
        msa.name = name.encode()

    # Cleanup
    temp_file.unlink()

    return msa

def make_msa(name:str, seqs:Dict[str,str]):
    
    dsb, names = prepare_seqs(seqs)
    msa = DigitalMSA(AB, name.encode(), sequences=dsb)
    
    return msa
    

def print_msa(msa:DigitalMSA, max_width=80, show_positions=True, file = None):
    
    # Get sequences as strings
    sequences = []
    names = []
    
    for i in range(len(msa.sequences)):
        seq = msa.alignment[i]
        seq_str = AB.decode(bytes(seq))
        sequences.append(seq_str)
        names.append(msa.sequences[i].name.decode() if msa.sequences[i].name else f"seq_{i}")

    # Print in chunks
    max_name_len = max(len(n) for n in names)
    seq_len = len(sequences[0])

    for start in range(0, seq_len, max_width):
        end = min(start + max_width, seq_len)

        if show_positions:
            print(f"\nPosition {start+1}-{end}:", file = file)
            print("=" * (max_name_len + max_width + 5), file = file)

        for name, seq in zip(names, sequences):
            chunk = seq[start:end]
            print(f"{name:<{max_name_len}}  {chunk}", file = file)
        print(file = file)

def msa_stats(msa):
    """Print MSA statistics"""
    num_seqs = len(msa.sequences)
    algn_len = len(msa.sequences[0])
    print(f"Number of sequences: {len(msa.sequences)}")
    print(f"Alignment length: {len(msa.alignment[0]) if msa.alignment else 0}")
    print(f"MSA name: {msa.name.decode() if msa.name else 'unnamed'}")

    # Conservation at each position
    if len(msa.alignment) > 0:
        align_seqs = [[c for c in AB.decode(seq)]
                            for seq in msa.alignment]
        
        gap_cts = []
        cons = []
        for i in range(algn_len):
            gap_ct = 0
            for seq in align_seqs:
                gap_ct += seq[i] == '-'
            
            gap_cts.append(gap_ct)
            cons.append(gap_ct==0)
        
        print(f"Positions with no gaps: {sum(cons)}/{algn_len}")
        print(f"Average gaps per position: {np.mean(gap_cts):.2f}")


def build_single_hmm(seq):
    
    db, names = prepare_seqs([seq])
    dseq = db[0]
    
    builder = Builder(AB)
    new_hmm, _, _ = builder.build(dseq, BG)
    
    return new_hmm

def build_hmm(msa=None, seqs:Dict[str,str] = {}):
    
    if not msa:
        msa = make_msa(seqs)
    
    builder = Builder(AB)
    hmm, prof, optprof = builder.build_msa(msa, BG)
    
    return hmm, prof, optprof
    

def search_genome(gm, hmm, chromes = None, start = 1e6, stop = None, batch_size = 64, batch_len = 2**16, min_score = None, max_e_value = 0.05, limit = None, min_query_cov = None, do_rc = True):
    
    if not chromes:
        chromes = [chrome for chrome in gm.iter_chromes()]
    
    if not min_score:
        min_score = 0.95 * hmm.M
    
    
    print(f"using min_score {min_score:0.1f}")
    
    pos = start
    
    ppl = Pipeline(AB, BG, T = min_score, E = max_e_value)
    
    batch_num = 0
    num_hits = 0
    
    for chrom in chromes:
        
        chr_len = gm.get_chrom_len(chrom)
        
        print(f"searching chr{chrom}")
        
        full_seq = "hi"
        
        while full_seq:
            
            t0 = time.perf_counter()
            full_seq = gm.get_sequence(chrom, pos, pos+batch_size*batch_len)
            dt = time.perf_counter() - t0
            print(f"took {1000*dt:0.2f}ms to obtain sequences length {batch_size*batch_len}")
            
            seqs = []
            rcseqs =[]
            for nb in range(batch_size):
                
                seq = full_seq[nb*batch_len:(nb+1)*batch_len]
                seqs.append(seq)
                if do_rc:
                    rcseqs.append(reverse_complement(seq))
            
            seqs += rcseqs
            
            dsb, names = prepare_seqs(seqs, do_rc = do_rc)
            
            t0 = time.perf_counter()
            ths = ppl.search_hmm(hmm, dsb)
            dt = time.perf_counter() - t0
            print(f"took {dt:0.2f}s to search sequences with hmm")
            
            
            prog = pos / chr_len
            print(f"obtained {len(ths)} hits from batch {batch_num} at chr{chrom}:{pos} ({prog:0.2%})")
            
            for i, th in enumerate(ths):
                    
                strand = "+"
                bn = int(th.name)
                
                # if bn < 0:
                #     strand = "-"
                #     bn = abs(bn)
                    
                if bn > batch_size:
                    bn //= 2
                    strand = "-"
                
                sr = SearchResult.from_hit(th, chrom=chrom, offset = pos + bn*batch_len, strand=strand)
                
                if min_query_cov and sr.query_cov < min_query_cov:
                    continue
                
                tgt_seq = seqs[int(th.name)][sr.target_from:sr.target_to]
                sr._target_seq = tgt_seq
                
                num_hits += 1
                
                yield sr
                    
            pos += batch_size*batch_len
            batch_num += 1
            
            if stop and pos > stop:
                break
            if limit and num_hits > limit:
                break
            
        if limit and num_hits > limit:
            break
            


