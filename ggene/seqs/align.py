
from typing import Dict
import random


from pathlib import Path
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt

from ggene import draw
from ggene.display.colors import Colors
from ggene.seqs import process, bio, heal

from Bio import Align

class AlignmentResult:
    
    def __init__(self, alignment:Align.Alignment):
        
        self.length = alignment.length
        self.target = alignment.target
        self.query = alignment.query
        
        cts = alignment.counts()
        self.gaps, self.identities, self.mismatches = cts.gaps, cts.identities, cts.mismatches
        self.score = alignment.score
        self.ham = self.length - self.identities
        
        self.target_algn, self.query_algn = alignment[0], alignment[1]
        self.cons = bio.merge_with_gaps(self.target_algn, self.query_algn, default = '-')
        
        self.subs = alignment.substitutions
        self.t_algn, self.q_algn = alignment.aligned
        
        self._alignment = alignment
    
    @property
    def consensus(self):
        return self.cons
    
    def get_aligned_seqs(self):
        return self.target_algn, self.query_algn
    
    def get_mutation_spectrum(self):
        
        spec = {}
        
        for b, row in zip("ATGC", self.subs):
            rsum = sum(row)
            for bb, v in zip("ATGC", row):
                
                if b==bb:
                    continue
                spec[bb+b] = v/rsum
        
        return spec
        
    def print(self, chunksz = 128, emph_match = True, emph_subs = True, emph_indels = False, color_subs = False, show_middle = False, show_consensus = False, file=None):
        
        print(Colors.BOLD + f"Score: {self.score:0.1f} with {self.identities} identities, {self.gaps} gaps, {self.mismatches} mismatches" + Colors.RESET, file=file)
        print(file=file)
        
        tgt_algn, q_algn = self.get_aligned_seqs()
        print_aligned_seqs(tgt_algn, q_algn, chunksz = chunksz, consensus = self.consensus if show_consensus else "", emph_match=emph_match,  emph_subs=emph_subs,  emph_indels=emph_indels, color_subs = color_subs, show_middle=show_middle, file=file)
    
    def __len__(self):
        return self.length

import itertools
pair_cost_dna = {(a,b):float(a==b) for a, b in itertools.product('ATGC','ATGC')}

pair_cost = {
    ("A","T"):0.0,
    ("G","C"):0.0,
    ("T","A"):0.0,
    ("C","G"):0.0,
    
    # G-U, common
    ("T","G"):0.2,
    ("G","T"):0.2,
    
    # G-A, common
    ("A","G"):0.4,
    ("G","A"):0.4,
    
    ("A","A"):1.0,
    ("T","T"):1.0,
    ("G","G"):1.0,
    ("C","C"):1.0,
    
    ("A","C"):1.0,
    ("T","C"):1.0,
    ("C","A"):1.0,
    ("C","A"):1.0,
}

_all_scores = {
    
    "open_left_deletion_score":0,
    "open_internal_deletion_score":-1,
    "open_right_deletion_score":0,
    
    "open_left_insertion_score":0,
    "open_internal_insertion_score":-1,
    "open_right_insertion_score":0,
    
    "extend_left_deletion_score":0,
    "extend_internal_deletion_score":-1,
    "extend_right_deletion_score":0,
    
    "extend_left_insertion_score":0,
    "extend_internal_insertion_score":-1,
    "extend_right_insertion_score":0,
    
}

_default_scores = {
    "match_score":1,
    "mismatch_score":-0.9,
    
    "left_gap_score":-1,
    "right_gap_score":-1,
    "internal_gap_score":-1,
}

def get_scores(**scores):
    sc = _default_scores.copy()
    sc.update(scores)
    # return sc
    return {}

def get_consensus_scores():
    cscores = _all_scores.copy()
    cscores.update(
        open_internal_deletion_score=-2,        
        open_internal_insertion_score=-2,        
        extend_internal_deletion_score=-1,        
        extend_internal_insertion_score=-1,
    )
    return cscores

def score_sequences(seqa, seqb, **scores):
    scores = get_scores(**scores)
    score = Align.PairwiseAligner(**scores).score(seqa, seqb)
    return score

def score_consensus(cons, seqb):
    scores = get_consensus_scores()
    return score_sequences(cons, seqb, **scores)

def align_sequences(seqa, seqb, topk = 1, **scores):
    scores = get_scores(**scores)
    
    aligner = Align.PairwiseAligner(**scores)
    res = aligner.align(seqa, seqb)
    
    return [AlignmentResult(res[k]) for k in range(topk)]

def align_consensus(seqa, seqb, topk = 1):
    scores = get_consensus_scores()
    return align_sequences(seqa, seqb, topk = topk, **scores)

def align_multi(seqs, topk = 1):
    
    min_len = min([len(s) for s in seqs])
    seqs_min = [s[:min_len] for s in seqs]
    
    algnr = Align.MultipleSeqAlignment(seqs_min)
    
    return algnr


def align_sequences_muscle(seqs: Dict[str, str], name: str = "msa"):
    """Align sequences using MUSCLE"""

    # Write sequences to temp FASTA
    temp_in = Path(f"/tmp/{name}_input.fasta")
    temp_out = Path(f"/tmp/{name}_aligned.fasta")

    with open(temp_in, 'w') as f:
        for seq_name, seq in seqs.items():
            f.write(f">{seq_name}\n{seq}\n")

    # Run MUSCLE
    subprocess.run([
        "muscle",
        "-in", str(temp_in),
        "-out", str(temp_out)
    ], check=True, capture_output = True)

    # Read aligned FASTA back
    aligned_seqs = {}
    with open(temp_out) as f:
        current_name = None
        current_seq = []

        for line in f:
            if line.startswith('>'):
                if current_name:
                    aligned_seqs[current_name] = ''.join(current_seq)
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_name:
            aligned_seqs[current_name] = ''.join(current_seq)

    # Cleanup
    temp_in.unlink()
    temp_out.unlink()

    # return msa
    return aligned_seqs

def trim_msa(msa:Dict[str,str], min_occ = 2):
    
    seq_len = len(next(iter(msa.values())))
    
    names = list(msa.keys())
    
    start = 0
    end = seq_len
    
    for i in range(seq_len):
        num_bs = sum([1 for n in names if msa[n][i] != '-'])
        if num_bs < min_occ:
            start = i+1
        else:
            break
        
    for i in range(seq_len-1, -1, -1):
        
        num_bs = sum([1 for n in names if msa[n][i] != '-'])
        if num_bs < min_occ:
            end = i
        else:
            break
    
    print(f"chose start {start} and end {end}")
    
    trimmed_seqs = {n:seq[start:end] for n, seq in msa.items()}
    return trimmed_seqs
    

def print_aligned_seqs(algn_target, algn_query, 
                        chunksz = 128, 
                        consensus = "",
                        emph_match = False, 
                        emph_subs = False, 
                        emph_indels = False, 
                        color_subs = False,
                        show_middle = False, 
                        show_subs = False,
                        file = None,
                    ):
    
    line_frm = "{:<4}{:>4} {}"
    
    num_chunks = len(algn_target)//chunksz + 1
    for n in range(num_chunks):
        sub_tgt = algn_target[n*chunksz:(n+1)*chunksz]
        sub_q =  algn_query[n*chunksz:(n+1)*chunksz]
        ctgt, cquery = get_colored_aligned_seqs(sub_tgt, sub_q, emph_match=emph_match,  emph_subs=emph_subs,  emph_indels=emph_indels, color_subs = color_subs)
        print(line_frm.format("T", n*chunksz, ctgt), file=file)
        
        if show_middle:
            mid = get_middle(sub_tgt, sub_q, show_subs = show_subs)
            print(line_frm.format("", "", mid))
        
        print(line_frm.format("Q", "", cquery), file=file)
        
        if consensus:
            sub_cons = consensus[n*chunksz:(n+1)*chunksz]
            print(line_frm.format("C","", sub_cons), file=file)
        
        print(file=file)
    
def get_colored_aligned_seqs(tgt, query, emph_match = True, emph_subs = False, emph_indels = False, color_subs = False, allow_alias = False):
    
    curr_tcol = Colors.SUBTLE
    curr_qcol = Colors.SUBTLE
    
    tgt_col = [curr_tcol]
    q_col = [curr_qcol]
    
    cmap = get_colormap(emph_match=emph_match, emph_subs=emph_subs, emph_indels=emph_indels, color_subs = color_subs)
    
    for t, q in zip(tgt, query):
        ct, cq = get_aligned_color(t, q, colormap = cmap, allow_alias = allow_alias)
        
        if ct != curr_tcol:
            tgt_col.append(Colors.RESET + ct)
            curr_tcol = ct
        
        if cq != curr_qcol:
            q_col.append(Colors.RESET + cq)
            curr_qcol = cq
        
        tgt_col.append(t)
        q_col.append(q)
    
    tgt_col.append(Colors.RESET)
    q_col.append(Colors.RESET)
    
    return "".join(tgt_col), "".join(q_col)

css = "\x1b[38;5;232m"
cs = "\x1b[38;5;238m"
cm = "\x1b[38;5;250m"
csub = "\x1b[38;5;14m"

cry = "\x1b[38;5;140m"
csw = "\x1b[38;5;160m"
ckm = "\x1b[38;5;112m"

cins = "\x1b[38;5;94m"
cdel = "\x1b[38;5;88m"

mtn_syms = {
    "C":"≬",
    "I":"⌇",
    "V":"∫",
}

def get_middle(tgt, query, show_subs = False):
    
    mid = [css]
    for t, q in zip(tgt, query):
        if t==q:
            mid.append("|")
        elif t=='-' or q=='-':
            mid.append(" ")
        elif show_subs:
            m = heal.identify_mutation_type(t, q)
            sym = mtn_syms.get(m)
            mid.append(sym)
        else:
            mid.append(" ")
    mid.append(Colors.RESET)
    return "".join(mid)

def get_colormap(emph_match, emph_subs, emph_indels, color_subs):
    
    colormap = {k:cs for k in ["","m","ins","del","C","I","V"]}
    if emph_match:
        colormap["m"] = cm
    if emph_subs:
        if color_subs:
            colormap["C"] = csw
            colormap["I"] = cry
            colormap["V"] = ckm
        else:
            colormap["C"] = csub
            colormap["I"] = csub
            colormap["V"] = csub
            
    if emph_indels:
        colormap["ins"] = cins
        colormap["del"] = cdel
    
    return colormap

def get_aligned_color(tgt_base, q_base, colormap = {}, allow_alias = False):
    
    if allow_alias:
        tgt_ali = bio.ALIASES.get(tgt_base, '-')
        q_ali = bio.ALIASES.get(q_base, '-')
        same = bool(set(q_ali).intersection(tgt_ali))
    else:
        same = tgt_base == q_base
    
    if same:
        return colormap.get("m", css), colormap.get("m", css)
    elif tgt_base == '-':
        return colormap.get("", css), colormap.get("ins", css)
    elif q_base == '-':
        return colormap.get("del", css), colormap.get("", css)
    else:
        mt = heal.identify_mutation_type(tgt_base, q_base)
        return colormap.get(mt, css), colormap.get(mt, css)

def get_substitution_color(tgt_base, q_base):
    
    m = heal.identify_mutation(tgt_base, q_base)
    
    if m in "SW":
        return csw, csw
    elif m in "YR":
        return cry, cry
    elif m in "MK":
        return ckm, ckm
    else:
        return cm, cm
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
############## just some stuff #######################


class IndexMap:
    
    def __init__(self, seq_len):
        self.global_seq_len = seq_len
        self.seq_len = self.global_seq_len
        self.mut_seq_len = self.global_seq_len
        self.mapping = None
        self.pdls = []
        
    def add_delta(self, pos, delta):
        if abs(delta) < 1:
            return
        mapped_pos = pos
        self.pdls.append((mapped_pos, delta))
        
        self.mut_seq_len += delta
        # self.global_seq_len += abs(delta)
        self.global_seq_len += max(delta, 0)
    
    def global_to_init(self, i):
        
        ii = i
        for p, d in self.pdls:
            dd = max(0, d)
            if ii >= p and ii < p + dd:
                ii = p
            elif ii >= p + dd:
                ii -= dd
        
        return ii
    
    def init_to_global(self, i):
        ii = i
        for p, d in reversed(self.pdls):
            dd = max(0, d)
            if ii >= p:
                ii += dd
        
        return ii
    
    def global_to_mut(self, i):
        ii = i
        for p, d in self.pdls:
            dd = min(0, d)
            if ii >= p and ii < p - dd:
                ii = p
            elif ii >= p - dd:
                ii += dd
        
        return ii
    
    def mut_to_global(self, i):
        ii = i
        for p, d in reversed(self.pdls):
            dd = min(0, d)
            if ii >= p:
                ii -= dd
        
        return ii
    
    # def init_to_mut(self, i):
    #     pass

    def unmutate(self, mutated_seq):
        """Convert mutated sequence to global representation with gaps for deletions."""
        if not self.pdls:
            return mutated_seq

        result = []
        mut_pos = 0

        # Sort mutations by position for proper ordering
        sorted_pdls = sorted(self.pdls)

        # Track cumulative effect of mutations
        cumulative_deletion = 0

        for orig_pos, delta in sorted_pdls:
            # Adjust position for previous deletions
            adjusted_pos = orig_pos - cumulative_deletion

            # Add characters before this mutation
            if mut_pos < adjusted_pos and mut_pos < len(mutated_seq):
                result.extend(mutated_seq[mut_pos:min(adjusted_pos, len(mutated_seq))])
                mut_pos = min(adjusted_pos, len(mutated_seq))

            if delta > 0:
                # Insertion: copy the inserted characters
                if mut_pos < len(mutated_seq):
                    result.extend(mutated_seq[mut_pos:min(mut_pos + delta, len(mutated_seq))])
                    mut_pos = min(mut_pos + delta, len(mutated_seq))
            else:
                # Deletion: add gaps for deleted characters
                result.extend(['-'] * (-delta))
                cumulative_deletion += (-delta)

        # Add remaining characters
        if mut_pos < len(mutated_seq):
            result.extend(mutated_seq[mut_pos:])

        return ''.join(result)

    def mutate(self, init_seq):
        """Convert initial sequence to global representation with gaps for insertions."""
        if not self.pdls:
            return init_seq

        result = []
        init_pos = 0

        # Sort mutations by position for proper ordering
        sorted_pdls = sorted(self.pdls)

        for mut_pos, delta in sorted_pdls:
            # Add characters before this mutation
            if init_pos < mut_pos and init_pos < len(init_seq):
                result.extend(init_seq[init_pos:min(mut_pos, len(init_seq))])

            if delta > 0:
                # Insertion: add gaps in the global representation
                result.extend(['-'] * delta)
                init_pos = mut_pos
            else:
                # Deletion: include deleted characters in global representation
                if mut_pos < len(init_seq):
                    result.extend(init_seq[mut_pos:min(mut_pos - delta, len(init_seq))])
                init_pos = mut_pos - delta

        # Add remaining characters
        if init_pos < len(init_seq):
            result.extend(init_seq[init_pos:])

        return ''.join(result)
    
    def plot(self, xdata = "g", ydata = ["i","m"], fname = ""):
        
        
        iglb = range(self.global_seq_len)
        iinit = [self.global_to_init(i) for i in iglb]
        imut = [self.global_to_mut(i) for i in iglb]
        
        ddict = {"g":iglb, "i":iinit, "m":imut}
        legdict = {"g":"Global","i":"WT","m":"Mutant"}
        
        xd = ddict.get(xdata)
        xlabel = f"{legdict.get(xdata)} Coordinate"
        
        leg = []
        
        f, ax = plt.subplots()
        for yk in ydata:
            ax.plot(xd, ddict.get(yk))
            leg.append(legdict.get(yk))
            
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Mapped coordinate")
        ax.legend(leg)
        
        if fname:
            fpath = os.path.join("./data", fname)
            f.savefig(fpath)
        
        return f, ax


def initialize_align(seqa, seqb, length = 3):
    
    prba = seqa[:length]
    prbb = seqb[:length]
    abstart =  len(seqa)
    bastart = len(seqb)
    
    for i in range(len(seqb) - length):
        bsub = seqb[i:i+length]
        if prba == bsub:
            abstart = i
            break
    
    for i in range(len(seqa) - length):
        asub = seqa[i:i+length]
        if prbb == asub:
            bastart = i
            break
    
    return abstart, bastart

def get_align_checkpoints(seqa, seqb, err_ratio = 0.2, num_checkpoints = 2, max_cons_entr = 0.07):
    
    seq_len = len(seqa)
    half_len = seq_len//2
    max_err = int(err_ratio*np.sqrt(seq_len))
    
    runs, inds, shifts, errs = process.correlate_longest_subseq_err(seqa, seqb, fill = None, max_err = max_err)
    
    tops, topdatas = process.extract_max_runs(seqa, seqb, runs, inds, shifts, 2*num_checkpoints)
    init_checkpts = []
    entrs = []
    
    for i in range(len(tops)):
        ssa, ssb = tops[i]
        cons = bio.merge(ssa, ssb)
        entr = bio.consensus_entropy(cons)
        
        if entr > max_cons_entr:
            break
        
        entrs.append(entr)
        
        r, ind, sh = topdatas[i]
        init_checkpts.append((ind, ind + r, ind - sh, ind - sh + r, entr))
    
    checkpts = []
    
    for i in range(len(init_checkpts)):
        
        mina, maxa, minb, maxb, entr = init_checkpts[i]
        keep = True
        for j in range(i+1, len(init_checkpts)):
            minaa, maxaa, minbb, maxbb, entrr = init_checkpts[j]
            overlap = min(maxa, maxaa) - max(minaa, mina)
            
            if overlap < 0:
                continue
            
            if entrr*(maxaa-minaa) < entr * (mina - maxa):
                print(f"hmm idk. {entr:0.3f} vs {entrr:0.3f}")
                keep = False
        
        if keep:
            checkpts.append(init_checkpts[i])
        
        if len(checkpts) == num_checkpoints:
            break
    
    # checkpts = list(sorted(checkpts, key = lambda k:k[4]*(k[1]-k[0]))) # sort by entropy or by start ind?
    checkpts = list(sorted(checkpts, key = lambda k:k[0])) # sort by entropy or by start ind?
    
    return checkpts
    

def align_sequences_bad(seqa, seqb):
    # brute force
    
    alen = len(seqa)
    blen = len(seqb)
    i = j = 0
    diffs = []
    
    iab, iba = initialize_align(seqa, seqb, length=3)
    if iab < blen and iab < iba:
        j = iab
        diffs.append((0, iab))
    elif iba < alen and iba < iab:
        i = iba
        diffs.append((iba, 0))
    
    algn = True
    while True:
        
        ba = seqa[i]
        bb = seqb[j]
        
        if ba == bb:
            if not algn:
                algn = True
                diffs.append((i, j))
            i+=1
            j+=1
        else:
            j += 1
            algn = False
        
        if i >= alen or j >= blen:
            break
    diffs.append((alen, blen))
    
    return diffs

def align_sequences_chkpt(seqa, seqb, checkpoints):
    
    alen = len(seqa)
    blen = len(seqb)
    nchk = 0
    curr_chk = checkpoints[nchk]
    i = curr_chk[0]
    j = curr_chk[2]
    diffs = []
    
    algn = True
    while True:
        
        ba = seqa[i]
        bb = seqb[j]
        
        if nchk +1 < len(checkpoints):
            next_chk = checkpoints[nchk + 1]
        else:
            next_chk = (alen, blen, alen, blen, 0)
        
        if ba == bb:
            if not algn:
                algn = True
                diffs.append((i, j))
            i+=1
            j+=1
        else:
            shift = next_chk[2] - next_chk[0]
            curr_shift = j - i
            if shift < curr_shift:
                i += 1
            elif shift > curr_shift:
                j += 1
            else:
                r = random.randint(0, 1)
                i += r
                j += 1-r
            
            algn = False
        
        if i > next_chk[0] and j > next_chk[2]:
            nchk += 1
        
        if i >= alen or j >= blen:
            break
    
    diffs.append((alen, blen))
    return diffs

def fill_aligned(seqa, seqb, diffs, fill = '-'):
    afill = []
    bfill = []
    i=j=0
    for da, db in diffs:
        afill.append(seqa[i:da])
        bfill.append(seqb[j:db])
        if da == len(seqa) or db == len(seqb):
            break
        diff = da - db
        if diff < 0:
            afill.append(fill*abs(diff))
        elif diff > 0:
            bfill.append(fill*abs(diff))
        else:
            print("weird")
        i, j = da, db
    
    return "".join(afill), "".join(bfill)

def print_aligned(seqa, seqb, diffs):
    print(f"sequence a: {seqa}")
    print(f"sequence b: {seqb}")
    print(f"{len(diffs)} differences")
    
    i = j = 0
    for a, b in diffs:
        print(seqa[i:a])
        print(seqb[j:b])
        i = a
        j = b
    print(seqa[i:])
    print(seqb[j:])
    
    

