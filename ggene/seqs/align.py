
import random

import numpy as np
import os
import matplotlib.pyplot as plt

from ggene.seqs import process
from ggene.seqs import bio

from Bio import Align

class AlignmentResult:
    
    def __init__(self, alignment:Align.Alignment, aligner:Align.PairwiseAligner = None):
        
        self.length = alignment.length
        self.target = alignment.target
        self.query = alignment.query
        
        cts = alignment.counts()
        self.gaps, self.identities, self.mismatches = cts.gaps, cts.identities, cts.mismatches
        self.score = alignment.score
        self.ham = self.length - self.identities
        
        self.target_algn, self.query_algn = alignment[0], alignment[1]
        
        self.subs = alignment.substitutions
        self.t_algn, self.q_algn = alignment.aligned
        
        # if aligner:
        #     self.match_score = aligner.match_score
        #     self.mismatch_score = aligner.mismatch_score
        #     self.ogap_score = aligner.open_gap_score
        #     self.egap_score = aligner.extend_gap_score
        # else:
        #     self.match_score=self.mismatch_score=self.ogap_score=self.egap_score=0
        
        self._alignment = alignment
    
    def get_aligned_seqs(self):
        return self.target_algn, self.query_algn
        
    def print(self, chunksz = 256):
        
        print(f"score: {self.score:0.1f}")
        # print(f"score: {self.score}, hamming: {self.ham}")
        print(f"identities: {self.identities}, gaps: {self.gaps}, mismatches: {self.mismatches}")
        # print(f"match score={self.match_score}, mismatch score={self.mismatch_score}, open gap={self.ogap_score}, extend gap={self.egap_score}")
        
        tgt_algn, q_algn = self.get_aligned_seqs()
        num_chunks = len(tgt_algn)//chunksz + 1
        for n in range(num_chunks):
            print(tgt_algn[n*chunksz:(n+1)*chunksz])
            print(q_algn[n*chunksz:(n+1)*chunksz])
            print()
        
    def __len__(self):
        return self.length

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

_default_scores = {
    "match_score":1,
    "mismatch_score":-0.9,
    
    # "open_left_deletion_score":0,
    # "open_internal_deletion_score":-1,
    # "open_right_deletion_score":0,
    
    # "open_left_insertion_score":0,
    # "open_internal_insertion_score":-1,
    # "open_right_insertion_score":0,
    
    # "extend_left_deletion_score":0,
    # "extend_internal_deletion_score":-1,
    # "extend_right_deletion_score":0,
    
    # "extend_left_insertion_score":0,
    # "extend_internal_insertion_score":-1,
    # "extend_right_insertion_score":0,
    
    "left_gap_score":-1,
    "right_gap_score":-1,
    "internal_gap_score":-1,
}

def get_scores(scores):
    sc = _default_scores.copy()
    sc.update(scores)
    return sc

def score_sequences(seqa, seqb, **scores):
    scores = get_scores(scores)
    score = Align.PairwiseAligner(**scores).score(seqa, seqb)
    return score

def align_sequences(seqa, seqb, topk = 1, **scores):
    scores = get_scores(scores)
    
    aligner = Align.PairwiseAligner(**scores)
    res = aligner.align(seqa, seqb)
    
    return [AlignmentResult(res[k], aligner=aligner) for k in range(topk)]

# def print_alignment(algn:Align.Alignment):
    
#     print(algn.target)
#     print(algn.query)
    
    # algn.aligned

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
    
    

