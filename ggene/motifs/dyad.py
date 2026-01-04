

import numpy as np
import string

from typing import Optional, Tuple, List

from ggene.seqs.bio import reverse_complement, COMPLEMENT_MAP, ALIASES, ORDER
from ggene.seqs.find import compare_sequences

vocab = ORDER[:4]

class Dyad:

    def __init__(self,
                 stem_length: int,
                 stem_start: int,
                 loop_length: int,
                 reverse_stem_start: int,
                 ref_pos = 0,
                 sequence: Optional[str] = None):
        self.stem_length = stem_length
        self._stem_start = stem_start
        self.loop_length = loop_length
        self._reverse_stem_start = reverse_stem_start
        self.ref_pos = ref_pos
        self.sequence = sequence
        self._dyad_seq = None  # Cache for extracted sequence

        # Validate consistency
        if self._reverse_stem_start != self._stem_start + self.stem_length + self.loop_length:
            raise ValueError(f"Inconsistent dyad positions: reverse_stem should be at "
                           f"{self.stem_start + self.stem_length + self.loop_length}, "
                           f"but is at {self.reverse_stem_start}")

    @classmethod
    def from_tuple(cls, dyad_tuple: Tuple[int, int, int, int], sequence: Optional[str] = None):
        """Create a Dyad from the old tuple format (stem_length, stem_start, loop_length, reverse_stem_start)."""
        return cls(*dyad_tuple, sequence=sequence)

    @classmethod
    def from_parameters(cls, seq_len, stem_length = None, stem_start = None, loop_start=None, loop_length = None, reverse_stem_start = None, reverse_stem_end = None, center_pos = None):
        
        return None

    def to_tuple(self) -> Tuple[int, int, int, int]:
        """Convert back to tuple format for compatibility."""
        ss = self.stem_start
        rss = self.reverse_stem_start
        return (self.stem_length, ss, self.loop_length, rss)

    def copy(self):
        return Dyad(*self.to_tuple(), ref_pos = self.ref_pos, sequence = self.sequence)

    @property
    def stem_start(self):
        return self._stem_start + self.ref_pos
    
    @property
    def stem_end(self):
        return self.loop_start
    
    @property
    def reverse_stem_start(self):
        return self._reverse_stem_start + self.ref_pos
    
    @property
    def reverse_stem_end(self):
        return self.end_position
    
    @property
    def end_position(self) -> int:
        """End position of the entire dyad in the sequence."""
        return self._reverse_stem_start + self.stem_length + self.ref_pos

    @property
    def total_length(self) -> int:
        """Total length of the dyad structure."""
        return 2 * self.stem_length + self.loop_length

    @property
    def loop_start(self) -> int:
        """Start position of the loop region."""
        return self._stem_start + self.stem_length + self.ref_pos

    @property
    def loop_end(self):
        return self.reverse_stem_start
    
    @property
    def center_position(self) -> float:
        """Center position of the dyad."""
        return self._stem_start + self.total_length / 2 + self.ref_pos

    @property
    def RC(self):
        tup = self.to_tuple()
        newtup = (-tup[0], tup[3], -tup[1], tup[2])
        return Dyad.from_tuple(newtup, sequence = self.sequence)

    def extract_sequence(self, seq: Optional[str] = None) -> str:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        # Cache the result if using stored sequence
        if seq is self.sequence and self._dyad_seq is not None:
            return self._dyad_seq

        # dyad_seq = seq[self.stem_start:self.end_position]
        dyad_seq = self.extract_sequence_cls(seq, *self.to_tuple())

        if seq is self.sequence:
            self._dyad_seq = dyad_seq

        return dyad_seq

    @classmethod
    def extract_sequence_cls(cls, seq, stem_length, stem_start, loop_length, reverse_stem_start):
        return seq[stem_start:stem_start+2*stem_length+loop_length]

    def extract_stems(self, seq: Optional[str] = None) -> Tuple[str, str]:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        first_stem = seq[self.stem_start:self.loop_start]
        second_stem = seq[self.reverse_stem_start:self.end_position]
        return first_stem, second_stem

    def extract_loop(self, seq: Optional[str] = None) -> str:
        """Extract the loop sequence."""
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        return seq[self.loop_start:self.reverse_stem_start]

    def test_dyad(self, seq: Optional[str] = None,
                  err_tol: Optional[int] = None,
                  wobble_tol = None) -> int:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")
        
        if self.stem_start < 0 or self.reverse_stem_end > len(self.sequence):
            return False
        
        first_stem, second_stem = self.extract_stems(seq)

        # Reverse complement the first stem
        first_rc = reverse_complement(first_stem)

        return compare_sequences(first_rc, second_stem, err_tol, wobble_tol =wobble_tol)

    def refer_to(self, position):
        self.ref_pos = position

    def shift(self, offset):
        self.ref_pos += offset

    def dereference(self):
        ref = self.ref_pos
        self._stem_start -= ref
        self._reverse_stem_start -= ref
        self.ref_pos = 0

    def check_subdyad(self, other: 'Dyad') -> bool:
        # Ensure self is the larger dyad
        if other.stem_length > self.stem_length:
            return other.check_subdyad(self)
        elif other.stem_length == self.stem_length:
            return self == other
        
        return self.check_subdyad_cls(self, other, seq = self.sequence)

    @classmethod
    def check_subdyad_cls(cls, da:'Dyad', db:'Dyad', seq=None):
        
        if not seq:
            seq = da.sequence
            
        # Check containment
        return (db.stem_start >= da._stem_start and
                db.stem_start + db.stem_length <= da.stem_start + da.stem_length and
                db.reverse_stem_start >= da.reverse_stem_start and
                db.reverse_stem_start + db.stem_length <= da.reverse_stem_start + da.stem_length)
        
    def check_mutual_subdyad(self, other: 'Dyad') -> Tuple[bool, Optional['Dyad']]:
        # Ensure self is the larger dyad
        if other.stem_length > self.stem_length:
            return other.check_mutual_subdyad(self)
        
        return self.check_mutual_subdyad_cls(self, other, self.sequence)

    @classmethod
    def check_mutual_subdyad_cls(cls, da:'Dyad', db:'Dyad', seq=None):
        if not seq:
            seq = da.sequence
        
        if da.center_position == db.center_position:
            # Create the minimal common superdyad
            stem_start = min(da.stem_start, db.stem_start)
            loop_length = min(da.loop_length, db.loop_length)
            reverse_stem_start = min(da.reverse_stem_start, db.reverse_stem_start)
            stem_length = reverse_stem_start - loop_length - stem_start

            superdyad = cls(stem_length, stem_start, loop_length, reverse_stem_start,
                           sequence=seq)
            return True, superdyad

        return False, None

    def expand(self, seq: Optional[str] = None,
               complement_map: Optional[dict] = None) -> Tuple[Optional['Dyad'], int]:
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")

        # Check bounds
        if self._stem_start == 0 or self.end_position >= len(seq):
            return None, 1

        expanded = Dyad(
            self.stem_length + 1,
            self._stem_start - 1,
            self.loop_length,
            self._reverse_stem_start,
            ref_pos = self.ref_pos,
            sequence=seq
        )

        err = expanded.test_dyad(seq, complement_map)
        return expanded, err

    def contract(self, seq: Optional[str] = None) -> Tuple[Optional['Dyad'], int]:
        if seq is None:
            seq = self.sequence

        # Check minimum size and loop length
        if self.stem_length <= 1 or self.loop_length <= 2:
            return None, 1

        contracted = Dyad(
            self.stem_length + 1,
            self._stem_start,
            self.loop_length - 2,
            self._reverse_stem_start - 1,
            ref_pos = self.ref_pos,
            sequence=seq
        )

        err = contracted.test_dyad(seq)
        return contracted, err

    def maximize(self, seq=None, nmax = None)->'Dyad':
        if seq is None:
            seq = self.sequence
        if seq is None:
            raise ValueError("No sequence provided or stored")
        dinit = self
        
        if not nmax:
            nmax_ctr = self.loop_length // 2
            nmax_exp = min(self.stem_start, len(seq) - self.reverse_stem_start)
        else:
            nmax_ctr = nmax
            nmax_exp = nmax
        
        n = 0
        while n < nmax_ctr:
            dctr, err = dinit.contract()
            if err == 0:
                dinit = dctr
            else:
                break
            n+=1
        
        n=0
        while n < nmax_exp:
            dexp, err = dinit.expand()
            if err == 0:
                dinit = dexp
            else:
                break
            n+=1
            
        return dinit

    def shift_copy(self, offset: int) -> 'Dyad':
        return Dyad(
            self.stem_length,
            self._stem_start + offset,
            self.loop_length,
            self._reverse_stem_start + offset,
            ref_pos = self.ref_pos,
            sequence=self.sequence
        )

    # @staticmethod
    # def _compare_sequences(s1: str, s2: str, err_tol: Optional[int] = None, allow_wobble = True) -> int:
    #     """Helper to compare two sequences and count mismatches."""
    #     if len(s1) != len(s2):
    #         return len(s1)

    #     if err_tol is None:
    #         err_tol = len(s1)

    #     nerr = 0
    #     for a, b in zip(s1, s2):
    #         if a != b:
    #             nerr += 1
    #         if nerr > err_tol:
    #             return nerr
    #     return nerr

    def print(self, total_len = -1, colors = {}):
        self.print_dyad(self.sequence, *self.to_tuple(), total_len = total_len, colors = colors)
    
    @classmethod
    def print_dyad(cls, seq, stem_length, stem_start, loop_length, reverse_stem_start, total_len = -1, colors = {}):
        cr ="\x1b[0m"
        cs = colors.get("seq","\x1b[38;5;8m") # bright black = gray
        cd = colors.get("stem","\x1b[38;5;3m") # yellow
        cl = colors.get("loop","\x1b[38;5;6m") # cyan
        
        d = cls(stem_length, stem_start, loop_length, reverse_stem_start, sequence = seq)
        st, rst = d.extract_stems(seq)
        lp = d.extract_loop()
        full = d.extract_sequence(seq)
        raw = seq[d.stem_start:d.stem_start+2*d.stem_length+d.loop_length]
        
        head_start = 0
        tail_end = len(seq)
        if total_len:
            crop = total_len - 2*stem_length - loop_length
            head_start = stem_start - crop//2
            tail_end = d.end_position + (crop+1)//2
        
        if head_start < 0:
            tail_end += abs(head_start)
            head_start = 0
        if tail_end > len(seq):
            head_start -= (tail_end - len(seq))
            tail_end = len(seq)
            
        head = seq[head_start:d.stem_start]
        tail = seq[d.end_position:tail_end]
        
        parts = []
        for p, c in zip([head,st, lp, rst, tail], [cs, cd, cl, cd, cs]):
            parts.append(f"{c}{p}{cr}")
        
        print("".join(parts))

    def __eq__(self, other: 'Dyad') -> bool:
        """Check equality based on all position attributes."""
        if not isinstance(other, Dyad):
            return False
        return (self.stem_length == other.stem_length and
                self._stem_start == other._stem_start and
                self.loop_length == other.loop_length and
                self._reverse_stem_start == other._reverse_stem_start)

    def __repr__(self) -> str:
        """String representation of the dyad."""
        seq_info = f", has_seq={self.sequence is not None}"
        return (f"Dyad(stem={self.stem_length}, start={self._stem_start}, "
                f"loop={self.loop_length}, rev_start={self._reverse_stem_start}{seq_info})")

    def __str__(self) -> str:
        """Human-readable string representation."""
        if self.sequence:
            stems = self.extract_stems()
            loop = self.extract_loop()
            return f"Dyad: {stems[0]}-{loop}-{stems[1]} @ pos {self._stem_start}"
        return f"Dyad: {self.stem_length}bp stems, {self.loop_length}bp loop @ pos {self._stem_start}"

    def __len__(self):
        return 2*self.stem_length + self.loop_length
    
    def __hash__(self):
        return hash(self.to_tuple())
    
    @classmethod
    def build_superdyads(cls, dyads:List['Dyad'])->List['Dyad']:
                
        dyads = cls.sort_dyads(dyads, sortby = "center_position")
        seq = dyads[0].sequence
        seq_len = len(seq)
        
        sups = []
        ncomps = 0
        
        n = 0
        for cc in range(int(2*dyads[0].center_position), int(2*dyads[-1].center_position)+1):
            c = cc/2
            
            minss = seq_len+1
            minll = seq_len+1
            minrss = seq_len+1
            nn=0
            for d in dyads[n:]:
                
                if d.center_position>c:
                    break
                
                ncomps += 1
                nn+=1
                
                if d.stem_start < minss:
                    minss = d.stem_start
                
                if d.loop_length < minll:
                    minll = d.loop_length
                
                if d.reverse_stem_start < minrss:
                    minrss = d.reverse_stem_start
                
            if nn == 0:
                continue
            elif nn==1:
                sups.append(dyads[n])
                n += 1
                continue
            n += nn
            
            sl = minrss - minll - minss
            sup = Dyad(sl, minss, minll, minrss, sequence = seq)
            if sup.end_position >= seq_len:
                continue
            
            err= sup.test_dyad()
            if err == 0:
                supmax = sup.maximize()
                sups.append(supmax)
            else:
                # print(f"dyad failed at {c}: {sup.to_tuple()}, sequence length {seq_len}")
                pass
        # print(f"compared {ncomps} dyads")
        
        return sups
    
    @classmethod
    def sort_dyads(cls, dyads:List['Dyad'], sortby = "stem_start", reverse = False):
        
        if len(dyads) < 2:
            return dyads
        
        if not hasattr(dyads[0], sortby):
            raise Exception
        
        m = -1 if reverse else 1
        
        return list(sorted(dyads, key = lambda d:m*getattr(d, sortby)))
    
    @classmethod
    def remove_subdyads(cls, dyads:List['Dyad']):
        
        dyads = cls.sort_dyads(dyads, sortby = "center_position")
        
        n = 0
        for c in range(dyads[0].center_position, dyads[-1].center_position):
            nn=n
            for d in dyads[n:]:
                if d.center_position>c:
                    break
                nn += 1
            
            cdyads = cls.sort_dyads(dyads[n:nn], sortby = "stem_length", reverse = True)



# Lookup table: base_i XOR complement_j == 3 for valid pairs
_BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
_COMP_XOR = 3  # A(0)^T(3)=3, C(1)^G(2)=3

def find_dyads_numpy(seq, arm_length, spacer_range=(0, 20)):
    """
    Fast dyad detection using numpy.
    
    For a dyad: seq[i:i+k] == revcomp(seq[j:j+k]) where j = i+k+spacer
    Which means: seq[i+m] complements seq[j+k-1-m] for all m in [0,k)
    """
    seq = seq.upper()
    n = len(seq)
    arr = np.array([_BASE_MAP.get(b, 4) for b in seq], dtype=np.uint8)

    min_sp, max_sp = spacer_range
    dyads = []

    for spacer in range(min_sp, max_sp + 1):
        total_len = 2 * arm_length + spacer
        if total_len > n:
            continue

        # For position i, compare arr[i:i+arm] with arr[i+arm+spacer:i+total][::-1]
        # Using XOR trick: complementary bases XOR to 3

        for i in range(n - total_len + 1):
            left = arr[i:i + arm_length]
            right = arr[i + arm_length + spacer:i + total_len][::-1]

            # Check all positions complement (XOR to 3, ignoring N's)
            mask = (left < 4) & (right < 4)  # valid bases
            if mask.all() and ((left ^ right) == _COMP_XOR).all():
                dyads.append((i, i + total_len, spacer))

    return dyads


def find_all_dyads_zalgo(seq, min_arm=4):
    """
    Find all reverse-complement palindromes using Z-algorithm.
    O(n) time complexity.
    """
    seq = seq.upper()
    rc = reverse_complement(seq)

    print(seq)
    print(rc)

    # Concatenate: seq + '$' + reversed(rc)
    # Matches indicate palindrome centers
    combined = seq + '$' + rc[::-1]

    # Compute Z-array
    z = compute_z_array(combined)

    n = len(seq)
    dyads = []

    # Z[n+1+i] tells us match length starting at position i in reversed(rc)
    # which corresponds to match at position (n-1-i) in rc
    # A match of length L at position i means seq[i:i+L] == rc[n-1-i-L+1:n-1-i+1]

    for i in range(n + 1, len(combined)):
        match_len = z[i]
        print(match_len, min_arm)
        if match_len >= min_arm:
            seq_pos = i - n - 1
            # This indicates a palindrome centered around seq_pos
            dyads.append((seq_pos, match_len))

    return dyads

def compute_z_array(s):
    """Standard Z-algorithm."""
    n = len(s)
    z = [0] * n
    z[0] = n
    l, r = 0, 0

    for i in range(1, n):
        if i < r:
            z[i] = min(r - i, z[i - l])
        while i + z[i] < n and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] > r:
            l, r = i, i + z[i]

    return z

def reverse_complement(seq):
    comp = str.maketrans('ACGT', 'TGCA')
    return seq.translate(comp)[::-1]









def get_next_symbols(num_symbols):
    
    global vocab, COMPLEMENT_MAP
    if num_symbols %2 == 1:
        num_symbols -= 1
    
    syms = string.ascii_uppercase
    out =[]
    s = ''
    for sym in syms:
        if sym in vocab:
            continue
        out.append(sym)
        if s:
            COMPLEMENT_MAP[s] = sym
            COMPLEMENT_MAP[sym] = s
            
        else:
            s = sym
        
        if len(out) == num_symbols:
            break
    vocab = vocab + "".join(out)
    
    return out

def encode_symbols(seq):
    global aliases
    
    nsyms = 2*len(seq)
    seq = seq
    cs = "".join(reverse_complement(seq))
    sym, csym = get_next_symbols(2)
    
    ALIASES[sym] = seq
    ALIASES[seq] = sym
    ALIASES[csym] = cs
    ALIASES[cs] = csym
    return sym, seq, csym, cs

def encode_dyad(seq, dyad:Dyad):
    
    # dyad_seq = extract_dyad_seqs(seq, [dyad])[0]
    dyad_seq = dyad.extract_sequence()
    
    sym, subseq, csym, csubseq = encode_symbols(dyad_seq)
    
    newseq, shift = substitute_symbols(seq, sym, subseq, csym, csubseq)
    
    return newseq, sym

def substitute_symbols(seq, sym, subseq, csym, csubseq):
    
    len_seq = len(seq)
    len_ss = len(subseq)
    out = []
    shift = 0
    i = 0
    while i < len_seq:
        testseq = seq[i:i+len_ss]
        if testseq == subseq:
            out.append(sym)
            shift += len_ss - 1
            i += len_ss - 1
            # print(f"subbed symbol {sym} at {i}")
        elif testseq == csubseq:
            out.append(csym)
            shift += len_ss - 1
            i += len_ss - 1
            # print(f"subbed comp symbol {csym} at {i}")
        else:
            out.append(seq[i])
        i+=1
            
    return "".join(out), shift

def combine_templates(reference, templates):
    ref = [r for r in reference]
    
    for i in range(len(ref)):
        alltmp = {b:0 for b in vocab}
        for t in templates:
            alltmp[t[i]] += 1
        ref[i] = alltmp
    
    for d in ref:
        for b in vocab:
            if d[b] == 0:
                del d[b]

    return ref
    
def make_crosscorr(reference, templates):
    cclist = []
    for t in templates:
        tccs = []
        for i, (a, b) in enumerate(zip(reference, t)):
            
            if a!=b:
                tccs.append(1)
            else:
                tccs.append(0)
            
        cclist.append(tccs)
    cc = np.array(cclist)
    idk1 = cc @ cc.T
    idk2 = np.linalg.pinv(cc)
    return idk1, idk2

def frequency_rank_dyads(seq, dyads:List[Dyad]):
    
    max_n = 0
    max_seq = ""
    
    freqs = {}
    for d in dyads:
        # ds = d.extract_sequence(seq)
        ds, cds = d.extract_stems()
        n = freqs.get(ds, 0) + 1
        freqs[ds] = n
        if n > max_n:
            max_n = n
            max_seq = ds
    
    return freqs, max_seq

def search_dyad(seq, stem_len, max_loop, min_loop, err_tol = 0):
    """
    dumb algorithm to find palindromes/dyads given their length
    a dyad should be unique given a sequence.. right?
    """
    
    dyad_start = dyad_rc_start = loop_len = -1
    n_loop_iters = max_loop - min_loop
    n_dyad_iters = len(seq) - 2*stem_len - min_loop
    
    done = False
    allres = []
    
    for nd in range(n_dyad_iters):
        
        test_seq = seq[nd:nd + stem_len]
        test_seq_rc = "".join(reverse_complement(test_seq))
        
        for nl in range(min_loop, max_loop+1):
        
            comp_seq = seq[nd + stem_len + nl:nd + 2*stem_len + nl]
            err = compare_sequences(test_seq_rc, comp_seq, err_tol = 0)
            
            # if comp_seq == test_seq_rc:
            if err < 1:
                done = True
                
            if done:
                dyad_start = nd
                dyad_rc_start = nd + stem_len + nl
                loop_len = dyad_rc_start - dyad_start - stem_len
                break
        
        if done:
            break
    
    return stem_len, dyad_start, loop_len, dyad_rc_start

def search_dyad_inout(seq, stem_len, max_loop, min_loop):
    """
    dumb algorithm to find palindromes/dyads given their length
    a dyad should be unique given a sequence.. right?
    """
    
    dyad_start = dyad_rc_start = loop_len = -1
    n_loop_iters = max_loop - min_loop
    n_dyad_iters = len(seq) - 2*stem_len - min_loop
    
    done = False
    allres = []
    
    for nl in range(min_loop, max_loop+1):
        
        for nd in range(n_dyad_iters):
            test_seq = seq[nd:nd + stem_len]
            test_seq_rc = "".join(reverse_complement(test_seq))
            comp_seq = seq[nd + stem_len + nl:nd + 2*stem_len + nl]
            
            if comp_seq == test_seq_rc:
                dyad_start = nd
                dyad_rc_start = nd + stem_len + nl
                loop_len = dyad_rc_start - dyad_start - stem_len
                allres.append((dyad_start, loop_len, dyad_rc_start))
                done = True
        
        if done:
            break
    
    return allres
    
def seek_dyad(seq, max_stem, min_stem, max_loop, min_loop):
    """
    really dumb algorithm to find palindromes/dyads
    """
    allres= {}
    
    for dl in range(min_stem, max_stem+1):
        res = search_dyad(seq, dl, max_loop, min_loop)
        if res[0] > -1:
            dyad_len, dyad_start, loop_len, dyad_rc_start = res
            k = (dyad_start, dyad_rc_start)
            if k in allres:
                _dl, _loop_len = allres[k]
                if dl > _dl:
                    allres[k] = (dl, loop_len)
            else:
                allres[k] = (dl, loop_len)
    
    return allres
    

def seek_dyad_rev(seq, max_stem, min_stem, max_loop, min_loop):
    """
    start big and early exit?
    """
    for dl in range(min_stem, max_stem+1, -1):
        res = search_dyad(seq, dl, max_loop, min_loop)
        if res[0] > -1:
            return dl, *res
    return -1, -1, -1, -1


def find_all_dyads(seq, min_stem, max_stem = -1, max_loop = -1, min_loop = -1, err_tol = None, allow_wobble = True)->List[Dyad]:
    
    len_seq = len(seq)
    if max_stem < 0:
        max_stem = len_seq // 2
    
    if max_loop < 0:
        max_loop = len_seq - 2*min_stem
    if min_loop < 0:
        min_loop = 3
    
    if err_tol is None:
        err_tol = 0
    
    ncomps = 0
    
    allres = []
    
    for sl in range(min_stem, max_stem):
        
        for ss in range(0, len_seq - 2*sl):
            
            test_seq = seq[ss:ss+sl]
            test_seq_rc = reverse_complement(test_seq)
            
            for rss in range(ss + sl + min_loop, ss + sl + max_loop + 1):
                
                if allow_wobble:
                    wobble_tol = len(test_seq) // 5
                else:
                    wobble_tol = None
                
                comp_seq = seq[rss:rss+sl]
                err = compare_sequences(test_seq_rc, comp_seq, err_tol = err_tol, wobble_tol = wobble_tol)
                if err <= err_tol:
                    loop_len = rss - sl - ss
                    nd = Dyad(sl, ss, loop_len, rss, sequence = seq)
                    
                    if nd.test_dyad(err_tol = err_tol, wobble_tol = wobble_tol):
                        allres.append(nd)
                ncomps += 1
    
    # print(f"tested {ncomps} sequences")
    
    return allres

def find_all_dyads_window(seq, windowsz, min_stem, max_stem = -1, max_loop = -1, min_loop = -1, merge = False, err_tol = None):
    
    nchunks = len(seq) // windowsz
    
    allds = []
    
    for n in range(nchunks):
        start = n*windowsz
        subseq = seq[start:start+windowsz]
        subds = find_all_dyads(subseq, min_stem, max_stem=max_stem, max_loop=max_loop, min_loop = min_loop, err_tol = err_tol)
        
        subds = [d.shift(start) for d in subds]
        
        allds.extend(subds)
        
        if merge:
            allds = Dyad.build_superdyads(allds)
        
    return allds

def find_all_dyads_expand(seq, min_stem, min_loop = 3, max_loop = 8):
    
    len_seq = len(seq)
    max_stem = min_stem + 1
    
    if max_loop < 0:
        max_loop = len_seq - 2*min_stem
    if min_loop < 0:
        min_loop = 3
    
    min_semiloop = (min_loop+1)//2
    max_semiloop = max_loop // 2
    
    ncomps = 0
    
    allres = []
    
    for dc in range(min_stem, len_seq - min_loop):
        _max_stem = min(max_stem, dc - max_semiloop)
        
        for sl in range(min_stem, _max_stem):
            
            for sll in range(min_semiloop, max_semiloop):
                # d = (dl, dc - dsl - dl, 2*dsl, dc + dsl)
                d = Dyad()
                
                # test_seq = extract_dyad_seqs()
        # for dl in range(min_stem, max_stem):
            
        #     for ds in range(0, len_seq - 2*dl):
                
        #         test_seq = seq[ds:ds+dl]
        #         test_seq_rc = "".join(reverse_complement(test_seq))
                
        #         for dds in range(ds + dl + min_loop, ds + dl + max_loop + 1):
                    
        #             comp_seq = seq[dds:dds+dl]
        #             err = compare_seqs(test_seq_rc, comp_seq)
        #             if err < 1:
        #                 loop_len = dds - dl - ds
                        
        #                 d_new = (dl, ds, loop_len, dds)
                        
        #                 d_new = maximize_dyad(seq, d_new) # produces maximal dyad given a center position
                        
        #                 allres.append(d_new)
        #             ncomps += 1
    
    # print(f"tested {ncomps} sequences")
    
    return allres
    
    pass

def maximize_dyads(seq, dyads:List[Dyad], nmax = None):
    
    dyads = Dyad.sort_dyads(dyads, sortby = "stem_length")
    out = []
    for d in dyads:
        dmax = d.maximize(nmax = nmax)
        out.append(dmax)
    
    out = Dyad.sort_dyads(out, sortby = "stem_start")
    return out

def try_substitute(seq, dyad, search_window = 20):
    
    dyad_seq = seq[dyad[1]:dyad[0]+dyad[1]]
    newseq, newsym = encode_dyad(seq, dyad)
    
    subseq = newseq[dyad[1] - search_window: dyad[1] + search_window]
    
    ds = find_all_dyads(subseq, 3, 4, 8, 4)
    
    return newseq, ds, newsym


def find_all_dyads_build(seq, min_stem, max_loop = -1, min_loop = -1, err_tol = None):
    
    max_stem = min_stem + 1
    
    dyads = find_all_dyads(seq, min_stem, max_stem = max_stem, max_loop = max_loop, min_loop = min_loop, err_tol = err_tol)
    
    supers = Dyad.build_superdyads(dyads)
    
    return dyads, supers

def find_all_dyads_chunk(seq, min_stem, chunksz):
    
    seq_len = len(seq)
    nchunks = seq_len // chunksz
    all_supers = []
    
    for nc in range(2*nchunks - 1):
        start = nc * chunksz // 2
        subseq = seq[start:start+chunksz]
        
        ds, ss = find_all_dyads_build(subseq, min_stem)
        
        # print(f"chunk {nc}, {len(ds)} naive, {len(ss)} superdyads")
        
        all_supers.extend(ss)
        all_supers = Dyad.build_superdyads(all_supers)
        # print(f"chunk {nc}, {len(ss)} combined superdyads")
    
    return all_supers

def remove_subdyads(dyads:List[Dyad]):
    
    dyads = list(sorted(dyads, key = lambda a:-a[0]))
    out = [dyads.pop(0)]
    
    while len(dyads) > 0:
        slay = True
        db = dyads.pop(0)
        for da in out:
            if da.check_subdyad(db):
                slay = False
                break
            else:
                pass
        if slay:
            out.append(db)
    
    return out
