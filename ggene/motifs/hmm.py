
from typing import Any, Optional, List, Dict
from dataclasses import dataclass

import pyhmmer

from ggene import DATA_DIR
from ggene import draw
from .motif import BaseMotif

import h5py

hmm_path = DATA_DIR / "dfam" / "human_repeats.hmm"

AB = pyhmmer.easel.Alphabet.dna()

# class HMMLibrary:
    
#     _families = ['Alu']
    
#     def __init__(self):
#         pass

# @dataclass
# class SearchResult:
#     target:str
#     query:str
#     score:float
    
#     _best_domain:Any=None
#     _alignment:Any = None
    
#     @classmethod
#     def from_hit(cls, hit:pyhmmer.plan7.Hit):
        
#         bd = hit.best_domain
#         algn = bd.alignment
        
#         return SearchResult(
#             algn.target_sequence,
#             algn.hmm_sequence,
#             bd.score,
            
#             _best_domain = bd,
#             _alignment=  algn
            
#         )
        
#         pass
    
class HMMMotif(BaseMotif):
    
    def __init__(self, name, hmm, scoring_function, allow_rc = True, motif_class = ""):
        super().__init__(name)
        
        self.hmm = hmm
        self.score_func = scoring_function
        self.allow_rc = allow_rc
        self.motif_class = motif_class
    
    def __call__(self, seq):
        return 0
    
    def score(self, seq, return_positions=False):
        return None
    
    def count_instances(self, seq):
        return 0
    
    def find_instances(self, seq, threshold=None):
        
        insts = []
        
        res = search_hmm(self.hmm, [seq])
        
        for hit in res:
            
            if threshold and hit.score < threshold:
                continue
            algn = hit.best_domain.alignment
            insts.append((algn.target_from, algn.target_to, hit.score))
        
        return insts

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
        
        new_motif = HMMMotif(name, hmm, None)
        motifs.append(new_motif)
    
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
    pass

def load_line1s():
    return load_hmms(name_filt = "L1")
    
    pass
    
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

def display_group(te_group, start=None, end=None, display_length = 256, print_info = True, show_ruler = True):
    

    
    fwds = [f for f in te_group if f.strand=="+"]
    revs = [f for f in te_group if f.strand=="-"]
    
    if not start:
        start = min([f.start-f.attributes.get("start_hmm",0) for f in fwds] + [f.start - f.attributes.get("end_hmm", 0) for f in revs])
    if not end:
        end = max([f.end+f.attributes.get("end_hmm",0) for f in fwds] + [f.end + f.attributes.get("start_hmm", 0) for f in revs])
    
    pos_str = f"{start}-{end} ({end-start}nt)"
    if len(te_group) < 1:
        print(f"Nothing found at {pos_str} ... ")
        return
    
    fwd_lines = build_display_lines(fwds, start, end, display_length=display_length)
    rev_lines = build_display_lines(revs, start, end, display_length=display_length, rev = True)
    
    print(f"Repeat Cluster at {te_group[0].chrom}:{pos_str}")
    for i in range(4):
        print()
    
    for row in rev_lines:
        print(row)
    for row in fwd_lines:
        print(row)
    
    if show_ruler:
        ruler, _ = draw.make_ruler(0, end-start, num_cols = display_length, num_labels = 5, ticks = 0, minor_ticks = 2*display_length, formatter = "genomic")
        print(ruler)
    print()
    
    if print_info:
        print_repeats(te_group, start_offset = start)
        print()
    
def build_display_lines(repeat_insts, start, end, display_length = 256, rev = False):
    """
    Build display lines for repeat instances with dynamic row allocation.

    Each repeat shows:
    - Full HMM consensus extent in gray dashes
    - Actual matched region (f.start to f.end) in color with arrowheads
    """

    scale = display_length / (end - start)

    # Gray color for HMM consensus extent
    gray_color = "\x1b[38;5;240m"

    # Track which columns are occupied in each row (for conflict detection)
    # row_occupancy[row_idx] = list of (start_col, end_col) tuples
    row_occupancy = []

    # Store drawing instructions for each repeat: (row_idx, repeat_obj)
    row_assignments = []

    # Assign each repeat to the topmost available row
    for f in repeat_insts:
        # Calculate full HMM extent in display coordinates
        hmm_start = f.attributes.get("start_hmm", 0)
        hmm_end = f.attributes.get("end_hmm", 0)
        motif_len = f.attributes.get("motif_len", hmm_end)

        # Full HMM span in genomic coordinates
        hmm_full_start = f.start - hmm_start
        hmm_full_end = f.end + (motif_len - hmm_end)

        # Convert to display coordinates
        hmm_full_start_d = max(0, int(scale * (hmm_full_start - start)))
        hmm_full_end_d = min(display_length - 1, int(scale * (hmm_full_end - start)))

        # Find the first row where this repeat fits
        assigned_row = None
        for row_idx, occupancy_list in enumerate(row_occupancy):
            # Check if this row has space
            has_conflict = False
            for occupied_start, occupied_end in occupancy_list:
                # Check for overlap (with 1 char buffer for readability)
                if not (hmm_full_end_d < occupied_start - 1 or hmm_full_start_d > occupied_end + 1):
                    has_conflict = True
                    break

            if not has_conflict:
                assigned_row = row_idx
                break

        # If no existing row fits, create a new one
        if assigned_row is None:
            assigned_row = len(row_occupancy)
            row_occupancy.append([])

        # Mark this region as occupied
        row_occupancy[assigned_row].append((hmm_full_start_d, hmm_full_end_d))
        row_assignments.append((assigned_row, f))

    # Create rows based on how many we need
    num_rows = len(row_occupancy)
    rows = [[" "] * display_length for _ in range(num_rows)]

    # Draw each repeat
    for row_idx, f in row_assignments:
        # Get repeat attributes
        hmm_start = f.attributes.get("hmm-start", 0)
        hmm_end = f.attributes.get("hmm-end", 0)
        motif_len = f.attributes.get("motif_len", hmm_end)

        # Calculate positions
        # Full HMM extent
        hmm_full_start = f.start - hmm_start
        hmm_full_end = f.end + (motif_len - hmm_end)
        hmm_full_start_d = max(0, int(scale * (hmm_full_start - start)))
        hmm_full_end_d = min(display_length - 1, int(scale * (hmm_full_end - start)))

        # Actual match extent
        match_start_d = max(0, int(scale * (f.start - start)))
        match_end_d = min(display_length - 1, int(scale * (f.end - start)))

        # Get color for this repeat
        col = get_display_color(f)

        # Draw full HMM extent in gray dashes
        for i in range(hmm_full_start_d, hmm_full_end_d + 1):
            if rows[row_idx][i] == " ":
                rows[row_idx][i] = gray_color + "·"
        
        head = "<>"
        tail = "┤├"
        body = "─"
        
        # Draw actual match region in color
        # Start arrowhead
        if match_start_d < display_length:
            rows[row_idx][match_start_d] = col + (head[0] if rev else tail[1])

        # End arrowhead
        if match_end_d < display_length:
            rows[row_idx][match_end_d] = (tail[0] if rev else head[1]) + "\x1b[0m"

        # Fill in the middle with dashes
        for i in range(match_start_d + 1, match_end_d):
            if i < display_length:
                rows[row_idx][i] = body

        # Overlay repeat name in the middle
        name_len = len(f.name)
        name_start = match_start_d + int((match_end_d - match_start_d) / 2 - name_len / 2)

        for i, char in enumerate(f.name):
            name_pos = name_start + i
            if match_start_d <= name_pos <= match_end_d and name_pos < display_length:
                rows[row_idx][name_pos] = char

    # Add color reset at end of each row
    result_rows = []
    for row in rows:
        row.append("\x1b[0m")
        result_rows.append("".join(row))

    return result_rows

def get_display_color(f):
    
    if "Alu" in f.name:
        if "AluY" in f.name:
            return "\x1b[38;5;142m"
        elif "AluJ" in f.name:
            return "\x1b[38;5;136m"
        elif "AluS" in f.name:
            return "\x1b[38;5;130m"
        else:
            return "\x1b[38;5;124m"
    elif "L1" in f.name:
        return "\x1b[38;5;121m"
    elif "SVA" in f.name:
        if "SVA_A" in f.name:
            return "\x1b[38;5;93m"
        elif "SVA_B" in f.name:
            return "\x1b[38;5;92m"
        elif "SVA_C" in f.name:
            return "\x1b[38;5;91m"
        elif "SVA_D" in f.name:
            return "\x1b[38;5;90m"
        else:
            return "\x1b[38;5;89m"
    else:
        return "\x1b[0m"
    