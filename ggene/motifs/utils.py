

from typing import List
import random

from ggene.seqs.find import find_subsequence, find_subsequences
from ggene.draw import Color, RESET, CC, CB

def print_complements(seq, vocab, color_complement = True):
    
    if not color_complement:
        print(seq)
        return
    
    noncomp = vocab[::2]
    comp = vocab[1::2]
    
    iscomp = -1
    colored_seq = []
    for s in seq:
        scomp = s in comp
        sym = s
        
        if scomp == iscomp:
            pass
        else:
            colored_seq.append(RESET)
            colored_seq.append(str(CC) if iscomp else str(CB))
            iscomp = scomp
        
        colored_seq.append(sym)
    
    colored_seq.pop(0)
    colored_seq.append(RESET)
    printseq = "".join(colored_seq)
    print(printseq)

def make_start_ends(feature_name, feature_positions, feature_length, starts = {}, ends = {}):
    
    for p in feature_positions:
        if not p in starts:
            starts[p] = []
        starts[p].append(feature_name)
        
        end_pos = p + feature_length
        if not end_pos in ends:
            ends[end_pos] = []
        ends[end_pos].append(feature_name)
    
    return starts, ends

def make_key(features, colors):
    
    parts = ["key:"]
    for f in features:
        cf = colors.get(f)
        cstr = f"\x1b[38;5;{cf}m"
        fstr = "".join([cstr, f, RESET])
        parts.append(fstr)
    
    return " ".join(parts)
    
def highlight_features(seq, features, feature_starts, feature_ends, colors = {}):
    """
    feature_starts: Dict:feature -> List[feature_start_pos]
    feature_ends: Dict:feature -> List[feature_end_pos]
    """
    # Build the colored string
    
    baseline_color = 240
    bc = f"\x1b[38;5;{baseline_color}m"
    
    result = []
    active_seqs = []  # Currently active sequences
    current_color = 0
    last_ansi = bc
    
    for s in features:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    
    key = make_key(features, colors)
    
    result.append(bc)  # Start with baseline color

    for i in range(len(seq)):
        if i in feature_ends:
            for s in feature_ends[i]:
                if s in active_seqs:
                    active_seqs.remove(s)

        if i in feature_starts:
            for s in feature_starts[i]:
                if s not in active_seqs:
                    active_seqs.append(s)

        if not active_seqs:
            current_color = baseline_color
        elif len(active_seqs) == 1:
            current_color = colors[active_seqs[0]]
        elif len(active_seqs) == 2:
            color_sum = sum(colors[s] for s in active_seqs)
            current_color = 20 + (color_sum % 211)
        else:
            current_color = 255

        ansi = f"\x1b[38;5;{current_color}m"

        if ansi != last_ansi:
            result.append(ansi)
            last_ansi = ansi

        result.append(seq[i])

    result.append(RESET)

    print("".join(result))
    print(key)

def highlight_sequence(seq, subseq, colors = {}):
    
    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there
    
    if not subseq in colors:
        colors[subseq] = random.randint(20, 230)

    start_pos = find_subsequence(seq, subseq)
    starts, ends = make_start_ends(subseq, start_pos, len(subseq), starts=starts, ends = ends)

    highlight_features(seq, [subseq], starts, ends)
    pass

def highlight_sequences(seq:str, subseqs:List[str], do_rc = False, min_len = 5, colors={}):

    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there
    
    for s in subseqs:
        if not s in colors:
            colors[s] = random.randint(20, 230)

    start_pos = find_subsequences(seq, subseqs, do_rc = do_rc)
    for s in subseqs:
        if len(s) < min_len:
            continue

        starts, ends = make_start_ends(s, start_pos[s], len(s), starts=starts, ends = ends)

    highlight_features(seq, subseqs, starts, ends)

    
def highlight_sequences_in_frame(seq:str, subseqs:List[str], frame_start, min_len = 5):
    
    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there

    # Baseline color - a visible gray (color 240 is a nice medium gray)
    baseline_color = 240
    bc = f"\x1b[38;5;{baseline_color}m"
    
    colors = {}
    for s in subseqs:
        colors[s] = random.randint(20, 230)

    # Find all occurrences of each subsequence
    for s in subseqs:
        if len(s) < min_len:
            continue
        
        seq_pos = find_subsequence(seq, s, frame_start = frame_start)
        for p in seq_pos:
            if not p in starts:
                starts[p] = []
            starts[p].append(s)
        
            end_pos = p + len(s)
            if not end_pos in ends:
                ends[end_pos] = []
            ends[end_pos].append(s)

    # Build the colored string
    result = []
    active_seqs = []  # Currently active sequences
    current_color = 0
    last_ansi = bc

    result.append(bc)  # Start with baseline color

    for i in range(len(seq)):
        if i in ends:
            for s in ends[i]:
                if s in active_seqs:
                    active_seqs.remove(s)

        if i in starts:
            for s in starts[i]:
                if s not in active_seqs:
                    active_seqs.append(s)

        if not active_seqs:
            current_color = baseline_color
        elif len(active_seqs) == 1:
            current_color = colors[active_seqs[0]]
        elif len(active_seqs) == 2:
            color_sum = sum(colors[s] for s in active_seqs)
            current_color = 20 + (color_sum % 211)
        else:
            current_color = 255

        ansi = f"\x1b[38;5;{current_color}m"

        if ansi != last_ansi:
            result.append(ansi)
            last_ansi = ansi

        result.append(seq[i])

    result.append(RESET)

    print("".join(result))
