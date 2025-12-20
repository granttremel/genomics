
from typing import List, Tuple, Dict, Any, Optional

import random
import regex


from ggene.seqs.bio import complement
from ggene.seqs.find import find_subsequence, find_subsequences
from .colors import Colors

def highlight_subsequences(subseqs, colors, delimit = ""):
    
    if isinstance(colors[0], int):
        colors = [f"\x1b[38;5;{c}m" for c in colors]
    
    lastcolor = colors[0]
    out = [colors[0]]
    for ss, c in zip(subseqs, colors):
        if c == lastcolor:
            pass
        else:
            out.append(Colors.RESET)
            if delimit:
                out.append(delimit)
            out.append(c)
            lastcolor = c
        out.append(ss)
    out.append(Colors.RESET)
    return "".join(out)

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

def make_spans(starts, ends):
    
    feats = set()
    istarts = {}
    iends = {}
    spans = {}
    
    for s, fs in starts.items():
        
        for f in fs:
            if not f in istarts:
                feats.add(f)
                istarts[f] = []
            istarts[f].append(s)
    
    for e, fs in ends.items():
        for f in fs:
            if not f in iends:
                feats.add(f)
                iends[f] = []
            iends[f].append(e)
    
    for f in feats:
        st = istarts[f]
        en = iends[f]
        spans[f] = [(s, e) for s, e in zip(sorted(st), sorted(en))]
    
    return spans

def make_key(features, colors):
    
    parts = ["key:"]
    for f in features:
        cf = colors.get(f)
        cstr = f"\x1b[38;5;{cf}m"
        fstr = "".join([cstr, f, Colors.RESET])
        parts.append(fstr)
    
    return " ".join(parts)

def highlight_matching(seqa, seqb, colors = None, do_rc = False, do_both = False, suppress = False, color_bg = False, chunksz = 256):
    
    do_fwd = (not do_rc)
    if do_both:
        do_rc = True
        do_fwd = True
    
    bg_frm = "\x1b[48;5;{}m"
    fg_frm = "\x1b[38;5;{}m"
    
    if color_bg:
        cfrm = bg_frm
    else:
        cfrm = fg_frm
    
    base_color_fg = fg_frm.format(248)
    base_color_bg = cfrm.format(234)
    base_color = Colors.Colors.RESET + base_color_fg + base_color_bg
    if not colors:
        color = cfrm.format(142)
        rcolor = Colors.MOTIF
    else:
        c, cr = colors
        color = cfrm.format(c)
        rcolor = cfrm.format(cr)
    
    color = color + Colors.BOLD
    rcolor = rcolor + Colors.BOLD
    
    seqah = [base_color]
    seqbh = [base_color]
    rcseqbh = [base_color]
    
    rseqb = reversed(seqb)
    
    n = 0
    for sa, sb, rsb in zip(seqa, seqb, rseqb):
        
        apre = ""
        bpre = ""
        rcbpre = ""
        post = ""
        
        if do_rc and sa == complement(rsb):
            apre = rcbpre = rcolor
            post = base_color
        if do_fwd and sa == sb:
            apre = bpre = color
            post = base_color
        
        seqah.append(f"{apre}{sa}{post}")
        seqbh.append(f"{bpre}{sb}{post}")
        rcseqbh.append(f"{rcbpre}{rsb}{post}")
        n+=1
    
    seqah.append(Colors.Colors.RESET)
    seqbh.append(Colors.Colors.RESET)
    rcseqbh.append(Colors.Colors.RESET)
    
    seqaout = "".join(seqah)
    seqbout = "".join(seqbh)
    rcseqbout = "".join(rcseqbh)
    
    if not suppress:
        for n in range(len(seqa)//chunksz + 1):
            print(seqaout[n*chunksz:(n+1)*chunksz])
            print(seqbout[n*chunksz:(n+1)*chunksz])
            print(rcseqbout[n*chunksz:(n+1)*chunksz])
        print(Colors.Colors.RESET)
    
    return seqaout, seqbout, rcseqbout

def highlight_correlated(full_seq, shift, colors = None, suppress = False):
    
    if not colors:
        ca = Colors.RCMOTIF
        cab= Colors.HIGHLIGHT
        cb = Colors.MOTIF
    
    seq_len = len(full_seq)
    
    seqah = []
    seqbh = []
    
    sas = min(0, shift)
    sbs = max(0, -shift)
    
    for n in range(len(full_seq)):
        
        sf = full_seq[n]
        sa = sb = ""
        
        if sas + n < seq_len:
            sa = full_seq[sas + n]        
        if sbs + n < seq_len:
            sb = full_seq[sbs + n]
        
        if sa == sb:
            pre = ca
            post = Colors.Colors.RESET
        else:
            pre = ""
            post = ""
        
        if sa and sb:
            seqah.append(f"{cab}{sf}{post}")
        elif sa:
            seqah.append(f"{ca}{sf}{post}")
        elif sb:
            seqah.append(f"{cb}{sf}{post}")
        
    if not suppress:
        print("".join(seqah))
        print("".join(seqbh))
    return "".join(seqah+seqbh)

def highlight_features(seq, features, feature_spans = {}, feature_starts = {}, feature_ends = {}, colors = {}, show_key = True, suppress = True, break_features = False):
    """
    feature_starts: Dict:feature -> List[feature_start_pos]
    feature_ends: Dict:feature -> List[feature_end_pos]
    """
    # Build the colored string
    
    if not feature_spans:
        # feature_spans = {f:[(s, e) for s, e in zip(feature_starts[f], feature_ends[f])]for f in features}
        feature_spans = make_spans(feature_starts, feature_ends)
        # print(feature_spans)
    
    baseline_color = 240
    bc = f"\x1b[38;5;{baseline_color}m"
    
    result = []
    current_color = 0
    last_ansi = bc
    
    for s in features:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    
    key = make_key(features, colors)
    
    result.append(bc)  # Start with baseline color

    for i in range(len(seq)):
        
        current_color = baseline_color
        
        pre_break = False
        done = False
        for f in features:
            for s, e in feature_spans[f]:
                if i < s:
                    continue
                elif i >= e:
                    continue
                else:
                    current_color = colors.get(f)
                    pre_break = i == s 
                    done = True
                if done:
                    break
            if done:
                break

        ansi = f"\x1b[38;5;{current_color}m"

        if ansi != last_ansi:
            if break_features:
                result.append(" ")
            result.append(ansi)
            last_ansi = ansi
        elif break_features and pre_break:
            result.append(" ")
        
        result.append(seq[i])        

    result.append(Colors.RESET)
    if not suppress:
        print("".join(result))
        if show_key:
            print(key)
    return result, key

def highlight_sequence(seq, subseq, color = None, suppress = True, show_key = True):
    
    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there
    
    if not color:
        color = random.randint(20, 230)
    colors = {subseq:color}

    start_pos = find_subsequence(seq, subseq)
    starts, ends = make_start_ends(subseq, start_pos, len(subseq), starts=starts, ends = ends)

    return highlight_features(seq, [subseq], feature_starts = starts, feature_ends = ends, colors=colors, suppress = suppress, show_key = show_key)

def highlight_sequence_fuzzy(seq, subseq, colors = None, max_err = 1, suppress = True, show_key = True, break_features = False):
    
    if not colors:
        colors = []
    if len(colors) <= max_err:
        colors.extend([random.randint(20, 230) for i in range(len(colors), max_err+1)])
    
    feats = [subseq]
    spans = {subseq:[]}
    clrs = {subseq:colors[0]}
    
    key_info = {}
    
    ptrn_str = "(%s){e<=%s}" % (subseq, str(max_err))
    ptrn = regex.compile(ptrn_str, regex.BESTMATCH)
    
    matches = regex.finditer(ptrn, seq)
    
    for m in matches:
        if not m:
            continue
        
        fsubseq = m.groups()[0]
        
        err = sum(m.fuzzy_counts)
        start, end = m.spans()[0]
        
        if not fsubseq in feats:
            feats.append(fsubseq)
            spans[fsubseq] = []
            clrs[fsubseq] = colors[err]
        if not err in key_info:
            key_info[err] = set()
            
        key_info[err].add(fsubseq)
        spans[fsubseq].append((start, end))
    
    out, _ = highlight_features(seq, feats, feature_spans = spans, colors = clrs, show_key = False, suppress = suppress, break_features = break_features)
    key = ", ".join([f"\x1b[38;5;{colors[err]}merr={err}:[{",".join(subseqs)}]" for err, subseqs in key_info.items()]) + "\x1b[0m"
    if show_key and not suppress:
        print(key)
    
    return out, key

def highlight_sequences(seq:str, subseqs:List[str], start_pos = None, do_rc = False, min_len = 5, colors={}, show_key = True, suppress = True):

    starts = {}  # position -> list of sequences starting there
    ends = {}    # position -> list of sequences ending there
    
    for s in subseqs:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    
    if not start_pos:
        start_pos = find_subsequences(seq, subseqs, do_rc = do_rc)
    for s in subseqs:
        if len(s) < min_len:
            continue

        starts, ends = make_start_ends(s, start_pos[s], len(s), starts=starts, ends = ends)

    highlight_features(seq, subseqs, feature_starts = starts, feature_ends = ends, colors = colors, show_key = show_key, suppress = suppress)

def _highlight_runs_auto(seqa, top_seqs, top_datas, **kwargs):
    all_seqs = []
    start_pos = {}
    
    for ns in range(len(top_seqs)):
        run, ind, shift = top_datas[ns]
        start_pos.update({top_seqs[ns][0]: [ind], top_seqs[ns][1]: [ind - shift]})
        all_seqs.extend(top_seqs[ns])
        
    highlight_sequences(seqa, all_seqs, start_pos = start_pos, **kwargs)

def _highlight_runs_diff(seqa, seqb, top_seqs, top_datas, **kwargs):
    
    seqs_a = []
    start_pos_a = {}
    seqs_b = []  
    start_pos_b = {}
    for ns in range(len(top_seqs)):
        run, ind, shift = top_datas[ns]
        start_pos_a.update({top_seqs[ns][0]: [ind]})
        start_pos_b.update({top_seqs[ns][1]: [ind - shift]})
        seqs_a.append(top_seqs[ns][0])
        seqs_b.append(top_seqs[ns][1])
        
    highlight_sequences(seqa, seqs_a, start_pos = start_pos_a, **kwargs)
    highlight_sequences(seqb, seqs_b, start_pos = start_pos_b, **kwargs)
    

def highlight_run(seqa, seqb, top_seqs, top_datas, suppress = False, show_key = True, **kwargs):
    
    colors = kwargs.pop("colors", {})
    
    for s in top_seqs:
        if not s in colors:
            colors[s] = random.randint(20, 230)
    kwargs["colors"] = colors
    
    if seqa == seqb:
        _highlight_runs_auto(seqa, top_seqs, top_datas, suppress = suppress, show_key = show_key, **kwargs)
    else:
        _highlight_runs_diff(seqa, seqb, top_seqs, top_datas, suppress = suppress, show_key = show_key, **kwargs)

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

    result.append(Colors.RESET)

    print("".join(result))

def highlight_sequence_by_span(seq, span_colors = {}, default_color = '\033[97m'):
    # span is (start, stop):color
    
    spans = sorted(span_colors.keys(), key=lambda k:k[0])
    ccurr = default_color
    colored_seq = []
    for i, b in enumerate(seq):
        
        in_span = False
        for st, sp in spans:
            c = span_colors[(st, sp)]
            if i>sp:
                continue
            elif i<st:
                break
            
            if i>=st and i<sp:
                in_span = True
                new_c = c
                if new_c == ccurr:
                    continue
                else:
                    colored_seq.append(new_c)
                    ccurr = new_c
        
        if not in_span:
            colored_seq.append(default_color)
        
        colored_seq.append(b)
    colored_seq.append(Colors.RESET)
    return "".join(colored_seq)

def highlight_dyads(seq, dyads):
    
    sec = {(d.stem_start, d.end_position): random.randint(20, 230) for d in sorted(dyads, key = lambda d:d.stem_start)}
    
    outseq = []
    
    for i in range(len(seq)):
        
        pre = ""
        post = ""
        for st, en in sec:
            if i < st:
                break
            if i > en:
                continue
            c = sec[(st, en)]
            pre = f"\x1b[38;5;{c}m"
            post = Colors.RESET
            
        outseq.append(f"{pre}{seq[i]}{post}")
    
    return "".join(outseq)
 