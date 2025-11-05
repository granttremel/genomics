
import numpy as np

from ggene import seqs, draw
from ggene.seqs import bio, find, process, heal, align
from ggene.seqs.bio import reverse_complement, reverse, complement
from ggene.seqs.heal import Healer

from ggene.genomemanager import GenomeManager

def load_genome():
    return GenomeManager()


def get_test_seq():
    return "AATGAATGAATGACTGAATGAATG"

def get_natural_seq(gm, chr, start, seq_len):
    return gm.get_sequence(chr, start, start+seq_len)

def show_autocorr(data, rcdata, fg, bg, **kwargs):
    
    print(f"data after bin {len(data)}")
    if len(data) > 256:
        rbin = len(data)//256
        data = process.bin_data(data, bin_size = rbin)
        rcdata = process.bin_data(rcdata, bin_size = rbin)
    
    print(f"data before bin {len(data)}")
    res = draw.scalar_to_text_nb(data, fg_color = fg, bg_color = bg, **kwargs)
    rcres = draw.scalar_to_text_nb(rcdata, fg_color = fg-12, bg_color = bg, flip=True, **kwargs)
    for r in res:
        print(r)
        
    for r in rcres:
        print(r)
    print()
    
def render_stats_table(data, label1, data2=None, label2=""):
    outlines = []
    marg = ""*10
    outlines.append("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}".format(*[marg,"Mean", "SD", "CV", "Min", "Max"]))
    outlines.append(f"{label1:<10}{np.mean(data):<10.2g}{np.std(data):<10.2g}{np.std(data)/np.mean(data):<10.2g}{min(data):<10.2g}{max(data):<10.2g}")
    if data2 is not None:
        marg2 = label2
        outlines.append(f"{marg2:<10}{np.mean(data2):<10.2g}{np.std(data2):<10.2g}{np.std(data2)/np.mean(data2):<10.2g}{min(data2):<10.2g}{max(data2):<10.2g}")
    
    for l in outlines:
        print(l)
    print()

def get_rep_templates(hlr):

    highs, lows = hlr.locate_repeats(zscore = 2)
    if not highs:
        print("weak autocorrelation")
        return None, None
    
    tempseqs = [t[2] for t in highs]
    temps_grp = find.group_templates(tempseqs)
    alltemps = {}
    # print("forward")
    for grp, temps in temps_grp.items():
        temps = find.compile_templates(temps)
        tempstr = find.make_template_str(temps)
        alltemps[grp] = temps
        
        # print(f"group template: {grp}")
        # print(f"forward consensus: {tempstr}")
    
    rchighs, rclows = hlr.locate_repeats(zscore = 1.5, do_rc = True, do_low = False)
    rctempseqs = [t[2] for t in rchighs]
    rctemps_grp = find.group_templates(rctempseqs)
    
    rcalltemps = {}
    # print("revcomp:")
    for rcgrp, rctemps in rctemps_grp.items():
        if not rctemps:
            continue
        rctemps = find.compile_templates(rctemps)
        rcalltemps[rcgrp] = rctemps
        rctempstr = find.make_template_str(rctemps)
        
        # print(f"group template: {rcgrp}")
        # print(f"rc consensus: {rctempstr}")

    
    return alltemps, rcalltemps
    
def disp_temps(temps):
    
    tempstr = find.make_template_str(temps)
    print("template:")
    print(tempstr)
    entr = find.get_template_entropy(temps)
    for i in range(len(entr)):
        if entr[i] > 0:
            pos_cons = find.make_template_str([temps[i]])
            print(f"entropy at position {i} with cons. {pos_cons}: {entr[i]:0.3f}")
    print()

def disp_seq(seq, max_disp = 256):
    seq_len = len(seq)
    ndisp = seq_len // max_disp
    for n in range(ndisp):
        print(n, seq[n*max_disp:(n+1)*max_disp])
    print(n+1, seq[(n+1)*max_disp:])
    print()

def search_for_seq(gm, chr, start, seq_len):
    bg, fg = draw.get_color_scheme("test")
    while True:
        print(f"displaying {chr}:{start+shift}-{start+shift+seq_len}")
        seq = get_natural_seq(gm, chr, start+shift, seq_len)

        hlr = Healer(seq)
        
        show_autocorr(hlr.runs, hlr.rcorr, fg, bg, bit_depth = 16)
        render_stats_table(hlr.runs, "Corr", hlr.rcorr, "RCCorr")
    
        res  = input("keep going?")
        if 'n' in res.lower():
            break
        shift += seq_len
    return seq

def test_align(seqa, seqb):
    
    # ts1 = "ATGCATGCATGC"
    # ts2 = "ATATATGCATGCATGC"
    # ts3 = "ATGCACCCCCTGCATGC"
    # ts4 = "CCCCATGCATGCCGTACGTAATGCATGCGGGG"
    # print(ts1)
    
    # test_seqs = [ts2, ts3, ts4]
    # test_seqs = [ts2]
    
    diffs = align.align_sequences3(seqa, seqb)
    print(diffs)
    # ts1f, tsf = align.fill_aligned(seqa, seqb, diffs)
    # print(ts1f)
    # print(tsf)
        # align.print_aligned(ts1, ts, diffs)
    

def main():
    
    bg, fg = draw.get_color_scheme("test")
    
    gm = load_genome()
    
    spots = [
        (1, 1064880, 0),
        (12, 1000043),
        (12, 10000043),
        (12, 1073590),
        (12, 1016000, 0),
        (12, 2255880),
        (12, 2280290),
        (12, 2356020),
        # (12, 1016300, 0),
        (12, 2446120), # deletion in microsatellite in intron (CACNA1C)
        (13, 86039060), # deletion in another msat TGCC in intergenic 
    ]
    
    # seq = get_test_seq()
    
    # print(am)
    chr, start, startb = spots[-1]
    seq_len = 512
    shift = 0
    seq = get_natural_seq(gm, chr, start, seq_len)
    seq = seq[:seq_len]
        # seqb = get_natural_seq(gm, chr, startb, seq_len)
    print(f"len seq: {len(seq)}")

    hlr = Healer(gm, chr, start, seq_len)
    print(f"len hlr seq: {len(hlr.base_seq)}")
    show_autocorr(hlr.runs, hlr.rcruns, fg, bg, bit_depth = 16)
    render_stats_table(hlr.runs, "Corr", hlr.rcruns, "RCCorr")
    show_autocorr(hlr.corr8, hlr.rcorr8, fg, bg, bit_depth = 8)
    show_autocorr(hlr.corr16, hlr.rcorr16, fg, bg, bit_depth = 8)
    show_autocorr(hlr.corr32, hlr.rcorr32, fg, bg, bit_depth = 8)
    
    disp_seq(seq, max_disp = 93)

    return
    
    rpts, rptinds = find.locate_repeats(seq, seq, scale = 512)
    for rpt, (i, rplen) in zip(rpts, rptinds):
        print(rpt, i, rplen)
    print()
    rpts_trim, rptinds_trim = find.trim_repeats(rpts, rptinds, seq_len)
    for rpt, (i, rplen) in zip(rpts_trim, rptinds_trim):
        print(rpt, i, rplen)
        print(seq[i:i+rplen])
    
    # tmps, rctmps = get_rep_templates(hlr)
    # print(tmps)
    # print("forward templates:")
    # for grp, tmp in tmps.items():
    #     print(find.make_template_str(tmp))
    # print()
    # print("rc templates:")
    # for grp,tmp in rctmps.items():
    #     print(find.make_template_str(tmp))
        
    # disp_temps(tmps)
    # disp_temps(rctmps)
    
    # print("forward:")
    # disp_temps(tmps)
    
    # print("revcomp:")
    # disp_temps(rctmps)
    


if __name__=="__main__":
    main()



