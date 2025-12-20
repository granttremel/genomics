import re
from dataclasses import dataclass, field, replace
from typing import Optional, List, Dict, Any, Tuple
import numpy as np

from ggene import DATA_DIR, other_paths
from ggene.database.genome_manager import GenomeManager
from ggene import draw
from ggene.display.colors import Colors
from ggene.draw.scalar_plot import ScalarPlot

from ggene.seqs import process, lambdas, bio, combi, vocab, find, align
    
@dataclass
class CorrelateConfig:
    scale:int = None
    step:int = 1
    keep:int = None
    shift_step:int = 1

    @property
    def options(self):
        return self.__dict__

@dataclass
class SearchConfig:
    method: str
    seq: str
    rcseq: str
    correlate_cfg:CorrelateConfig = field(default_factory=CorrelateConfig)
    ptrn: Optional[re.Pattern] = None
    rcptrn: Optional[re.Pattern] = None
    ali_dict: Optional[Dict[str, str]] = None
    ali_seq: Optional[str] = None
    rc_ali_seq: Optional[str] = None
    convert_cns: str = ""

@dataclass
class SearchResult:
    score: float
    rcscore: float
    argmax: Optional[int] = None
    rc_argmax: Optional[int] = None
    corr: Optional[np.ndarray] = None
    rccorr: Optional[np.ndarray] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    @property
    def max_score(self):
        return max(self.score, self.rcscore)

    @property
    def best_strand(self):
        return "+" if self.score >= self.rcscore else "-"

@dataclass
class SearchStep:
    method: str = "correlate"
    scale_ratio: float = 0.0
    collect_chunks:int = 3
    convert_cns: str = ""
    subseq_len:int = None
    correlate_cfg:CorrelateConfig = field(default_factory=CorrelateConfig)
    num_chunks: Optional[int] = 240
    chunksz: float = None
    topk: int = 3
    min_chunksz: int = 64

    def with_updates(self, **kwargs) -> "SearchStep":
        return replace(self, **kwargs)
    
    def with_corr_config(self, **kwargs):
        ccfg = CorrelateConfig(**kwargs)
        return replace(self, correlate_cfg = ccfg)
    
    def next_pos(self, center, length, last_chunksz):
        
        if self.scale_ratio:
            # print(f"using scale ratio {self.scale_ratio}")
            new_length = length // self.scale_ratio
        else:
            new_length = int(self.collect_chunks * self.num_chunks * length / last_chunksz)
            # print(f"using chunks {self.collect_chunks} from prev length {length}, chunk size {last_chunksz} to new length {new_length}")
        
        # if new_length // self.num_chunks < self.min_chunksz:
        #     new_length = self.num_chunks * self.min_chunksz
        
        new_start = center - new_length // 2
        return new_start, new_length

@dataclass
class TopResult:
    ind:int
    pos:int
    score:float
    is_rc:bool = False
    
    @property
    def scaled_score(self):
        return self.score*(4/3 if self.is_rc else 1.0)

@dataclass
class ChromosomeSearchResult:
    chr: str
    scores: List[float]
    rcscores: List[float]
    tops: List[TopResult]
    chunksz: int
    start: int
    length: int

    @property
    def argmax(self):
        if not self.scores:
            return None
        return int(np.argmax(self.scores))

    @property
    def rc_argmax(self):
        if not self.rcscores:
            return None
        return int(np.argmax(self.rcscores))

    @property
    def best_position(self):
        if self.tops:
            return self.tops[0][0]
        return None
    
    def iter_tops(self):
        for i, res in enumerate(self.tops):
            yield i, res

@dataclass
class HomologyResult:
    query: str
    start: int
    end: int
    center: int
    max_score: float
    chr: str
    seq: str
    features: List[Any]
    data: List[float]
    rcdata: List[float]
    topdata: List[Tuple[int, float]]
    iterations: List[ChromosomeSearchResult] = field(default_factory=list)
    steps_used: List[SearchStep] = field(default_factory=list)
    is_rc:bool = False

    @property
    def length(self):
        return self.end - self.start
    
    def print_summary(self):
        print("Homology search result summary:")
        print(f"Initial query: {self.query}")
        print(f"Final region: {get_position_string(self.chr, self.start, length=self.length)}" + (" (RC)" if self.is_rc else ""))
        print()
        print(f"Sequence with length {self.length}:")
        print(self.seq[:2*256] + "...")
        print("..." + self.seq[-2*256:])
        print()
        
        print("Features:")
        for f in self.features:
            print(f)
        
        pass

def search_region(gm:GenomeManager, chr, start, length, config: SearchConfig) -> SearchResult:

    test_seq = gm.get_sequence(chr, start = start, end = start+length)

    if config:
        method = config.method
    else:
        method = kwargs.get("method")

    if method=="correlate":
        return _search_correlate(test_seq, config=config)
    elif method=="pattern":
        return _search_pattern(test_seq, config=config)
    elif method=="align":
        return _search_align(test_seq, config=config)

    else:
        print("no method provided!")
        return SearchResult(score=-1, rcscore=-1)

def _search_correlate(test_seq, config: SearchConfig) -> SearchResult:

    seq = config.seq
    ali_dict = config.ali_dict

    cf = sf = None

    if ali_dict:
        cf = lambda a, b: ali_dict.get(a, -1) == ali_dict.get(b, -2)
        sf = lambda a, b, _sc: int(ali_dict.get(a, -1) == ali_dict.get(b, -2))
    else:
        cf = lambda a, b: a==b
        sf = lambda a, b, _sc: int(a==b)/_sc

    corr, rccorr = process.correlate(seq, test_seq, comparison_func=cf, score_func = sf, fill = None, **config.correlate_cfg.options)

    # score = np.percentile(corr, [0.9])[0]
    # rcscore = np.percentile(rccorr, [0.9])[0]
    score = max(corr)
    rcscore = max(rccorr)

    return SearchResult(
        score=score,
        rcscore=rcscore,
        argmax=int(np.argmax(corr)),
        rc_argmax=int(np.argmax(rccorr)),
        corr=corr,
        rccorr=rccorr
    )

def _search_pattern(test_seq, config: SearchConfig) -> SearchResult:

    ptrn = config.ptrn
    rcptrn = config.rcptrn

    matches = list(re.finditer(ptrn, test_seq))
    rc_matches = list(re.finditer(rcptrn, test_seq))

    score = len(matches)
    rcscore = len(rc_matches)

    argmax = matches[0].start() if matches else None
    rc_argmax = rc_matches[0].start() if rc_matches else None

    return SearchResult(
        score=score,
        rcscore=rcscore,
        argmax=argmax,
        rc_argmax=rc_argmax,
        extra={"matches": matches, "rc_matches": rc_matches}
    )

def _search_align(test_seq, config: SearchConfig) -> SearchResult:

    seq = config.seq
    rcseq = config.rcseq
    min_seq_len = min(len(seq), len(test_seq))
    max_seq_len = max(len(seq), len(test_seq))
    chk_stride = 2
    chunk_sz = int(min_seq_len / chk_stride)
    nchks = int(max_seq_len/chunk_sz)
    
    # print(f"align with min seq len {min_seq_len}, max {max_seq_len}, chunk size {chunk_sz}, {nchks} chunks")
    
    smax = 0
    argmax =-1
    rcsmax = 0
    rc_argmax = -1
    for nn in range(nchks):
        subseq =test_seq[nn*chunk_sz:(nn+chk_stride)*chunk_sz]
        # print(f"chunk {nn}, subseq len {len(subseq)}")
        score = align.score_sequences(seq, subseq) / chunk_sz
        if score > smax:
            smax = score
            argmax = nn*chunk_sz + chunk_sz//2
        rcscore = align.score_sequences(rcseq, subseq) / chunk_sz
        if rcscore > rcsmax:
            rcsmax = rcscore
            rc_argmax = nn*chunk_sz + chunk_sz//2

    return SearchResult(
        score=smax,
        rcscore=rcsmax,
        argmax = argmax,
        rc_argmax = rc_argmax,
    )

def _search_hmm(test_seq, config: SearchConfig):
    
    
    pass

def search_chr(gm, chr, num_chunks, chunksz, topk, chr_start = None, chr_length = None, buffer = 1e6, config: SearchConfig = None, **search_kwargs) -> ChromosomeSearchResult:

    chr_scores = []
    chr_rcscores = []
    chr_argmaxes = []
    chr_rc_argmaxes = []
    chr_tops = []

    max_ind = gm.gene_map.max_indices.get(chr) - buffer

    if num_chunks:
        chunksz = int((max_ind - buffer) / num_chunks)

    if not chr_start:
        chr_start = int(buffer)
    if not chr_length:
        chr_length = int(max_ind)

    chunksz = int(chunksz)
    nchks = int(chr_length / chunksz)

    print(f"starting {get_position_string(chr, chr_start, length=chr_length)} with chunk size {chunksz}, {nchks} chunks")

    for n in range(nchks):
        start = chr_start + n*chunksz
        length = chunksz

        result = search_region(gm, chr, start, length, config=config, **search_kwargs)

        chr_scores.append(result.score)
        chr_rcscores.append(result.rcscore)
        chr_argmaxes.append(result.argmax if result.argmax is not None else chunksz // 2)
        chr_rc_argmaxes.append(result.rc_argmax if result.rc_argmax is not None else chunksz // 2)

    if topk:
        top_scores = sorted(enumerate(chr_scores), key=lambda k: -k[1])[:topk]
        top_rcscores = sorted(enumerate(chr_rcscores), key=lambda k: -k[1])[:topk]
        # Use actual argmax position within chunk instead of chunk center
        chr_tops = [TopResult(i, chr_start + i*chunksz + chr_argmaxes[i], s, is_rc=False) for i, s in top_scores]
        chr_tops += [TopResult(i, chr_start + i*chunksz + chr_rc_argmaxes[i], s, is_rc=True) for i, s in top_rcscores]
        chr_tops = sorted(chr_tops, key=lambda res: -res.scaled_score)[:topk]

    return ChromosomeSearchResult(
        chr=chr,
        scores=chr_scores,
        rcscores=chr_rcscores,
        tops=chr_tops,
        chunksz=chunksz,
        start=chr_start,
        length=chr_length
    )

def setup_search_options(seq, step:SearchStep,
                        #  method, convert_cns = "", subseq_len = None
                         ) -> SearchConfig:

    print(seq)
    rcseq = bio.reverse_complement(seq)
    print(rcseq)
    
    if step.subseq_len:
        hl = len(seq)//2
        hssl = step.subseq_len//2
        seq = seq[hl-hssl:hl+hssl]
        rcseq = seq[hl-hssl:hl+hssl]
    
    ali_dict = None
    ali_seq = seq
    rc_ali_seq = rcseq
    ptrn = None
    rcptrn = None

    if step.convert_cns:
        ali_dict = {}
        ali_seq = bio.consensus_to_alias(seq, step.convert_cns)
        rc_ali_seq = bio.consensus_to_alias(rcseq, step.convert_cns)
        for b in step.convert_cns:
            bbs = bio.ALIASES.get(b)
            for bb in bbs:
                ali_dict[bb] = b

    if step.method == "pattern":
        ptrn = re.compile(find.consensus_to_re(ali_seq))
        rcptrn = re.compile(find.consensus_to_re(rc_ali_seq))

    return SearchConfig(
        method=step.method,
        seq=seq,
        rcseq=rcseq,
        correlate_cfg = step.correlate_cfg,
        ptrn=ptrn,
        rcptrn=rcptrn,
        ali_dict=ali_dict,
        ali_seq=ali_seq,
        rc_ali_seq=rc_ali_seq,
        convert_cns=step.convert_cns
    )

def do_homology_search(gm:GenomeManager, seq, chrs = None, convert_cns = "SW", method = "correlate", chunksz = 1e4, num_chunks = None, show_plots = True, topk = None) -> Dict[str, ChromosomeSearchResult]:

    if not chrs:
        chrs = [str(i) for i in range(1, 23)] + ["X","Y"]

    config = setup_search_options(seq, method, convert_cns=convert_cns)
    # config = setup_search_options(seq, step)

    buffer = int(1e6)

    results = {}

    for chr in chrs:

        max_ind = int(gm.gene_map.max_indices.get(chr)) - buffer

        if method == "correlate":
            seq = process.pad_sequence(seq, chunksz)

        chr_result = search_chr(gm, chr, num_chunks, chunksz, topk, buffer=buffer, config=config)

        results[chr] = chr_result

        if show_plots:
            print(f"result for chr{chr}")
            show_homology_plots(chr_result.scores, chr_result.rcscores, xmin = buffer, xmax = max_ind, ruler=True)

    return results

def show_homology_plots(chr_scores, chr_rcscores, **kwargs):

    minval = kwargs.get("minval", 0)
    xmin = kwargs.get("xmin")
    xmax = kwargs.get("xmax")
    ruler = kwargs.get("ruler", False) and xmin and xmax

    marks = kwargs.get("marks", [])
    xlabel = None
    if marks:
        xlabel = [" "] * len(chr_scores)
        for m, s in marks:
            if s=="":
                s = "*"
            xlabel[m] = s
        xlabel = Colors.HIGHLIGHT + "".join(xlabel) + Colors.RESET

    max_row = 256
    max_disp = 3*max_row

    if len(chr_scores) > max_disp:
        subdata = resample_data(chr_scores, max_disp)
        rcsubdata = resample_data(chr_rcscores, max_disp)
    else:
        subdata = list(chr_scores)
        rcsubdata = list(chr_rcscores)
    
    rx = xmax-xmin
    
    u, f, fstr = select_unit(rx)
    ux, fx, _ = select_unit((xmin+xmax)/2)
    if f < fx:
        fmtr = lambda v: format_position(v - xmin, unit=u)
    else:
        fmtr = lambda v: format_position(v, unit=ux)

    # print("")

    try:
        sc1 = ScalarPlot(subdata, add_range = True, minval = minval)
        sc2 = ScalarPlot(rcsubdata, add_range = True, minval=minval, xmin = xmin, xmax = xmax, ruler = ruler, ruler_formatter=fmtr, num_labels = 11)
        ScalarPlot.show_paired(sc1, sc2, chunksz = max_row, xlabel = xlabel, center_xlabel = True)
    except Exception as e:
        print(f"error occurred during scalar plots: {str(e)}")
        print(subdata)

def resample_data(data, new_length):
    
    window_sz = len(data) // new_length + 1
    
    subdata = []
    
    for nn in range(new_length):
        
        start = window_sz//2 + nn*window_sz
        stop = start + window_sz
        
        submean = np.mean(data[start:stop])
        
        subdata.append(submean)
    
    return subdata

def search_homology_iter(
    gm: GenomeManager,
    query: str,
    chr: str,
    steps: Optional[List[SearchStep]] = None,
    start:int = None,
    end:int = None,
    show_plots: bool = True,
    **kwargs
) -> HomologyResult:

    buffer = int(1e6)
    if not start:
        start = buffer
    if not end:
        end = int(gm.gene_map.max_indices.get(chr) - buffer)

    length = end - start
    center = (start + end) // 2

    iterations = []
    steps_used = []

    for ni, step in enumerate(steps):

        # config = setup_search_options(seq, step.method, convert_cns=step.convert_cns)
        config = setup_search_options(query, step)

        # iter_chunksz = step.chunksz
        # if step.num_chunks:
        iter_chunksz = length // step.num_chunks
        
        print(f"starting iteration {ni} with method {step.method}, consensus {step.convert_cns}, chunk size {iter_chunksz} x {step.num_chunks}")
        
        chr_result = search_chr(
            gm, chr,
            num_chunks=None,
            chunksz=iter_chunksz,
            topk=step.topk,
            chr_start=start,
            chr_length=length,
            config=config
        )

        iterations.append(chr_result)
        steps_used.append(step)

        if show_plots:
            marks = [(res.ind, str(i)) for i,res in enumerate(chr_result.tops) if res.score > 0.0]
            show_homology_plots(chr_result.scores, chr_result.rcscores, xmin=start, xmax=start+length, ruler=True, marks = marks)

            
        nz = [(i, res) for i, res in chr_result.iter_tops() if res.score > 0.0]
        if step.topk==1:
            nz = [nz[0]]
        
        if len(nz) > 1:
            ind = choose_next_region(chr_result.tops)
        elif len(nz) == 1:
            ind, _ = nz[0]
        else:
            print("no more score..")
            ind = 0
        
        res = chr_result.tops[ind]
        center = res.pos
        max_hom = res.score
        
        start, length = step.next_pos(center, length, iter_chunksz)
        
        if ni < len(steps)-1:    
            pos_str = get_position_string(chr, start, length=length)
            print(f"Finished iteration {ni} ({step.method}), adjusting to {pos_str} about {format_position(center)} where score = {max_hom}\n")
        else:
            print(f"before update: start={start}, length={length}")
            clc_length = int(step.collect_chunks * length / step.num_chunks)
            start = center - length//2
            result_is_rc = res.is_rc
            print(f"collecting sequence length {clc_length} at start {start}, final step collect_chunks = {step.collect_chunks}, is_rc {result_is_rc}")
            
    print(f"Located maximum homology about {format_position(center)} with a score of {max_hom:0.2f}\n")
    
    region_info = [f for f in gm.annotations.stream_by_types(["gene", "repeats"], chr, start=start, end=start+clc_length)][:100]
    region_seq = gm.get_sequence(chr, start=start, end=start+clc_length)
    
    if result_is_rc:
        region_seq = bio.reverse_complement(region_seq)

    return HomologyResult(
        query=query,
        start=start,
        end=start + clc_length,
        center=center,
        max_score=max_hom,
        chr=chr,
        seq=region_seq,
        features=region_info,
        data=chr_result.scores,
        rcdata=chr_result.rcscores,
        topdata=chr_result.tops,
        iterations=iterations,
        steps_used=steps_used,
        is_rc = result_is_rc
    )

def choose_next_region(top_regions):
    
    for i, res in enumerate(top_regions):
        if res.score <= 0.0:
            continue
        print(f"{i}: score = {res.score:0.1f} at {res.pos//1e6:0.2f}M")
        
    res = input("\nselect a region:\n")
    
    v = 0
    if res:
        try:
            v = int(res.strip())
        except:
            print(f"failed to parse input {res} ")
    
    return v

def get_position_string(chr, start, end=None, length=None):
    
    if length:
        end = start + length
    else:
        length = end-start    
    
    su, sf, _ = select_unit(start)
    eu, ef, _ = select_unit(end)
    if sf > ef:
        unit = su
    else:
        unit = eu
    
    startstr = format_position(start, unit=unit)
    endstr = format_position(end, unit=unit)
    
    lu, lf, lfstr = select_unit(length)
    
    if lf<max(sf,ef)/2:
        lengthstr = format_position(length)
        return f"chr{chr}:{startstr}+{lengthstr}"
    else:
        return f"chr{chr}:{startstr}-{endstr}"

def select_unit(pos, unit = None):
    if pos > 0.5e6 or unit=="M":
        unit = "M"
        fact = 1e6
        fstr = "0.2f"
    elif pos > 0.5e3 or unit=="k":
        unit = "k"
        fact = 1e3
        fstr = "0.1f"
    else:
        unit = ""
        fact=1
        fstr = ""
    
    return unit, fact, fstr

def format_position(pos, unit = None):
    unit, fact, fstr = select_unit(pos, unit=unit)
    return format(pos/fact, fstr) + unit
    