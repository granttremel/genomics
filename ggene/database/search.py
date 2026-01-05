
from typing import List, Dict, Tuple, Any, Optional, Callable, TYPE_CHECKING
from dataclasses import dataclass, field
import logging
import re
import time

from ggene import motifs
from ggene.database.ufeature import UFeature
from ggene.seqs import bio
from ggene.seqs.find import consensus_to_re
from ggene.motifs import hmm as ghmm

if TYPE_CHECKING:
    from ggene.database.genome_manager import GenomeManager

logger = logging.getLogger(__name__)
logger.setLevel("CRITICAL")

@dataclass
class SearchResult:
    chrom: str
    start: int
    end: int
    strand: str
    
    seq:str
    full_seq:str
    features:List[UFeature] = field(default_factory = list)
    attributes:Dict[str,Any] = field(default_factory = dict)
    
    
@dataclass
class SearchConfig:
    chromes:List[str] = field(default_factory = list)
    start:int = 1
    end:int = None
    chunksz:int = None
    num_chunks:int = None
    batch_len:int = 32
    personal:bool = False
    collect_sequences:bool = True
    collect_features:bool = True
    do_rc:bool = True
    
    context_sz:int = 32
    hit_predicate:Callable = None
    feature_predicate:Callable = None
    
    limit:int = None
    timeout:float = None # s i guess
    
    _t0 = 0
    
    def test_hit(self, sr):
        if not self.hit_predicate:
            return False
        else:
            return self.hit_predicate(sr)
    
    def test_feature(self, f):
        if not self.feature_predicate:
            return True
        else:
            return self.feature_predicate(f)
        
    def start_timer(self):
        self._t0 = time.perf_counter()
    
    def restart_timer(self):
        self.start_timer()
    
    def get_elapsed(self):
        return time.perf_counter() - self._t0
    
    def check_timer(self):
        dt = self.get_elapsed()
        if self.timeout and dt > self.timeout:
            return False
        else:
            return True

def get_search_config(chromes = [], start = 1, end = None, chunksz = 256, batch_len = 32, personal = False, context_sz = 32, hit_predicate = None, feature_predicate = None, limit = None, timeout = None, do_rc = True ):
    
    return SearchConfig(
        chromes=chromes,
        start=start,
        end=end,
        chunksz = chunksz,
        batch_len=batch_len,
        personal=personal,
        context_sz = context_sz,
        hit_predicate = hit_predicate,
        feature_predicate=feature_predicate,
        limit=limit,
        timeout=timeout,
        do_rc = do_rc
    )
    
class GenomeSearch:
    
    def __init__(self, genome_manager:'GenomeManager'):
        self.gm = genome_manager

    def simple_sequence_search(self, sseq, search_config = None, **kwargs):
        
        if not search_config:
            cfg = SearchConfig(**kwargs)
        else:
            cfg = search_config
        
        bothseqs = [sseq]
        if search_config.do_rc:
            rcsseq = bio.reverse_complement(sseq)
            bothseqs.append(rcsseq)
        done = False
        
        ni = 0
        
        cfg.start_timer()
        for chrom, poss, seqs in self.iter_chromes(cfg.chromes, cfg.start, end=cfg.end, chunksz = cfg.chunksz, batch_len = cfg.batch_len, personal = cfg.personal):
            
            sseq_len = len(sseq)
            
            for pos, seq in zip(poss, seqs):
                
                strand = "+"
                for frseq in bothseqs:
                    insts = find_all_instances(frseq, seq)
                    if insts:
                        logger.info(f"insts: {len(insts)} ({strand})")
                        for ind in insts:
                            hit_start = int(pos+ind)
                            hit_end = int(hit_start+sseq_len)
                            
                            feats = self.collect_features(chrom, hit_start, hit_end, cfg) if cfg.collect_features else []
                            
                            hitseq = seq[ind:ind+sseq_len]
                            fullseq = seq[ind - cfg.context_sz//2:ind+sseq_len+cfg.context_sz//2]
                            
                            res = SearchResult(chrom, hit_start, hit_end, strand, hitseq, fullseq, feats)
                            yield res
                            ni += 1
                    
                    strand = '-'
                    
                    if cfg.limit and ni >= cfg.limit:
                        logger.info(f"collected requested results {cfg.limit}, ending search")
                        done = True
                    if not cfg.check_timer():
                        logger.info(f"Config with timout {cfg.timeout}s timed out, ending search")
                        done = True
            
            if done:
                break
    
    def consensus_search(self, cons, search_config = None, **kwargs):
        
        ptrn = re.compile(consensus_to_re(cons))
        
        yield from self.pattern_search(ptrn, search_config=search_config, **kwargs)
        
    
    def pattern_search(self, ptrn, search_config = None, **kwargs):
        
        if not search_config:
            cfg = SearchConfig(**kwargs)
        else:
            cfg = search_config
        
        if isinstance(ptrn, str):
            from ggene.seqs import find
            ptrn = re.compile(ptrn)
        
        done = False
        
        ni = 0
                
        cfg.start_timer()
        for chrom, poss, seqs in self.iter_chromes(cfg.chromes, cfg.start, end=cfg.end, chunksz = cfg.chunksz, batch_len = cfg.batch_len, personal = cfg.personal):
            
            strand = "+"
            for pos, seq in zip(poss, seqs):
                for inst in re.finditer(ptrn, seq):
                    subpos, subend = inst.regs[0]
                    hit_start = pos+subpos
                    hit_end = pos+subend              
                    
                    feats = self.collect_features(chrom, hit_start, hit_end, cfg) if cfg.collect_features else []
                    
                    hitseq = inst.group()
                    fullseq = seq[subpos - cfg.context_sz//2:subend+cfg.context_sz//2]
                    res = SearchResult(chrom, hit_start, hit_end, strand, hitseq, fullseq, feats)
                    yield res
                    ni += 1
                
                if cfg.limit and ni >= cfg.limit:
                    logger.info(f"collected requested results {cfg.limit}, ending search")
                    done = True
                if not cfg.check_timer():
                    logger.info(f"Config with timout {cfg.timeout}s timed out, ending search")
                    done = True
                    
                if done:
                    break
            if done:
                break
    
    def hmm_search(self, hmm, search_config = None, **kwargs):
                
        if not search_config:
            cfg = SearchConfig(**kwargs)
        else:
            cfg = search_config
        
        ppl = kwargs.get("pipeline", None)
        if not ppl:
            ppl = ghmm.Pipeline(ghmm.AB, ghmm.BG, Z = 10, E = 0.1)
                
        done = False
        
        ni = 0
                
        cfg.start_timer()
        for chrom, poss, seqs in self.iter_chromes(cfg.chromes, cfg.start, end=cfg.end, chunksz = cfg.chunksz, batch_len = cfg.batch_len, personal = cfg.personal):
            
            strand = "+"
            
            poss = list(poss)
            dseqs, names = ghmm.prepare_seqs(seqs)
            
            hits = ppl.search_hmm(hmm, dseqs)
            
            for i, h in enumerate(hits):
                
                strand = "+"
                bn = int(h.name)
                
                if bn > cfg.batch_len:
                    bn //= 2
                    strand = "-"
                
                pos = poss[bn]
                
                offset = pos + bn*cfg.chunksz
                sr = ghmm.SearchResult.from_hit(h, chrom=chrom, offset = offset, strand=strand)
                
                if cfg.test_hit(sr):
                    
                    hit_start = sr.target_from
                    hit_end = sr.target_to
                    
                    try:
                        seq = ghmm.AB.decode(dseqs[bn].sequence)
                    except Exception as e:
                        print(str(e))

                    feats = self.collect_features(chrom, hit_start, hit_end, cfg) if cfg.collect_features else []
                    
                    # indices might not quite be the tea...
                    hitseq = seq[hit_start:hit_end]
                    fullseq = seq[hit_start - cfg.context_sz//2:hit_end+ cfg.context_sz//2]
                    
                    full_sr = SearchResult(chrom, pos + hit_start, pos+hit_end, strand, hitseq, fullseq, features=feats, attributes = sr.__dict__)
                    
                    yield full_sr
                    ni += 1
                else:
                    logger.info(f"rejected hit {sr}")
                
                if cfg.limit and ni >= cfg.limit:
                    logger.info(f"collected requested results {cfg.limit}, ending search")
                    done = True
                if not cfg.check_timer():
                    logger.info(f"Config with timout {cfg.timeout}s timed out, ending search")
                    done = True
                    
                if done:
                    break
            if done:
                break
        
    
    def feature_search(self, search_config = None, **kwargs):
        
        if not search_config:
            cfg = SearchConfig(**kwargs)
        else:
            cfg = search_config
        
        done = False
        
        nf = 0
        
        cfg.start_timer()
        for chrom in self.gm.iter_chromes():
            
            if cfg.chromes and not chrom in cfg.chromes:
                continue
            
            for f in self.gm.annotations.stream_all(chrom, cfg.start, end=cfg.end):
                
                if cfg.feature_predicate(f):
                    if cfg.collect_sequences:
                        seq, full_seq = self.collect_sequences(chrom, f.start, f.end, cfg.context_sz, personal = cfg.personal)
                    
                    sr = SearchResult(chrom, f.start, f.end, f.strand, seq, full_seq, features = [f])
                    
                    yield sr
                    nf += 1
                
                if cfg.limit and nf >= cfg.limit:
                    logger.info(f"collected requested results {cfg.limit}, ending search")
                    done = True
                if not cfg.check_timer():
                    logger.info(f"Config with timout {cfg.timeout}s timed out, ending search")
                    done = True
            
                if done:
                    break
        
            if done:
                break

    def collect_sequences(self, chrom, start, end, contextsz, personal = False):
        r = contextsz//2
        if personal:
            full_seq = self.gm.annotations.get_personal_sequence(chrom, start-r, end+r)
        else:
            full_seq = self.gm.annotations.get_sequence(chrom, start-r, end+r)
        seq = full_seq[r:-r-1]
        
        return "".join(seq), "".join(full_seq)
    
    def collect_features(self, chrom, start, end, cfg):
        return [f for f in self.gm.annotations.query_range(chrom, start, end) if cfg.test_feature(f)]
    
    def iter_chromes(self, chromes = [], start=1e6, end=None, chunksz = None, num_chunks = None, batch_len = 1, personal = False):
        
        for chrom in self.gm.iter_chromes():
            
            if chromes and not chrom in chromes:
                continue
            
            logger.info(f"GenomeSearch beginning chrom {chrom}")
            
            for poss, seqs in self.iter_seq(chrom, start, end, chunksz=chunksz, num_chunks = num_chunks, batch_len = batch_len, personal = personal):
                yield chrom, poss, seqs
            
    
    def iter_seq(self, chrom, start, end = None, chunksz = None, num_chunks = None, batch_len = 1, personal = False):
        
        if not chunksz:
            chunksz = (start - end) / num_chunks
        start = int(start)
        end = int(end) if end else None
        chunksz = int(chunksz)
        ntake = batch_len * chunksz
        
        if not end:
            end = self.gm.get_chrom_len(chrom)
        
        # logger.info(f"enter iter_seq with {chrom}, {start}, {end}, {chunksz}, {batch_len}")
        
        while start < end:
            
            take_end = min(end, start+ntake)
            if personal:
                seq = self.gm.annotations.get_personal_sequence(chrom, start, take_end)
            else:
                seq = self.gm.annotations.get_sequence(chrom, start, take_end)
            
            qseqs = []
            starts = []
            for nb in range(batch_len):
                qseq = seq[nb*chunksz:(nb+1)*chunksz]
                qseqs.append(qseq)
                starts.append(start + nb*chunksz)
                
            yield starts, qseqs
        
            start += ntake

def find_all_instances(kseq, qseq):
    
    rmind = 0
    klen = len(kseq)
    
    insts = []
    
    done = False
    while not done:
        
        try:
            ssind = qseq.index(kseq)
        except:
            done = True
            break
        
        insts.append(ssind+rmind)
        
        rmind += ssind+klen
        qseq = qseq[ssind+klen:]
        
        if not qseq:
            break
    
    return insts
        