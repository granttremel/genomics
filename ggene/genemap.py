
import gzip
import pysam

import logging
logger = logging.getLogger(__name__)
from .features import Feature, Gene

class GeneMap:
    
    def __init__(self, gtf_path = ''):

        if not gtf_path:
            gtf_path = './data/GRCh38_sorted.gtf.gz'
        self.gtf_path = gtf_path
        
        self.tabix = pysam.TabixFile(self.gtf_path)
        
        contigs = self.tabix.contigs
        self.chromes = [a for a in contigs if len(a) < 4]
        self.other = [a for a in contigs if len(a) >= 4]
        
        # self.skip_nonsense = True
        self.only_coding = True
        # SLOW!!
        # print('loading gtf as pyrange')
        # self.gr = pr.read_gtf(gtf_path)
        # print('done loading gtf as pyrange')
    
    def filter_nonsense(self, features):
        
        fs = []
        for f in features:
            if 'info' in f:
                if 'transcript_biotype' in f['info']:
                    if f['info']['transcript_biotype'] == 'nonsense_mediated_decay':
                        continue
            fs.append(f)
        
        return fs
    
    def list_genes(self, chrom, start=0, end=None):
        gs = [g for g in self.fetch(chrom,start, end=end,features = ['gene'])]
        
        for g in gs:
            g['nfeatures'] = self.count(chrom,g['start'],g['end'])
            
        genelist = [
            [g['info'].get('gene_name','?'),g['chrom'],g['start'],g['end'],
             g['end']-g['start'],g['strand'],g['nfeatures']] 
            for g in gs]
        # genestrs = [','.join([str(a) for a in g]) for g in genelist]
        return genelist
    
    def summarize_features(self, chrom, start, end):
        
        featcount = {}
        
        for g in self.fetch(chrom, start, end=end, features = tuple()):
            
            fd1 = featcount.get(g['info']['gene_biotype'],{})
            fc = fd1.get(g['feature'],0)
            featcount[g['feature']] = fc+1
        
        return featcount
    
    def get_gene(self, name, chrom):
        if not chrom:
            # chrom = self.find_gene_chromosome(name)
            chrom, genes = self.find_gene(name)
            if not chrom:
                return False
            
        gene_data = {}
        for feat in Gene.FEATURE_TYPES + Gene.SUB_FEATURE_TYPES:
            infquer = [v for v in self.by_gene_name(chrom, name, features = (feat,))]
            gene_data[feat] = infquer
        return gene_data
    
    def make_gene(self, name, chrom):
        gene_data = self.get_gene(name, chrom)
        newgene = Gene(gene_data)
        return newgene
    
    def find_gene(self, name):
        name = name.upper()
        
        found = False
        genes = []
        
        for chrom in self.chromes:
            logger.info(f'checking chromosome {chrom}')
            r={}
            try:
                gen = self.fetch(chrom, start = 0, info_field = ('gene_name',name))
                for r in gen:
                    if r:
                        found = True
                        logger.info(f'found something: {r}!')
                        break
            except Exception as e:
                logger.info(f'exception occurred searching chromosome {chrom}: {str(e)}')
                pass
            
            if found:
                genes = [r] + [g for g in gen]
                break
        
        if genes == []:
            return None, None
        else:
            return chrom, genes
    
    def by_gene_name(self, chrom, gene_name, features = ('gene')):
        
        inftup = ('gene_name',gene_name)
        res = [v for v in self.fetch(chrom, 0, features=features, info_field = inftup)]
        
        return res
    
    def get_feature(self, chrom, pos):
        fs = self.fetch(chrom, pos, pos + 1, features = tuple())
        flist = [f for f in fs]
        return flist
    
    def get_feature_at_variant(self, var):
        chrom = var.CHROM.lstrip('chr')
        fs = self.fetch(chrom, var.POS, var.POS + 1, features = tuple())
        flist = [f for f in fs]
        return flist
    
    def get_features_at_variant(self, vars):
        feature_list = []
        
        for v in vars:
            chrom = v.CHROM.lstrip('chr')
            fs = self.fetch(chrom, v.POS, v.POS + 1, features = tuple())
            feature_list.append([f['feature'] for f in fs])
        
        return feature_list
    
    def count(self, chrom, start, end=None, features=tuple(), info_field=tuple()):
        chrom = str(chrom)
        if not isinstance(features, tuple):
            features = tuple(features)
        
        if info_field:
            k, v = info_field
            info_query = f"{k} \"{v}\""
        else:
            info_query = ""
        
        if self.only_coding:
            k,v = ('gene_biotype','protein_coding')
            cond1 = f"{k} \"{v}\""
        else:
            cond1=""
        ct = 0
        try:
            for line in self.tabix.fetch(chrom, start, end):
                fields = line.split('\t')
                if features and fields[2] not in features:
                    continue
                if info_query and not info_query in fields[8]:
                    continue
                if cond1 and not cond1 in fields[8]:
                    continue
                ct += 1
        
        except Exception as e:
            logger.error(f'Exception occurred: {str(e)}')
            return  # Region not in file
        return ct
        
    
    def fetch(self, chrom, start, end=None, features=('gene',), info_field = tuple()):
        chrom = str(chrom)
        if not isinstance(features, tuple):
            features = tuple(features)
        
        if info_field:
            k, v = info_field
            info_cond = f"{k} \"{v}\""
        else:
            info_cond = ""
        
        if self.only_coding:
            k,v = ('gene_biotype','protein_coding')
            cond1 = (k,f"{k} \"{v}\"")
            k,v = ('transcript_biotype','protein_coding')
            cond2 = (k,f"{k} \"{v}\"")
        else:
            cond1=cond2=""
        
        try:
            for line in self.tabix.fetch(chrom, start, end):
                fields = line.split('\t')
                if features and fields[2] not in features:
                    continue
                if info_cond and not info_cond in fields[8]:
                    continue
                yield {
                    'chrom': fields[0],
                    'feature': fields[2],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'info': self.parse_info(fields[8])
                }
        except Exception as e:
            logger.error(f'Exception occurred: {str(e)}')
            return  # Region not in file
    
    def fetch_raw(self, chrom, start, end=None, features='gene'):
                
        chrom = str(chrom)
        if not isinstance(features, tuple):
            features = tuple(features)
        try:
            for line in self.tabix.fetch(chrom, start, end):
                fields = line.split('\t')
                if features and fields[2] not in features:
                    continue
                yield fields
        except:
            return
    
    def fetch_all(self, chrom, start, end, features = ['gene']):
        
        gen = self.fetch(chrom, start, end, features = features)
        out = []
        try:
            for g in gen:
                out.append(g)
        except StopIteration:
            pass
        return out
        
    
    def query(self, chrom, start, end, feature = ['gene'], take = 5):
        
        gen = self.fetch(chrom, start, end, features=feature)
        igen = iter(gen)
        out = [next(igen) for i in range(take)]
            
        return out
    
    def neighbors(self, chrom, pos, features = ['gene']):
        if not isinstance(features, tuple):
            features = tuple(features)
        current = None
        next = None
        for g in self.fetch(chrom, pos, features = features):
            if g['end'] < pos:
                # last = g
                pass
            elif g['start'] > pos:
                next = g
                break
            else:
                current = g
        
        step = pos // 10
        cg = current if current else next
        cgstart = cg['start'] if cg else pos
        
        #backtrack to find prior gene
        last = None
        while not last:
            newpos = max(pos - step, 0)
            for g in self.fetch(chrom, newpos, cgstart - 1, features = features):
                if g['end'] < pos:
                    last = g
            step *= 2
            if newpos == 0:
                break
        
        return last, current, next
        
    
    def stream(self, chromosome, features=['gene']):
        """Parse GTF without loading entire file"""
        chromosome = str(chromosome)
        open_func = gzip.open if self.gtf_path.endswith('.gz') else open
        with open_func(self.gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chrom = fields[0]
                feature = fields[2]
                
                # Filter early
                if chromosome and not chrom==chromosome:
                    continue
                if features and feature not in features:
                    continue
                
                # Parse attributes
                attrs = {}
                for attr in fields[8].split(';'):
                    if attr.strip():
                        key_value = attr.strip().split(' ', 1)
                        if len(key_value) == 2:
                            key = key_value[0]
                            value = key_value[1].strip('"')
                            attrs[key] = value
                
                yield {
                    'chrom': chrom,
                    'source': fields[1],
                    'feature': feature,
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'strand': fields[6],
                    'attributes': attrs
                }
    
    def parse_info(self, infostr):
        attrs = {}
        for attr in infostr.split(';'):
            if attr.strip():
                key_value = attr.strip().split(' ', 1)
                if len(key_value) == 2:
                    key = key_value[0]
                    value = key_value[1].strip('"')
                    attrs[key] = value
        return attrs
        
    