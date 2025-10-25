


from typing import Dict, List, Optional, Tuple, Union, Iterator, Any, Callable
from collections import OrderedDict
import os
import traceback
from cyvcf2 import Variant

import sys

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
from . import utils


MIN_ATTRIBUTES = [
    'type', 'sfid', 'start', 'end', 'start_relative', 'end_relative',
    'transcript_name', 'exon_number', 'protein_version', 'alt', 'ref', 'qual'
]
def variant_to_dict(variant, gene_name, chrom, strand):
    delta = len(variant.ALT[0]) - len(variant.REF)
    variant_dict = {
        "feature": "variant",
        "gene_name": gene_name,
        "chrom": chrom,
        "start": variant.start,
        "end": variant.end,
        "strand": strand,
        "pos": variant.POS,
        "alt": variant.ALT[0] if variant.ALT else None,
        "ref": variant.REF,
        "qual": variant.QUAL,
        "delta": delta,
        "var_type": getattr(variant, 'var_subtype', ''),
        "genotype": variant.gt_bases[0]
        # "genotype": str(variant.genotype) if hasattr(variant, 'genotype') else ''
    }
    if variant.num_het > 0:
        zyg = "heterozygous"
    elif variant.num_hom_alt > 0:
        zyg = "homozygous alt"
    elif variant.num_hom_ref > 0:
        zyg = "homozygous ref"
    else:
        zyg = ""
    variant_dict["zygosity"] = zyg
    
    return variant_dict

def order(feature_type):
    return _order.get(feature_type,_order['default'])

_order = {
    'gene':0,
    'transcript':1,
    'exon':2,
    'intron':2,
    'CDS':3,
    'variant':4,
    'five_prime_utr':4,
    'three_prime_utr':4,
    'start_codon':4,
    'stop_codon':4,
    'default':4
}
_orders = tuple(range(max(_order.values())))

class Feature:
    """Represents a genomic feature like genes, exons, transcripts, etc."""
    
    FEATURE_TYPES = [
        'gene', 'CDS', 'transcript', 'five_prime_utr', 'start_codon',
        'exon', 'stop_codon', 'three_prime_utr', 'variant'
    ]
    
    MIN_ATTRIBUTES = [
        'type', 'sfid', 'start', 'end', 'start_relative', 'end_relative',
        'transcript_name', 'exon_number', 'protein_id', 'alt', 'ref', 'qual','delta','var_type'
    ]
    
    SHOW_DIGITS = 5
    
    def __init__(self, feature_data: Dict[str, Any]) -> None:
        """Initialize a Feature from a dictionary.
        
        Args:
            feature_dict: Dictionary containing feature information
        """
        if isinstance(feature_data, Variant):
            feature_data = variant_to_dict(feature_data)
        
        self._attributes: List[str] = ['type']
        self.type: str = feature_data.get('feature', 'unknown')
        self.sfid: str = ''
        self.subfeatures: List['Feature'] = []
        self.num_subfeatures = 0
        
        # Parent hierarchy tracking
        self.parent_gene: Optional['Gene'] = None
        self.parents: Optional[List['Feature']] = []
        
        self._hashes = set()
        # Set basic attributes
        for key, value in feature_data.items():
            if key in ('info', 'feature'):
                continue
            setattr(self, key, utils.try_parse_number(value))
            self._attributes.append(key)
        
        # Set relative positions
        self.start_relative: int = 0
        if hasattr(self, 'end') and hasattr(self, 'start'):
            self.end_relative: int = self.end - self.start
        else:
            self.end_relative: int = 0
        self._attributes.extend(['start_relative', 'end_relative'])
          
        # Process info section
        if 'info' in feature_data:
            for key, value in feature_data['info'].items():
                setattr(self, key, utils.try_parse_number(value))
                self._attributes.append(key)
        
    def set_feature_id(self, feature_id: Union[str, int]) -> None:
        """Set the subfeature ID.
        
        Args:
            feature_id: ID to assign to the subfeature
        """
        self.sfid = f'{self.type}-{str(feature_id)}'
    
    def refer_to(self, reference: 'Feature') -> None:
        """Set relative positions based on a reference feature.
        
        Args:
            reference: Reference feature to calculate relative positions from
        """
        if hasattr(self, 'start') and hasattr(reference, 'start'):
            self.start_relative = self.start - reference.start
            self.end_relative = self.end - reference.start
    
    def insert_feature(self, subfeature: 'Feature', set_parent = True) -> None:        
        """Add a subfeature to this feature.
        
        Args:
            subfeature: Feature to add as a subfeature
        """
        self.subfeatures.append(subfeature)
        self._hashes.add(hash(subfeature))
        self.num_subfeatures += 1
        
        if set_parent:
            subfeature.add_parent(self)
    
    def try_insert_feature(self, subfeature, set_parent = True):
        
        o1 = order(self.type)
        o2 = order(subfeature.type)
        if o2 <= o1:
            return False
        h = hash(subfeature)
        if subfeature in self.subfeatures or h in self._hashes:
            return False

        if self.type == "gene" and hasattr(subfeature, "gene_name"):
            if self.name == subfeature.gene_name:
                return True
        
        elif self.type == "transcript" and hasattr(subfeature, "transcript_name"):
            if subfeature.transcript_name == self.transcript_name:
                self.insert_feature(subfeature, set_parent = set_parent)
                return True
        
        elif self.type == "exon" and hasattr(subfeature, "exon_number"):
            if self.exon_number == subfeature.exon_number:
                self.insert_feature(subfeature)
                return True
        
        if subfeature in self:
            self.insert_feature(subfeature)
            return True
        
        return False
    
    def order(self):
        return order(self.type)
    
    def add_parent(self, parent: 'Feature') -> None:
        """Set the direct parent of this feature.
        
        Args:
            parent: The direct parent feature
        """
        o1 = self.order()
        o2 = parent.order()
        if o2 >= o1:
            return False
        h = hash(parent)
        if h in self._hashes:
            return False
        # if not self.type == "variant":
        #     logger.debug(f"adding parent {parent} to {self}")
        self.parents.append(parent)
        self._hashes.add(h)
        
        return True
        
    def set_parent_gene(self, gene: 'Gene') -> None:
        """Set the parent gene for this feature.
        
        Args:
            gene: The parent gene
        """
        self.parent_gene = gene
        self.refer_to(gene)
    
    def get_hierarchy_path(self) -> List['Feature']:
        """Get the full hierarchy path from gene to this feature.
        
        Returns:
            List of features from gene down to this feature
        """
        path = []
        if self.parent_gene:
            path.append(self.parent_gene)
        if self.parent_transcript:
            path.append(self.parent_transcript)
        if self.parent_exon:
            path.append(self.parent_exon)
        if self.parent_cds:
            path.append(self.parent_cds)
        path.append(self)
        return path
    
    def show_hierarchy(self):
        
        hier = [str(self)]
        if self.parents:
            heir += self.parents[0].show_hierarchy()
        
        return '>'.join(reversed(hier))
        
    def five_prime(self):
        if self.strand == '+':
            return self.start
        else:
            return self.end
        
    def three_prime(self):
        if self.strand == '+':
            return self.end
        else:
            return self.start
        
    def __contains__(self, other: 'Feature') -> bool:
        """Check if this feature contains another feature.
        
        Args:
            other: Feature to check containment for
            
        Returns:
            True if this feature contains the other feature
        """
        return (hasattr(self, 'start') and hasattr(self, 'end') and
                hasattr(other, 'start') and hasattr(other, 'end') and
                self.start <= other.start and self.end >= other.end)
    
    def count_overlap(self, other):
        #needs work
        return max(min(self.end - other.start, other.end - self.start),0)
    
    def overlaps(self, other: 'Feature') -> bool:
        """Check if this feature overlaps with another feature.
        
        Args:
            other: Feature to check overlap with
            
        Returns:
            True if features overlap
        """
        if not all(hasattr(f, attr) for f in (self, other) for attr in ('start', 'end')):
            return False
        discontinuous = (self.start > other.end or other.start > self.end)
        return not discontinuous
    
    def sort(self) -> None:
        """Sort subfeatures based on strand direction."""
        if not hasattr(self, 'strand') or not self.subfeatures:
            return
            
        self.subfeatures = self._sort_features(self.subfeatures)
        
    def _get_sort_key(self, force_forward = False):
        if not hasattr(self, 'strand') or self.strand == '+' or force_forward:
            return lambda f: getattr(f, 'start', 0)
        else:
            return lambda f: -getattr(f, 'end', 0)
    
    def _sort_features(self, feature_list: List[Any], force_forward: bool = False) -> List[Any]:
        """Sort features based on strand direction.
        
        Args:
            feature_list: List of features to sort
            force_forward: Force forward sorting regardless of strand
            
        Returns:
            Sorted list of features
        """
        if not feature_list:
            return []
        
        sort_key = self._get_sort_key()
        
        return sorted(feature_list, key=sort_key)
    
    def to_dict(self, abbreviate: bool = True, 
                parent_keys: Optional[List[str]] = None,
                include_parents: bool = True) -> Dict[str, Any]:
        """Convert feature to dictionary representation.
        
        Args:
            abbreviate: Whether to abbreviate subfeatures
            parent_keys: Keys from parent object (for filtering)
            include_parents: Whether to include parent information
            
        Returns:
            Dictionary representation of the feature
        """
        
        if parent_keys is None:
            parent_keys = []
            
        keys = [attr for attr in self.MIN_ATTRIBUTES if hasattr(self, attr)]
        result = {attr: getattr(self, attr) for attr in keys}
        
        # Add parent information
        if include_parents:
            if self.parents:
                result['parents'] = [p.to_str_abbr() for p in self.parents]
        
        if self.subfeatures:
            if abbreviate:
                result['subfeatures'] = [sf.to_str_abbr() for sf in self.subfeatures]
            else:
                result['subfeatures'] = [sf.to_dict(parent_keys = list(result.keys()), include_parents=include_parents) for sf in self.subfeatures]
        
        return result
    
    def __str__(self) -> str:
        """String representation of the feature."""
        
        if hasattr(self,'sfid') and self.sfid:
            _id = self.sfid
        else:
            _id = self.type
        
        if not hasattr(self, 'start') or not hasattr(self, 'end'):
            return f'{_id}(no position data)'
            
        if self.SHOW_DIGITS > 0:
            start_str = '..' + str(self.start)[-self.SHOW_DIGITS:]
            end_str = '..' + str(self.end)[-self.SHOW_DIGITS:]
        else:
            start_str = str(self.start)
            end_str = str(self.end)
        
        parts = []
        if hasattr(self, 'gene_name'):
            chrom = getattr(self, 'chrom', '?')
            strand = getattr(self, 'strand', '?')
            parts.append(f"gene={self.gene_name} ({chrom}:{start_str}-{end_str},{strand})")
            
        if hasattr(self, "transcript_name"):
            parts.append(f"transcript={self.transcript_name}")
        if hasattr(self, "exon_number"):
            parts.append(f"exon={self.exon_number}")
        if self.subfeatures:
            parts.append(f"{len(self.subfeatures)} subfeatures")
            
        content = ','.join(parts) if parts else f"{start_str}-{end_str}"
        return f'{_id}({content})'
    
    def to_str_abbr(self):
        if self.type == 'variant':
            return f'{self.sfid} ({self.var_type} {self.delta})'
        elif self.type == "gene":
            return f'{self.type} {self.name}'
        elif not self.sfid:
            return f'unknown {self.type}'
        else:
            return self.sfid
    
    def __repr__(self) -> str:
        return str(self)
    
    def __len__(self) -> int:
        """Return the length of the feature."""
        if not hasattr(self, 'start') or not hasattr(self, 'end'):
            return 0
        return abs(self.end - self.start)
    
    def __hash__(self) -> int:
        """Return hash based on type and position."""
        start = getattr(self, 'start', 0)
        end = getattr(self, 'end', 0)
        return hash((self.type, start, end))
    
    @classmethod
    def from_dict(cls, feature_dict: Dict[str, Any]) -> 'Feature':
        """Reconstruct a Feature from a dictionary.
        
        Args:
            feature_dict: Dictionary containing feature data
            
        Returns:
            Reconstructed Feature object
        """
        # Prepare the feature dictionary
        prepared_dict = {}
        
        # Copy basic attributes
        for key, value in feature_dict.items():
            if key not in ['subfeatures', 'parents']:
                prepared_dict[key] = value
        
        # Ensure feature type is set
        if 'type' in feature_dict and 'feature' not in prepared_dict:
            prepared_dict['feature'] = feature_dict['type']
        elif 'feature' not in prepared_dict:
            prepared_dict['feature'] = 'unknown'
        
        # Create the feature
        feature = cls(prepared_dict)
        
        # Note: subfeatures and parents will be restored later by the Gene class
        # Store references for later linking
        if 'subfeatures' in feature_dict:
            feature._subfeature_refs = feature_dict['subfeatures']
        if 'parents' in feature_dict:
            feature._parent_refs = feature_dict['parents']
        
        return feature

class Gene(Feature):
    """Represents a gene with its transcripts, exons, and other features."""
    
    FEATURE_TYPES = ['gene', 'transcript', 'exon']
    SUB_FEATURE_TYPES = [
        'CDS', 'intron', 'five_prime_utr', 'start_codon', 
        'stop_codon', 'three_prime_utr'
    ]
    
    def __init__(self, feature_map: Dict[str, List[Dict]]) -> None:
        """Initialize a Gene from a feature map.
        
        Args:
            feature_map: Dictionary mapping feature types to feature data
        """
        if not feature_map.get('gene'):
            raise ValueError("Gene feature map must contain 'gene' entry")
            
        super().__init__(feature_map['gene'][0])
        self.name: str = getattr(self, 'gene_name', '')
        self._attributes.insert(0, 'name')
        
        # ID tracking sets
        self.ids: Dict[str,set] = {}
        
        # Counters
        self.counts: Dict[str,int] = {}
        
        # Feature collections
        self._all_features: List[Feature] = [] # place to put things before organization
        self.transcripts: OrderedDict = OrderedDict()
        self.exons: OrderedDict = OrderedDict()
        self.proteins: OrderedDict = OrderedDict()
        self.cds: List[Feature] = []
        self.variants: List[Feature] = []
        self.subfeatures: List[Feature] = []
        
        self._import_features(feature_map)
        self._organize()
        self._establish_hierarchy()
        self._sort_all()
    
    def get_transcript(self, 
                      transcript_name: Optional[str] = None,
                      transcript_key: Optional[str] = None, 
                      transcript_value: Optional[Any] = None, 
                      pos: Optional[int] = None) -> Optional[Feature]:
        """Get transcript by name, attribute, or position.
        
        Args:
            transcript_name: Name of the transcript
            transcript_key: Attribute name to search by
            transcript_value: Value to match for the attribute
            pos: Position to find transcript containing it
            
        Returns:
            Transcript feature or None if not found
        """
        if transcript_name:
            return self.transcripts.get(transcript_name)
        
        if transcript_key and transcript_value:
            for transcript in self.transcripts.values():
                if hasattr(transcript, transcript_key):
                    if getattr(transcript, transcript_key) == transcript_value:
                        return transcript
        
        if pos is not None:
            for transcript in self.transcripts.values():
                if (hasattr(transcript, 'start') and hasattr(transcript, 'end') and
                    transcript.start <= pos <= transcript.end):
                    return transcript
        
        return None
    
    def get_exon(self, exon_number: Optional[int] = None, 
                pos: Optional[int] = None) -> Optional[Feature]:
        """Get exon by number or position.
        
        Args:
            exon_number: Exon number to find
            pos: Position to find exon containing it
            
        Returns:
            Exon feature or None if not found
        """
        if exon_number is not None:
            return self.exons.get(exon_number)
        
        if pos is not None:
            for exon in self.exons.values():
                if (hasattr(exon, 'start') and hasattr(exon, 'end') and
                    exon.start <= pos <= exon.end):
                    return exon
        
        return None 
    
    def get_feature(self, type, pf = None) -> Feature:
        
        if type == "gene":
            return self
        
        if type == "transcript":
            checklist = self.transcripts.values()
        elif type == "exon":
            checklist = self.exons.values()
        elif type == "CDS":
            checklist = self.cds
        elif type == "variant":
            checklist = self.variants
        else:
            checklist = self.subfeatures
            
        out = list()
        
        if pf is None:
            pf = lambda i,v: True
        
        for i,sf in enumerate(checklist):
            if not sf.type == type:
                continue
            
            if not pf(i, sf):
                continue
                
            out.append(sf)
        
        return out

    def add_new_protein(self, protein_id) -> Feature:
        
        newp = Feature(dict(
            feature="protein",
            protein_id=protein_id,
            gene_name=self.name,
        ))
        self.proteins[protein_id] = newp
        pids = self.ids.get("proteins",set())
        pids.add(protein_id)
        self.ids["proteins"] = pids
        newp.set_feature_id(protein_id)
        
        logger.debug(f'added new protein {newp}')
        
        ct = self.counts.get("proteins", 0)+1
        self.counts["proteins"]=ct
        return newp
        
    def _import_features(self, feature_map) -> None:
        
        for k in feature_map:
            if k == 'gene':
                continue
            for eachf in feature_map[k]:
                if isinstance(eachf, Variant):
                    vdict = variant_to_dict(eachf, gene_name=self.name,chrom=self.chrom,strand=self.strand)
                    f = Feature(vdict)
                else:
                    f = Feature(eachf)
                self.try_insert_feature(f)
        
    def _organize(self) -> None:
        
        for sf in self._all_features:
            ct = self.counts.get(sf.type,0)+1
            if sf.type == "transcript":
                _id = sf.transcript_name
                self.transcripts[_id] = sf
            elif sf.type == "exon":
                _id = sf.exon_number
                k = hash(sf)
                self.exons[k] = sf
            elif sf.type == "CDS":
                _id = ct
                self.cds.append(sf)
            elif sf.type == "variant":
                _id = ct
                self.variants.append(sf)
            else:
                _id = ct
                self.subfeatures.append(sf)
                
            if hasattr(sf, "protein_id") and not sf.protein_id in self.proteins:
                _ = self.add_new_protein(sf.protein_id)
            
            
            sf.set_feature_id(_id)
            sf.set_parent_gene(self)
            
            ids = self.ids.get(sf.type,set())
            ids.add(sf.sfid)
            self.ids[sf.type] = ids 
            
            self.counts[sf.type] = ct
        
        
    def _establish_hierarchy(self):
        
        #adds transcripts to gene
        for ts in self.transcripts.values():
            # self.try_insert_feature(ts)
            # res = self.try_insert_feature(ts)
            res = False
            pass
            if res:
                pass
        
        #add exons to transcripts
        for ex in self.exons.values():
            for _sf in self._all_features:
                res = _sf.try_insert_feature(ex)
                if res:
                    pass
        
        #add cds to exons
        for cds in self.cds:
            for _sf in self._all_features:
                res = _sf.try_insert_feature(cds)
                if res:
                    pass
        
        self._infer_introns()
        
        #add subfeats to feats
        for sf in self.subfeatures:
            for _sf in self._all_features:
                res = _sf.try_insert_feature(sf)
                if res:
                    pass
        
        #add variants to feats
        for var in self.variants:
            for _sf in self._all_features:
                res = _sf.try_insert_feature(var)
                if res:
                    pass
        
    def insert_feature(self, subfeature: 'Feature', set_parent = True) -> None:        
        """Add a feature to this gene.
        
        Args:
            subfeature: Feature to add as a feature
        """
        self._all_features.append(subfeature)
        self._hashes.add(hash(subfeature))
    
    def try_insert_feature(self, subfeature: 'Feature', set_parent = False) -> bool:
        h = hash(subfeature)
        if subfeature in self._all_features or h in self._hashes:
            # print(f'cannot insert {subfeature} into {self}: already present')
            return False
        # if subfeature.type == "exon":
        #     for ex in self.
        
        self.insert_feature(subfeature, set_parent = set_parent)
        return True
    
    def sort_all_features(self):
        self._all_features = self._sort_features(self._all_features)
    
    def _sort_all(self):
        sort_key = self._get_sort_key()
        dict_sort_key = lambda item : sort_key(item[1])
        self.transcripts = OrderedDict(sorted(self.transcripts.items(), key = dict_sort_key))
        self.exons = OrderedDict(sorted(self.exons.items(), key = dict_sort_key))
        self.proteins = OrderedDict(sorted(self.proteins.items(), key = dict_sort_key))
        
        self.cds = self._sort_features(self.cds)
        self.subfeatures = self._sort_features(self.subfeatures)
        self.variants = self._sort_features(self.variants)
        
        for f in self._all_features:
            f.sort()
        
    
    def _infer_introns(self):
        # logger.debug(f'Inferring introns from {len(self.exons)} exons')

        istart = -1
        iend = -1
        
        introns = list()
        n = 1
        for exid, ex in self.exons.items():
            iend = ex.start
            
            if istart >= 0:
                #make new intron
                intr = self.make_intron(len(introns)+1, istart, iend)
                intr.set_parent_gene(self)
                res = self.try_insert_feature(intr)
                if res:
                    introns.append(intr)
                    intr.set_feature_id(n)
                    self.subfeatures.append(intr)
                    n+=1
                    
            
            istart = ex.end
        
        self.sort_all_features()
        logger.debug(f'Found {len(introns)} introns')
        
    def make_intron(self, number, start, end):
        
        newint = Feature(
            {
                "feature":"intron",
                "chrom":self.chrom,
                "start":start,
                "end":end,
                "strand":self.strand,
                "intron_number":number,
                "gene_name": self.name,
            })
        
        return newint
    
    def sort(self, featlist, force_forward = False):
        
        if not featlist:
            return []
        
        if self.strand == '+' or force_forward:
            sk = lambda f:f.start
        else:
            sk = lambda f:-f.end
        
        fd = list(sorted(featlist, key = sk))
        return fd
    
    def print_hierarchy(self):
        
        for sf in self._all_features:
            logger.debug(sf.show_hierarchy)
        
        
    def to_dict(self) -> Dict[str, Any]:
        """Convert gene to dictionary representation.
        
        Args:
            include_ordered_features: Whether to include 5' to 3' ordered features per transcript
            
        Returns:
            Dictionary representation of the gene
        """
        dict_attrs = self._attributes + [
            'transcript_ids', 'exon_ids', 'protein_ids'
        ]
        result = {attr: getattr(self, attr) for attr in dict_attrs if hasattr(self, attr)}
        
        # Basic transcript information
        result['transcripts'] = {
            transcript_name: transcript.to_dict(abbreviate=True, parent_keys=list(result.keys())) 
            for transcript_name, transcript in self.transcripts.items() if transcript
        }
        
        # Exon information
        result['exons'] = {
            exon_num: exon.to_dict(abbreviate=True, parent_keys=list(result.keys())) 
            for exon_num, exon in self.exons.items() if exon
        }
        
        # Protein information
        result['proteins'] = dict()
        for pk in self.proteins:
            prote = self.proteins[pk]
            result['proteins'][pk] = prote.to_dict()
        
        # CDS features
        result['cds'] = [
            subfeature.to_dict(parent_keys=dict_attrs, include_parents=True) 
            for subfeature in self.cds if subfeature
        ]
        
        # Other subfeatures
        result['subfeatures'] = [
            subfeature.to_dict(parent_keys=dict_attrs, include_parents=True) 
            for subfeature in self.subfeatures if subfeature
        ]
        
        # Variants
        result['variants'] = [
            variant.to_dict(include_parents=True) 
            for variant in self.variants if variant
        ]
        
        return result
    
    @classmethod
    def from_dict(cls, gene_dict: Dict[str, Any]) -> 'Gene':
        """Reconstruct a Gene object from a dictionary (e.g., from JSON).
        
        Args:
            gene_dict: Dictionary containing gene data
            
        Returns:
            Reconstructed Gene object with all relationships
        """
        # First, create a feature map structure expected by __init__
        feature_map = {'gene': [{}]}
        
        # Copy basic gene attributes
        gene_attrs = {}
        for key, value in gene_dict.items():
            if key not in ['transcripts', 'exons', 'proteins', 'cds', 'subfeatures', 'variants',
                          'transcript_features_ordered', 'transcript_ids', 'exon_ids', 
                          'protein_ids', 'num_transcripts', 'num_exons', 'num_proteins',
                          'ids', 'counts']:
                gene_attrs[key] = value
        
        gene_attrs['feature'] = 'gene'  # Ensure feature type is set
        feature_map['gene'][0] = gene_attrs
        
        # Create the gene object with just the basic attributes
        gene = cls(feature_map)
        
        # Now rebuild all the features and relationships
        try:
            gene._rebuild_from_dict(gene_dict)
        except Exception as e:
            logger.error(f"Exception rebuilding gene: {str(e)}")
            import traceback
            traceback.print_exc()
        
        return gene
    
    def _rebuild_from_dict_old(self, gene_dict: Dict[str, Any]) -> None:
        """Rebuild all features and relationships from dictionary.
        
        Args:
            gene_dict: Dictionary containing complete gene data
        """
        # Clear existing collections
        self.transcripts.clear()
        self.exons.clear()
        self.proteins.clear()
        self.cds.clear()
        self.subfeatures.clear()
        self.variants.clear()
        self._all_features.clear()
        
        # self.transcript_ids.clear()
        # self.exon_ids.clear()
        # self.protein_ids.clear()
        # self.cds_ids.clear()
        # self.subfeature_ids.clear()
        # self.variant_ids.clear()
        
        # Store references to all features by sfid for linking
        feature_registry = {}
        
        # First pass: Create all feature objects
        
        # Create transcripts
        if 'transcripts' in gene_dict:
            for transcript_name, transcript_data in gene_dict['transcripts'].items():
                transcript = Feature(self._prepare_feature_dict(transcript_data, 'transcript'))
                transcript.transcript_name = transcript_name
                self.transcripts[transcript_name] = transcript
                # self.transcript_ids.add(transcript_name)
                feature_registry[transcript.sfid] = transcript
                transcript.set_parent_gene(self)
        
        # Create exons
        if 'exons' in gene_dict:
            for exon_num, exon_data in gene_dict['exons'].items():
                exon = Feature(self._prepare_feature_dict(exon_data, 'exon'))
                exon_number = int(exon_num) if isinstance(exon_num, str) else exon_num
                self.exons[exon_number] = exon
                self.exon_ids.add(exon_number)
                feature_registry[exon.sfid] = exon
                exon.set_parent_gene(self)
        
        # Create CDS features
        if 'cds' in gene_dict:
            for cds_data in gene_dict['cds']:
                cds = Feature(self._prepare_feature_dict(cds_data, 'CDS'))
                self.cds.append(cds)
                feature_registry[cds.sfid] = cds
                cds.set_parent_gene(self)
        
        # Create other subfeatures
        if 'subfeatures' in gene_dict:
            for subfeature_data in gene_dict['subfeatures']:
                feature_type = subfeature_data.get('type', 'unknown')
                subfeature = Feature(self._prepare_feature_dict(subfeature_data, feature_type))
                self.subfeatures.append(subfeature)
                # self.subfeature_ids.add(subfeature.sfid)
                feature_registry[subfeature.sfid] = subfeature
                subfeature.set_parent_gene(self)
        
        # Create variants
        if 'variants' in gene_dict:
            for variant_data in gene_dict['variants']:
                variant = Feature(self._prepare_feature_dict(variant_data, 'variant'))
                self.variants.append(variant)
                feature_registry[variant.sfid] = variant
                variant.set_parent_gene(self)
        
        # Second pass: Rebuild parent-child relationships
        
        # Link transcripts to their exons and subfeatures
        for transcript_name, transcript_data in gene_dict.get('transcripts', {}).items():
            transcript = self.transcripts[transcript_name]
            
            # Add subfeatures to transcript
            for sf_ref in transcript_data.get('subfeatures', []):
                sf_id = self._extract_sfid(sf_ref)
                if sf_id in feature_registry:
                    subfeature = feature_registry[sf_id]
                    transcript.subfeatures.append(subfeature)
                    subfeature.set_parent(transcript)
                    subfeature.set_parent_transcript(transcript)
        
        # Link exons to their subfeatures
        for exon_num, exon_data in gene_dict.get('exons', {}).items():
            exon = self.exons[int(exon_num) if isinstance(exon_num, str) else exon_num]
            
            # Set transcript parent
            if hasattr(exon, 'transcript_name') and exon.transcript_name in self.transcripts:
                exon.set_parent_transcript(self.transcripts[exon.transcript_name])
            
            # Add subfeatures to exon
            for sf_ref in exon_data.get('subfeatures', []):
                sf_id = self._extract_sfid(sf_ref)
                if sf_id in feature_registry:
                    subfeature = feature_registry[sf_id]
                    exon.subfeatures.append(subfeature)
                    subfeature.set_parent(exon)
                    subfeature.set_parent_exon(exon)
        
        # Rebuild parent relationships based on stored parent info
        for sfid, feature in feature_registry.items():
            # Check if parent info was stored in the data
            parent_info = None
            
            # Find the original data for this feature
            for exon_data in gene_dict.get('exons', {}).values():
                if exon_data.get('sfid') == sfid and 'parents' in exon_data:
                    parent_info = exon_data['parents']
                    break
                    
            if not parent_info:
                for sf_data in gene_dict.get('subfeatures', []):
                    if sf_data.get('sfid') == sfid and 'parents' in sf_data:
                        parent_info = sf_data['parents']
                        break
                        
            if not parent_info:
                for v_data in gene_dict.get('variants', []):
                    if v_data.get('sfid') == sfid and 'parents' in v_data:
                        parent_info = v_data['parents']
                        break
            
            if parent_info:
                # Set transcript parent
                if 'transcript' in parent_info and parent_info['transcript'] in self.transcripts:
                    feature.set_parent_transcript(self.transcripts[parent_info['transcript']])
                
                # Set exon parent
                if 'exon' in parent_info:
                    try:
                        exon_num = int(parent_info['exon'])
                        if exon_num in self.exons:
                            feature.set_parent_exon(self.exons[exon_num])
                    except ValueError:
                        pass
                
                # Set CDS parent
                if 'cds' in parent_info and parent_info['cds'] in feature_registry:
                    feature.set_parent_cds(feature_registry[parent_info['cds']])
        
        # Update counters
        self.num_transcripts = len(self.transcripts)
        self.num_exons = len(self.exons)
        self.num_proteins = len(self.proteins)
        self.num_cds = len(self.cds)
        self.num_subfeatures = len(self.subfeatures)
        self.num_variants = len(self.variants)
        
        # Re-establish any remaining hierarchy
        self._establish_hierarchy()
    
    def _rebuild_from_dict(self, gene_dict: Dict[str, Any]) -> None:
        """Rebuild all features and relationships from dictionary.
        
        Args:
            gene_dict: Dictionary containing complete gene data
        """
        # Clear existing collections
        self.transcripts.clear()
        self.exons.clear()
        self.proteins.clear()
        self.cds.clear()
        self.subfeatures.clear()
        self.variants.clear()
        self._all_features.clear()
        
        # Restore ids and counts if present
        if 'ids' in gene_dict:
            self.ids = gene_dict['ids']
        if 'counts' in gene_dict:
            self.counts = gene_dict['counts']
        
        # Store references to all features by sfid for linking
        feature_registry = {}
        
        # First pass: Create all feature objects using Feature.from_dict
        
        # Create transcripts
        if 'transcripts' in gene_dict:
            for transcript_name, transcript_data in gene_dict['transcripts'].items():
                transcript = Feature.from_dict(self._prepare_feature_dict(transcript_data, 'transcript'))
                transcript.transcript_name = transcript_name
                self.transcripts[transcript_name] = transcript
                feature_registry[transcript.sfid] = transcript
                transcript.set_parent_gene(self)
                self._all_features.append(transcript)
        
        # Create exons  
        if 'exons' in gene_dict:
            for exon_key, exon_data in gene_dict['exons'].items():
                exon = Feature.from_dict(self._prepare_feature_dict(exon_data, 'exon'))
                # Use hash as key for dictionary
                key = hash(exon)
                self.exons[key] = exon
                feature_registry[exon.sfid] = exon
                exon.set_parent_gene(self)
                self._all_features.append(exon)
        
        # Create proteins
        if 'proteins' in gene_dict:
            for protein_id, protein_data in gene_dict['proteins'].items():
                protein = Feature.from_dict(self._prepare_feature_dict(protein_data, 'protein'))
                self.proteins[protein_id] = protein
                feature_registry[protein.sfid] = protein
                protein.set_parent_gene(self)
                self._all_features.append(protein)
        
        # Create CDS features
        if 'cds' in gene_dict:
            for cds_data in gene_dict['cds']:
                cds = Feature.from_dict(self._prepare_feature_dict(cds_data, 'CDS'))
                self.cds.append(cds)
                feature_registry[cds.sfid] = cds
                cds.set_parent_gene(self)
                self._all_features.append(cds)
        
        # Create other subfeatures
        if 'subfeatures' in gene_dict:
            for subfeature_data in gene_dict['subfeatures']:
                feature_type = subfeature_data.get('type', subfeature_data.get('feature', 'unknown'))
                subfeature = Feature.from_dict(self._prepare_feature_dict(subfeature_data, feature_type))
                self.subfeatures.append(subfeature)
                feature_registry[subfeature.sfid] = subfeature
                subfeature.set_parent_gene(self)
                self._all_features.append(subfeature)
        
        # Create variants
        if 'variants' in gene_dict:
            for variant_data in gene_dict['variants']:
                variant = Feature.from_dict(self._prepare_feature_dict(variant_data, 'variant'))
                self.variants.append(variant)
                feature_registry[variant.sfid] = variant
                variant.set_parent_gene(self)
                self._all_features.append(variant)
        
        # Second pass: Rebuild relationships using the stored references
        
        # Process each dictionary that might have abbreviated references
        for dict_name in ['transcripts', 'exons', 'proteins']:
            if dict_name in gene_dict:
                for key, data in gene_dict[dict_name].items():
                    # Get the actual feature object
                    if dict_name == 'transcripts':
                        feature = self.transcripts.get(key)
                    elif dict_name == 'exons':
                        # Find exon by sfid
                        feature = None
                        for exon in self.exons.values():
                            if exon.sfid == data.get('sfid'):
                                feature = exon
                                break
                    elif dict_name == 'proteins':
                        feature = self.proteins.get(key)
                    
                    if not feature:
                        continue
                    
                    # Restore subfeature references
                    if 'subfeatures' in data:
                        for sf_ref in data['subfeatures']:
                            sf_id = self._extract_sfid(sf_ref)
                            if sf_id in feature_registry:
                                subfeature = feature_registry[sf_id]
                                if subfeature not in feature.subfeatures:
                                    feature.subfeatures.append(subfeature)
                                    subfeature.add_parent(feature)
        
        # Process lists that might have parent references
        for list_name in ['cds', 'subfeatures', 'variants']:
            if list_name in gene_dict:
                source_list = gene_dict[list_name]
                if list_name == 'cds':
                    target_list = self.cds
                elif list_name == 'subfeatures':
                    target_list = self.subfeatures
                else:
                    target_list = self.variants
                
                for i, data in enumerate(source_list):
                    if i < len(target_list):
                        feature = target_list[i]
                        
                        # Restore parent references
                        if 'parents' in data:
                            for parent_ref in data['parents']:
                                parent_id = self._extract_sfid(parent_ref)
                                if parent_id in feature_registry:
                                    parent = feature_registry[parent_id]
                                    if parent not in feature.parents:
                                        feature.add_parent(parent)
                                    # Also add to parent's subfeatures if not there
                                    if feature not in parent.subfeatures:
                                        parent.subfeatures.append(feature)
        
        # Third pass: Ensure all transcript-exon-CDS relationships are properly established
        for transcript in self.transcripts.values():
            # Find all exons that belong to this transcript
            for exon in self.exons.values():
                if hasattr(exon, 'transcript_name') and exon.transcript_name == transcript.transcript_name:
                    if exon not in transcript.subfeatures:
                        transcript.subfeatures.append(exon)
                    if transcript not in exon.parents:
                        exon.add_parent(transcript)
            
            # Find all CDS that belong to this transcript
            for cds in self.cds:
                if hasattr(cds, 'transcript_name') and cds.transcript_name == transcript.transcript_name:
                    # Check if CDS overlaps with any exon of this transcript
                    for exon in transcript.subfeatures:
                        if exon.type == 'exon' and cds.overlaps(exon):
                            if cds not in exon.subfeatures:
                                exon.subfeatures.append(cds)
                            if exon not in cds.parents:
                                cds.add_parent(exon)
        
        # Re-establish hierarchy for any remaining connections
        self._establish_hierarchy()
        
        # Sort all features
        self._sort_all()
    
    def _prepare_feature_dict(self, data: Dict[str, Any], feature_type: str) -> Dict[str, Any]:
        """Prepare feature dictionary for Feature constructor.
        
        Args:
            data: Raw feature data
            feature_type: Type of feature
            
        Returns:
            Dictionary ready for Feature constructor
        """
        feature_dict = data.copy()
        feature_dict['feature'] = feature_type
        
        # Propagate gene-level attributes if missing
        if not feature_dict.get('gene_name') and hasattr(self, 'gene_name'):
            feature_dict['gene_name'] = self.gene_name
        if not feature_dict.get('chrom') and hasattr(self, 'chrom'):
            feature_dict['chrom'] = self.chrom
        if not feature_dict.get('strand') and hasattr(self, 'strand'):
            feature_dict['strand'] = self.strand
            
        return feature_dict
    
    def _extract_sfid(self, sf_ref: str) -> str:
        """Extract sfid from a subfeature reference string.
        
        Args:
            sf_ref: Subfeature reference (e.g., "variant-1 (ts 0)" or "exon-1")
            
        Returns:
            The sfid part of the reference
        """
        if isinstance(sf_ref, str):
            # Handle cases like "variant-1 (ts 0)"
            if ' (' in sf_ref:
                return sf_ref.split(' (')[0]
            return sf_ref
        return str(sf_ref)
    
    def __iter__(self):
        """Custom iteration for genes: yield all unique features"""
        yield self
        
        # Yield all transcripts
        for transcript in self.transcripts.values():
            yield transcript
        
        # Yield all unique exons (even if shared)
        seen_exons = set()
        for transcript in self.transcripts.values():
            for exon in transcript.subfeatures:
                if exon.id not in seen_exons:
                    seen_exons.add(exon.id)
                    yield exon
                    
                    # Yield exon subfeatures (CDS, UTR)
                    for subfeature in exon.subfeatures:
                        yield subfeature
                        
