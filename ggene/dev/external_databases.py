"""Integration with external genomic databases.

Provides easy access to JASPAR, ENCODE, Rfam, RepeatMasker, and other
public genomic annotation databases.
"""

import requests
import gzip
import json
from pathlib import Path
from typing import Dict, List, Optional, Iterator
import logging
from dataclasses import dataclass
import numpy as np

from ggene.database.unified_stream import UFeature, AnnotationStream

logger = logging.getLogger(__name__)


class DatabaseDownloader:
    """Download and cache external database files."""
    
    def __init__(self, cache_dir: str = "~/.ggene/cache"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
    def download_file(self, url: str, filename: str, force: bool = False) -> Path:
        """Download file if not cached."""
        filepath = self.cache_dir / filename
        
        if filepath.exists() and not force:
            logger.info(f"Using cached file: {filepath}")
            return filepath
            
        logger.info(f"Downloading {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                
        logger.info(f"Downloaded to {filepath}")
        return filepath


@dataclass
class JASPARMotif:
    """JASPAR motif with PWM."""
    matrix_id: str
    name: str
    tf_class: str
    species: str
    pwm: np.ndarray  # Position Weight Matrix
    consensus: str
    
    def scan(self, sequence: str, threshold: float = 0.8) -> List[tuple]:
        """Scan sequence for motif matches."""
        from motifs.pwm import PWM
        
        pwm_obj = PWM(self.pwm)
        return pwm_obj.scan_sequence(sequence, threshold)


class JASPARDatabase:
    """Interface to JASPAR transcription factor binding motifs."""
    
    BASE_URL = "https://jaspar.genereg.net/api/v1"
    
    def __init__(self, cache_dir: str = "~/.ggene/cache/jaspar"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.motifs = {}
        
    def download_collection(self, collection: str = "CORE", species: str = "Homo sapiens"):
        """Download a JASPAR collection."""
        # Download matrix list
        url = f"{self.BASE_URL}/matrix/?collection={collection}&species={species}&format=json"
        response = requests.get(url)
        matrices = response.json()
        
        logger.info(f"Found {len(matrices)} motifs in {collection} collection")
        
        # Download each matrix
        for matrix_info in matrices:
            matrix_id = matrix_info['matrix_id']
            self.download_motif(matrix_id)
            
    def download_motif(self, matrix_id: str) -> JASPARMotif:
        """Download a specific motif."""
        cache_file = self.cache_dir / f"{matrix_id}.json"
        
        if cache_file.exists():
            with open(cache_file, 'r') as f:
                data = json.load(f)
        else:
            # Download motif data
            url = f"{self.BASE_URL}/matrix/{matrix_id}/"
            response = requests.get(url)
            data = response.json()
            
            # Cache it
            with open(cache_file, 'w') as f:
                json.dump(data, f)
        
        # Parse PWM
        pfm = data['pfm']  # Position Frequency Matrix
        pwm = np.array([
            pfm['A'],
            pfm['C'],
            pfm['G'],
            pfm['T']
        ])
        
        motif = JASPARMotif(
            matrix_id=matrix_id,
            name=data['name'],
            tf_class=data.get('class', ''),
            species=data.get('species', ''),
            pwm=pwm,
            consensus=data.get('sequence', '')
        )
        
        self.motifs[matrix_id] = motif
        return motif
    
    def scan_sequence(self, sequence: str, threshold: float = 0.8) -> List[UFeature]:
        """Scan sequence with all loaded motifs."""
        features = []
        
        for matrix_id, motif in self.motifs.items():
            hits = motif.scan(sequence, threshold)
            for start, score, matched_seq, _ in hits:
                features.append(UFeature(
                    chrom='',  # Will be set by caller
                    start=start,
                    end=start + len(matched_seq) - 1,
                    feature_type='tf_binding',
                    source='JASPAR',
                    score=score,
                    name=motif.name,
                    id=matrix_id,
                    attributes={
                        'tf_class': motif.tf_class,
                        'consensus': motif.consensus,
                        'matched_seq': matched_seq
                    }
                ))
        
        return features


class ENCODEDatabase:
    """Interface to ENCODE project data."""
    
    BASE_URL = "https://www.encodeproject.org"
    
    def __init__(self, cache_dir: str = "~/.ggene/cache/encode"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.downloader = DatabaseDownloader(self.cache_dir)
        
    def search_experiments(self, biosample: str = "K562", 
                          assay: str = "ChIP-seq",
                          target: Optional[str] = None) -> List[Dict]:
        """Search for ENCODE experiments."""
        params = {
            'type': 'Experiment',
            'biosample_ontology.term_name': biosample,
            'assay_title': assay,
            'format': 'json',
            'limit': 'all'
        }
        
        if target:
            params['target.label'] = target
            
        response = requests.get(f"{self.BASE_URL}/search/", params=params)
        data = response.json()
        
        return data['@graph']
    
    def download_peaks(self, experiment_id: str, assembly: str = "GRCh38") -> Path:
        """Download peak calls from an experiment."""
        # Get experiment metadata
        url = f"{self.BASE_URL}/experiments/{experiment_id}/?format=json"
        response = requests.get(url)
        experiment = response.json()
        
        # Find peak file for the right assembly
        for file_info in experiment.get('files', []):
            if (file_info.get('assembly') == assembly and
                file_info.get('output_type') == 'peaks' and
                file_info.get('file_format') == 'bed'):
                
                file_url = f"{self.BASE_URL}{file_info['href']}"
                filename = f"{experiment_id}_{assembly}_peaks.bed.gz"
                return self.downloader.download_file(file_url, filename)
        
        raise ValueError(f"No peak file found for {experiment_id} on {assembly}")
    
    def load_chip_peaks(self, target: str, cell_line: str = "K562") -> List[UFeature]:
        """Load ChIP-seq peaks for a target."""
        experiments = self.search_experiments(
            biosample=cell_line,
            assay="ChIP-seq",
            target=target
        )
        
        features = []
        for exp in experiments[:1]:  # Just use first experiment
            try:
                peak_file = self.download_peaks(exp['accession'])
                # Parse BED file
                with gzip.open(peak_file, 'rt') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        parts = line.strip().split('\t')
                        features.append(UFeature(
                            chrom=parts[0],
                            start=int(parts[1]) + 1,
                            end=int(parts[2]),
                            feature_type='chip_peak',
                            source='ENCODE',
                            score=float(parts[4]) if len(parts) > 4 else None,
                            name=f"{target}_peak",
                            attributes={
                                'target': target,
                                'cell_line': cell_line,
                                'experiment': exp['accession']
                            }
                        ))
            except Exception as e:
                logger.warning(f"Failed to load peaks from {exp['accession']}: {e}")
        
        return features


class RfamDatabase:
    """Interface to Rfam RNA families database."""
    
    BASE_URL = "https://rfam.org"
    
    def __init__(self, cache_dir: str = "~/.ggene/cache/rfam"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.downloader = DatabaseDownloader(self.cache_dir)
        
    def download_families(self):
        """Download Rfam family list."""
        url = f"{self.BASE_URL}/families?format=json"
        response = requests.get(url)
        return response.json()
    
    def download_cm(self, rfam_acc: str) -> Path:
        """Download covariance model for an RNA family."""
        url = f"{self.BASE_URL}/family/{rfam_acc}/cm"
        filename = f"{rfam_acc}.cm"
        return self.downloader.download_file(url, filename)
    
    def search_sequence(self, sequence: str, rfam_acc: str) -> List[UFeature]:
        """Search sequence for RNA family matches using Infernal."""
        # This would use Infernal's cmsearch
        # For now, return empty list
        return []


class RepeatMaskerDatabase:
    """Interface to RepeatMasker repeat annotations."""
    
    UCSC_URL = "https://hgdownload.soe.ucsc.edu/goldenPath"
    
    def __init__(self, cache_dir: str = "~/.ggene/cache/repeatmasker"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.downloader = DatabaseDownloader(self.cache_dir)
        
    def download_annotations(self, assembly: str = "hg38") -> Path:
        """Download RepeatMasker annotations for an assembly."""
        url = f"{self.UCSC_URL}/{assembly}/database/rmsk.txt.gz"
        filename = f"{assembly}_repeatmasker.txt.gz"
        return self.downloader.download_file(url, filename)
    
    def parse_annotations(self, filepath: Path) -> Iterator[UFeature]:
        """Parse RepeatMasker annotations."""
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 17:
                    continue
                    
                yield UFeature(
                    chrom=parts[5],
                    start=int(parts[6]) + 1,
                    end=int(parts[7]),
                    feature_type='repeat',
                    source='RepeatMasker',
                    strand=parts[9],
                    name=parts[10],
                    attributes={
                        'repeat_class': parts[11],
                        'repeat_family': parts[12],
                        'sw_score': int(parts[1]),
                        'percent_div': float(parts[2]),
                        'percent_del': float(parts[3]),
                        'percent_ins': float(parts[4])
                    }
                )


class ClinVarDatabase:
    """Interface to ClinVar clinical variant database."""
    
    FTP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
    
    def __init__(self, cache_dir: str = "~/.ggene/cache/clinvar"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.downloader = DatabaseDownloader(self.cache_dir)
        
    def download_vcf(self, assembly: str = "GRCh38") -> Path:
        """Download ClinVar VCF."""
        url = f"{self.FTP_URL}/vcf_{assembly}/clinvar.vcf.gz"
        filename = f"clinvar_{assembly}.vcf.gz"
        return self.downloader.download_file(url, filename)
    
    def get_pathogenic_variants(self, chrom: str, start: int, end: int) -> List[UFeature]:
        """Get pathogenic variants in a region."""
        # Would parse VCF and filter for pathogenic
        return []


class dbSNPDatabase:
    """Interface to dbSNP variant database."""
    
    FTP_URL = "https://ftp.ncbi.nih.gov/snp"
    
    def __init__(self, cache_dir: str = "~/.ggene/cache/dbsnp"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.downloader = DatabaseDownloader(self.cache_dir)
        
    def download_common_variants(self, assembly: str = "GRCh38") -> Path:
        """Download common variants (MAF > 0.01)."""
        # Simplified URL - real one is more complex
        url = f"{self.FTP_URL}/organisms/human_{assembly}/VCF/common_all.vcf.gz"
        filename = f"dbsnp_common_{assembly}.vcf.gz"
        return self.downloader.download_file(url, filename)


class IntegratedDatabaseManager:
    """Manage all external databases in one place."""
    
    def __init__(self, cache_dir: str = "~/.ggene/cache"):
        self.cache_dir = Path(cache_dir).expanduser()
        self.jaspar = JASPARDatabase(self.cache_dir / "jaspar")
        self.encode = ENCODEDatabase(self.cache_dir / "encode")
        self.rfam = RfamDatabase(self.cache_dir / "rfam")
        self.repeatmasker = RepeatMaskerDatabase(self.cache_dir / "repeatmasker")
        self.clinvar = ClinVarDatabase(self.cache_dir / "clinvar")
        self.dbsnp = dbSNPDatabase(self.cache_dir / "dbsnp")
        
        logger.info("Initialized integrated database manager")
        
    def initialize_core_databases(self, assembly: str = "hg38"):
        """Download core annotation databases."""
        logger.info(f"Initializing core databases for {assembly}...")
        
        # Download key datasets
        tasks = [
            ("JASPAR motifs", lambda: self.jaspar.download_collection("CORE")),
            ("RepeatMasker", lambda: self.repeatmasker.download_annotations(assembly)),
            ("ClinVar", lambda: self.clinvar.download_vcf("GRCh38" if assembly == "hg38" else assembly)),
        ]
        
        for name, task in tasks:
            try:
                logger.info(f"Downloading {name}...")
                task()
            except Exception as e:
                logger.warning(f"Failed to download {name}: {e}")
    
    def annotate_region(self, chrom: str, start: int, end: int, 
                       sequence: Optional[str] = None) -> Dict[str, List[UFeature]]:
        """Get all annotations for a region."""
        annotations = {}
        
        # TF binding sites (requires sequence)
        if sequence:
            annotations['tf_binding'] = self.jaspar.scan_sequence(sequence)
            # Adjust positions
            for feature in annotations['tf_binding']:
                feature.chrom = chrom
                feature.start += start - 1
                feature.end += start - 1
        
        # Other annotations would be added here
        # annotations['repeats'] = self.repeatmasker.query_region(chrom, start, end)
        # annotations['clinical'] = self.clinvar.get_pathogenic_variants(chrom, start, end)
        
        return annotations


# Example usage
def example_external_databases():
    """Demonstrate external database integration."""
    
    # Initialize manager
    db_manager = IntegratedDatabaseManager()
    
    # Download some motifs
    print("Downloading JASPAR motifs...")
    db_manager.jaspar.download_motif("MA0139.1")  # CTCF
    db_manager.jaspar.download_motif("MA0099.3")  # AP1
    
    # Scan a sequence
    test_sequence = "CCCTCCCAGGTCAGCTGCCCTCCCGGGATCAGAT"
    print(f"\nScanning sequence: {test_sequence}")
    
    hits = db_manager.jaspar.scan_sequence(test_sequence, threshold=0.8)
    for hit in hits:
        print(f"  Found {hit.name} at position {hit.start}-{hit.end} (score: {hit.score:.2f})")
    
    # Search ENCODE experiments
    print("\nSearching ENCODE for CTCF ChIP-seq...")
    experiments = db_manager.encode.search_experiments(
        biosample="K562",
        assay="ChIP-seq", 
        target="CTCF"
    )
    print(f"Found {len(experiments)} CTCF ChIP-seq experiments")
    
    # Get all annotations for a region
    print("\nGetting all annotations for chr1:1000000-1001000...")
    annotations = db_manager.annotate_region("chr1", 1000000, 1001000, "A" * 1000)
    
    for source, features in annotations.items():
        print(f"  {source}: {len(features)} features")


if __name__ == "__main__":
    example_external_databases()