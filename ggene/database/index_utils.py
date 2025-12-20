"""Utilities for creating and managing tabix indexes for genomic files."""

import os
import subprocess
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def ensure_indexed(filepath: str, file_type: str = "auto") -> bool:
    """Ensure a genomic file is bgzipped and tabix-indexed.
    
    Args:
        filepath: Path to the genomic file
        file_type: Type of file (vcf, bed, gff, gtf, or auto to detect)
        
    Returns:
        True if file is indexed or indexing succeeded
    """
    filepath = Path(filepath)
    
    # Detect file type if auto
    if file_type == "auto":
        if "vcf" in filepath.suffix.lower():
            file_type = "vcf"
        elif "bed" in filepath.suffix.lower():
            file_type = "bed"
        elif "gff" in filepath.suffix.lower() or "gtf" in filepath.suffix.lower():
            file_type = "gff"
        else:
            logger.warning(f"Cannot auto-detect file type for {filepath}")
            return False
    
    # Check if already indexed
    if filepath.suffix == '.gz':
        index_files = [
            Path(str(filepath) + '.tbi'),
            Path(str(filepath) + '.csi')
        ]
        if any(idx.exists() for idx in index_files):
            logger.info(f"{filepath} is already indexed")
            return True
        
        # Try to create index
        try:
            logger.info(f"Creating tabix index for {filepath}...")
            cmd = ['tabix', '-p', file_type, str(filepath)]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info(f"Successfully indexed {filepath}")
                return True
            else:
                logger.error(f"Failed to index {filepath}: {result.stderr}")
                return False
                
        except FileNotFoundError:
            logger.error("tabix not found. Install with: conda install -c bioconda tabix")
            return False
            
    else:
        # Need to bgzip first
        bgz_file = Path(str(filepath) + '.gz')
        
        if bgz_file.exists():
            # Already bgzipped, recurse to check index
            return ensure_indexed(bgz_file, file_type)
            
        try:
            logger.info(f"Compressing {filepath} with bgzip...")
            cmd = ['bgzip', '-c', str(filepath)]
            
            with open(bgz_file, 'wb') as out:
                result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
                
            if result.returncode == 0:
                logger.info(f"Created {bgz_file}")
                # Now index it
                return ensure_indexed(bgz_file, file_type)
            else:
                logger.error(f"Failed to bgzip {filepath}: {result.stderr}")
                return False
                
        except FileNotFoundError:
            logger.error("bgzip not found. Install with: conda install -c bioconda htslib")
            return False


def check_all_indexes(file_list: list) -> dict:
    """Check indexing status of multiple files.
    
    Args:
        file_list: List of file paths
        
    Returns:
        Dictionary mapping files to their index status
    """
    status = {}
    
    for filepath in file_list:
        filepath = Path(filepath)
        if not filepath.exists():
            status[str(filepath)] = "missing"
            continue
            
        if filepath.suffix == '.gz':
            index_files = [
                Path(str(filepath) + '.tbi'),
                Path(str(filepath) + '.csi')
            ]
            if any(idx.exists() for idx in index_files):
                status[str(filepath)] = "indexed"
            else:
                status[str(filepath)] = "not_indexed"
        else:
            bgz_file = Path(str(filepath) + '.gz')
            if bgz_file.exists():
                status[str(filepath)] = "bgzipped"
                # Check if bgzipped version is indexed
                index_files = [
                    Path(str(bgz_file) + '.tbi'),
                    Path(str(bgz_file) + '.csi')
                ]
                if any(idx.exists() for idx in index_files):
                    status[str(filepath)] = "bgzipped_indexed"
            else:
                status[str(filepath)] = "not_compressed"
                
    return status


def create_all_indexes(file_list: list) -> bool:
    """Create indexes for all files in list.
    
    Args:
        file_list: List of file paths
        
    Returns:
        True if all files were successfully indexed
    """
    success = True
    
    for filepath in file_list:
        if not ensure_indexed(filepath):
            success = False
            
    return success


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) > 1:
        for filepath in sys.argv[1:]:
            print(f"\nChecking {filepath}...")
            if ensure_indexed(filepath):
                print(f"✓ {filepath} is indexed")
            else:
                print(f"✗ Failed to index {filepath}")
    else:
        print("Usage: python index_utils.py <file1> [file2] ...")
        print("\nChecking common genomic files...")
        
        common_files = [
            "./data/GRCh38_sorted.gtf.gz",
            "./genome/vcf/HG00096.vcf.gz"
        ]
        
        status = check_all_indexes(common_files)
        for filepath, stat in status.items():
            print(f"  {filepath}: {stat}")