
import numpy as np
from pathlib import Path
import os
from typing import Optional, Dict, List

from ggene.config import DEFAULT_CHRDATA_CACHE_PATH, DEFAULT_LIBRARY
from ggene.database.genome_manager import GenomeManager
from ggene.processing.chrome_mapper import ChromeMapper
from ggene.database.annotations import chr_lens
from ggene.seqs import lambdas


class CacheSampler:
    """
    Sampler interface for cached genomic quantities.

    Wraps cached numpy arrays and provides convenient sampling with automatic
    chromosome switching and resolution adjustment.

    Usage:
        cache = GenomeCache(gm)
        sampler = cache.get_sampler("gc_content")

        # Sample a region
        data = sampler.sample("1", start=1000000, length=100000, chunksz=1000)

        # Or use in an artist
        artist = MapArtist(sampler=sampler)
    """

    def __init__(
        self,
        cache: "GenomeCache",
        seq_spec: str,
        base_resolution: int = 4096
    ):
        """
        Initialize sampler for a specific quantity.

        Args:
            cache: GenomeCache instance
            seq_spec: The quantity specification (e.g., "gc_content")
            base_resolution: Resolution the cache was computed at
        """
        self.cache = cache
        self.seq_spec = seq_spec
        self.base_resolution = base_resolution

        # Lazy-loaded chromosome data
        self._chr_data: Dict[str, np.ndarray] = {}
        self._current_chrom: Optional[str] = None

    def _load_chrom(self, chrom: str) -> Optional[np.ndarray]:
        """Load chromosome data if not already cached."""
        chrom = str(chrom).removeprefix("chr")

        if chrom not in self._chr_data:
            data = self.cache.load_quantity_cache(chrom, self.seq_spec)
            if data is not None:
                self._chr_data[chrom] = data
            else:
                return None

        self._current_chrom = chrom
        return self._chr_data[chrom]

    def sample(
        self,
        chrom: str,
        start: int,
        length: int,
        chunksz: int,
        agg_func=None,
        num_chunks: int = None
    ) -> Optional[np.ndarray]:
        """
        Sample data for a genomic region.

        Args:
            chrom: Chromosome name
            start: Start position (1-based)
            length: Length of region in bp
            chunksz: Desired chunk size (bp per data point)
            agg_func: Aggregation function for resampling (default: np.mean)
            num_chunks: Exact number of output samples desired (if provided,
                       guarantees output length matches this value)

        Returns:
            Numpy array of sampled values, or None if no data
        """
        if agg_func is None:
            agg_func = np.mean

        chrom = str(chrom).removeprefix("chr")
        data = self._load_chrom(chrom)

        if data is None:
            return None

        # Convert genomic coordinates to array indices
        # Data was computed with base_resolution bp per element
        start_idx = max(0, (start - 1) // self.base_resolution)
        end_idx = min(len(data), (start - 1 + length) // self.base_resolution + 1)

        if start_idx >= end_idx:
            return None

        region_data = data[start_idx:end_idx]

        # Calculate target number of samples if not provided
        if num_chunks is None:
            num_chunks = int(length / chunksz)

        # Resample to exact target length
        if len(region_data) != num_chunks:
            region_data = self._resample(region_data, num_chunks, agg_func)

        return region_data

    def _resample(
        self,
        data: np.ndarray,
        target_length: int,
        agg_func
    ) -> np.ndarray:
        """
        Resample data to an exact target length.

        Args:
            data: Input array
            target_length: Exact number of output samples desired
            agg_func: Aggregation function for downsampling

        Returns:
            Resampled array with exactly target_length elements
        """
        if len(data) == target_length:
            return data

        if target_length <= 0:
            return np.array([agg_func(data)]) if len(data) > 0 else np.array([0])

        # Use interpolation to get exact target length
        # This works for both upsampling and downsampling
        x_old = np.linspace(0, 1, len(data))
        x_new = np.linspace(0, 1, target_length)
        return np.interp(x_new, x_old, data)

    def get_full_chrom(self, chrom: str) -> Optional[np.ndarray]:
        """Get full chromosome data without resampling."""
        return self._load_chrom(chrom)

    def available_chroms(self) -> List[str]:
        """List chromosomes that have cached data."""
        chroms = []
        for chrom in chr_lens.keys():
            fname = self.cache.get_cache_fname(chrom, self.seq_spec)
            if fname.exists():
                chroms.append(chrom)
        return chroms

    def preload(self, chroms: Optional[List[str]] = None):
        """
        Preload chromosome data into memory.

        Args:
            chroms: List of chromosomes to load (default: all available)
        """
        if chroms is None:
            chroms = self.available_chroms()

        for chrom in chroms:
            self._load_chrom(chrom)

    def clear_cache(self):
        """Clear loaded chromosome data from memory."""
        self._chr_data.clear()
        self._current_chrom = None

    def __repr__(self):
        loaded = list(self._chr_data.keys())
        return f"CacheSampler('{self.seq_spec}', base_res={self.base_resolution}, loaded={loaded})"


class GenomeCache:
    
    def __init__(self, gm:GenomeManager = None, base_resolution = 4096, dtype = 'float32'):
        self.gm = gm
        self.cache_dir:Path = DEFAULT_CHRDATA_CACHE_PATH 
        self.base_resolution = base_resolution
        self.dtype=  dtype
    
    def bind_genome_manager(self, genome_manager:GenomeManager):
        self.gm = genome_manager
    
    def get_cache_fname(self, chrom, seq_spec):
        return self.cache_dir / f"chr{chrom}" / f"chr{chrom}_{seq_spec}.npz"
    
    def cache_all_chromes(self, seq_specs,overwrite = False):
        
        if not self.gm:
            print("genome manager not set! cannot compute chromosomal quantities")
            return
        
        chromes =[chrom for chrom in  self.gm.iter_chromes()]
        self.cache_chromes(chromes, seq_specs, overwrite=overwrite)
    
    def cache_chromes(self, chromes, seq_specs, overwrite = False):
        
        if not self.gm:
            print("genome manager not set! cannot compute chromosomal quantities")
            return
        
        for chrom in chromes:
            self.cache_quantities(chrom, seq_specs, overwrite=overwrite)
    
    def cache_quantities(self, chrom, seq_specs, overwrite = False):
        
        if not self.gm:
            print("genome manager not set! cannot compute chromosomal quantities")
            return
        
        datas = self.precompute_quantities(chrom, seq_specs)
        
        for sp, data in zip(seq_specs, datas):
            print(f"caching quantity {sp} at chr{chrom}")
            self.write_quantity_cache(data, chrom, sp, overwrite=overwrite)
    
    def cache_quantity(self, chrom, seq_spec, overwrite = False):
        
        if not self.gm:
            print("genome manager not set! cannot compute chromosomal quantities")
            return
        
        data = self.precompute_quantities(chrom, [seq_spec],)[0]
        self.write_quantity_cache(data, chrom, seq_spec, overwrite=overwrite)
    
    def precompute_quantities(self, chrom, seq_specs):
        
        if not self.gm:
            print("genome manager not set! cannot compute chromosomal quantities")
            return
        
        chr_len = chr_lens.get(chrom, -1)
        num_chunks = int(chr_len / self.base_resolution) + 1
        chunksz = int(chr_len / num_chunks)
        
        sg, fg = ChromeMapper.get_generators(self.gm, None)
        cm = ChromeMapper(sg, fg)
        
        needs = []
        for seq_spec in seq_specs:
            needs.extend(lambdas.needs_features(seq_spec))
        
        needs = list(set(needs))
        
        print(f"beginning precompute with chr{chrom}, seq_specs {", ".join(seq_specs)}")
        
        qts, _ = cm.get_chromosomal_quantities(chrom, seq_specs, chunksz = chunksz, start = 1, needs_feats = needs)
        
        return [np.array(qts[i], dtype = self.dtype) for i in range(len(seq_specs))]
        

    def write_quantity_cache(self, qt:np.ndarray, chrom, seq_spec, overwrite = False):
        
        fname = self.get_cache_fname(chrom, seq_spec)
        
        if not fname.parent.exists():
            fname.parent.mkdir(exist_ok = True)
        
        if fname.exists() and not overwrite:
            print(f"file {fname} for cache with chr{chrom} and seq_spec {seq_spec} could not be overwritten!")
            return
        
        if np.isclose(np.mean(qt),0.0):
            print("quantity is null, skipping write")
            return
        
        print(f"writing array with shape {qt.shape}, mean {np.mean(qt)}, range {np.min(qt)}-{np.max(qt)} to file {fname}")
        
        with open(fname, "wb") as f:
            np.save(f, qt, allow_pickle = False)
        
        sz = os.path.getsize(fname)
        print(f"cache written with size {sz}b")
        
    def load_quantity_cache(self, chrom, seq_spec):

        fname = self.get_cache_fname(chrom, seq_spec)

        if not fname.exists():
            return

        with open(fname, "rb") as f:
            data = np.load(f)

        return data

    def get_sampler(self, seq_spec: str, base_resolution: int = 4096) -> CacheSampler:
        """
        Get a sampler for a cached quantity.

        Args:
            seq_spec: The quantity specification (e.g., "gc_content")
            base_resolution: Resolution the cache was computed at

        Returns:
            CacheSampler instance for convenient data access

        Example:
            sampler = cache.get_sampler("gc_content")
            data = sampler.sample("chr1", start=1000000, length=100000, chunksz=1000)
        """
        return CacheSampler(self, seq_spec, base_resolution)

    def list_cached_specs(self, chrom: str = "1") -> List[str]:
        """
        List available cached seq_specs for a chromosome.

        Args:
            chrom: Chromosome to check (default: "1")

        Returns:
            List of seq_spec names that have cached data
        """
        chrom = str(chrom).removeprefix("chr")
        chrom_dir = self.cache_dir / f"chr{chrom}"

        if not chrom_dir.exists():
            return []

        specs = []
        for f in chrom_dir.glob(f"chr{chrom}_*.npz"):
            # Extract seq_spec from filename: chr1_gc_content.npy -> gc_content
            spec = f.stem.replace(f"chr{chrom}_", "")
            specs.append(spec)

        return specs


