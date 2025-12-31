"""JASPAR format I/O for Position Weight Matrices.

This module handles loading PWMs from JASPAR format files.
The core PWM functionality is in pwm.py - this module just handles I/O.

JASPAR format:
    >MA0001.1 AHR::ARNT
    A [ 3  0  0  0  0  0 ]
    C [ 8  0  23 10  0  0 ]
    G [ 2  23 0  0  23 0 ]
    T [ 10 0  0  13 0  23 ]

Classes:
    Jaspar: PWM loaded from JASPAR format (extends PWM)
    JasparLibrary: Collection loaded from directory (extends PWMLibrary)
"""

from typing import Dict, List, Optional
from tabulate import tabulate
from pathlib import Path
import numpy as np

from ggene.motifs.pwm import PWM, PWMLibrary


class Jaspar(PWM):
    """PWM loaded from JASPAR format file.

    Extends PWM with JASPAR-specific I/O and display methods.
    All scanning functionality is inherited from PWM.

    Attributes:
        path: Path to source .jaspar file
    """

    def __init__(self, name: str, counts: np.ndarray,
                 motif_id: str = "", path: str = "",
                 pseudocount: float = 0.1,
                 background: Dict[str, float] = None):
        """Initialize Jaspar PWM.

        Args:
            name: Motif name (e.g., "AHR::ARNT")
            counts: 4xN count matrix, rows are [A, C, G, T]
            motif_id: JASPAR ID (e.g., "MA0001.1")
            path: Path to source file
            pseudocount: Pseudocount for probability calculation
            background: Background frequencies
        """
        super().__init__(name, counts, pseudocount=pseudocount,
                        background=background, motif_id=motif_id)
        self.path = path

    @classmethod
    def from_file(cls, file_path: str, pseudocount: float = 0.1) -> 'Jaspar':
        """Load Jaspar PWM from a .jaspar format file.

        Args:
            file_path: Path to .jaspar file
            pseudocount: Pseudocount for probability calculation

        Returns:
            Jaspar instance

        Raises:
            ValueError: If file format is invalid
        """
        bases_order = []
        all_weights = None
        header = None

        with open(file_path, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line:
                    continue

                # First non-empty line is header
                if header is None:
                    header = line
                    continue

                # Parse weight line: "A [ 3  0  0  0  0  0 ]"
                parts = line.split()
                base = parts[0].strip()
                bases_order.append(base)

                # Remove brackets and parse weights
                weight_str = line.split('[')[1].split(']')[0]
                weights = np.array([int(v) for v in weight_str.split()])

                if all_weights is None:
                    all_weights = np.zeros((4, len(weights)))
                all_weights[len(bases_order) - 1] = weights

        if all_weights is None:
            raise ValueError(f"No weight data found in {file_path}")

        # Reorder rows to standard [A, C, G, T] order
        standard_order = ['A', 'C', 'G', 'T']
        base_map = [bases_order.index(b) for b in standard_order]
        counts = all_weights[base_map]

        # Parse header: ">MA0001.1 AHR::ARNT"
        header = header.lstrip('>')
        parts = header.split(maxsplit=1)
        motif_id = parts[0].strip()
        name = parts[1].strip() if len(parts) > 1 else motif_id

        return cls(name, counts, motif_id=motif_id, path=file_path,
                  pseudocount=pseudocount)

    @classmethod
    def from_path(cls, dir_path: str, max_motifs: Optional[int] = None,
                  pseudocount: float = 0.1) -> 'JasparLibrary':
        """Load all .jaspar files from a directory.

        Args:
            dir_path: Path to directory containing .jaspar files
            max_motifs: Maximum number of motifs to load (None for all)
            pseudocount: Pseudocount for probability calculation

        Returns:
            JasparLibrary containing all loaded PWMs
        """
        dir_path = Path(dir_path).absolute()
        library = JasparLibrary()

        for jasp_file in sorted(dir_path.iterdir()):
            if jasp_file.suffix != ".jaspar":
                continue

            try:
                jaspar = cls.from_file(str(jasp_file), pseudocount=pseudocount)
                library.add_jaspar(jaspar)
            except Exception as e:
                print(f"Error loading {jasp_file}: {e}")
                continue

            if max_motifs and len(library) >= max_motifs:
                break

        return library

    def to_consensus(self, min_p_ratio: float = 2.0) -> str:
        """Get consensus sequence (alias for inherited consensus property)."""
        return self._compute_consensus(min_p_ratio)

    def print(self):
        """Print formatted PWM information."""
        from ggene.seqs import bio

        tab = [
            ["Name", self.name],
            ["ID", self.id],
            ["Num pos.", self.length],
            ["Consensus", self.consensus],
            ["Total IC", f"{self.total_ic:.2f} bits"],
            ["Path", self.path or "N/A"],
        ]
        print(tabulate(tab))

        # Print weight matrix
        headers = ["Base"] + [f"P{i+1}" for i in range(self.length)]
        matrix_tab = []
        for bi, base in enumerate(self.bases):
            row = [base] + list(self.probs[bi, :])
            matrix_tab.append(row)

        print(tabulate(matrix_tab, headers=headers, floatfmt=".3f"))


class JasparLibrary(PWMLibrary):
    """Collection of Jaspar PWMs with JASPAR-specific I/O.

    Extends PWMLibrary with methods for loading JASPAR format files
    and display. All scanning functionality is inherited from PWMLibrary.
    """

    def add_jaspar(self, jaspar: Jaspar) -> None:
        """Add a Jaspar PWM to the library.

        Args:
            jaspar: Jaspar instance to add
        """
        self.add(jaspar)

    @classmethod
    def from_path(cls, dir_path: str, max_motifs: Optional[int] = None,
                  pseudocount: float = 0.1) -> 'JasparLibrary':
        """Load all .jaspar files from a directory.

        Args:
            dir_path: Path to directory containing .jaspar files
            max_motifs: Maximum number of motifs to load
            pseudocount: Pseudocount for probability calculation

        Returns:
            JasparLibrary instance
        """
        return Jaspar.from_path(dir_path, max_motifs, pseudocount)

    @property
    def num_jaspars(self) -> int:
        """Number of Jaspar PWMs (alias for num_motifs)."""
        return self.num_motifs

    @property
    def jaspars(self) -> Dict[str, Jaspar]:
        """Dict of Jaspar PWMs (alias for internal storage)."""
        return self._pwms

    @property
    def jaspar_ids(self) -> List[str]:
        """List of Jaspar IDs (alias for motif_ids)."""
        return self._pwm_ids

    def print(self):
        """Print formatted library summary."""
        tab = []

        for pwm_id in self._pwm_ids:
            pwm = self._pwms[pwm_id]
            tab.append([
                pwm.id,
                pwm.name,
                pwm.length,
                f"{pwm.entropy:.3f}",
                pwm.consensus,
                getattr(pwm, 'path', 'N/A'),
            ])

        print(tabulate(tab, headers=[
            "ID", "Name", "Len", "Entropy", "Consensus", "Path"
        ]))

    def __repr__(self) -> str:
        return f"JasparLibrary({len(self._pwms)} motifs)"
