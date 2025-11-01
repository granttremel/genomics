
from typing import Optional, Dict, List, Tuple, Any

def align_sequences_with_gaps(seqa: str, seqb: str, features: List[Dict[str, Any]], window_start = 0):
    """Align reference and personal sequences by inserting gaps at indel positions.
    
    Returns:
        Tuple of (aligned_ref, aligned_pers, variant_positions)
    """
    # Find variant features to identify indels
    variant_features = [f for f in features if f.feature_type == 'variant']
    
    # Create a list of indel events in the window
    indels = []
    
    for var in variant_features:
        var_start = var.start
        ref = var.get('ref', '')
        alt = var.get('alt', '')
        
        # Calculate position in our window (0-based)
        window_pos = var_start - window_start
        
        # Only process if the variant starts within our window
        if 0 <= window_pos < len(seqa):
            if len(ref) != len(alt):
                indels.append({
                    'pos': window_pos,
                    'ref_len': len(ref),
                    'alt_len': len(alt),
                    'ref': ref,
                    'alt': alt
                })
    
    # Sort indels by position
    indels.sort(key=lambda x: x['pos'])
    
    # Build aligned sequences with gaps
    aligned_a = []
    aligned_b = []
    variant_positions = []  # Track where variants are for coloring
    
    idx_a = 0
    idx_b = 0
    
    for i in range(len(seqa)):
        # Check if we're at an indel position
        current_indel = None
        for indel in indels:
            if indel['pos'] == i:
                current_indel = indel
                break
        
        if current_indel:
            ref_len = current_indel['ref_len']
            alt_len = current_indel['alt_len']
            
            if alt_len > ref_len:
                # Insertion - add gaps to reference
                for j in range(ref_len):
                    aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                    aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                    variant_positions.append(True)
                    idx_a += 1
                    idx_b += 1
                
                # Add the inserted bases with gaps in reference
                for j in range(alt_len - ref_len):
                    aligned_a.append('-')
                    aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                    variant_positions.append(True)
                    idx_b += 1
                    
            elif ref_len > alt_len:
                # Deletion - add gaps to personal  
                for j in range(alt_len):
                    aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                    aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                    variant_positions.append(True)
                    idx_a += 1
                    idx_b += 1
                
                # Add the deleted bases with gaps in personal
                for j in range(ref_len - alt_len):
                    aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                    aligned_b.append('-')
                    variant_positions.append(True)
                    idx_a += 1
                    
            else:
                # SNP - just add both
                aligned_a.append(seqa[idx_a] if idx_a < len(seqa) else 'N')
                aligned_b.append(seqb[idx_b] if idx_b < len(seqb) else 'N')
                variant_positions.append(seqa[idx_a] != seqb[idx_b])
                idx_a += 1
                idx_b += 1
        else:
            # No indel at this position
            if idx_a < len(seqa) and idx_b < len(seqb):
                aligned_a.append(seqa[idx_a])
                aligned_b.append(seqb[idx_b])
                variant_positions.append(seqa[idx_a] != seqb[idx_b])
                idx_a += 1
                idx_b += 1
            elif idx_a < len(seqa):
                aligned_a.append(seqa[idx_a])
                aligned_b.append('-')
                variant_positions.append(True)
                idx_a += 1
            elif idx_b < len(seqb):
                aligned_a.append('-')
                aligned_b.append(seqb[idx_b])
                variant_positions.append(True)
                idx_b += 1
    
    return ''.join(aligned_a), ''.join(aligned_b), variant_positions


