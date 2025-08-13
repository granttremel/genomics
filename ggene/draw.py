
import matplotlib.pyplot as plt
import matplotlib.patches as patches



def draw_gene_structure(gene_features, gene_name, strand):
    """Draw a simple gene structure cartoon"""
    fig, ax = plt.subplots(figsize=(12, 3))
    
    # Sort features by position
    features = sorted(gene_features, key=lambda x: x['start'])
    
    # Gene baseline
    gene_start = features[0]['start']
    gene_end = features[-1]['end']
    ax.plot([gene_start, gene_end], [0.5, 0.5], 'k-', linewidth=1)
    
    # Draw features
    for feat in features:
        if feat['type'] == 'exon':
            # Determine if CDS or UTR
            height = 0.3 if 'CDS' in feat.get('extra_types', []) else 0.15
            color = 'darkblue' if 'CDS' in feat.get('extra_types', []) else 'lightblue'
            
            rect = patches.Rectangle(
                (feat['start'], 0.5 - height/2), 
                feat['end'] - feat['start'], 
                height,
                facecolor=color,
                edgecolor='black'
            )
            ax.add_patch(rect)
    
    # Add arrow for strand
    arrow_y = 0.5
    if strand == '+':
        ax.arrow(gene_end, arrow_y, 1000, 0, head_width=0.05, head_length=500, fc='red')
    else:
        ax.arrow(gene_start, arrow_y, -1000, 0, head_width=0.05, head_length=500, fc='red')
    
    # Formatting
    ax.set_ylim(0, 1)
    ax.set_xlim(gene_start - 5000, gene_end + 5000)
    ax.set_title(f"{gene_name} ({strand} strand)")
    ax.set_xlabel("Genomic Position")
    ax.get_yaxis().set_visible(False)
    
    plt.tight_layout()
    return fig

