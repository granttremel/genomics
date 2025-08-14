
import json

from ggene.genomemanager import GenomeManager



geneman = GenomeManager()

chromdict = {}
for chrom in list(range(1,23)) + ['X','Y','MT']:
    chromdict[chrom]={}

    
    print(f"starting chrome {chrom}")
    i=0
    mystery=0
    for g in geneman.gene_map.stream(chrom):
        gatts = g['attributes']
        
        if not 'gene_name' in gatts:
            gene_name = f'unknown-gene{mystery}'
            mystery+=1
        else:
            gene_name = gatts['gene_name']
        
        gdict=dict(
            name=gene_name,
            start=g['start'],
            end=g['end'],
            strand=g['strand'],
            id = gatts['gene_id'],
            biotype=gatts['gene_biotype']
        )      
        chromdict[chrom][gene_name]=gdict
        
        if i%100==0:
            print(f'completed gene {i}, named {gene_name} on chrom {chrom}')
        
        i+=1

outfile = './data/chromosome_map.json'
with open(outfile,'w') as f:
    json.dump(chromdict,f)
            

