
from ggene.genomemanager import GenomeManager
from ggene.genemap import GeneMap
from ggene.features import variant_to_dict

import json
gm = GenomeManager()
gmap = gm.gene_map
vcf = gm.vcf
ref = gm.ref

min_len = 200

outpath_frm = "./data/{c}-var-mindelta{min_len}.json"

vardict = {}

max_var = -1
max_var_per_chrom={}

for c in gmap.chromes:
    vardict[c]=[]
    max_var_c =0
    chrome_len = ref.get_reference_length(c)
    for v in vcf(gm._make_index(c,0,chrome_len)):
        delta = len(v.ALT[0])-len(v.REF)
        if abs(delta)>min_len:
            vdict = variant_to_dict(v, gene_name="", chrom=c, strand='+')
            genes = gmap.fetch_all(c, vdict['start'],vdict['end'],features=['gene'])
            CDSs = gmap.fetch_all(c, vdict['start'],vdict['end'],features=['CDS'])
            gene_names = [g.get('gene_name',g.get('info',{}).get('gene_name','?')) for g in genes]
            CDSnames = [f"CDS({cds.get("chrom",'?')}:{cds.get('start',-1)}-{cds.get('end',-1)})" for cds in CDSs]
            vdict['genes'] = genes
            vdict['cds'] = CDSs
            vardict[c].append(vdict)
            max_var_c = max(max_var_c, abs(delta))
            print(delta, 
                v.var_type, 
                v.var_subtype, 
                ','.join([g.get('info',{}).get('gene_name',"") for g in genes]),
                ','.join([c['transcript_name'] for c in CDSs])
            )
    max_var = max(max_var_c, max_var)
    max_var_per_chrom[c] = max_var_c


    outpath=outpath_frm.format(c=c, min_len=min_len)
    with open(outpath,'w') as f:
        json.dump(vardict[c],f, indent = 2)
    
print(f"maximum variant size: {max_var}")
print(f"maximum variant size by chromosome: {max_var_per_chrom}")