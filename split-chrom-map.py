


import json


with open('./data/chromosome_map.json','r') as f:
    
    data = json.load(f)
    
for c in data:
    
    outfile = f'./data/chrom{c}_map.json'
    
    with open(outfile,'w') as ff:
        json.dump(data[c],ff)


