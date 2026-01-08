
from ggene.genomemanager import GenomeManager

gm = GenomeManager()
gmap = gm.gene_map
vcf = gm.vcf

# print(dir(vcf))
# print(vcf.seqlens)

# print(dir(gmap.tabix))
# print(dir(gm.ref))
# print(gm.ref.lengths)
# print(gm.ref.references)

for chrom in gmap.chromes:
    chromelen = gm.ref.lengths[gm.ref.references.index(chrom)]
    print(chrom, chromelen)





 