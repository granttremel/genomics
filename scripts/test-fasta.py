import pysam

# Open the reference
ref = pysam.FastaFile("./data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")

# 1. Check what chromosomes/contigs are present
print("Chromosomes in reference:")
print(f"Number of sequences: {ref.nreferences}")
print(f"Main chromosomes: {ref.references[:25]}")  # First 25
print(f"Lengths: {ref.lengths[:25]}")

# 2. Verify chromosome naming (important - some use 'chr1', others just '1')
if "chr1" in ref.references:
    chr_prefix = "chr"
else:
    chr_prefix = ""
print(f"\nChromosome naming style: '{chr_prefix}1' format")

# 3. Check a known sequence (like mitochondrial genome start)
mt_name = f"{chr_prefix}MT" if f"{chr_prefix}MT" in ref.references else f"{chr_prefix}M"
if mt_name in ref.references:
    mt_start = ref.fetch(mt_name, 0, 20)
    print(f"\nMitochondrial genome start: {mt_start}")
    # Should be something like "GATCACAGGTCTATCACCCT"

# 4. Check total genome size
total_size = sum(ref.lengths)
print(f"\nTotal genome size: {total_size:,} bp")
# Should be ~3.1 billion for human

# 5. Look for common features
# Check centromeric region (lots of Ns)
centro_region = ref.fetch(f"{chr_prefix}1", 121500000, 121500100)
print(f"\nCentromeric region: {centro_region}")
print(f"N content: {centro_region.count('N')}%")

# 6. Spot check some genes you know
# Beta globin region (classic genetics region)
hbb_region = ref.fetch(f"{chr_prefix}11", 5246000, 5246100)
print(f"\nBeta globin region: {hbb_region[:50]}...")

# 7. Check for lowercase (soft-masked repeats)
sample = ref.fetch(f"{chr_prefix}1", 10000, 11000)
lowercase_count = sum(1 for c in sample if c.islower())
print(f"\nSoft-masking: {lowercase_count/len(sample)*100:.1f}% lowercase in sample")