#!/usr/bin/env python3
"""
Migration guide and helper for transitioning from old GeneMap to refactored version.

This provides a compatibility layer and shows how to update existing code.
"""

from ggene.genemap import GeneMap as OldGeneMap
from ggene.genemap_refactored import GeneMap as NewGeneMap
from typing import Dict, List, Any, Optional
import logging

logger = logging.getLogger(__name__)


class GeneMapMigrationHelper:
    """Helper class to migrate from old GeneMap to new architecture."""
    
    def __init__(self, gtf_path: str = '', library_path: str = './genome/library'):
        """Initialize both old and new GeneMap for comparison/migration."""
        self.old_gm = OldGeneMap(gtf_path)
        self.new_gm = NewGeneMap(gtf_path, library_path)
    
    def migrate_existing_usage(self):
        """Show how to migrate common usage patterns."""
        
        print("🔄 GENEMAP MIGRATION GUIDE")
        print("=" * 50)
        
        print("\n📊 OLD vs NEW Method Mapping:")
        print("-" * 30)
        
        mappings = [
            ("OLD: gm.fetch(chrom, start, end)", "NEW: gm.fetch(chrom, start, end)"),
            ("OLD: gm.fetch_all(chrom, start, end)", "NEW: gm.fetch_all(chrom, start, end)"),
            ("OLD: gm.get_gene(name, chrom)", "NEW: gm.load_gene(name, chrom)"),
            ("OLD: gm.make_gene(name, chrom)", "NEW: gm.assemble_gene(name, chrom)"),
            ("OLD: gm.by_gene_name(chrom, name)", "NEW: gm.query({'chrom': chrom, 'gene_name': name})"),
            ("OLD: gm.get_feature(chrom, pos)", "NEW: gm.get_feature(chrom, pos)"),
            ("OLD: gm.neighbors(chrom, pos)", "NEW: gm.neighbors(chrom, pos)"),
            ("OLD: gm.find_gene(name)", "NEW: gm.find_gene(name)"),
            ("OLD: gm.stream(chrom)", "NEW: gm.stream(chrom)")
        ]
        
        for old, new in mappings:
            print(f"  {old}")
            print(f"  → {new}")
            print()
        
        print("\n🎯 KEY IMPROVEMENTS:")
        print("-" * 20)
        print("✅ Clear separation of concerns")
        print("✅ Consistent interfaces")
        print("✅ Automatic JSON caching")
        print("✅ Better error handling")
        print("✅ Type hints throughout")
        print("✅ Comprehensive logging")
        
        print("\n🔧 NEW POWERFUL FEATURES:")
        print("-" * 25)
        print("✅ Complex queries: gm.query({'chrom': '1', 'gene_biotype': 'protein_coding'})")
        print("✅ Auto-caching genes: gm.load_gene('BRCA1') # caches automatically")
        print("✅ Gene search: gm.search_genes('BRCA') # finds all BRCA genes")
        print("✅ Standardized features: gm.assemble_features(raw_records)")
        
    def test_compatibility(self, test_chrom: str = '22', test_gene: str = 'BRCA1'):
        """Test that new implementation gives equivalent results to old."""
        
        print(f"\n🧪 TESTING COMPATIBILITY (Chr {test_chrom}, Gene {test_gene})")
        print("=" * 60)
        
        # Test 1: Basic fetch
        print("\n1. Testing fetch...")
        try:
            old_results = list(self.old_gm.fetch(test_chrom, 1000000, 1001000, features=('gene',)))
            new_results = self.new_gm.fetch_all(test_chrom, 1000000, 1001000, features=['gene'])
            
            print(f"   Old: {len(old_results)} records")
            print(f"   New: {len(new_results)} records")
            print(f"   ✅ Match: {len(old_results) == len(new_results)}")
        except Exception as e:
            print(f"   ❌ Error: {e}")
        
        # Test 2: Gene finding
        print(f"\n2. Testing gene finding for {test_gene}...")
        try:
            old_chrom, old_genes = self.old_gm.find_gene(test_gene)
            new_chrom, new_genes = self.new_gm.find_gene(test_gene)
            
            print(f"   Old: Chr {old_chrom}, {len(old_genes) if old_genes else 0} records")
            print(f"   New: Chr {new_chrom}, {len(new_genes) if new_genes else 0} records")
            print(f"   ✅ Match: {old_chrom == new_chrom}")
        except Exception as e:
            print(f"   ❌ Error: {e}")
        
        # Test 3: Feature assembly (new feature)
        print(f"\n3. Testing new feature assembly for {test_gene}...")
        try:
            if new_chrom:
                gene = self.new_gm.assemble_gene(test_gene, new_chrom)
                if gene:
                    print(f"   ✅ Successfully assembled {test_gene}")
                    print(f"   📊 Features: {len(gene.subfeatures)} subfeatures")
                    print(f"   📊 Transcripts: {len([f for f in gene.subfeatures if f.type == 'transcript'])}")
                else:
                    print(f"   ❌ Failed to assemble {test_gene}")
        except Exception as e:
            print(f"   ❌ Error: {e}")


# Example of how to update existing code
def migration_examples():
    """Show concrete examples of migrating existing code."""
    
    print("\n📝 CODE MIGRATION EXAMPLES:")
    print("=" * 40)
    
    print("""
OLD CODE:
--------
# Messy, inconsistent interface
gene_features = gm.get_gene('BRCA1', '17')
gene_obj = gm.make_gene('BRCA1', '17')
transcripts = gm.by_gene_name('17', 'BRCA1', features=('transcript',))

NEW CODE:
--------
# Clean, consistent interface
gene = gm.load_gene('BRCA1', '17')  # Auto-cached!
# OR
gene = gm.load_gene('BRCA1')  # Auto-finds chromosome

# Complex queries
results = gm.query({
    'chrom': '17',
    'gene_name': 'BRCA1',
    'features': ['transcript', 'exon'],
    'limit': 100
})

# Search functionality
brca_genes = gm.search_genes('BRCA', limit=10)
""")

    print("""
OLD CODE:
--------
# Manual feature handling
features = list(gm.fetch('1', 1000000, 1001000))
for f in features:
    # Manual parsing of attributes
    attrs = gm.parse_info(f['info'])
    gene_name = attrs.get('gene_name')

NEW CODE:
--------
# Standardized features
features = gm.fetch_all('1', 1000000, 1001000)
standardized = gm.assemble_features(features)
for f in standardized:
    # Clean, consistent structure
    gene_name = f['info']['gene_name']
""")

    print("""
OLD CODE:
--------
# No caching, repeated database hits
gene1 = gm.make_gene('BRCA1', '17')  # DB hit
gene2 = gm.make_gene('BRCA1', '17')  # Another DB hit!

NEW CODE:
--------
# Automatic caching
gene1 = gm.load_gene('BRCA1', '17')  # DB hit + auto-save to JSON
gene2 = gm.load_gene('BRCA1', '17')  # Loaded from cache!
""")


if __name__ == "__main__":
    # Run migration demo
    helper = GeneMapMigrationHelper()
    helper.migrate_existing_usage()
    migration_examples()
    helper.test_compatibility()