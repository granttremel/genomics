"""
Fixed _rebuild_from_dict method for Gene class.
This will be integrated into features.py
"""

def _rebuild_from_dict(self, gene_dict: Dict[str, Any]) -> None:
    """Rebuild all features and relationships from dictionary.
    
    Args:
        gene_dict: Dictionary containing complete gene data
    """
    # Clear existing collections
    self.transcripts.clear()
    self.exons.clear()
    self.proteins.clear()
    self.cds.clear()
    self.subfeatures.clear()
    self.variants.clear()
    self._all_features.clear()
    
    # Restore ids and counts if present
    if 'ids' in gene_dict:
        self.ids = gene_dict['ids']
    if 'counts' in gene_dict:
        self.counts = gene_dict['counts']
    
    # Store references to all features by sfid for linking
    feature_registry = {}
    
    # First pass: Create all feature objects using Feature.from_dict
    
    # Create transcripts
    if 'transcripts' in gene_dict:
        for transcript_name, transcript_data in gene_dict['transcripts'].items():
            transcript = Feature.from_dict(self._prepare_feature_dict(transcript_data, 'transcript'))
            transcript.transcript_name = transcript_name
            self.transcripts[transcript_name] = transcript
            feature_registry[transcript.sfid] = transcript
            transcript.set_parent_gene(self)
            self._all_features.append(transcript)
    
    # Create exons  
    if 'exons' in gene_dict:
        for exon_key, exon_data in gene_dict['exons'].items():
            exon = Feature.from_dict(self._prepare_feature_dict(exon_data, 'exon'))
            # Use hash as key for dictionary
            key = hash(exon)
            self.exons[key] = exon
            feature_registry[exon.sfid] = exon
            exon.set_parent_gene(self)
            self._all_features.append(exon)
    
    # Create proteins
    if 'proteins' in gene_dict:
        for protein_id, protein_data in gene_dict['proteins'].items():
            protein = Feature.from_dict(self._prepare_feature_dict(protein_data, 'protein'))
            self.proteins[protein_id] = protein
            feature_registry[protein.sfid] = protein
            protein.set_parent_gene(self)
            self._all_features.append(protein)
    
    # Create CDS features
    if 'cds' in gene_dict:
        for cds_data in gene_dict['cds']:
            cds = Feature.from_dict(self._prepare_feature_dict(cds_data, 'CDS'))
            self.cds.append(cds)
            feature_registry[cds.sfid] = cds
            cds.set_parent_gene(self)
            self._all_features.append(cds)
    
    # Create other subfeatures
    if 'subfeatures' in gene_dict:
        for subfeature_data in gene_dict['subfeatures']:
            feature_type = subfeature_data.get('type', subfeature_data.get('feature', 'unknown'))
            subfeature = Feature.from_dict(self._prepare_feature_dict(subfeature_data, feature_type))
            self.subfeatures.append(subfeature)
            feature_registry[subfeature.sfid] = subfeature
            subfeature.set_parent_gene(self)
            self._all_features.append(subfeature)
    
    # Create variants
    if 'variants' in gene_dict:
        for variant_data in gene_dict['variants']:
            variant = Feature.from_dict(self._prepare_feature_dict(variant_data, 'variant'))
            self.variants.append(variant)
            feature_registry[variant.sfid] = variant
            variant.set_parent_gene(self)
            self._all_features.append(variant)
    
    # Second pass: Rebuild relationships using the stored references
    
    # Process each dictionary that might have abbreviated references
    for dict_name in ['transcripts', 'exons', 'proteins']:
        if dict_name in gene_dict:
            for key, data in gene_dict[dict_name].items():
                # Get the actual feature object
                if dict_name == 'transcripts':
                    feature = self.transcripts.get(key)
                elif dict_name == 'exons':
                    # Find exon by sfid
                    feature = None
                    for exon in self.exons.values():
                        if exon.sfid == data.get('sfid'):
                            feature = exon
                            break
                elif dict_name == 'proteins':
                    feature = self.proteins.get(key)
                
                if not feature:
                    continue
                
                # Restore subfeature references
                if 'subfeatures' in data:
                    for sf_ref in data['subfeatures']:
                        sf_id = self._extract_sfid(sf_ref)
                        if sf_id in feature_registry:
                            subfeature = feature_registry[sf_id]
                            if subfeature not in feature.subfeatures:
                                feature.subfeatures.append(subfeature)
                                subfeature.add_parent(feature)
    
    # Process lists that might have parent references
    for list_name in ['cds', 'subfeatures', 'variants']:
        if list_name in gene_dict:
            source_list = gene_dict[list_name]
            if list_name == 'cds':
                target_list = self.cds
            elif list_name == 'subfeatures':
                target_list = self.subfeatures
            else:
                target_list = self.variants
            
            for i, data in enumerate(source_list):
                if i < len(target_list):
                    feature = target_list[i]
                    
                    # Restore parent references
                    if 'parents' in data:
                        for parent_ref in data['parents']:
                            parent_id = self._extract_sfid(parent_ref)
                            if parent_id in feature_registry:
                                parent = feature_registry[parent_id]
                                if parent not in feature.parents:
                                    feature.add_parent(parent)
                                # Also add to parent's subfeatures if not there
                                if feature not in parent.subfeatures:
                                    parent.subfeatures.append(feature)
    
    # Third pass: Ensure all transcript-exon-CDS relationships are properly established
    for transcript in self.transcripts.values():
        # Find all exons that belong to this transcript
        for exon in self.exons.values():
            if hasattr(exon, 'transcript_name') and exon.transcript_name == transcript.transcript_name:
                if exon not in transcript.subfeatures:
                    transcript.subfeatures.append(exon)
                if transcript not in exon.parents:
                    exon.add_parent(transcript)
        
        # Find all CDS that belong to this transcript
        for cds in self.cds:
            if hasattr(cds, 'transcript_name') and cds.transcript_name == transcript.transcript_name:
                # Check if CDS overlaps with any exon of this transcript
                for exon in transcript.subfeatures:
                    if exon.type == 'exon' and cds.overlaps(exon):
                        if cds not in exon.subfeatures:
                            exon.subfeatures.append(cds)
                        if exon not in cds.parents:
                            cds.add_parent(exon)
    
    # Re-establish hierarchy for any remaining connections
    self._establish_hierarchy()
    
    # Sort all features
    self._sort_all()