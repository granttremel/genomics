


from ggene.database.annotations import UFeature, AnnotationStream, TabularStream, ColumnSpec, DerivedSpec
from ggene.database.sequences import FASTAStream


class NoncodeV6Stream(FASTAStream):
    
    def __init__(self, filepath):
        super().__init__(filepath)
        
    



