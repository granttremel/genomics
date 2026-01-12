
from typing import Dict, Any

from ggene.database.uobject.uobject import UObject

class UInfo(UObject):
    _core_fields = frozenset({
        'name', "description", 'id'
    })

    _defaults = {
        'name': "",
        "description":"",
        'id': "",
    }

    __slots__ = ('name', "description", 'id', 'attributes', '_canonical')

