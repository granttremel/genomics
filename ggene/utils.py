
import inspect as insp
from enum import Enum, StrEnum
import numpy as np
import json
import os
import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Union

import logging
logger = logging.getLogger(__name__)
 
def try_parse_number(num: Union[str, int, float]) -> Union[str, int, float]:
    """Attempt to parse a string as a number.
    
    Args:
        num: Input value to parse
        
    Returns:
        Parsed number (int or float) or original string if parsing fails
    """
    if not isinstance(num, str):
        return num
    
    # Check if it's a valid number format
    _num = num.replace('.', '')
    if not _num.isdecimal():
        return num
    
    # Try integer first
    try:
        return int(num)
    except ValueError:
        pass
        
    # Try float
    try:
        return float(num)
    except ValueError:
        pass
        
    return num


class locals:
    
    def __init__(self, fname = 'locals',fdir = './analysis', autowrite = True, overwrite = False):
        self.vars = {}
        self.tracking = {}
        self.meta = {}
        
        self.fpath = os.path.join(fdir,fname + '.json')
        self.overwrite = overwrite
        if autowrite:
            import atexit
            atexit.register(self.on_exit)

    def track(self, name, obj, regtype = ''):
        
        self.add(name, obj, regtype = regtype)
        
        return self.get(name)
    
    def get(self, name):
        return self.tracking.get(name)
    
    def add(self, name, obj, regtype = ''):
        
        self.tracking[name] = obj
        if not regtype:
            regtype = "object"
            if type(obj) is np.ndarray:
                regtype = "data"
        self.meta[name] = register_type(regtype)
    
    def convert(self, name):
        
        if not name in self.tracking:
            return
        
        obj = self.tracking[name]
        meta = self.meta[name]
        if meta == 'data':
            _obj = stat.from_data(obj).to_dict()
        else:
            if hasattr(obj, 'to_dict'):
                _obj = obj.to_dict()
                _obj = make_serializable(_obj)
            elif hasattr(obj,'__dict__'):
                _obj = dict(obj)
                _obj = make_serializable(_obj)
            else:
                try:
                    _ = json.dumps(obj)
                    _obj = obj
                except:
                    print(f'converting object {name}:{str(obj)} to dict the hard way!')
                    _obj = object_to_dict(obj)
        
        self.vars[name] = _obj
        del self.tracking[name]
        
    
    def convert_all(self):
        names = list(self.tracking.keys())
        for name in names:
            self.convert(name)
    
    def write(self, fname):
        
        fdir = os.getcwd()
        fpath = os.path.join(fdir, fname)
        try:
            write_object(self.vars, fpath, overwrite = self.overwrite)
        except Exception as e:
            print(e)
            return False
        return True
    
    @staticmethod
    def load():
        
        
        pass
    
    def on_exit(self):
        self.convert_all()
        res = self.write(self.fpath)

class register_type(StrEnum):
    OBJECT = "object"
    DATA = "data"

@dataclass
class stat:
    mean: float=0
    sd: float=0
    median: float=0
    min: float=0
    max: float=0
    
    @staticmethod
    def from_data(data):
        return stat(mean = np.mean(data),
                    sd = np.std(data),
                    median = np.median(data),
                    min = np.min(data),
                    max = np.max(data),
                    )
    
    def __str__(self):
        strs = []
        for f in ["mean","sd","median","min","max"]:
            v = getattr(self, f)
            abv = np.abs(v)
            fstr = '0.3f'
            if abv < 1e-4 or abv > 1e4:
                fstr = '0.3e'
            fieldstr = format(v, fstr)
            strs.append(fieldstr)
        return '\n'.join(strs)
    
    def to_dict(self):
        outdict = {}
        for f in ["mean","sd","median","min","max"]:
            v = getattr(self, f)
            abv = np.abs(v)
            fstr = '0.3f'
            if abv < 1e-4 or abv > 1e4:
                fstr = '0.3e'
            fieldstr = format(v, fstr)
            outdict[f] = fieldstr 
        return outdict

@dataclass
class dictstatus:
    seen: List[int]=field(default_factory=list)
    depth: float=0
    maxdepth: Optional[int]=5
    flag:str=''



def inspect(obj, unders = False, callables = False, deep = False, do_print = True):
    outstrs = []
    outstrs.append(f"type: {type(obj).__name__}")
    
    # atts = dir(obj)
    avs = insp.getmembers(obj)
    
    boringtypes = ['list','tuple','dict',
                   'str', 'int','float','bool','NoneType',
                   'ndarray']
    methods = []
    
    for a,vv in avs:
        if a[0] == '_' and not unders:
            continue
        # try:
        #     v = getattr(obj, a)
        # except Exception as e:
        #     v = str(e)
        v = vv
        
        if callable(v):
            methods.append(v.__name__)
        else:
            vstr = str(v)
            vtyp = type(v).__name__
        
        if '<' in vstr and '>' in vstr and deep:
            vstr = inspect(v, deep = False, do_print = False)
            vstr = vstr.replace('\n','\n\t')
        
        if len(vstr) > 200:
            vstr = vstr[:200] + f'...({len(vstr)} chars)'
        
        if vtyp in boringtypes:
            outstrs.append(f"{a}: {vstr}")
        else:
            outstrs.append(f"{a}: {vstr} (<{vtyp}>)")
            
        # outstrs.append(f"{a}: {vstr} (<{vtyp}>)")
    
    if callables:
        methlist = ', '.join(methods)
        methstr = f'Methods: [{methlist}]'
        outstrs.insert(1, methstr)
    
    outstr = '\n'.join(outstrs)
    
    if do_print:
        print(outstr)
        
    return outstr

types = set(('bool','str','int','float','ndarray','NoneType'))
tostr = set(('ndarray',))
clcs = set(('dict','list','set','tuple','OrderedDict'))

nonser = [set]

def make_serializable(objdict):
    
    newdict = dict()
    
    for k in objdict:
        v = objdict[k]
        
        if isinstance(v, str):
            newv = v
        elif isinstance(v,dict):
            newv = make_serializable(v)
        elif type(v) in nonser:
            newv = list(v)
        elif type(v).__name__ in tostr:
            newv = str(v)
        else:
            newv = v
            
        newdict[k] = newv
    
    return newdict


def object_to_dict(obj, status = None):
    
    if not status:
        status = dictstatus()
    
    objtypename = type(obj).__name__
    if objtypename in types:
        if objtypename in tostr:
            return str(obj)
        else:
            return obj
    
    objid = id(obj)
    if objid in status.seen:
        status.flag = 'circle'
        return objtypename
    
    status.depth += 1
    if status.depth > status.maxdepth:
        print(f'depth too deep! ({status.depth}) abort!')
        status.depth -= 1
        return objtypename
    
    if status.flag:
        if status.flag == 'circle':
            status.flag = ''
            return objtypename
    
    status.seen.append(id(obj))
    
    # if hasattr(obj, '__getitem__'):
    if type(obj).__name__ in clcs:
        if isinstance(obj, dict):
            outobj = dict()
            for k in obj:
                outobj[k] = object_to_dict(obj[k], status = status)
        else:
            outobj = list()
            for item in obj:
                newitem = object_to_dict(item, status = status)
                outobj.append(newitem)
        
        status.depth -= 1
        return outobj
    
    outdict = dict(__name__ = objtypename)
    mems = insp.getmembers(obj)
    
    for a, v in mems:
        if a[0] == '_':
            continue
        
        if callable(v):
            continue

        outdict[a] = object_to_dict(v, status=status)
    
    
    status.depth -= 1
    return outdict


def get_next_file(fpath):
    
    fn, fext = os.path.splitext(fpath)
    fdir, fname = os.path.split(fn)
    res = re.match(f'(.*)-(\d+)$', fname)
    
    if not res:
        nn = 1
    else:
        fname = res.groups()[0]
        n = res.groups()[1]
        nn = int(n) + 1
    
    while nn < 99999:
        
        fname_new = f'{fname}-{str(nn)}{fext}'
        fpath_new = os.path.join(fdir, fname_new)
        if not os.path.exists(fpath_new):
            return fpath_new
        
        nn = nn + 1
    
    return ''

def write_object(obj, fpath, overwrite = False):
    
    if not fpath.endswith('.json'):
        fpath = fpath + '.json'
    
    if os.path.exists(fpath):
        # print(f'file {fpath} exists')
        if not overwrite:
            fpnew = get_next_file(fpath)
            if not fpnew:
                return
            fpath = fpnew

        else:
            print('overwriting...')
    
    with open(fpath,'w') as f:
        json.dump(obj, f, indent = 3)
        
    if os.path.exists(fpath):
        return True
    return False
