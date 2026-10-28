"""GMP integer-only acceleration. BCOV equations remain in direct_jet.py."""
from pathlib import Path
import ctypes as C
from fractions import Fraction as Q
import bcov_yy as alg

def encode(m):return sum(v*8**i for i,v in enumerate(m))
def decode(code):return tuple((code//8**i)%8 for i in range(5))

class NativeSource:
    def __init__(self):
        self.lib=C.CDLL(str(Path(__file__).with_name('integer_source.so')))
        self.lib.add_vertex.argtypes=[C.c_int,C.c_int,C.POINTER(C.c_int),C.POINTER(C.c_char_p),C.c_char_p]
        self.lib.compute_source.argtypes=[C.c_int];self.lib.compute_source.restype=C.c_void_p
        self.lib.free_result.argtypes=[C.c_void_p]
        self.lib.reset_vertices();self.installed=set()
    def put(self,g,p):
        ints,den=p.integer_form();items=list(ints.items())
        codes=(C.c_int*len(items))(*(encode(m) for m,_ in items))
        numbers=(C.c_char_p*len(items))(*(str(v).encode() for _,v in items))
        self.lib.add_vertex(g,len(items),codes,numbers,str(den).encode())
    def __call__(self,g,vertices,previous_two):
        for h in range(1,g):
            if h not in self.installed:self.put(h,vertices[h]);self.installed.add(h)
        self.put(-1,previous_two)
        ptr=self.lib.compute_source(g)
        try:lines=C.string_at(ptr).decode().splitlines()
        finally:self.lib.free_result(ptr)
        assert lines[0]=='REDUCED';terms={}
        for line in lines[1:]:
            code,num,den=map(int,line.split())
            # GMP has already certified coprimality and positive denominator.
            # Older Python versions use the ordinary safe constructor.
            if num:terms[decode(code)]=Q._from_coprime_ints(num,den) if hasattr(Q,'_from_coprime_ints') else Q(num,den)
        return alg.Poly(terms)
