import sys, time; sys.path.insert(0,'pybuild')
import symfn
from sage.all import SymmetricFunctions, QQ

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("the Sage arm")
Sym = SymmetricFunctions(QQ); s = Sym.schur()
# Distinct shapes, each computed exactly once -> no cache reuse either side.
shapes=[[3,2,1],[4,2,1],[4,3,1],[4,3,2],[5,3,2],[5,4,2],[4,3,2,1],[5,3,2,1],[5,4,2,1],[5,4,3,2,1]]
print("single large products, each computed once (no cache reuse):")
for sh in shapes:
    st=time.perf_counter(); s(s[sh]*s[sh]).monomial_coefficients(); ts=time.perf_counter()-st
    st=time.perf_counter(); symfn.schur_multiply([(sh,1)],[(sh,1)]); tr=time.perf_counter()-st
    print("  s%-14s Sage %8.4fs | symfn %8.4fs -> %5.2fx"%(str(sh),ts,tr,ts/tr))
