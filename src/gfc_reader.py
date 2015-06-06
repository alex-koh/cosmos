import numpy as np
import os.path

def read_gfc(filepath):
    """
        Read params of gravity model from *.gfc. 
        Represent C and S parameters as complex number C - j * S 
        to simplify calculation.
    """
    assert os.path.exists(filepath)

    result = {}

    f = open(filepath, "rb")
    try:
        s = "_"
        while s != "":  s = f.readline().strip()
        while s == "":  s = f.readline().strip()
        while s != "":
            ss = s.split(None, 1)
            assert len(ss)==2, s
            if ss[0]=="earth_gravity_constant" or ss[0]=="radius":
                sss = ss[1].split(None,1)
                result[ss[0]] = np.float64(sss[0])
            elif ss[0] == "max_degree":
                result[ss[0]] = int(ss[1])
            else:
                result[ss[0]] = ss[1]
            s = f.readline().strip()
        while s == "":  s = f.readline().strip()
        s = f.readline().strip()
        ss = s.split()
        assert len(ss)==2 and "end_of_head" == ss[0] and ss[1].startswith("====="), s
        
        n = result["max_degree"]
        assert n > 0, result

        idx = np.arange(n+1)
        idx.cumsum(out=idx)

        data = np.ndarray(((n+1)*(n+2)/2, ), dtype=np.complex128)
        
        s = f.readline().strip()
        while s != "":
            ss = s.split(None, 6)
            assert len(ss) >= 5, s
            assert ss[0] == "gfc"
            n = int(ss[1])
            m = int(ss[2])
            data[idx[n] + m] = np.float64(ss[3])+1j*np.float64(ss[4])
            s = f.readline().strip()
        result["data"] = data
        return result
    finally:
        f.close()


