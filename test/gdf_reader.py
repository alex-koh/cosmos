import os.path
import numpy as np


def read_gdf(filepath):
    """
        Read and check *.gdf. Represent it as a key-value map.
    """
    assert os.path.exists(filepath), filepath
    result = {}
    floats = set([
        "gmrefpot",
        "radiusrefpot",
        "flatrefpot",
        "omegarefpot",
        "normal_potential",

        "weighted_mean",
        "maxvalue",
        "minvalue",
        "signal_wrms"
    ])
    ints = set([
        "number_of_gridpoints",
        "longitude_parallels",
        "latitude_parallels"
    ])
    f = open(filepath, "rb")
    try:
        s = f.readline().strip()
        while s != "":
            ss = s.split(None,1)
            assert len(ss)==2, s
            ss[1] == ss[1].strip()
            if ss[0] in ints:
                result[ss[0]] = int(ss[1])
            elif ss[0] == "height_over_ell":
                sss = ss[1].split(None, 1)
                assert len(sss) == 2 and sss[1].strip() == "m", s
                result[ss[0]] = float(sss[0])
            elif ss[0] in floats:
                sss = ss[1].split(None, 1)
                result[ss[0]] = float(sss[0])
                assert not ss[0]+"_desc" in result
                result[ss[0]+"_desc"] = sss[1]
            else:
                result[ss[0]] = ss[1]
            s = f.readline().strip()

        s = f.readline()
        s = f.readline()
        s = f.readline()
        ss = s.split()
        assert len(ss)==2 and "end_of_head" == ss[0] and ss[1].startswith("====="), s
        
        size = result["number_of_gridpoints"]
        lon_size = result["longitude_parallels"]
        lat_size = result["latitude_parallels"] 
        assert size == lat_size * lon_size, str(result)
        assert result["long_lat_unit"] == "degree", str(result)
        assert float(result["gridstep"]) == 1., str(result)
        
        data = np.ndarray((lon_size*lat_size, 3), dtype=np.float64)
        
        s, i = f.readline().strip(), 0
        while s != "":
            ss = s.split()
            assert len(ss)==3, s
            data[i] = [np.float64(ss[0]), np.float64(ss[1]), np.float64(ss[2])]
            s, i = f.readline().strip(), i + 1

        result["data"] = data
    finally:
        f.close()
    data = result["data"]
    def cmp(a,b, e):
        return 2*np.abs((a-b)/(a+b)) < e
    assert cmp(data[:,2].min(), result["minvalue"], 1e-7), "%e %e" % (data.min(), result["minvalue"])
    assert cmp(data[:,2].max(), result["maxvalue"], 1e-7), "%e %e" % (data.max(), result["maxvalue"])
    return result
