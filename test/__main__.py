from unittest.loader import TestLoader
from unittest.runner import TextTestRunner

import os.path
import sys

def test_path(path):
    if not os.path.exists(path): raise Exception("cannot find path "+path)

bdir = os.path.dirname(os.path.abspath(__file__))

sys.__proj__ = os.path.dirname(bdir)
sys.__resources__ = os.path.join(sys.__proj__,"resources")

test_path(sys.__proj__)
test_path(sys.__resources__)

del test_path

sys.path.append(os.path.join(sys.__proj__,"src"))

loader = TestLoader()
tests = loader.discover(bdir)
#print tests
runner = TextTestRunner()
runner.run(tests)
