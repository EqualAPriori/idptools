import os
path = os.path.dirname(os.path.realpath(__file__)) #__thismodule__.__path__ will give same thing
datadir = "/".join(path.split("/")[:-1]) + "/data/"
#sys.path.insert(1, modulepath)

