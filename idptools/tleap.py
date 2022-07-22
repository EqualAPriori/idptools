# Wrapper for AmberTools
# See tutorials:
#   https://ambermd.org/tutorials/basic/tutorial7/index.php
#   https://ambermd.org/tutorials/basic/tutorial3/section1.htm
#
import os

import utility.log as log
from . import config
from . import aa

templatedir = config.datadir + "/templates/"


def add_line(x,y):
    return x + "\n" + y

def setup(sequence=None,settings=None):
    """
    Notes:
        Default behavior: sequence overrides whatever is in settings
    """
    if settings == None:
        print("using default tleap settings")
        settings = log.json_load(templatedir+"tleap_defaults.json")

    script = ""

    # load ff
    if "ff" not in settings:
        raise ValueError("missing ff specification")
    else:
        ff = validate_ff(settings["ff"])
        script += "source leaprc.{}".format(ff)
    if "water" not in settings:
        raise ValueError("missing water specification")
    else:
        water = validate_water(settings["water"])
        if water.startswith("mbond"):
            script += "\nset default PBRadii {}".format(water)
        else:
            script += "\nsource leaprc.water.{}".format(water)

    # parse sequence
    if sequence is None:
        if "seq" in settings:
            sequence = settings["seq"]
        else:
            if "fullseq" in settings:
                sequence = settings["fullseq"]
            else:
                sequence = ""

    if sequence.upper() in ["","__SEQ__"]:
        writeseq = "__SEQ__"
        settings["seq"] = ""
        settings["fullseq"] = ""
    else:
        writeseq = " ".join(aa.expand(sequence))
        settings["seq"] = aa.abbreviate(writeseq)
        settings["fullseq"] = writeseq

    # save
    if "name" not in settings:
        print("using default molecule name `test`")
        molname = "test"
        settings["name"] = "test"
    molname = settings["name"]

    script += "\n"
    script += "\n{0} = sequence {{ {1} }}".format(molname,writeseq)
    script += "\nsaveoff {0} {0}.lib".format(molname)
    script += "\nsavepdb {0} {0}.pdb".format(molname)
    script += "\nsaveamberparm {0} {0}.prmtop {0}.rst7".format(molname)
    script += "\nsaveamberparm {0} {0}.parm7 {0}.rst7".format(molname)

    script += "\n\nquit"

    with open("tleap.in","w") as f:
        f.write(script)

    log.json_write("settings.json",settings)

def run():
    os.system("tleap -f tleap.in")

def minimize(name):
    """
    name: prefix of topology and forcefield files
    """
    os.system("sander -O -i {0}/min1.in -o min1.out -p {1}.prmtop -c {1}.rst7 -r min1.ncrst".format(templatedir,name))

    os.system("ambpdb -p {}.prmtop -c min1.ncrst > min1.pdb".format(name))

# ===== VALIDATION
def validate_ff(ff):
    if ff.lower() in ["99sb","f99sb","ff99sb"]:
        return "ff99SB"
    elif ff.lower() in ["protein.ff14sbonlysc","ff14sbonlysc"]:
        return "protein.ff14SBonlysc"
    elif ff.lower().startswith("ff"):
        return "ff"+ff[2:].upper()
    else:
        print("Unknown ff setting {}, using as is".format(ff))
        return ff
def validate_water(ff):
    if ff.lower() in ["tip3p"]:
        return "tip3p"
    elif ff.lower() in ["opc"]:
        return "opc"
    elif ff.lower().startswith("igb"):
        igb = int(ff[3]) #the igb number, see https://github.com/openmm/openmm/issues/3708#issuecomment-1192126842 
        igb_pdbradii_dict = {1:"mbondi",2:"mbondi2",5:"mbondi2",7:"mbondi",8:"mbondi3"}
        return igb_pdbradii_dict[igb]
    else:
        print("Unknown water ff {}, using as is".format(ff))
        return ff

# ===== MAIN
if __name__ == "__main__":
    cmdln()

def cmdln(unknown=None):
    """
    Note:
        after installing, call this script via "pyleap"
        -q "<AAVPGXG" -ff ff99SB -w tip3p -n myPeptide -r -m
    """
    import argparse
    parser = argparse.ArgumentParser("set up tleap for single chain")
    parser.add_argument("-r",action="store_true",help="run tleapflag")
    parser.add_argument("-m",action="store_true",help="run sander minimization")
    parser.add_argument("-s",default="settings.json",help="settings file")
    parser.add_argument("-q",default=None,type=str,help="sequence")
    parser.add_argument("-n",default=None,type=str,help="name")
    parser.add_argument("-ff",default=None,type=str,help="ff")
    parser.add_argument("-w",default=None,type=str,help="w")
    if unknown is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(unknown)
    sequence = args.q

    if os.path.exists(args.s):
        settings = log.json_load(args.s)
    else:
        print("settings file {} doesn't exist, loading from template".format(args.s))
        settings = log.json_load(templatedir+"tleap_defaults.json")
        log.json_write(args.s,settings)

    if args.n is not None:
        settings["name"] = args.n
    if args.ff is not None:
        settings["ff"] = args.ff
    if args.w is not None:
        settings["water"] = args.w

    setup(sequence, settings)
    if args.r:
        run()
    if args.m:
        minimize(settings["name"])

