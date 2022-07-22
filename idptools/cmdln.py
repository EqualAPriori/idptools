# Main script for centralized handling of commandline arguments
# instead of having to add more arguments to the setup.py
# requires each module to expose their own customized handling of commandline arguments
# another potential paradigm is to have each module have their own dispatchers.
#
import os, sys, argparse
import utility.log as log
from . import config
from . import aa
from . import tleap

def dispatch():
    parser = argparse.ArgumentParser("master cmdline dispatcher")
    parser.add_argument("cmd",type=str,help="command name")
    args, unknown = parser.parse_known_args()

    if args.cmd == "pyleap":
        tleap.cmdln(unknown)
    elif args.cmd == "implicit_relax":
        implicit_relax(unknown)
    elif args.cmd == "ELP_homo":
        ELP_homo(unknown)
    else:
        print("Unrecognized command: {}".format(args.cmd))

# ===== functions =====
# -----> create homo-ELP sequence
#pyleap -q VPGVGVPGVGVPGVGVPGVGVPGVG -ff ff14SBonlysc -r -n V5
def ELP_homo(cmdln_args):
    parser = argparse.ArgumentParser("Make ELP chain and relax using ff14SBonlysc + igb8 (defaults)")
    parser.add_argument("X",type=str,help="Guest residue")
    parser.add_argument("n",type=int,help="number of pentamer repeats")
    parser.add_argument("-centered",action="store_true",help="use centered repeat definition PGXGV")
    parser.add_argument("-reversed",action="store_true",help="use reversed repeat definition XGPVG following Chilkoti")
    args = parser.parse_args(cmdln_args)
    X,nrepeat = args.X,args.n

    if args.centered:
        pentamer = "PG{}GV".format(X)
        name = f"{X}{nrepeat}_centered"
    elif args.reversed:
        pentamer = "{}GPVG".format(X)
        name = f"{X}{nrepeat}_reversed"
    else:
        pentamer = "VPG{}G".format(X)
        name = f"{X}{nrepeat}"
    seq = pentamer*nrepeat

    cmd = f"pyleap -q {seq} -ff ff14SBonlysc -w igb8 -n {name} -r"
    print("executing: {}".format(cmd))
    os.system(cmd)


# -----> simple implicit solvent relaxation on whatever system is in settings.json
#python ~/lib/idp/idptools/scripts/md_implicit.py test -ps 1000 -stride 100 -solv 5 > z.out &
#implicit_relax(ps=None,stride=None,T=None,salt=None,**kwargs):
def implicit_relax(cmdln_args=None):
    parser = argparse.ArgumentParser("simple implicit solvent relaxation")
    parser.add_argument("-ps",type=int,default=None,help="ps")
    parser.add_argument("-stride",type=int,default=None,help="stride (ps)")
    parser.add_argument("-T",type=float,default=None,help="T (K)")
    parser.add_argument("-salt",type=float,default=None,help="salt (M)")
    parser.add_argument("-device",type=int,default=0,help="GPU device")
    if cmdln_args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(cmdln_args)
    ps,stride,T,salt = args.ps,args.stride,args.T,args.salt

    settings = log.json_load("settings.json")
    if "implicit_relax" not in settings:
        implicit_relax = {"ps":1000*500, "stride":100, "T":298.15, "salt":0.1}
    else:
        implicit_relax = settings["implicit_relax"]

    if ps is not None:
        implicit_relax["ps"] = ps
    if stride is not None:
        implicit_relax["stride"] = stride
    if T is not None:
        implicit_relax["T"] = T
    if salt is not None:
        implicit_relax["salt"] = salt
    
    cmd = "python " + config.path + "/scripts/md_implicit.py"
    cmd += " {}".format(settings["name"])
    if settings["water"].startswith("igb"):
       cmd += " -solv {}".format( int(settings["water"][3]) )
    cmd += " -ps {}".format(implicit_relax["ps"])
    cmd += " -stride {}".format(implicit_relax["stride"])
    cmd += " -T {}".format(implicit_relax["T"])
    cmd += " -salt {}".format(implicit_relax["salt"])
    cmd += " -device {}".format(args.device)

    print("executing: {}".format(cmd))
    settings["implicit_relax"] = implicit_relax
    log.json_write("settings.json",settings)

    os.system(cmd)



