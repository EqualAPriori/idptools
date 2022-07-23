# Main script for centralized handling of commandline arguments
# instead of having to add more arguments to the setup.py
# requires each module to expose their own customized handling of commandline arguments
# another potential paradigm is to have each module have their own dispatchers.
#
import os, sys, argparse, distutils
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
    elif args.cmd == "md_simple":
        md_simple(unknown)
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
#md_simple(ps=None,stride=None,T=None,salt=None,**kwargs):
def md_simple(cmdln_args=None):
    parser = argparse.ArgumentParser("simple implicit solvent relaxation")
    parser.add_argument("-ps",type=int,default=None,help="ps")
    parser.add_argument("-stride",type=int,default=None,help="stride (ps)")
    parser.add_argument("-T",type=float,default=None,help="T (K)")
    parser.add_argument("-p",type=float,default=None,help="pressure (bar)")
    parser.add_argument("-salt",type=float,default=None,help="salt (M)")
    parser.add_argument("-cutoff",type=float,default=None,help="cutoff (nm)")
    parser.add_argument("-init",default=None,help="initial structure (pdb)")
    parser.add_argument("-device",type=int,default=0,help="GPU device")
    parser.add_argument("-style",type=str,default="amber",choices=["amber","omm","gromacs"],help="GPU device")
    parser.add_argument("-dispcorr",default=None,type=lambda x:bool(distutils.util.strtobool(x)),help="toggle on dispersion correction")
    if cmdln_args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(cmdln_args)
    ps,stride,T,p,salt = args.ps,args.stride,args.T,args.p,args.salt
    cutoff,device = args.cutoff,args.device

    settings = log.json_load("settings.json")
    if "md_settings" not in settings: #default settings
        md_settings = {"ps":1000*500, 
                            "stride":100, 
                            "T":298.15, 
                            "p":None,
                            "salt":None, 
                            "device":0,
                            "dispcorr":False,
                            "cutoff":None,
                            "init":None}
    else:
        md_settings = settings["md_settings"]

    #handling arguments
    if ps is not None:
        md_settings["ps"] = ps
    if stride is not None:
        md_settings["stride"] = stride
    if T is not None:
        md_settings["T"] = T
    if p is not None:
        md_settings["p"] = p
    if salt is not None:
        md_settings["salt"] = salt
    if cutoff is not None:
        md_settings["cutoff"] = cutoff
    if args.init is not None:
        md_settings["init"] = args.init
    if args.device is not None:
        md_settings["device"] = args.device
    if args.dispcorr is not None:
        md_settings["dispcorr"] = args.dispcorr
    print(md_settings)
   
    #running -- potentially have the *simulation script* document this section...
    #cmd = "python " + config.path + "/scripts/md_implicit.py"
    if args.style == "amber":
        cmd = "python " + config.path + "/scripts/md_amber.py"
    elif args.style == "omm":
        cmd = "python " + config.path + "/scripts/md_omm.py"
    cmd += " {}".format(settings["name"])
    cmd += " -ps {}".format(md_settings["ps"])
    cmd += " -stride {}".format(md_settings["stride"])
    cmd += " -T {}".format(md_settings["T"])
    if settings["water"].startswith("igb"):
       cmd += " -solv {}".format( int(settings["water"][3]) )
    if md_settings["p"] not in [None,0.]:
        cmd += " -p {}".format(md_settings["p"])
    if md_settings["salt"] not in [None,0.]:
        cmd += " -salt {}".format(md_settings["salt"])
    if md_settings["cutoff"] not in [None,0.]:
        cmd += " -cutoff {}".format(md_settings["cutoff"])
    if md_settings["init"] not in [None]:
        cmd += " -init {}".format(md_settings["init"])
    if md_settings["dispcorr"]:
        cmd += " -dispcorr"
    cmd += " -device {}".format(md_settings["device"])

    print("executing: {}".format(cmd))
    sys.stdout.flush()

    settings["md_settings"] = md_settings
    log.json_write("settings.json",settings)

    os.system(cmd)



