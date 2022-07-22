# Main script for centralized handling of commandline arguments
# instead of having to add more arguments to the setup.py
# requires each module to expose their own customized handling of commandline arguments
# another potential paradigm is to have each module have their own dispatchers.
#
import sys, argparse
from . import config
from . import aa
from . import tleap

def dispatch():
    parser = argparse.ArgumentParser("master cmdline dispatcher")
    parser.add_argument("cmd",type=str,help="command name")
    args, unknown = parser.parse_known_args()

    if args.cmd == "pyleap":
        tleap.cmdln(unknown)
    else:
        print("Unrecognized command: {}".format(args.cmd))
