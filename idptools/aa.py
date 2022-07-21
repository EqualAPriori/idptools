import re

# Utilities for working with amino acids
aa_lib = {
            ">":"NME",
            "<":"ACE",
            "A":"ALA",
            "R":"ARG",
            "N":"ASN",
            "D":"ASP",
            "C":"CYS",
            "E":"GLU",
            "Q":"GLN",
            "G":"GLY",
            "H":"HIS",
            "I":"ILE",
            "L":"LEU",
            "K":"LYS",
            "M":"MET",
            "F":"PHE",
            "P":"PRO",
            "S":"SER",
            "T":"THR",
            "W":"TRP",
            "Y":"TYR",
            "V":"VAL",
        }
aa_vals = list(aa_lib.values())
aa_keys = list(aa_lib.keys())
aa_symbol = dict(zip(aa_vals,aa_keys))

def abbreviate(seq):
    """
    Note:
        should have space between each AA
    """
    #expanded = expand(seq)
    if isinstance(seq,str):
        expanded = seq.split()
    else:
        expanded = seq

    abbreviated = ""
    for key in expanded:
        key = key.upper()
        if len(key) == 4:
            key = key[1:]
        if len(key) == 3 and key in aa_symbol:
            val = aa_symbol[key]
        elif len(key) == 1 and key in aa_keys:
            val = key
        else:
            raise ValueError("Unrecognized AA: {}".format(key))
            
        abbreviated += val
    return abbreviated

    

def expand(seq):
    """
    Notes:
        at first, only accept a single chain, not multiple chains
        also, only single letter! otherwise can't distinguish between NME=endcap or NME = 3 AA.
    """
    if isinstance(seq,str):
        #re.split('(?=[A-Z])', 'theLongAndWindingRoad') #python 3.7, split on capital letters
        #matches = re.findall('[A-Z<>][^A-Z]*', seq)
        matches = re.findall('[A-Z<>][a-z]*', seq)
        #>>> seq = ">AlaThr<" 
        #['>', 'Ala', 'Thr', '<'] #after the RegExp step
    elif isinstance(seq,(list,tuple)):
        matches = list(seq)
    else:
        raise ValueError("unparseable type {} for sequence {}".format(type(seq),seq))
    matches = [match.strip() for match in matches]

    expanded = []
    for im,match in enumerate(matches):
        if match.upper() in aa_lib:
            #direct match
            val = aa_lib[match.upper()]
        elif match.upper() in aa_vals:
            #if given full name of AA
            val = match.upper()
        elif match[1:].upper() in aa_vals:
            #if C-terminated or N-terminated
            if im > 0 and im < len(matches)-1:
                raise ValueError("Can't have terminated AA {} in middle of chain".format(match))
            elif match[0].upper() in ["N","C"]:
                val = match[1:]
            else:
                raise ValueError("Unrecognized AA {}".format(match))
        else:
            raise ValueError("Unrecognized AA {}".format(match))

        if val not in ["ACE","NME"]:
            if im == 0:
                val = "N" + val
            elif im == len(matches)-1:
                val = "C" + val

        expanded.append(val)

    return expanded



