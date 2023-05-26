#!/usr/bin/python3
import argparse
from pathlib import Path
from argparse import RawDescriptionHelpFormatter
import glob

def source_regions(src_name, directory, pythoncall="/home/majestix/hdd/python/bin/python3.8", filename=None):

    sas_init = glob.glob(directory+"/sas_*.sh")
    if len(sas_init) != 1:
        raise Exception("Cannot find sas-init-file in "+directory)
    obsID = sas_init[0][-13:-3]
    print("# obsID form sas-init-file: ",obsID)
    
    conf_fn = glob.glob(directory+"/xmmpy*.conf")
    if len(conf_fn) == 1:
        print("# Using observation specific xmmpy-configuration ("+conf_fn[0]+")")
        conf_fn = conf_fn[0]
    elif len(conf_fn)>1:
        print("# More than one observations specific xmmpy-configuration ("+str(conf_fn)+"), ignoring...")
        conf_fn = None
    else:
        conf_fn = None
      
    r="import os\n"
    r+="import logging\n"
    r+="from xmmpy import Obs\n"
    r+="from xmmpy.etc import path4\n\n"
    if conf_fn:
        r+="o = Obs(conf_file = \""+conf_fn+"\")\n"
    else:
        r+="o = Obs(\""+obsID+"\")\n"
    r+="o.exposures_from_directory()\n"    
    
    r+="cfn = os.path.expanduser(str(path4(o.config, which=\"datadir\").joinpath(\""+src_name.replace(" ","_")+"_"+obsID+".conf\")))\n"
    r+="print(\"Writing config to \", cfn)\n"
    r+="o.regions4source(\""+src_name+"\", ofn=cfn)\n"
    r+="ll = logging.getLogger(\"xmmpy\")\n"
    r+="ll.info(\"New config-file for source=\\\""+src_name+"\\\" and obsID=\"+o.config[\"obsID\"]+\": \\n\"+cfn)\n"
    
    ret = pythoncall+" <<< '\n"+r+"'\n"  
    
    if filename is not None:
        with open(filename,"w") as oo:
            oo.write(ret)
     
    return ret

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog = 'xmm_source_regions',
                    description = 'Generate fits-region files for source and new source-specific conf-file.\nInformation for file locations with be read from config-file (script assumes xmmpy{obsID}.conf-file to be present).',
                    epilog = 'Use at your own discretion...', formatter_class = RawDescriptionHelpFormatter)
    
    parser.add_argument('directory', default='.', nargs=1, help="directory must contain xmmpy{obsID}.conf and sas_{obsID}.sh. ")
    parser.add_argument('source', nargs=1, help="Simbad findable name.")    
    parser.add_argument('--script', default=None)
    
    args = parser.parse_args()
    
    print(source_regions(args.source[0],str(Path(args.directory[0]).resolve()),filename=args.script))
    
