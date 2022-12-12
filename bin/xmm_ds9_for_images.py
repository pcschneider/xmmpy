#!/usr/bin/python3
from xmmpy.etc import path4, read_config, conffile_reader
import argparse
from argparse import RawDescriptionHelpFormatter
from pathlib import Path
import os
import glob
from difflib import SequenceMatcher


def best_matching_file(files, name):
    """
    Return filename best matching 'name'.
    """
    br, bc = 0, None
    for fn in files:
            r = SequenceMatcher(None, name, fn).ratio()
            #print(pc, " -> ", r)
            if r > br: 
                bc=fn
                br=r
    return bc

@conffile_reader(arg=0)
def ds9_call(conf):
    """
    Generate ds9-commandline call to display source and background regions
    
    Parameters
    ----------
    conf : conf-dictionary
    """
    #print("Using ",conf)
    r = "ds9"
    for det in conf["DATA"]["detectors"]:
        r+=" "
        evt_fn = os.path.abspath(path4(conf, which=det+"_evt"))
        src_reg_fn = os.path.abspath(path4(conf, which="src_reg"))
        bkg_reg_fn = os.path.abspath(path4(conf, which="bkg_"+det+"_reg"))
        r+=evt_fn+" -region load "+src_reg_fn+" -region load "+bkg_reg_fn
    print(r)
    return r

def get_conf_file(arg1, source_name=None):
    """
    
    """
    if source_name is not None and len(source_name)==0: source_name=None
    print("get_conf_file '%s' '%s'. " % (arg1, source_name))
    if os.path.isfile(arg1):
        try:
            return ds9_call(arg1)
        except Exception as EE:
            print("Cannot read ",arg1, " as config-file.")
            raise EE
    elif os.path.isdir(arg1):
        cnf_files = glob.glob(arg1+"/*.conf")
        print("Possible conf-files: ",cnf_files)
        if source_name is not None:
            bst, cnf_file = 0, None
            return ds9_call(best_matching_file(cnf_files, source_name))
        print("No source name provided, trying default region names.")
        return ds9_call(best_matching_file(cnf_files, "xmmpy"))
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                prog = 'xmm_ds9_for_images',
                description = 'Display ds9-commandline string for a source-specific conf-file or directory and source name.\n',
                epilog = 'Use at your own discretion...', formatter_class = RawDescriptionHelpFormatter)

    parser.add_argument('directory', default='.', nargs=1, help="directory or conf-file")
    parser.add_argument('source', nargs='?', help="Simbad findable name.")    
    parser.add_argument('--script', default=None)
    
    args = parser.parse_args()
    conf = get_conf_file(args.directory[0], args.source)
    
    
