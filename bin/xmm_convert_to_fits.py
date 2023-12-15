#!/home/majestix/hdd/python/bin/python3.8
import argparse
# import tempfile
from pathlib import Path
import glob
import os
import re

pythoncall='p38'

def correct_regions(filename, pythoncall=pythoncall, ofn=None, image_file=None, overwrite=True):
    """
      description
    """
    from xmmpy.etc import read_config, path4, fits_region_file_writer, fits_multi_region_writer
    from astropy.io import fits as pyfits
    from regions import CirclePixelRegion, PixCoord
    p = re.compile('([a-z]+\(.+\))')
    if ofn is None: ofn = filename

    if not os.path.exists(filename):
        raise Exception("No such file (" + str(filename)+")")
    if not os.path.isfile(filename):
        raise Exception("Is not a file (" + str(filename)+")")
    
    try:
        ff = pyfits.open(filename)
        ff.close()
    except:
        print(filename, " is not a true fits-file")
        infile = open(filename, "r")
        
        lines = []
        for l in infile:
            print("  Processing: ",l.strip())
            m = p.match(l.strip())
            if m is None: continue
            print("   ->  \'%s\'" % m.group())
            lines.append(l.strip())
        fits_multi_region_writer(lines, ofn=ofn, overwrite=overwrite)
        
    return ofn #updated
            
if __name__ == "__main__":
    # requires:
    #
    # export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools
    # export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/scripttools
    parser = argparse.ArgumentParser(
                    prog = 'xmm_convert_to_fits.py',
                    description = 'Convert supposedly fits-file regions to real fits-files',
                    epilog = 'Use at your own discretion...')
    
    # parser.add_argument('obsID')    
    parser.add_argument('file', help="File that is some region file, but with a \'.fits\' name.")
    parser.add_argument('ofn', nargs='?', help="Filename for the new, real fits-file.", default=None)
    #parser.add_argument('--script', default=None)
    
    args = parser.parse_args()
    #print(args.ofn)
    if args.ofn is not None:
        tmp = args.ofn
        ofn = str(Path(tmp).resolve())
    else: ofn=None
    print("Updated regions: ",correct_regions(str(Path(args.file).resolve()),  ofn=ofn))
    
