#!/home/majestix/hdd/python/bin/python3.8
import argparse
# import tempfile
from pathlib import Path
import glob
import os

pythoncall='p38'

def correct_regions(directory, pythoncall=pythoncall, filename=None):
    """
      Check for region in directory and convert supposedly fits-file regions to real fits-files
    """
    from xmmpy.etc import read_config, path4, fits_region_file_writer
    from astropy.io import fits as pyfits
    from regions import CirclePixelRegion, PixCoord

    if os.path.exists(directory) and os.path.isfile(directory):
        try:
            cfg = read_config(directory)
        except:
            print("Cannot read %s a config-file." % directory)
            raise Exception("No config file (" + str(directory)+")")
    
    if os.path.exists(directory) and os.path.isdir(directory):
        cfg_fns = glob.glob(directory+"/*.conf")
        if len(cfg_fns) > 1:
            raise Exception("Too many config-files in directory ("+str(directory)+"):\n   "+ str(cfg_fns))
        elif len(cfg_fns)==0:
            raise Exception("No config-file in directory ("+str(directory))
        cfg = read_config(cfg_fns[0])
    
    fnames = [os.path.expanduser(path4(cfg, "src_reg"))]
    for det in ["pn", "m1", "m2"]:
        if det in cfg["DATA"]["detectors"]:
            reg_fn = path4(cfg, "bkg_"+det+"_reg")
            fnames.append(os.path.expanduser(reg_fn))
            print(det," -> ", reg_fn)
    # print(fnames)

    updated = []
    for fn in fnames:
        if not os.path.exists(fn):
            print(fn, " does not exis.")
            continue
        
        try:
            ff = pyfits.open(fn)
            ff.close()
        except:
            print(fn, " is not a true fits-file")
            infile = open(fn, "r")
            for l in infile:
                # print(l)
                if "circle" in l:
                    x = l.strip()[7:-1]
                    # print("x",x)
                    a,b,c = x.split(",")
                    a,b,c = float(a), float(b), float(c)
                    # print(a,b,c)
                    #Regions.PixRegion
                    regions = [CirclePixelRegion(PixCoord(x=a, y=b), radius=c)]
                # else:
                #     print(fn, " ----> Region is not a circle, aborting...")
            # print("regions:",regions)
            if len(regions)>1:
                raise Exception("More than one region in "+fn+"("+str(len(regions))+")")
            elif len(regions)==0:
                raise Exception("No circular region in "+fn)
            fits_region_file_writer(regions[0], fn)
            updated.append(fn)
    return updated
            
if __name__ == "__main__":
    # requires:
    #
    # export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools
    # export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/scripttools
    parser = argparse.ArgumentParser(
                    prog = 'xmm_correct_region_file_format.py',
                    description = 'Convert supposedly fits-file regions to real fits-files',
                    epilog = 'Use at your own discretion...')
    
    # parser.add_argument('obsID')    
    parser.add_argument('directory', default='.', help="directory must contain xmmpy{obsID}.conf and sas_{obsID}.sh.")
    parser.add_argument('--script', default=None)
    
    args = parser.parse_args()
    
    print("Updated regions: ",correct_regions(str(Path(args.directory).resolve()),filename=args.script))
    
