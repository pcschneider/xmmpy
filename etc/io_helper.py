import functools
from astropy.io import fits as pyfits
import logging
from ..obstools import ds9_to_physical
#from xmmpy.obstools import yyy


def region_io_support(n, verbose=1):    
    """
    """
    def region_io_func(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            args_tmp = list(args)
            if len(args)<n:
                raise TypeError("region_io_support: Expecting at least n (="+str(n)+") arguments, but too few provided.")
            #regs =[]
            for i in range(n):
                if type(args[i]) == type("xxx"):
                    if verbose>0: print("region_io_support: Assuming ", args[i], " is a region file.")
                    if "evt_file" in kwargs: 
                        reg = reg2physical(args[i], evt_fn=kwargs["evt_file"])
                    else:    
                        reg = reg2physical(args[i])
                    #regs.append(reg)s
                    args_tmp[i] = reg
                    r = func(*tuple(args_tmp), **kwargs)
            return r    
        return wrapper
    return region_io_func


def reg2physical(fn, evt_fn=None):
        
    src = ds9_to_physical(fn, evt_fn)
    #bkg = ds9_to_physical(bkg_fn, evt_fn)
    
    if src is None:
        print("Cannot transfer source file regions to physical units, provide event file for conversion. Exiting...")
       
    try:
        import pToolsUtils as ptu
        import pToolsRegion as ptr
    except:
        raise ImportError("Cannot import pTools")
    
    
    circ = ptr.circle((src.center.x, src.center.y, src.radius))
    ll.debug("circ: %s" % circ)
    x, y = ptu.realCentroid(evt_fn, circ, eLo=cen_lo, eHi=cen_hi)
    ll.debug(" -> centroid: %f, %f" % (x, y))
    ll.info(" -> Centroid is offset by dx, dy = %f, %f (physical, region - centrod)." % (src.center.x-x, src.center.y-y))
    
    if centroid:
        src.center.x = x
        src.center.y = y
        ll.info("Using centroid position as source region center (x, y=%f, %f)." % (src.center.x, src.center.y))
    return src    
  
  
def fits_region_file_writer(pix_region, ofn, overwrite=True):
    """
    Currently only accepts "Circles"
    
    """
    cols = [pyfits.Column(name="SHAPE", array=["CIRCLE" ], format='16A')]
    
    nc = pyfits.Column(name="X",coord_type='RA---TAN', coord_unit = 'deg', array=[pix_region.center.x], format='E')
    cols.append(nc)
    
    nc = pyfits.Column(name="Y",coord_type='DEC--TAN', coord_unit = 'deg', array=[pix_region.center.y], format='E')
    cols.append(nc)
    
    nc = pyfits.Column(name="R", array=[pix_region.radius], format='E')
    cols.append(nc)

    hdu = pyfits.PrimaryHDU()

    hd = pyfits.BinTableHDU.from_columns(cols)
    hd.header["EXTNAME"] = "REGION"
    hd.header["HDUCLAS1"] = "REGION"
    hd.header["MFORM1"] = "X,Y"
    hd.header["HDUCLASS"] = "ASC"
    hd.header["MTYPE1"] = "pos"
    hdul = pyfits.HDUList([hdu, hd])
    for c in hdul[1].columns:
        
        print(c, type(c), hdul[1].data[c.name])
    hdul.writeto(ofn, overwrite=True)        
    
  
def ofn_support(func, verbose=1):
    """
    Adds the ofn keyword-option to a function, i.e., the return-value can be written to a file specified by 'ofn=x.x'
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if "ofn" in kwargs:
            ofn = kwargs["ofn"]
            del kwargs["ofn"]
        else:
            ofn = None
        r = func(*args, **kwargs)
        if ofn is not None:
            #ll = logging.get_logger()
            with open(ofn,"w") as oo:
                logging.info("Writing to \'"+str(ofn)+"\'")
                for line in r:
                    oo.write(line)
        return r
    return wrapper    
  
