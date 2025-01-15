import functools
from astropy.io import fits as pyfits
import logging
from ..obstools import ds9_to_physical
#from xmmpy.obstools import yyy
from collections.abc import Iterable

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
  

def fits_single_region_writer(line, ofn, overwrite=True, verbose=1):
    """
    Write a \'line\' as a fits file
    """
    import re
    p = re.compile('([a-z]*\()')
    m = p.match(line.strip())
    region_shape = m.group()[0:-1]
    if verbose>1: print("   ->  \'%s\'" % region_shape)
        
            
    if region_shape == "circle":
        pp = re.compile("[a-z]*\(((.+),(.+),(.+)\))")
        mm = pp.match(line.strip())
        x = mm.group(2)
        y = mm.group(3)
        r = mm.group(4)
        if verbose>1: print("x=<%s>, y=<%s>, r=<%s>" % (x, y, r))
        cols = [pyfits.Column(name="SHAPE", array=["CIRCLE" ], format='16A')]
        cx = pyfits.Column(name="X",coord_type='RA---TAN', coord_unit = 'deg', array=[x], format='E')
        cy = pyfits.Column(name="Y",coord_type='DEC--TAN', coord_unit = 'deg', array=[y], format='E')
        cr = pyfits.Column(name="R", array=[r], format='E')
        cols.append(cx)
        cols.append(cy)        
        cols.append(cr)
    elif region_shape == "panda":
        pp = re.compile("[a-z]*\(((.+),(.+),(.+),(.+),(.+),(.+),(.+),(.+)\))")
        mm = pp.match(line.strip())
        x = mm.group(2)
        y = mm.group(3)
        angle0 = mm.group(4)
        angle1 = mm.group(5)
        r0 = mm.group(7)
        r1 = mm.group(8)
        if verbose>1: print("  x=<%s>, y=<%s>, r0=<%s>, r1=<%s>, angle0=<%s>, angle1=<%s>" % (x, y, r0, r1, angle0, angle1))

        cols = [pyfits.Column(name="SHAPE", array=["Pie","Annulus"], format='16A')]        
        cx = pyfits.Column(name="X",coord_type='RA---TAN', coord_unit = 'pix', array=[x,x], format='E')
        cy = pyfits.Column(name="Y",coord_type='DEC--TAN', coord_unit = 'pix', array=[y,y], format='E')
        cr = pyfits.Column(name="R", array=[[0,0],[r0,r1]], format='2E')
        #cr = pyfits.Column(name="R", array=[r1,r0], format='E')
        ca = pyfits.Column(name="ROTANG", array=[[angle0,angle1],[angle0,angle1]], format='2E')
        cols.append(cx)
        cols.append(cy)           
        cols.append(cr)
        cols.append(ca)


    hdu = pyfits.PrimaryHDU()

    hd = pyfits.BinTableHDU.from_columns(cols)
    hd.header["EXTNAME"] = "REGION"
    hd.header["HDUCLAS1"] = "REGION"
    hd.header["MFORM1"] = "X,Y"
    hd.header["HDUCLASS"] = "ASC"
    hd.header["MTYPE1"] = "pos"
    hdul = pyfits.HDUList([hdu, hd])
    #if verbose>2:
        ##for c in hdul[1].columns:
            ##print(c, type(c), hdul[1].data[c.name])
    hdul.writeto(ofn, overwrite=True)        
    return ofn



def fits_multi_region_writer(lines, ofn, overwrite=True, verbose=1, physical=True):
    """
    Write \'lines\' as a fits file
    """
    import re
    p = re.compile('([a-z]*\()')
    col_vals = {"SHAPE":[], "X":[], "Y":[], "R":[], "ROTANG":[], "COMPONENT":[]}    
            
            
    for i, line in enumerate(lines):
        m = p.match(line.strip())
        region_shape = m.group()[0:-1]
        print("   ->  \'%s\'" % region_shape)            
        if region_shape == "circle":
            pp = re.compile("[a-z]*\(((.+),(.+),(.+)\))")
            mm = pp.match(line.strip())
            x = float(mm.group(2))
            y = float(mm.group(3))
            r = float(mm.group(4))
            print("x=<%s>, y=<%s>, r=<%s>" % (x, y, r))
            col_vals["SHAPE"].append("CIRCLE")
            col_vals["X"].append(x)
            col_vals["Y"].append(y)
            col_vals["R"].append(r)
            col_vals["ROTANG"].append([None])
            col_vals["COMPONENT"].append(i+1)
            
        elif region_shape == "panda":
            pp = re.compile("[a-z]*\(((.+),(.+),(.+),(.+),(.+),(.+),(.+),(.+)\))")
            mm = pp.match(line.strip())
            x = float(mm.group(2))
            y = float(mm.group(3))
            angle0 = float(mm.group(4))
            angle1 = float(mm.group(5))
            r0 = float(mm.group(7))
            r1 = float(mm.group(8))
            print("  x=<%s>, y=<%s>, r0=<%s>, r1=<%s>, angle0=<%s>, angle1=<%s>" % (x, y, r0, r1, angle0, angle1))

            col_vals["SHAPE"].append("PIE")
            col_vals["X"].append(x)
            col_vals["Y"].append(y)
            col_vals["R"].append([r0, r1])
            col_vals["ROTANG"].append([angle0, angle1])
            col_vals["SHAPE"].append("Annulus")
            col_vals["COMPONENT"].append(i+1)
            col_vals["X"].append(x)
            col_vals["Y"].append(y)
            col_vals["R"].append([r0, r1])
            col_vals["ROTANG"].append([angle0, angle1])
            col_vals["COMPONENT"].append(i+1)
            
        elif region_shape == "Ellipse":
            pp = re.compile("[a-z]*\(((.+),(.+),(.+),(.+),(.+)\))")
            mm = pp.match(line.strip())
            x = float(mm.group(2))
            y = float(mm.group(3))
            r0 = float(mm.group(4))
            r1 = float(mm.group(5))
            a = float(mm.group(6))
            print("x=<%s>, y=<%s>, r=<%s>" % (x, y, r))
            col_vals["SHAPE"].append("Circle")
            col_vals["X"].append(x)
            col_vals["Y"].append(y)
            col_vals["R"].append([r0, r1])
            col_vals["ROTANG"].append([a])
            col_vals["COMPONENT"].append(i+1)
            
        elif region_shape == "Annulus":
            pp = re.compile("[a-z]*\(((.+),(.+),(.+),(.+)\))")
            mm = pp.match(line.strip())
            x = float(mm.group(2))
            y = float(mm.group(3))
            r0 = float(mm.group(4))
            r1 = float(mm.group(5))
            print("x=<%s>, y=<%s>, r=<%s>" % (x, y, r))
            col_vals["SHAPE"].append("Annulus")
            col_vals["X"].append(x)
            col_vals["Y"].append(y)
            col_vals["R"].append([r0, r1])
            col_vals["ROTANG"].append([None])
            col_vals["COMPONENT"].append(i+1)
        
        else:
            raise Exception("Region shape "+str(region_shape)+" not implemented yet. Aborting...")
        
    
    cols = [pyfits.Column(name="SHAPE", array=col_vals["SHAPE"], format='16A')]        
    cx = pyfits.Column(name="X",coord_type='RA---TAN', coord_unit = 'pix', array=col_vals["X"], format='E')
    cy = pyfits.Column(name="Y",coord_type='DEC--TAN', coord_unit = 'pix', array=col_vals["Y"], format='E')
    max_entries = max([len(a) if isinstance(a, Iterable) else 1 for a in col_vals["R"]]) 
    print(col_vals["R"], "max entries: ",max_entries)
    if max_entries >2 or max_entries==1:
        cr = pyfits.Column(name="R", array=col_vals["R"], format='E')
    else:
        cr = pyfits.Column(name="R", array=col_vals["R"], format=str(max_entries)+'E')
    max_entries = max([len(a) if isinstance(a, Iterable) else 1  for a in col_vals["ROTANG"]])
    if max_entries >2 or max_entries==1:
        ca = pyfits.Column(name="ROTANG", array=col_vals["ROTANG"], format='E')
    else:
        ca = pyfits.Column(name="ROTANG", array=col_vals["ROTANG"], format=str(max_entries)+'E')
    cc = pyfits.Column(name="COMPONENT", array=col_vals["COMPONENT"], format='J')
    #print(cr) 
    cols.append(cx)
    cols.append(cy)           
    cols.append(cr)
    cols.append(ca)
    cols.append(cc)

    hdu = pyfits.PrimaryHDU()

    hd = pyfits.BinTableHDU.from_columns(cols)
    hd.header["EXTNAME"] = "REGION"
    hd.header["HDUCLAS1"] = "REGION"
    hd.header["MFORM1"] = "X,Y"
    hd.header["HDUCLASS"] = "ASC"
    hdu.header["CREATOR"] = "pcs"
    hdu.header["DATE"] = "today"
    hd.header["MTYPE1"] = "pos"
    hdul = pyfits.HDUList([hdu, hd])
    #if verbose>2:
        ##for c in hdul[1].columns:
            ##print(c, type(c), hdul[1].data[c.name])
    hdul.writeto(ofn, overwrite=True)        
    return ofn


  
def fits_region_file_writer(pix_region, ofn, overwrite=True, verbose=1):
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
    if verbose>2:
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
  
