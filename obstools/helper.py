import logging
import os
import glob
from regions import PixelRegion
import regions
import astropy.io.fits as pyfits
from astropy import wcs

def identify_longest_exposure(dr, verbose=1):
    """
    Get the longest exposure for all detectors in one directory

    Parameters
    ----------
    dr : string
      Directory

    Returns
    --------
    dictionary with det:fn-pairs
    """
    dct = {}
    for det in ["EPN", "EMOS1", "EMOS2"]:
        best_fn = None
        longest = None
        gstr = dr+"/*"+det+"*ImagingEvts.ds*"
        fnames = glob.glob(gstr)
        if len(fnames) == 0:
          if verbose>0: print("Cannot find evt-file for ",det, " in ", dr)
          continue
        elif len(fnames) == 1:
           dct[det] = fnames[0]
        else: # More than one 'ImagingEvts' file
            for fn in fnames:
                with pyfits.open(fn) as ff:
                    if ff[0].header['FILTER'].lower() not in ["closed", "calclosed"]:
                      if longest is None:
                         best_fn = fn
                         longest = ff[1].header['ONTIME']
                      elif longest is not None and longest < ff[1].header['ONTIME']:
                         best_fn = fn
                         longest = ff[1].header['ONTIME']
            dct[det] = best_fn
    ret = {}
    if 'EPN' in dct:
       ret["pn"] = dct["EPN"]
    if "EMOS1" in dct:
       ret["m1"] = dct["EMOS1"]
    if "EMOS2" in dct:
       ret["m2"] = dct["EMOS2"]       
    return ret
              
                                          
def he_lc_from_dr(dr, det='pn'):
    """det must be in [pn, m1,m2]"""
    mp = {"pn":"EPN", "m1":"EMOS1", "m2":"EMOS2"}
    fn = dr+"/"+det+".fits"
    if os.path.exists(fn):
       with pyfits.open(fn) as ff:
          expidstr = ff[0].header["EXPIDSTR"]
       he_fnames = glob.glob(dr+"/*"+mp[det]+"_"+expidstr+"*ImagingEvt_he_lc.fits")
       if len(he_fnames) == 1: return he_fnames[0]
    return None
   

def ds9_to_physical(fname, evt_fn=None):
  """
    Regions from ds9-region file
  """
  regs = regions.Regions.read(fname)
  phys = True
  for r in regs:
    if  not isinstance(r, PixelRegion):
      phys = False
  if phys == False:
    if evt_fn is None:
      print("file ",fname, " not in physical coordinates and no event-files provided.")
      return None
    w = wcs4xmm(evt_fn)
    rr = []
    for r in regs:
      rr.append(r.to_pixel(w))
    regs = rr  
  return regs    


def sky_to_physical(skycoord, evts):
    """
      Convert sky coordinate to physical (detector) coordinates
    """
    if isinstance(evts, str):
        w = wcs4xmm(evts)
    elif isinstance(evts, wcs.WCS):
        w = evts
    else:
        raise Exception("Parameter \'evts\' must be either filename or wcs-instance.")
    tmp = w.world_to_pixel(skycoord)    
    return (tmp[0]+1, tmp[1]+1) # Due to different numbering in ds9 and python (empirically somewhat verified)



def physical_to_sky(coord, evts):
    """
      Convert physical (detector) coordinate to sky coordinates

      Parameters
      ----------
      coord : (x, y)
    """
    if isinstance(evts, str):
        w = wcs4xmm(evts)
    elif isinstance(evts, wcs.WCS):
        w = evts
    else:
        raise Exception("Parameter \'evts\' must be either filename or wcs-instance.")
    tmp = w.pixel_to_world([coord[0]-1], [coord[1]-1])[0]   
    return (tmp.ra.degree, tmp.dec.degree) # Due to different numbering in ds9 and python (empirically somewhat verified)


def reg4det(det, config, which="src"):
    idir = config["FILES"]["basedir"]

    if which=="src":
        reg_prefix = config["FILES"]["src_reg_prefix"]
        gstr = idir+"/"+config["FILES"]["data_subdir"]+"/"+reg_prefix+".reg"
    elif which=="bkg":
        reg_prefix = config["FILES"]["bkg_reg_prefix"]
        gstr = idir+"/"+config["FILES"]["data_subdir"]+"/"+reg_prefix+"_"+det+".reg"
    else:
        raise BaseException("Don't know anything about region-type %s" % which)
    
    fn = glob.glob(gstr)
    if len(fn)==0:
        raise Exception(which+" region file for "+det+" not found (looking for "+gstr+")")
    elif len(fn)!=1:
        raise Exception("More than one matching "+which+" region file for "+det+" found ("+gstr+")")
    
    return fn[0]


def wcs4xmm(fn):
  """
  WCS-instance from XMM-Newton event file, or image
  """  
  ff = pyfits.open(fn)
  if "HDUCLAS1" in ff[0].header and ff[0].header["HDUCLAS1"].strip() == "IMAGE":
      return wcs.WCS(ff[0].header)
  
  w = wcs.WCS(naxis=2)
  w.wcs.crpix = [ff[1].header["REFXCRPX"], ff[1].header["REFYCRPX"]]
  w.wcs.cdelt = [ff[1].header["REFXCDLT"], ff[1].header["REFYCDLT"]]
  w.wcs.crval = [ff[1].header["REFXCRVL"], ff[1].header["REFYCRVL"]]
  w.wcs.ctype = [ff[1].header["REFXCTYP"], ff[1].header["REFYCTYP"]]
  return w
      
