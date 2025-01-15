from .io_helper import ofn_support
from .my_configs import cnf_support
from .my_paths import path4
import os

@ofn_support
def source_region(source_name, observation=None, radius=15):
    """
    Source region from Simbad

    Returns
    -------
    region : CircleSkyRegion in ds9-formatting
    """
    from astroquery.simbad import Simbad
    from astropy.coordinates import SkyCoord
    import regions
    import astropy.units as u

    customSimbad = Simbad()
    customSimbad.get_votable_fields()
    result_table = customSimbad.query_object(source_name)
    #print(len(result_table))
    if len(result_table) != 1:
        raise Exception("More than one entry for "+str(source_name)+" in Simbad")
    result_table = result_table[0]
    center_sky = SkyCoord(result_table["RA"], result_table["DEC"], unit=(u.hourangle, u.deg), frame='fk5')

    region_sky = regions.CircleSkyRegion(center=center_sky, radius=radius * u.arcsec)
    return region_sky.serialize(format='ds9')

def source_region_from_sourcelist(cnf_fn):
    """
    Nearest source from source list

    Returns
    -------
    coord : SkyCoord
        of nearest source
    """
    pass

@cnf_support(0)
def correct_region_format_to_fits(cnf, verbose=10):
    """
    Check if all regions used in the config are in the correct (fits) format
      and correct the format if necessary (and possible)
    
    Parameters
    ----------
    cnf : str, dict
        Configuration (filename or dictionary)
    
    Returns
    -------
    updated_files : list, [] if no files were updated
    """
    from astropy.io import fits as pyfits
    from .io_helper import fits_region_file_writer
    from regions import CirclePixelRegion, PixCoord
    def check_one_region(fn):
        if verbose>1: print("Checking ",fn)
        fn = os.path.abspath((os.path.expanduser(fn)))
        if not os.path.exists(fn):
            return None
        try:
            ff = pyfits.open(fn)
            ff.close()
        except:
            print(fn, " is not a fits-file")
            infile = open(fn, "r")
            for l in infile:
                print(l)
                if "circle" in l:
                    x = l.strip()[7:-1]
                    print("x",x)
                    a,b,c = x.split(",")
                    a,b,c = float(a), float(b), float(c)
                    print(a,b,c)
                    #Regions.PixRegion
                    regions = [CirclePixelRegion(PixCoord(x=a, y=b), radius=c)]
            updated.append(fn)
            if verbose>1: print("regions:",regions)
            fits_region_file_writer(regions[0], fn)

    updated = []
    check_one_region(path4(cnf, "src_reg"))
    for det in cnf["DATA"]["detectors"]:
        fn = path4(cnf, "bkg_"+det+"_reg")
        check_one_region(fn)
    return updated

def source_center(fn, region_index=0, wcs_ref_file=None, verbose=1):
    """
    RA, Dec from filename
    """
    import os
    from astropy.io import fits as pyfits
    from ..obstools import physical_to_sky

    fnextension = os.path.splitext(fn)[-1]
    if verbose>1: print("etc.source_center::  fnextension: ", fnextension)
    if fnextension == ".fits":
        with pyfits.open(fn) as rff:
            if len(rff[1].data["SHAPE"]) >= region_index:
                if rff[1].data["SHAPE"][region_index].upper() in ["CIRCLE"]:
                    sx, sy = rff[1].data["X"][region_index], rff[1].data["Y"][region_index]
    elif fnextension == ".reg":
        raise Exception("Not implemented yet")
    
    if wcs_ref_file is not None:
        return physical_to_sky((sx, sy), wcs_ref_file)
    return (sx, sy)

# def convert_physical_to_degree(sx, sy, wcs_ref_file=None):
#     """
#     Convert phsyical position to WCS

#     Returns
#     -------
#     RA, Dec : in degree (hopefully)
#     """
#     if wcs_ref_file is None: return sx, sy
#     from astropy.wcs import WCS
#     from astropy.coordinates import SkyCoord
#     import astropy.units as u
#     from astropy.io import fits as pyfits
#     from .helper import wcs4xmm
#     wcs = wcs4xmm(wcs_ref_file)
#     # with pyfits.open(wcs_ref_file) as ff:
#         # wcs = WCS(ff[0].header)
#         # try:
#             # wcsL = WCS(ff[0].header, key="L")
#         # except Exception as EE:
#             # print("No L-key WCS in ", wcs_fn)
#             # raise EE
#         # tmp = wcsL.wcs_world2pix(sx, sy, 1)
#         # print(tmp, sx, sy)
#     sky = wcs.wcs_pix2world(sx, sy, 1)
#     print(sky)
#     return sky[0], sky[1]

def aperture_from_fits(fn, wcs_ref_file=None, region_index=0):
   """
   Get photutils-aperture from fits-file

   Parameters
   ------------
   fn : str
       Region file
   wcs_reg_file : str
       File containing the WCS-information for converting physical to sky coordinates
   region_index : int
       Index of region in fits-file (if more than one region in file)  
   """
   from photutils.aperture import SkyCircularAperture, CircularAperture
   from astropy.io import fits as pyfits
   
   with pyfits.open(fn) as rff:
        if rff[1].data["SHAPE"].upper()=="CIRCLE":
            sx, sy, sr = rff[1].data["X"][region_index], rff[1].data["Y"][region_index], rff[1].data["R"][region_index]
        else:
            raise Exception("Only circle-regions allowed in xmmpy.etc.aperture_from_fits")

   if wcs_ref_file is not None:
        from astropy.wcs import WCS
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        with pyfits.open(wcs_ref_file) as ff:
            wcs = WCS(ff[0].header)
            try:
                wcsL = WCS(ff[0].header, key="L")
            except Exception as EE:
                print("No L-key WCS in ", wcs_fn)
                raise EE
        tmp = wcsL.wcs_world2pix(sx, sy, 1)
        sky = wcs.wcs_pix2world(tmp[0], tmp[1], 1)
    
        coord = SkyCoord(sky[0], sky[1], unit=(u.degree, u.degree))
        return SkyCircularAperture(coord, r=sr/20*u.arcsec)  
   ap = CircularAperture(sx, sy, sr)
   return ap
  


def fits_region(fn):
   """
   Extract region information from fits-file

   Parameters
   ----------
   fn : str, Path

   Returns
   -------
   regions : list
       ds9-style region string(s)
   """
   from astropy.io import fits as pyfits
   ff = pyfits.open(fn)
   ret = []
   for l in range(len(ff[1].data)):
       if ff[1].data["SHAPE"][l].upper() == "CIRCLE":
          rr = str("circle(%f, %f, %f)" % (ff[1].data["x"][l], ff[1].data["y"][l],ff[1].data["r"][l]))
          ret.append(rr)
       else:
          raise Exception("xmmpy.etc.fits_region - Only circle-regions allows")
   return ret
     
if __name__ == "__main__":
  #x = source_region("Arcturus",ofn="test.reg")    
  #print(x)
  pass
