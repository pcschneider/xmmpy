from astropy.io import fits as pyfits
import glob
import logging
from xmmpy.etc import default_config, read_config, path4, source_center
from astropy.io import fits as pyfits
import os
from astropy.coordinates import SkyCoord,search_around_sky
import astropy.units as u
import numpy as np
from astroquery.simbad import Simbad
from astropy.time import Time

mySimbad = Simbad()
mySimbad.add_votable_fields('distance_result')


def nearest_src(conf_fn, Simbadquery=True):
    ll = logging.getLogger("xmmpy")
    ll.debug("Reading conf-file: "+str(conf_fn))
    dr = os.path.dirname(conf_fn)
    cnf = read_config(conf_fn)
    if "pn" in cnf["DATA"]["detectors"]:
        src_list_fname_fix = dr+"/odata/*EPN*eboxlist_local.fits"
        evt_fn=path4(cnf, which="pn_evt")
    elif "m1" in cnf["DATA"]["detectors"]:
        src_list_fname_fix = dr+"/odata/*EMOS1*eboxlist_local.fits"
        evt_fn=path4(cnf, which="m1_evt")
    elif "m2" in cnf["DATA"]["detectors"]:
        src_list_fname_fix = dr+"/odata/*EMOS2*eboxlist_local.fits"
        evt_fn=path4(cnf, which="m2_evt")
    else:
        raise Exception("No detector found in "+str(conf_fn)+" to use for getting source list.")
    
    src_reg_fn = path4(cnf,"src_reg")
    ra_src, dec_src = source_center(src_reg_fn, wcs_ref_file=str(evt_fn))
    print(ra_src, dec_src, "from",src_reg_fn, 'using',evt_fn)
    co_src = SkyCoord(ra_src, dec_src, unit="degree, degree")

    src_list_fname = glob.glob(src_list_fname_fix) 
    if len(src_list_fname) != 1: 
        raise Exception("Don't know which source list to use: "+str(src_list_fname))
    src_list_fname=src_list_fname[0]
    
    ll.debug("Using source list filename: \'"+str(src_list_fname)+"\'")
    fsrc = pyfits.open(src_list_fname)
    do = fsrc[0].header["DATE-OBS"]
    epoch=Time(do, format='fits').decimalyear
    # print(fsrc[0].header["XPROC0"])
    # exit()
    gi = np.where((fsrc[1].data["ID_INST"]==1) & (fsrc[1].data["LIKE"]>6))[0]
    rate = fsrc[1].data["FLUX"][gi]
    si = np.argsort(rate)
    rate = rate[si]
    poserr = np.sqrt(fsrc[1].data["RADEC_ERR"][gi[si]]**2+1)
    ra, dec = fsrc[1].data["RA"][gi[si]], fsrc[1].data["DEC"][gi[si]]
    
    print(len(np.unique(ra)))
    cat_coords = SkyCoord(ra, dec, unit=("degree,degree"))
    ll.debug("     which cotains "+str(len(cat_coords))+" source entries")

    sep2cat = co_src.separation(cat_coords)
    dsi = np.argsort(sep2cat)
    # print("distances: ",sep2cat[dsi].arcsec)
    ret = {"sources":[], \
           "source_list_fn":str(os.path.abspath(src_list_fname)),\
           "conf_fn":str(os.path.abspath(conf_fn)),\
           "src_region_fn":str(os.path.abspath(src_reg_fn)),\
           "src_ra_used":ra_src, "src_dec_used":dec_src
           }
    for j, co, r, pe in zip(dsi, cat_coords[dsi], rate[dsi], poserr[dsi]):
        dct = {}
        sep = co_src.separation(co).arcsec
        if sep>30: continue
        print("sep:",sep)
        dct= {"ra":co.ra.degree, "dec":co.dec.degree, "separation":sep, "rate":r, "poserr":pe, "idx":j}
        print(co, r, pe)
        print("separation (arcsec)", sep)
        if Simbadquery:
            r = mySimbad.query_region(co, radius=30*u.arcsec, epoch="J"+str(epoch))
            if r is None: continue
            if len(r.errors)>0: print("Simbad Query Errors:",r.errors)
            if len(r)==0: continue
            mi = np.where(r["DISTANCE_RESULT"]/pe<5)
            print(r[mi])
        ret["sources"].append(dct)    
        print()
    print(epoch)
    return ret
    # src_list_fnames = glob.glob(directory+"/odata/"+src_list_postfix)
    # if len(src_list_fnames)!=1:
    #     raise Exception("Don't know which source list ot use:"+str(src_list_fnames))
