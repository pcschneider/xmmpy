from astropy.io import fits as pyfits
import glob
import logging
from xmmpy.etc import default_config, read_config, path4, source_center
from xmmpy.obstools import sky_to_physical, wcs4xmm
from astropy.io import fits as pyfits
import os
from regions import Regions, CirclePixelRegion, CircleSkyRegion, RegionMeta
from astropy.coordinates import SkyCoord,search_around_sky
import astropy.units as u
import numpy as np
from astroquery.simbad import Simbad
from astropy.time import Time

class Surrounding():
    """
    Coordinates are stored as
      1) physical coordinates in attributes 'px, py' (pyhsical coordinates as displayed in ds9 and as used by SAS, ie, with counting starting at 1)
      2) celestial coordinates as 'coord' in degree
    """
    def __init__(self, conf, max_sky_separation=30, source_list=True, simbad=True):
        """
        Parameters
        ----------
        max_sky_separation : float
            maximum distance between nominal source position and nearby sources (from source list and Simbad)
        """
        self.ll = logging.getLogger("xmmpy")
        self.conf_fn = conf
        self.max_sky_separation = max_sky_separation
        self.ll.debug("Reading conf-file: "+str(conf))
        dr = os.path.dirname(conf)
        self.conf = read_config(conf)
        self.ref_directory = dr
        self.init_source_pos()
        self.simbad_sources = {}
        self.source_list_sources = {}
        if source_list: self.populate_source_list_sources()
        if simbad: self.populate_simbad_sources()

    def init_source_pos(self):
        "Read source position from region file"
        if "pn" in self.conf["DATA"]["detectors"]:
            self.src_list_fname_fix = self.ref_directory+"/odata/*EPN*eboxlist_local.fits"
            self.evt_fn=path4(self.conf, which="pn_evt")
        elif "m1" in self.conf["DATA"]["detectors"]:
            self.src_list_fname_fix = self.ref_directory+"/odata/*EMOS1*eboxlist_local.fits"
            self.evt_fn=path4(self.conf, which="m1_evt")
        elif "m2" in self.conf["DATA"]["detectors"]:
            self.src_list_fname_fix = self.ref_directory+"/odata/*EMOS2*eboxlist_local.fits"
            self.evt_fn=path4(self.conf, which="m2_evt")
        else:
            raise Exception("No detector found in "+str(self.conf_fn)+" to use for getting source list.")
        src_reg_fn = path4(self.conf,"src_reg")
        ra_src, dec_src = source_center(src_reg_fn, wcs_ref_file=str(self.evt_fn))
        self.source_px, self.source_py = source_center(src_reg_fn)
        
        self.ll.debug(str(ra_src)+", "+str(dec_src)+ " from " + str(src_reg_fn)+ ' using '+str(self.evt_fn))
        self.source_coord = SkyCoord(ra_src, dec_src, unit="degree, degree")
        # print(self.source_px, self.source_py, " or ", sky_to_physical(self.source_coord, str(self.evt_fn)))
        with pyfits.open(self.evt_fn) as ff:
            do = ff[0].header["DATE-OBS"]
            self.epoch=Time(do, format='fits').decimalyear

    def populate_source_list_sources(self):
        src_list_fname = glob.glob(self.src_list_fname_fix) 
        if len(src_list_fname) != 1: 
            raise Exception("Don't know which source list to use: "+str(src_list_fname))
        src_list_fname=src_list_fname[0]
        self.ll.debug("Using source list filename: \'"+str(src_list_fname)+"\'")
        fsrc = pyfits.open(src_list_fname)
        gi = np.where((fsrc[1].data["ID_INST"]==1) & (fsrc[1].data["LIKE"]>6))[0]
        rate = fsrc[1].data["FLUX"][gi]
        si = np.argsort(rate)
        rate = rate[si]
        poserr = np.sqrt(fsrc[1].data["RADEC_ERR"][gi[si]]**2+1)
        ra, dec = fsrc[1].data["RA"][gi[si]], fsrc[1].data["DEC"][gi[si]]
        cat_coords = SkyCoord(ra, dec, unit=("degree,degree"))
        self.ll.debug("     which cotains "+str(len(cat_coords))+" source entries")
        sep2cat = self.source_coord.separation(cat_coords)
        dsi = np.argsort(sep2cat)
        # print("distances: ",sep2cat[dsi].arcsec)
        for j, co, r, pe in zip(dsi, cat_coords[dsi], rate[dsi], poserr[dsi]):
            if sep2cat[j].arcsec > self.max_sky_separation: 
                continue
            else:
                px, py = sky_to_physical(co, str(self.evt_fn))
                self.source_list_sources[fsrc[1].data["BOX_ID_SRC"][j]] = {"coord":co, "rate":r, "poserr":pe, "sky_separation":sep2cat[j], "px":px, "py":py}
    
    def populate_simbad_sources(self):
        tmp = nearby_Simbad_sources(self.source_coord, radius=self.max_sky_separation)
        for k, i in tmp.items():
            x, y = sky_to_physical(i["coord"], str(self.evt_fn))
            i["px"], i["py"] = x, y
            self.simbad_sources[k] = i
    
    def regions(self, coord_system="sky", radius=10.):
        """
        if 'ofn'== None: return str

        Parameters
        -----------
        ofn : str
        coord_system : str
            Must be in ["sky", "physical"]
        radius : float (in arcsec)
            Radius for nearby sources; the source region will have the size from the region-file
            Make cross(es) otherwise
        """        
        def physical_regions():
            # raise Exception("not implemented yet")
            src_reg_fn = path4(self.conf,"src_reg")
            source = Regions.read(src_reg_fn)
            tmp = []
            for t in source:
                t.meta = {"text":"xmmpy-source"}
                t.visual={"linewidth":2}
                t.origin='file'
                t.file=src_reg_fn
                tmp.append(t)
            tmp = Regions(tmp)
            for k, i in self.simbad_sources.items():
                if radius is None:
                    r =  1
                else:
                    r = CircleSkyRegion(center=i["coord"], radius=radius * u.arcsec, meta={"text":k}, visual={"color":'cyan', "linestyle":":"})
                r = r.to_pixel(wcs)
                r.origin='Simbad'
                # print(r.serialize(format='ds9'))
                tmp.append(r)
            for k, i in self.source_list_sources.items():
                r = CircleSkyRegion(center=i["coord"], radius=radius * u.arcsec, meta={"text":"sls"+str(k)}, visual={"color":"yellow", "linestyle":"--"})
                r = r.to_pixel(wcs)
                r.origin="sourcelist"
                # print(r.serialize(format='ds9'))
                tmp.append(r)
            return tmp
        
        def sky_regions():
            source = Regions.read(src_reg_fn)
            tmp = []
            for t in source:
                t.meta = {"text":"xmmpy-source"}
                t.visual={"linewidth":2}
                tmp.append(t.to_sky(wcs))
            tmp = Regions(tmp)
            for k, i in self.simbad_sources.items():
                r = CircleSkyRegion(center=i["coord"], radius=radius * u.arcsec, meta=RegionMeta({"text":k}), visual={"color":'cyan', "linestyle":"--"})
                # print(r.serialize(format='ds9'))
                tmp.append(r)
            for k, i in self.source_list_sources.items():
                r = CircleSkyRegion(center=i["coord"], radius=radius * u.arcsec, meta={"text":"sls"+str(k)}, visual={"color":"yellow", "linestyle":"--"})
                # print(r.serialize(format='ds9'))
                tmp.append(r)
            return tmp
        src_reg_fn = path4(self.conf,"src_reg")
        wcs = wcs4xmm(self.evt_fn)       
        if coord_system == "sky":
            regs = sky_regions()
        elif coord_system=="physical":
            regs = physical_regions()
        return Regions(regions=regs)

def nearby_Simbad_sources(coord, radius=30.*u.arcsec):
    """
    Dictionary of nearby simbad sources

    Returns
    """
    mySimbad = Simbad()
    mySimbad.add_votable_fields('distance_result')
    mySimbad.add_votable_fields("pm")
    mySimbad.add_votable_fields("fluxdata(V)")
    mySimbad.add_votable_fields("fluxdata(G)")

    result_table = mySimbad.query_region(coord, radius=radius*u.arcsec)
    if result_table is None or len(result_table) == 0:
        return None
    rc = {}
    for row in result_table:
        # print(row)
        ra, dec = row["RA"], row["DEC"]
        ra, dec = ra.replace(" ", ":"), dec.replace(" ",":")
        c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), pm_ra_cosdec = row["PMRA"]*u.mas/u.yr, pm_dec = row["PMDEC"]*u.mas/u.yr, obstime="J2000")
        rc[row["MAIN_ID"]] = {"coord":c, "sky_seperation":row["DISTANCE_RESULT"], "Gmag":row["FLUX_G"]}
    return rc

def nearest_src(conf_fn, Simbadquery=True):
    """
    Check source list and Simbad for sources that are near the current extraction region on the sky
    """
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
    
    nss = nearby_Simbad_sources(co_src)

    mySimbad = Simbad()
    mySimbad.add_votable_fields('distance_result')
    mySimbad.add_votable_fields("pm")
    mySimbad.add_votable_fields("fluxdata(V)")
    mySimbad.add_votable_fields("fluxdata(G)")

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
