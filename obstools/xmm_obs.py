import os
import glob
from pathlib import Path
import copy
#import configparser
#from configparser import ExtendedInterpolation
import logging
import astropy.io.fits as pyfits
import xmmpy.obstools.xmm_exp as xe
from xmmpy.etc import read_config, addFileHandler, setup_logging, default_config, ofn_support, update_source_in_config
#from xmmpy.etc import shell_scripts as xmm_scripts
from yaml import dump 
from ..etc import path4

def discover_file(directory, typ):
    """
    """
    directory = str(directory)
    if typ=="pn_exp_map":
        gstr = directory+"/*EPN*expmap*.ds"
    elif typ=="m1_exp_map":
        gstr = directory+"/*EMOS1*expmap*.ds"
    elif typ=="m2_exp_map":
        gstr = directory+"/*EMOS2*expmap*.ds"
    
    fnames = glob.glob(gstr)
    
    if len(fnames)==1:
        return fnames[0]
    
    ll = logging.getLogger("xmmpy")
    if len(fnames)==0:
        ll.error("No file for "+str(typ)+ " in "+str(directory)+ "("+gstr+")")
    else:
        ll.error("More than one file for "+str(typ)+ " in "+str(directory)+ " ("+gstr+")")
    return None

class Obs():
    """
    Deals with the observation
    
      
    :ivar exposures: dictionary holding the exposures, keys are the expIDs
    """

    def __init__(self, obsID=None, conf_file = None, directory = None, populate=True):
        """
        Parameters
        ----------
        conf_file : str
          config file containing the description for the names of pn-files etc.
        
        """
        
        self.exposures = {} 
        self.initialized = False
        yy = ""
        if conf_file is None and obsID is None:
            raise ValueError("You need to provide at least an obsID, or a suitable config-file.")
        elif conf_file is None:
            yy += "No config file provided, using default values.\n"
            self.config = default_config(obsID)
        else:
            self.config = read_config(conf_file)
            if obsID is not None:
                if self.config["obsID"] != obsID:
                    raise ValueError("obsID in provided config-file ("+str(self.config["obsID"])+") and as argument do not match ("+str(obsID)+").")
        
        ll = logging.getLogger("xmmpy")
        if len(ll.handlers)==1:
          addFileHandler(self.config["LOGGING"]["log_file"])
        
        if yy !="": ll.info(yy)
        
        if directory: 
            bd = Path(directory).resolve()
            bd = str(bd).replace(str(Path.home()), "~")
            self.config["DATA"].update({"basedir":str(bd)})
            ll.info("Using basedir="+bd)

        rr = dump(self.config)
        if rr:
            ll.debug(36*"=" + " CONFIG " + 36*"=")
            ll.debug("\n"+rr)
            ll.debug(80*"=")
        
        if populate:
            ne = self.populate_exposures()
            if ne == 0: ll.info("No exposures found for "+str(self) + " in "+self.config["DATA"]["basedir"])
            else: self.initialized = True
        
    def __str__(self):
        """
        ObsID plus some text.
        """
        return str("XMM Obs ("+self.config["obsID"]+")")
            
    def __getitem__(self, k):
        """
        Simplifies the access to common observation attributes.
        
        Parameters
        ----------
        k : str
            Must be within ["decimalyear"]
        """
        if k == "decimalyear":
            e = next(iter(self.exposures.values()))
            return e["decimalyear"]
        if k in ["pn_exp_map", "m1_exp_map", "m2_exp_map"]:
            return discover_file(path4(self.config, "odata"), k)
        if k == "obsID":
            return self.config["obsID"]
    
    def write_config(self, fn=None):
        """
          Write config to file. If no filename is provided, get path for config-file from current config
        """
        if fn is None:
            fn = path4(self.config, which="conf-file")
        rr = dump(self.config)
        with open(fn, "w") as oo:
            oo.write(rr)
            
    def download_odf(self):
        """
        Download ODF from archive
        
        Returns
        --------
        filename : str
            The filename of the odf-tar.gz
        """ 
        import urllib.request

        req_str = self.config["XMM"]["archive_request_string"]+self.config["obsID"]+"&level=ODF"
        logging.debug("req_str=%s." % req_str)
        
        ofn = path4(self.config, which="odf_file")
        logging.debug("ofn=%s" % ofn)
        local_filename, headers = urllib.request.urlretrieve(req_str, ofn)       
        urllib.request.urlcleanup()
        return local_filename
    
    def odf_reduction_script(self, fn=None, ofn=None):
        """
        Generate reduction script.
        """
        if fn is None: fn = path4(self.config, which="odf_file")
        if ofn is None: ofn = path4(self.config, which="odf_reduction_script_fn")
        
        script_dir = self.config["XMM"]["scripts_directory"]
        ifn = os.path.expanduser(Path(script_dir).joinpath("make_xmm.sh"))
        xmmfn = os.path.expanduser(Path(self.config["DATA"]["basedir"]).joinpath("make_xmm.sh"))
        logging.debug("using script_fn="+str(ifn)+ " and ofn= "+str(xmmfn))
        #import shutil
        #shutil.copyfile(ifn, ofn)
        
        r = "# Generated by 'xmmpy'\n\n"
        r+="pwd=$(pwd -P)\n"
        r+= self.config["XMM"]["tools"] + "\n"
        r+="cp "+os.path.abspath(ifn)+" "+os.path.abspath(xmmfn)+"\n"
        r+="cd " + os.path.abspath(self.config["DATA"]["basedir"]) + "\n"
        r+=". make_xmm.sh "+ os.path.basename(str(fn))
        r+="\ncd ${pwd}"
        
        with open(ofn, "w") as oo:
            oo.write(r)
        
        logging.info("ODF reduction script: "+str(ofn))
        return r

    @ofn_support
    def regions4source(self, source, coord=None):
        """
          Generate source and background sources
        """        
        from astroquery.simbad import Simbad
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        if coord is None:
            customSimbad = Simbad()
            customSimbad.add_votable_fields("pm")
            result_table = customSimbad.query_object(source)
            if result_table is None or len(result_table) == 0:
                raise Exception("Cannot retrieve coordinates from Simbad for "+str(source))
            ra, dec = result_table["RA"][0], result_table["DEC"][0]
            ra, dec = ra.replace(" ", ":"), dec.replace(" ",":")
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), pm_ra_cosdec = result_table["PMRA"][0]*u.mas/u.yr, pm_dec = result_table["PMDEC"][0]*u.mas/u.yr, obstime="J2000")
        else:
            c = coord            
        
        dt = self["decimalyear"]-2000.
        cc = c.apply_space_motion(dt=dt * u.year)
        
        update_source_in_config(self.config, source)
        
        for e in self.exposures.values():
            regs = e.regions4coordinates(cc, source_name=source)
            
        return dump(self.config)
        
    def exposures_from_directory(self):
        """
          Populate exposures using files found in particular directory
        """
        import glob
        detectors = ["pn", "m1", "m2"]
        for det in detectors:
            glb_str = os.path.expanduser(str(path4(self.config, which=det+"_evt")))
            fn = glob.glob(glb_str)
            if len(fn)==1:
                exp = xe.Exposure(fn[0], self.config)            
                self.exposures[exp.exp_id] = exp
        
    def populate_exposures(self):
        """
        Looks for the event-files at the expected location (from the current  config). 
        
        Returns
        -------
        cnt : int 
            The number of event-files, or rather exposures, that were successfully created.
        """
        detectors = self.config["DATA"]["detectors"]
        ll = logging.getLogger("xmmpy")
        #ll.debug(80*"=")
        ll.info("Populating exposures:")
        cnt = 0
        for det in detectors:
            filename = path4(self.config, det+"_evt")
            ll.info("  Detector %s -> %s" % (det, filename))
            try:
                exp = xe.Exposure(filename, self.config)
                if exp.exp_id in self.exposures:
                    det0 = self.exposures[exp.exp_id].det
                    det1 = exp.det
                    if det0 != det1: # Not the same exposure
                        ii = int(exp.exp_id)+1
                        nexpid = str(ii).zfill(len(exp.exp_id))
                        exp.exp_id = nexpid
                self.exposures[exp.exp_id] = exp
                
                cnt+=1
                ll.info("        with EXP_ID=%s" % exp.exp_id)
                 
            except Exception as EE:
                ll.info(str(EE))
        return cnt    
    
    def gen_spec_shell_scripts(self, sas_init=False):
        from ..scripttools import spec_script
        r = ""
        if sas_init:
            sfn = str(path4(self.config, "SAS_init_script"))
            r+="source "+sfn+"\n\n"
    
        ll = logging.getLogger("xmmpy")
        for e in self.exposures.values(): 
            ll.info("Generating spectral extraction script for "+str(e))
            ofn = path4(self.config, which = e.det+"_spec_script")
            x = spec_script(e, ofn = ofn)
            r+=x
        ofn = path4(self.config, which = "spec_script")
        ll.info("Writing spec script to "+str(ofn))
        with open(ofn, "w") as oo:
            oo.write(r)
        return r
    
    def gen_lc_shell_scripts(self, sas_init=False):
        from ..scripttools import lc_script
        ll = logging.getLogger("xmmpy")
        r = ""
        if sas_init:
            sfn = str(path4(self.config, "SAS_init_script"))
            r+="source "+sfn+"\n\n"
    
        for e in self.exposures.values():
            ll.info("Generating light curve extraction script for "+str(e))
            ofn = path4(self.config, which = e.det+"_lc_script")
            x = lc_script(e, ofn=ofn)
            r+=x
        ofn = path4(self.config, which = "lc_script")
        ll.info("Writing light curve script to "+str(ofn))
        with open(ofn, "w") as oo:
            oo.write(r)
        return r

    def gen_evt_shell_scripts(self, sas_init=False):
        from ..scripttools import evt_script
        ll = logging.getLogger("xmmpy")
        r = ""
        if sas_init:
            sfn = str(path4(self.config, "SAS_init_script"))
            r+="source "+sfn+"\n\n"
    
        for e in self.exposures.values():
            ll.info("Generating event extraction script for "+str(e))
            ofn = path4(self.config, which = e.det+"_event_script")
            x = evt_script(e, ofn=ofn)
            r+=x
        ofn = path4(self.config, which = "event_script")
        ll.info("Writing light curve script to "+str(ofn))
        with open(ofn, "w") as oo:
            oo.write(r)
        return r
        
        
    @ofn_support    
    def shell_scripts(self, spec=None, lc=None, evt=None, sas_init=True):
        """
        Generate shell scripts for spectra (if True) and light curves (if lc==True)
        """
        if spec is None:
            spec = self.config["SOURCE PRODUCTS"]["spectra"]
        if lc is None:
            lc = self.config["SOURCE PRODUCTS"]["light curves"]
        if evt is None:
            evt = self.config["SOURCE PRODUCTS"]["events"]
        ll = logging.getLogger("xmmpy")
        ll.info("Generating source products (lc="+str(lc)+", spectra="+str(spec)+", evts="+str(evt)+").")
        r = "# xmmpy ana script\n\n"
        
        if sas_init:
            sfn = str(path4(self.config, "SAS_init_script"))
            r+="source "+sfn+"\n\n"
        
        if spec:
            r+=self.gen_spec_shell_scripts()
        if lc:
            r+=self.gen_lc_shell_scripts()
        if evt:
            r+=self.gen_evt_shell_scripts()
        return r
    
        
        
        
if __name__ == "__main__":
    o = Obs("0892000201", populate=False)
    o.download_odf()
    print(o)
    
