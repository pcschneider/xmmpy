import astropy.io.fits as pyfits
#import xmmpy.scripttools as xss
import logging
import glob
import xmmpy.obstools.helper as helper
from xmmpy.etc import path4, fits_region_file_writer
import os 

class Exposure():
  def __init__(self, evt_filename, config):
      """
      config : configparser.ConfigParser instance
      """
      mp = {"EPN":"pn", "EMOS2":"m2", "EMOS1":"m1"}
      self.evt_filename = evt_filename
      with pyfits.open(self.evt_filename) as ff:
        self.exp_id = ff[0].header["EXP_ID"]
        self.det = mp[ff[0].header["INSTRUME"]]        
      self.config = config
      ll = logging.getLogger("xmmpy")
      ll.debug("Initialized image_file(oi, e.det)exposure with %s and evt-files \'%s\'." % (self.exp_id, self.evt_filename))
      
  def __str__(self):
      return "Exp: "+str(self.exp_id)+" ("+str(self.det)+" -> "+str(self.evt_filename)+")"
  
  def __getitem__(self, k):
      if isinstance(k, str) and k.lower() == "decimalyear":
          from astropy.time import Time
          with pyfits.open(self.evt_filename) as ff:
              x = ff[0].header["DATE-OBS"]
          d = Time(x, format='isot', scale='utc').decimalyear
          return d
      elif isinstance(k, str) and k.lower() == "filter":
          with pyfits.open(self.evt_filename) as ff:
              x = ff[0].header["Filter"]
      elif isinstance(k, str) and k.lower() == "exposure":
          print(self.evt_filename)
          with pyfits.open(self.evt_filename) as ff:
              x = ff[0].header["Duration"]         
          return x
        
  def regions4coordinates(self, coord, src_radius=15, bkg_radius=35, write=True, source_name="src"):
   """
     Parameters
     ----------
     coord - SkyCoord
     
     Returns
     -------
     regions : dict
         values in dict are CircleSkyRegions in celestial coordinates
   """
   #wcs = helper.wcs4xmm(self.evt_filename)
   #print(wcs)
   from regions import CircleSkyRegion, Regions
   from astropy.coordinates import Angle
   import numpy as np
   
   src = CircleSkyRegion(coord, Angle(src_radius, 'arcsec'), meta={"name":source_name})
   
   ll = logging.getLogger("xmmpy")
   
   with pyfits.open(self.evt_filename) as ff:
       pnt_angle = ff[0].header["PA_PNT"]
   
   rot = 5
   bkg_coord1 = coord.directional_offset_by(Angle(60-pnt_angle+rot, "degree"), Angle(80, "arcsec"))
   bkg_coord2 = coord.directional_offset_by(Angle(210-pnt_angle+rot, "degree"), Angle(80, "arcsec"))
   bkg_coord3 = coord.directional_offset_by(Angle(135-pnt_angle+rot, "degree"), Angle(80, "arcsec"))
   
   bkg1 = CircleSkyRegion(bkg_coord1, Angle(bkg_radius, 'arcsec'))
   bkg2 = CircleSkyRegion(bkg_coord2, Angle(bkg_radius, 'arcsec'))
   bkg3 = CircleSkyRegion(bkg_coord3, Angle(bkg_radius, 'arcsec'))
   
   dct = {"src":src, "bkg1":bkg1, "bkg2":bkg2, "bkg3":bkg3, "bkg": bkg3}
   
   wcs = helper.wcs4xmm(self.evt_filename)
   if write:
       ofn = str(path4(self.config, which="src_reg"))
       ll.info("Writing source region to "+ofn)    
       src_px = src.to_pixel(wcs)
       fits_region_file_writer(src_px, ofn, overwrite=True)

       #src_px.write(ofn, overwrite=True, format='fits')
       ofn = str(path4(self.config, which="bkg_"+self.det+"_reg"))
       ll.info("Writing bkg region to "+ofn)    
       bkg3_px = bkg3.to_pixel(wcs)
       fits_region_file_writer(bkg3_px, ofn, overwrite=True)
       
           
   #regions = Regions([src, bkg1, bkg2, bkg3])
            #regions.write("tmp.reg", overwrite=True)
   
   
   return dct
    
  def regions(self):
    """
    Generate src and bkg regions
    """
    idir = self.config["FILES"]["basedir"]
    src_reg_prefix =  self.config["FILES"]["src_reg_prefix"]
    bkg_reg_prefix =  self.config["FILES"]["bkg_reg_prefix"]
    det = self.det
    
    ll = logging.getLogger()

    src_reg = helper.reg4det(det, self.config, which="src")
    bkg_reg = helper.reg4det(det, self.config, which="bkg")
    ll.info("files: %s, %s" % (src_reg, bkg_reg))
      
    
    src = helper.ds9_to_physical(src_reg, evt_fn=self.evt_filename)
    if src is None:
        ll.error("Cannot transfer source file regions to physical units for reg-file=\'%s\' and evt-file=\'%s\'"  % (src_reg, self.evt_filename))
        raise BaseException("Conversion from sky to pixel coordinates failed for reg-file=\'%s\' and evt-file=\'%s\'" % (src_reg, self.evt_filename))
    if len(src)!=1:
        #ll.error("Not exactly one region in reg-file=\'%s\'"  % src_reg)
        raise BaseException("Not exactly one region in reg-file=\'%s\'"  % src_reg)
    else:
        src = src[0]
        
    bkg = helper.ds9_to_physical(bkg_reg, evt_fn=self.evt_filename)
    if bkg is None:
        ll.error("Cannot transfer source file regions to physical units for reg-file=\'%s\' and evt-file=\'%s\'"  % (bkg_reg, self.evt_filename))
        raise BaseException("Conversion from sky to pixel coordinates failed for reg-file=\'%s\' and evt-file=\'%s\'" % (bkg_reg, self.evt_filename))
    if len(bkg)!=1:
        #ll.error("Not exactly one region in reg-file=\'%s\'"  % src_reg)
        raise BaseException("Not exactly one region in reg-file=\'%s\'"  % bkg_reg)
    else:
        bkg = bkg[0]
    
    return src, bkg
  
  #def script_filename(self, which):
    #if which not in ["spec","lc"]:
        #raise BaseException("Exposure::scipt_filename called with unknown option (%s)." % which)
        
    #dd = self.config.getint('SCRIPTS', 'direcotry') if self.config.has_option('SCRIPTS', 'direcotry') else None
    #if dd == None:
        #odir = self.config["FILES"]["basedir"]+"/"+self.config["FILES"]["data_subdir"]+"/"
    #else:
        #odir = self.config["SCRIPTS"]["directory"]+"/"
        
    #ofn = odir+self.det+self.config["SCRIPTS"][which]
    #return ofn  
  
  #def lc_shell_script(self, binning=None):
    ##eLo, eHi = int(self.config["ENERGIES"]["lo"]), int(self.config["ENERGIES"]["hi"])
    #idir = self.config["FILES"]["basedir"]
    #elo, ehi = float(self.config["ENERGIES"]["lo"]), float(self.config["ENERGIES"]["hi"])
    #if binning is None: binning = float(self.config["TIMING"]["binning"])
    #ll = logging.getLogger()
    #ll.debug(80*"=")
    #ll.info("lc shell script for exp_id=%s (%s)" % (self.exp_id, self.det))

    
    #ofn = self.script_filename(which='lc')
    #src, bkg = self.regions()
    #lcf = xss.gen_new_lc_file(src, bkg, binning=binning, lo=elo, hi=ehi, ofn = ofn)
    #return ofn, lcf

  #def spec_shell_script(self):
    #idir = self.config["FILES"]["basedir"]
    #det = self.det
    #binning = int(self.config["SPECTRA"]["binning"])
    #ll = logging.getLogger()
    #ll.debug(80*"=")
    #ll.info("spec shell script for exp_id=%s (%s)" % (self.exp_id, det))

    #ofn = self.script_filename("spec")    
    #src, bkg = self.regions()
    #filtered = bool(self.config["SPECTRA"]["filtered"])
    #if filtered:
        #evt = os.path.relpath(self.filtered_evt_filename, os.path.dirname(ofn))
    #else:
        #evt = os.path.relpath(self.evt_filename, os.path.dirname(ofn))
    #ll.debug("Using evt-file \'%s\'" % evt)
    #sfile = xss.gen_new_spec_file(src, bkg, ofn=ofn, instr=det, binning=binning, verbose=1, evt_file=evt)
    #ll.info("Generated: \'%s\'." % ofn)  
    #return ofn, sfile 
        
    def gen_shell_scripts(self):
        #pass
        f0 = self.spec_shell_script()
        f1 = self.lc_shell_script()
            
        #ll.info("shell scripts: ", f0, f1)
        return (f0, f1)
