import yaml
from pathlib import Path
from .my_logger import *

def default_config(obsID=None):
    """
      Generate a default configuration file for an obsID
    """
    rr = "obsID : \""+str(obsID)+"\""
    rr += """
    
DATA:
  basedir : "." # All pathes in this config are relative to this directory
  specreldir : "specs"
  lcreldir : "lcs"
  detectors : ['pn', 'm1', 'm2']
  
LOGGING:
  log_file : "xmmpy.log"
  
XMM:
  scripts_directory : "~/hdd/scripts/XMM"
  tools : "heainit ; sasinit"
  archive_request_string : "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno="
  
LIGHT CURVES:
  binning : 1000
  energies : [300:1000, 1000:2000, 200:2000]
  script_pn : lc_ana_pn.sh
  script_m1 : lc_ana_m1.sh
  script_m2 : lc_ana_m2.sh
  script : lc_ana.sh
  use_filtered : FALSE
  
SPECTRA:
  binning_expression : "group min 3"
  script_pn : spec_ana_pn.sh
  script_m1 : spec_ana_m1.sh
  script_m2 : spec_ana_m2.sh
  script : spec_ana.sh
  use_filtered : FALSE

SOURCE PRODUCTS:
  spectra : True
  light curves : TRUE
  
REGIONS:
    src : reg_src.fits 
    bkg_pn : reg_bkg_pn.fits
    bkg_m1 : reg_bkg_m1_bkg.fits
    bkg_m2 : reg_bkg_m2_bkg.fits
    
FILENAMES:
    pn_evt_file : "pn.fits"
    m1_evt_file : "m1.fits"
    m2_evt_file : "m2.fits"    
    
    pn_evt_file_filt : "pn_filt.fits"
    m1_evt_file_filt : "m1_filt.fits"
    m2_evt_file_filt : "m2_filt.fits"    
    
    pn_src_spec_prefix : "spec_src_pn"
    m1_src_spec_prefix : "spec_src_m1"
    m2_src_spec_prefix : "spec_src_m2"
    pn_bkg_spec_prefix : "spec_bkg_pn"
    m1_bkg_spec_prefix : "spec_bkg_m1"
    m2_bkg_spec_prefix : "spec_bkg_m2"
    pn_bin_spec_prefix : "spec_mn3_pn"
    m1_bin_spec_prefix : "spec_mn3_m1"
    m2_bin_spec_prefix : "spec_mn3_m2"
    pn_rmf_spec_prefix : "spec_rmf_pn"
    m1_rmf_spec_prefix : "spec_rmf_m1"
    m2_rmf_spec_prefix : "spec_rmf_m2"
    pn_arf_spec_prefix : "spec_arf_pn"
    m1_arf_spec_prefix : "spec_arf_m1"
    m2_arf_spec_prefix : "spec_arf_m2"
    
    pn_src_lc_prefix : "lc_src_pn"
    m1_src_lc_prefix : "lc_src_m1"
    m2_src_lc_prefix : "lc_src_m2"
    pn_bkg_lc_prefix : "lc_bkg_pn"
    m1_bkg_lc_prefix : "lc_bkg_m1"
    m2_bkg_lc_prefix : "lc_bkg_m2"
    pn_crr_lc_prefix : "lc_crr_pn"
    m1_crr_lc_prefix : "lc_crr_m1"
    m2_crr_lc_prefix : "lc_crr_m2"
    
    ana_script : EPIC_ana.sh
    
    odf_reduction_fn : reduce_odf_
    config-file-postfix : .xmmpy
    """

    config = yaml.safe_load(rr)
    config["XMM"].update({"SAS_init_script":"sas"+str(obsID)+".sh"})
    bd = Path(config["DATA"]["basedir"]).resolve()
    bd = str(bd).replace(str(Path.home()), "~")
    config["DATA"].update({"basedir":str(bd)})
    return config
  
def update_config(dct1, dct2, verbose=1):
    """
      Update every item in dct1 with a possibly existint entry in dct2
    """
    def update_section(s1, s2):
        for k in s1.keys():
            if k in s2:
                s1[k] = s2[k]
                if verbose>1: print("    k",k," (",s1[k]," -> ",s2[k],")")
            else:
                if verbose>1: print("    k ",k, "(using default: ",s1[k],")")
                
    for s in dct1.keys():
        if verbose>1: print("config-file section",s)
        if s in dct2:
            if verbose>1: print(s, type(dct1[s]))
            if isinstance(dct1[s], dict):
                update_section(dct1[s], dct2[s])
            else:
                if verbose>1: print(" no dict ",type(s), "  ",  dct1[s]," -> ",dct2[s])
                dct1[s] = dct2[s]
        else:
            if verbose>1: print("Using default: ", dct1[s])
        if verbose>1: print()
    if verbose>2: print("XXX",dct1)
    return dct1

def update_source_in_config(config, source):
    """
      Update all relevant entries 
    """
    def update_key_value(leaf, key, value):
        old = Path(leaf[key])
        new =  str(old.parent.joinpath(value))
        leaf[key] = new
        
    def update_region(typ):
        src_old = Path(config["REGIONS"][typ])        
        src_new = str(src_old.parent.joinpath(src_on+"_reg_"+typ+".fits"))
        config["REGIONS"][typ] = src_new
        ll.debug("REGIONS:"+typ+" -> "+src_new)
        
    def update_script_name(typ, sect):
        scrpt_old = Path(config[sect]["script_"+typ])
        if sect == "SPECTRA": osect = "spec"
        elif sect == "LIGHT CURVES": osect = "lc"
        else: raise Exception("Only SPECTRA and LIGHT CURVES allowed as sections names.")
        scrpt_new =  str(scrpt_old.parent.joinpath(src_on+"_"+osect+"_script_"+typ+".sh"))
        config[sect]["script_"+typ] = scrpt_new
        
        ll.debug(sect+":script_"+typ+" -> "+scrpt_new)
        
        
    def update_source_product_names(d, typ):
        if typ == "spec_prefix": 
            lst = ["src","bkg","bin","rmf","arf"]
            mdl = "_spec_"
        elif typ == "lc_prefix": 
            lst = ["src","bkg", "crr"]
            mdl = "_lc_"
        else:
            raise Exception("Only \"spec_prefix\" and \"lc_prefix\" allowed for \'typ\' in <update_source_product_names!")
        
        for pdct in lst:
            k = d+"_"+pdct+"_"+typ
            old_pfx = Path(config["FILENAMES"][k])
            #print(d, typ, " -> ", old_fn)
            new_pfx =  str(old_pfx.parent.joinpath(src_on+mdl+pdct+"_"+d))
            config["FILENAMES"][k] = new_pfx
            ll.debug("FILENAMES:"+k+" -> "+new_pfx)
        
        
        
    src_on = source.strip().replace(" ","_")
    
    ll = logging.getLogger("xmmpy")
    ll.info("Updating configuration for source="+str(source)+"(using outname="+src_on+")")
    
    update_region("src")
    for det in config["DATA"]["detectors"]:
        update_region("bkg_"+det)
        update_script_name(det, "SPECTRA")
        update_script_name(det, "LIGHT CURVES")
        update_source_product_names(det, "spec_prefix")
        update_source_product_names(det, "lc_prefix")
    update_key_value(config["FILENAMES"], "ana_script", src_on+"_EPIC_ana.sh")   
    update_key_value(config["SPECTRA"], "script", src_on+"_spec_ana.sh")   
    update_key_value(config["LIGHT CURVES"], "script", src_on+"_lc_ana.sh")   
        
def read_config(filename):
    """
    Make sure that the configuration parameters are available
    
    Parameters
    ----------
    config - str or dict
    """
    df = default_config()
    
    import logging
    ll = logging.getLogger("xmmpy")
    
    if filename is None:
        ll.debug("No config-file provided...")
        return dict()
    
    if isinstance(filename, str):
        import logging
        ll.info(str("Reading configuration file \'%s\'." % filename))
        with open(filename, mode="r") as fp:
            config = yaml.safe_load(fp)
    elif isinstance(filename, dict):
        config = filename
    else:
        #import logging
        logging.error("Don't know what to do with %s as a configuration." % str(config))
        raise ValueError("read_config - Cannot use object provided as \'config\'.")
    
    return update_config(df, config)
    
    # Some checks:
      #TBI
    #return config        
    
