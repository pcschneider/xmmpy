import yaml
from pathlib import Path
from .my_logger import *
import os
import functools

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
  evtreldir : "evts"
  rgsreldir : "rgs"
  detectors : ['pn', 'm1', 'm2']
  source_name : None
  
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
  time_bins : None

EVENTS:
  script: event_ana.sh
  script_pn : event_ana_pn.sh
  script_m1 : event_ana_m1.sh
  script_m2 : event_ana_m2.sh

RGS:
  script : RGS_ana.sh
  gti : FALSE

SOURCE PRODUCTS:
  spectra : TRUE
  light curves : TRUE
  events : TRUE
  rgs : TRUE
  
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
    
    src_evt_prefix : "evt_src"
    bkg_evt_prefix : "evt_bkg"
    
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
    
    ana_script : ana.sh
    
    odf_reduction_fn : reduce_odf_
    config-file-postfix : .xmmpy
    """

    config = yaml.safe_load(rr)
    config["XMM"].update({"SAS_init_script":"sas_"+str(obsID)+".sh"})
    bd = Path(config["DATA"]["basedir"]).resolve()
    bd = str(bd).replace(str(Path.home()), "~")
    config["DATA"].update({"basedir":str(bd)})
    return config


def update_value_in_config(ifn, updater=None, ofn=None):
    """
    Change one parameter in config-file
    
    Example
    --------
    .. code-block::

        update_value_in_config()
    """
    from yaml import dump 
    cnf = read_config(ifn)
    
    if updater:
        updater(cnf)
    
    rr = dump(cnf)
    if ofn is None:
        ofn = ifn
    with open(ofn, "w") as oo:
        oo.write(rr)
        
  
def rewrite_config(ifn, base_config_func=None, config_manipulator_func=None, ofn=None):
    """
    Rewrite config-file
    
    Parameters
    ----------
    ifn : str
        Filename of config-file (ATTENTION: Will be overwritten if ofn is None)
    base_config_func : function
        Function that generates the base-config filename, `base_config_func` takes a config-object as argument
    config_manipulator_func : function
        Manipulates config-object, takes config-object as argument
    ofn : str
        Outfile
        
    Example
    -------
    
    .. code-block::
    
       rewrite_config('x.conf', lambda c: str(path4(c, "datadir"))+"/xmmpy"+c["obsID"]+".conf", lambda p:       p["XMM"].update({"SAS_init_script":"sas_"+str(p["obsID"])+".sh"}), "test.conf")
    
    """
    from yaml import dump 

    c = read_config(ifn)
    pure_config_fn = base_config_func(c)
    pconfig = read_config(pure_config_fn)
    
    update_source_in_config(pconfig, c["DATA"]["source_name"])
    
    if config_manipulator_func:
        config_manipulator_func(pconfig)
    
    rr = dump(pconfig)
    if ofn is None:
        ofn = ifn
    with open(ofn, "w") as oo:
        oo.write(rr)
        
  
def update_config2(dct1, dct2, verbose=1):
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


def update_config(dct1, dct2, verbose=1):
    """
      Update every item in dct1 with a possibly existint entry in dct2
      (dct1=default, dct2=new, return: dct1
    """
    def update_section(s1, s2):
        """
        Make sure every item in s2 is also present in s1
        (s1=default, s2=new)
        """
        for k in s2.keys():
            s1[k] = s2[k]
            if verbose>1: print("    k",k," (",s1[k]," -> ",s2[k],")")
                
    for s in dct2.keys():
        if verbose>1: print("config-file section",s)
        if s in dct1:
            if verbose>1: print(s, type(dct1[s]))
            if isinstance(dct1[s], dict):
                update_section(dct1[s], dct2[s])
            else: 
                if verbose>1: print(" no dict ",type(s), "  ",  dct1[s]," -> ",dct2[s])
                dct1[s] = dct2[s]
        else: # s not in dct1
            if verbose>1: print(s, "not in default, setting to: ", str(dct2[s]))
            dct1[s] = dct2[s]
        if verbose>1: print()
    if verbose>2: print("XXX",dct1)
    return dct1

def config_file_updater(ifn, dct, ofn=None, verbose=1, overwrite=True):
    """
    """
    cnf = read_config(ifn)
    n = update_config(cnf, dct)
    if ofn is not None:
        write_config(n, ofn=ofn, overwrite=overwrite)
    return n

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
        elif sect == "EVENTS": osect = "event"
        else: raise Exception("Only SPECTRA, LIGHT CURVES, and EVENTS allowed as sections names.")
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
    
    config["DATA"]["source_name"] = source.strip()
    
    update_region("src")
    for det in config["DATA"]["detectors"]:
        update_region("bkg_"+det)
        update_script_name(det, "SPECTRA")
        update_script_name(det, "LIGHT CURVES")
        update_script_name(det, "EVENTS")

        update_source_product_names(det, "spec_prefix")
        update_source_product_names(det, "lc_prefix")
        
    update_key_value(config["FILENAMES"], "src_evt_prefix", src_on+"_evt_src")
    update_key_value(config["FILENAMES"], "bkg_evt_prefix", src_on+"_evt_bkg")    
    update_key_value(config["FILENAMES"], "ana_script", src_on+"_EPIC_ana.sh")   
    update_key_value(config["SPECTRA"], "script", src_on+"_spec_ana.sh")   
    update_key_value(config["LIGHT CURVES"], "script", src_on+"_lc_ana.sh")   
    update_key_value(config["EVENTS"], "script", src_on+"_evt_ana.sh")   


def conffile_reader(arg=0, verbose=1):
    """
    Can be used to decorate functions that require a config-dictionary as input to also accept filenames with a config.
    
    The decorator assumes that the conf-argument is the arg-th argument in the function call. 
    """
    def creader(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if isinstance(arg, int):
                if isinstance(args[arg], dict):
                    return func(*args, **kwargs)
                elif isinstance(args[arg], str):
                    if os.path.isfile(args[arg]):
                        cnf = read_config(args[arg])
                        nargs=list(args)
                        nargs[arg] = cnf
                        args=tuple(nargs)
                        return func(*args, **kwargs)
                    else:
                        raise FileNotFoundError("Cannot find "+str(args[arg]))
                else:
                    raise TypeError("Argument for conf not compatible. ("+str(args[arg])+")")
                
            elif isinstance(kwargs[arg], str):
                if isinstance(kwargs[arg], dict):
                    return func(*args, **kwargs)
                elif isinstance(kwargs[arg], str):
                    if os.path.isfile(kwargs[arg]):
                        cnf = read_config(kwargs[arg])
                        kwargs[arg] = cnf
                        return func(*args, **kwargs)
                    else:
                        raise FileNotFoundError("Cannot find "+str(args[arg]))
                else:
                    raise TypeError("Argument for conf not compatible. ("+str(args[arg])+")")
            
            else:
                raise TypeError("Argument for 'arg' not compatible, must be str or int. ("+str(arg)+","+str(type(arg))+")")
        return wrapper    
    return creader


def make_abs_path(cnf_fn, ofn=None):
    """
    Replace relative pathes with pathes relative to '~'
    """
    cnf = read_config(cnf_fn)
    bd = cnf["DATA"]["basedir"]
    p = os.path.expanduser(Path(bd))
    if os.path.exists(p) and os.path.isdir(p):
        new_p = os.path.abspath(p)
    else:
        raise Exception("The config-file is invalid as 'exists="+str(os.path.exists(p))+"' and 'dir="+str(os.path.isdir(p))+"' for 'path="+str(p)+"'.")
    relp = os.path.relpath(new_p, start=os.path.expanduser("~"))
    update_p = "~/"+str(relp)
    print(new_p, relp, update_p)
    # cnf["DATA"]["basedir"] = update_p
    # cnf.write(ofn)
    update_value_in_config(cnf_fn, lambda c: c["DATA"].update({"basedir":str(update_p)}), ofn=ofn)


def read_config(filename, verbose=1):
    """
    Make sure that the configuration parameters are available
    
    Parameters
    ----------
    config - str or dict
    verbose - int
        Higher increaeses vebosity
    """
    df = default_config()
    
    import logging
    ll = logging.getLogger("xmmpy")
    
    if filename is None:
        ll.debug("No config-file provided...")
        return dict()
    
    if isinstance(filename, str):
        import logging
        if verbose>0: ll.info(str("Reading configuration file \'%s\'." % filename))
        with open(filename, mode="r") as fp:
            config = yaml.safe_load(fp)
    elif isinstance(filename, dict):
        config = filename
    else:
        logging.error("Don't know what to do with %s as a configuration." % str(config))
        raise ValueError("read_config - Cannot use object provided as \'config\'.")
    
    return update_config(df, config)


def write_config(cnf, ofn='test.conf', overwrite=False):
    from yaml import dump 

    rr = dump(cnf)
    if ofn is None:
        raise Exception("Need to provide 'ofn'.")
    if overwrite==False and os.path.exists(ofn):
        raise Exception("Filename ("+str(ofn)+") exists and 'overwrite==False'.")
    with open(ofn, "w") as oo:
        oo.write(rr)
        
def cnf_support(kw):
    """
    Adds the  keyword-option to a function, i.e., the return-value can be written to a file specified by 'cnf'
    """
    def deco(func, verbose=1):
        # def cnf_support(func, verbose=1):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            nargs = args
            if isinstance(kw, int):
                if len(args)<kw:
                    raise Exception("Not enough arguments provided, need at least "+str(kw)+" args.")
                cf = args[kw]
                # if cf isinstance()
                if not isinstance(cf, dict):
                    cf = read_config(cf)
                    nargs = *args[0:kw], cf, *args[kw+1:]
                    nargs = tuple([x for x in nargs if x!=() ])
                    if verbose>1: print("int arg, new args: ",nargs)
            elif kw in kwargs:
                cf = kwargs[kw]
                if not isinstance(cf, dict):
                    cf = read_config(cf)
                    kwargs[kw] = cf
                    if verbose>1: print("kw arg, new kwargs: ",kwargs)
            r = func(*nargs, **kwargs)
            return r
        return wrapper    
    return deco
    
  
def cnf_ofn_support(func, verbose=1):
    """
    Adds the ofn keyword-option to a function, i.e., a returned config can be written to a file specified by 'ofn=x.conf'
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
            write_config(r, ofn, overwrite=True)
        return r
    return wrapper    
  
