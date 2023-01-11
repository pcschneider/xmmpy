from pathlib import Path
import glob
from .my_configs import conffile_reader

@conffile_reader()
def path4(config, which="datadir", energy_range=None):
    """
    """
    #print("which=",which)
    if which == "datadir":
            return Path(config["DATA"]["basedir"]).joinpath(config["obsID"])
    elif which == "odata":
        return Path(config["DATA"]["basedir"]).joinpath(config["obsID"], "odata")
    elif which == "specdir":
        return path4(config,"datadir").joinpath(config["DATA"]["specreldir"])
    elif which == "lcdir":
        return path4(config,"datadir").joinpath(config["DATA"]["lcreldir"])
    elif which == "evtdir":
        return path4(config,"datadir").joinpath(config["DATA"]["evtreldir"])
    
    elif which == "pn_evt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["pn_evt_file"])
    elif which == "m1_evt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["m1_evt_file"])
    elif which == "m2_evt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["m2_evt_file"])
    
    elif which == "pn_evt_filt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["pn_evt_file_filt"])
    elif which == "m1_evt_filt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["m1_evt_file_filt"])
    elif which == "m2_evt_filt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["m2_evt_file_filt"])

# EVENTS
    elif which == "pn_src_evt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["pn_evt_file"])
    elif which == "m1_src_evt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["m1_evt_file"])
    elif which == "m2_src_evt":
        return path4(config, which="datadir").joinpath("odata", config["FILENAMES"]["m2_evt_file"])
    
# REGIONS    
    elif which == "src_reg":
        return path4(config, which="odata").joinpath( config["REGIONS"]["src"])
    elif which == "bkg_pn_reg":
        return path4(config, which="odata").joinpath(config["REGIONS"]["bkg_pn"])
    elif which == "bkg_m1_reg":
        return path4(config, which="odata").joinpath(config["REGIONS"]["bkg_m1"])
    elif which == "bkg_m2_reg":
        return path4(config, which="odata").joinpath(config["REGIONS"]["bkg_m2"])

# SPECTRA   
    elif which == "pn_src_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_src_spec_prefix"]+".fits")
        #return path4(config, which="datadir").joinpath(config["DATA"]["specreldir"], config["FILENAMES"]["pn_src_spec_prefix"]+".fits")
    elif which == "m1_src_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_src_spec_prefix"]+".fits")
    elif which == "m2_src_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_src_spec_prefix"]+".fits")
    
    elif which == "pn_bkg_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_bkg_spec_prefix"]+".fits")
    elif which == "m1_bkg_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_bkg_spec_prefix"]+".fits")
    elif which == "m2_bkg_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_bkg_spec_prefix"]+".fits")
    
    elif which == "pn_rmf":
         return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_rmf_spec_prefix"]+".fits")
    elif which == "m1_rmf":
         return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_rmf_spec_prefix"]+".fits")
    elif which == "m2_rmf":
          return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_rmf_spec_prefix"]+".fits")
   
    elif which == "pn_arf":
         return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_arf_spec_prefix"]+".fits")
    elif which == "m1_arf":
         return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_arf_spec_prefix"]+".fits")
    elif which == "m2_arf":
          return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_arf_spec_prefix"]+".fits")
   
    elif which == "pn_bin":
         return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_bin_spec_prefix"]+".fits")
    elif which == "m1_bin":
         return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_bin_spec_prefix"]+".fits")
    elif which == "m2_bin":
          return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_bin_spec_prefix"]+".fits")

# EVENTS
    elif which == "pn_src_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["src_evt_prefix"]+"_pn.fits")
    elif which == "m1_src_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["src_evt_prefix"]+"_m1.fits")
    elif which == "m2_src_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["src_evt_prefix"]+"_m2.fits")
    elif which == "pn_bkg_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["bkg_evt_prefix"]+"_pn.fits")
    elif which == "m1_bkg_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["bkg_evt_prefix"]+"_m1.fits")
    elif which == "m2_bkg_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["bkg_evt_prefix"]+"_m2.fits")
    

# LIGHT CURVES   
    elif which == "pn_src_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
    elif which == "m1_src_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
    elif which == "m2_src_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
        
    elif which == "pn_bkg_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
    elif which == "m1_bkg_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
    elif which == "m2_bkg_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
        
    elif which == "pn_crr_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
    elif which == "m1_crr_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
    elif which == "m2_crr_lc":
        if energy_range is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+energy_range+".fits")
    
    elif which == "odf_file":
        #if obsid == None: raise Exception("Need an observation for which='odf_file'")
        return Path(config["DATA"]["basedir"]).joinpath(config["obsID"]+".tar.gz")
    elif which == "odf_reduction_script_fn":
        #if obsid == None: raise Exception("Need an observation for which='odf_reduction_script_fn'")
        return Path(config["DATA"]["basedir"]).joinpath(config["FILENAMES"]["odf_reduction_fn"]+config["obsID"]+".sh")
    elif which == "SAS_init_script":
        return path4(config, which="datadir").joinpath(config["XMM"]["SAS_init_script"])

    elif which == "conf-file":
        return path4(config, which='datadir').joinpath(config["obsID"]).joinpath("xmmpy.conf")
    
    elif which == "pn_spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script_pn"])
    elif which == "m1_spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script_m1"])
    elif which == "m2_spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script_m2"])
    elif which == "spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script"])

    
    elif which == "pn_event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script_pn"])
    elif which == "m1_event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script_m1"])
    elif which == "m2_event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script_m2"])
    elif which == "event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script"])
    
    
    elif which == "pn_lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script_pn"])
    elif which == "m1_lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script_m1"])
    elif which == "m2_lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script_m2"])
    elif which == "lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script"])
    
    elif which == "ana_script":
        return path4(config, which="datadir").joinpath(config["FILENAMES"]["ana_script"])
    
    #elif which == "src_reg":
        ##if obsid == None: raise Exception("Need an observation for which='src_reg'")
        #if src == None: raise Exception("Need a source name for which='src_reg'")
        #return path4(config, which="datadir", obsid=obsid).joinpath(config["Sources"][src][obsid]["src_reg"])
    elif True:
        return None



def path44(config, which="datadir", obsid=None, src=None):
    """
    """
    if which == "datadir":
        if obsid is not None:
            return Path(config["DATA"]["basedir"]).joinpath(config["obsID"])
        return config["Data"]["basedir"]    
    elif which == "specdir":
        return Path(config["DATA"]["basedir"]).joinpath(config["DATA"]["specreldir"])
    elif which == "lcdir":
        return Path(config["DATA"]["basedir"]).joinpath(config["DATA"]["lcreldir"])
    elif which == "pn_evt":
        return path4(config, which="datadir",obsid=obsid).joinpath(config["FILENAMES"]["pn_evt_file"])
    elif which == "m1_evt":
        return path4(config, which="datadir").joinpath(config["Observations"][obsid]["m1_evt_file"])
    elif which == "m2_evt":
        return path4(config, which="datadir").joinpath(config["Observations"]["relpath"]).joinpath(config["Observations"][obsid]["m2_evt_file"])
    elif which == "odf_file":
        if obsid == None: raise Exception("Need an observation for which='odf_file'")
        return Path(config["DATA"]["basedir"]).joinpath(config["obsID"]+".tar.gz")
    elif which == "odf_reduction_script_fn":
        if obsid == None: raise Exception("Need an observation for which='odf_reduction_script_fn'")
        return Path(config["DATA"]["basedir"]).joinpath(config["FILENAMES"]["odf_reduction_fn"]+config["obsID"]+".sh")
    elif which == "src_reg":
        if obsid == None: raise Exception("Need an observation for which='src_reg'")
        if src == None: raise Exception("Need a source name for which='src_reg'")
        return path4(config, which="datadir", obsid=obsid).joinpath(config["Sources"][src][obsid]["src_reg"])
    elif True:
        return 
