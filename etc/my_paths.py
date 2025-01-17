from pathlib import Path
import glob
import os
from .my_configs import conffile_reader


def path4(config, which="datadir",postfix=None):
    tmp = path4_ll(config, which=which, postfix=postfix)
    if str(tmp).lower().strip()=="none": return None
    return Path(os.path.expanduser(str(tmp)) )

@conffile_reader()
def path4_ll(config, which="datadir", postfix=None):
    """
    Parameters
    -----------
    config : str or dict (config-instance)
    which : str
    postfix : str
      
    """
    
    which = which.lower()
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
    elif which == "rgsdir":
        return path4(config,"datadir").joinpath(config["DATA"]["rgsreldir"])
    elif which == "imgdir":
        return path4(config,"datadir").joinpath(config["DATA"]["imgreldir"])
    
    
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
        p1 = path4(config, which="odata")
        p2 = config["REGIONS"]["bkg_pn"]
        if p2 is None: return None
        return p1.joinpath(p2)
    elif which == "bkg_m1_reg":
        p1 = path4(config, which="odata")
        p2 = config["REGIONS"]["bkg_m1"]
        if p2 is None: return None
        return p1.joinpath(p2)
    elif which == "bkg_m2_reg":
        p1 = path4(config, which="odata")
        p2 = config["REGIONS"]["bkg_m2"]
        if p2 is None: return None
        return p1.joinpath(p2)


# IMAGES
    elif which == "pn_image":
        if postfix is None:
            return path4(config, which="imgdir").joinpath(config["FILENAMES"]["pn_image_prefix"]+".fits")
        else:
            return path4(config, which="imgdir").joinpath(config["FILENAMES"]["pn_image_prefix"]+"_"+postfix+".fits")
    elif which == "m1_image":
        if postfix is None:
            return path4(config, which="imgdir").joinpath(config["FILENAMES"]["m1_image_prefix"]+".fits")
        else:
            return path4(config, which="imgdir").joinpath(config["FILENAMES"]["m1_image_prefix"]+"_"+postfix+".fits")
    elif which == "m2_image":
        if postfix is None:
            return path4(config, which="imgdir").joinpath(config["FILENAMES"]["m2_image_prefix"]+".fits")
        else:
            return path4(config, which="imgdir").joinpath(config["FILENAMES"]["m2_image_prefix"]+"_"+postfix+".fits")

        
# LIGHT CURVES   
    elif which == "pn_src_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
    elif which == "m1_src_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
    elif which == "m2_src_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_src_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
        
    elif which == "pn_bkg_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
    elif which == "m1_bkg_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
    elif which == "m2_bkg_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_bkg_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
    elif which == "pn_crr_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["pn_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
    elif which == "m1_crr_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m1_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")
    elif which == "m2_crr_lc":
        if postfix is None:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s.fits")
        else:
            return path4(config, which="lcdir").joinpath(config["FILENAMES"]["m2_crr_lc_prefix"]+"_"+str(config["LIGHT CURVES"]["binning"])+"s_"+postfix+".fits")

    # SPECTRA   
    if postfix is None: postfix=""
    if which == "pn_src_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_src_spec_prefix"]+postfix+".fits")
        #return path4(config, which="datadir").joinpath(config["DATA"]["specreldir"], config["FILENAMES"]["pn_src_spec_prefix"]+".fits")
    if which == "m1_src_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_src_spec_prefix"]+postfix+".fits")
    if which == "m2_src_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_src_spec_prefix"]+postfix+".fits")
    
    if which == "pn_bkg_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_bkg_spec_prefix"]+postfix+".fits")
    if which == "m1_bkg_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_bkg_spec_prefix"]+postfix+".fits")
    if which == "m2_bkg_spec_file":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_bkg_spec_prefix"]+postfix+".fits")
    
    if which == "pn_rmf":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_rmf_spec_prefix"]+postfix+".fits")
    if which == "m1_rmf":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_rmf_spec_prefix"]+postfix+".fits")
    if which == "m2_rmf":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_rmf_spec_prefix"]+postfix+".fits")

    if which == "pn_arf":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_arf_spec_prefix"]+postfix+".fits")
    if which == "m1_arf":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_arf_spec_prefix"]+postfix+".fits")
    if which == "m2_arf":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_arf_spec_prefix"]+postfix+".fits")

    if which == "pn_bin":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["pn_bin_spec_prefix"]+postfix+".fits")
    if which == "m1_bin":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m1_bin_spec_prefix"]+postfix+".fits")
    if which == "m2_bin":
        return path4(config, which="specdir").joinpath(config["FILENAMES"]["m2_bin_spec_prefix"]+postfix+".fits")
    
# EVENTS
    if which == "pn_src_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["src_evt_prefix"]+"_pn.fits")
    if which == "m1_src_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["src_evt_prefix"]+"_m1.fits")
    if which == "m2_src_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["src_evt_prefix"]+"_m2.fits")
    if which == "pn_bkg_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["bkg_evt_prefix"]+"_pn.fits")
    if which == "m1_bkg_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["bkg_evt_prefix"]+"_m1.fits")
    if which == "m2_bkg_evt_file":
        return path4(config, which="evtdir").joinpath(config["FILENAMES"]["bkg_evt_prefix"]+"_m2.fits")
    
# EXPOSURE MAP
    if which == "pn_expmap":
        odata = str(path4(config, which="odata"))
        print(odata)
        dd = odata+"/*EPN*ImagingEvt_expmap.ds"
        fname = glob.glob(dd)
        if len(fname)==1:
            return fname[0]
        elif len(fname)==0: 
            raise Exception("No exposure map for pn found (",dd,")")
        else:
            raise Exception("More than one exposure map found for pn (",dd,")")
    if which == "m1_expmap":
        dd = str(path4(config, which="odata"))+"/*EMOS1*ImagingEvt_expmap.ds"
        fname = glob.glob(dd)
        if len(fname)==1:
            return fname[0]
        elif len(fname)==0: 
            raise Exception("No exposure map for m1 found (",dd,")")
        else:
            raise Exception("More than one exposure map found for m1 (",dd,")")
    if which == "m2_expmap":
        dd = str(path4(config, which="odata"))+"/*EMOS2*ImagingEvt_expmap.ds"
        fname = glob.glob(dd)
        if len(fname)==1:
            return fname[0]
        elif len(fname)==0: 
            raise Exception("No exposure map for m2 found (",dd,")")
        else:
            raise Exception("More than one exposure map found for m2 (",dd,")")

    if which == "odf_file":
        #if obsid == None: raise Exception("Need an observation for which='odf_file'")
        return Path(config["DATA"]["basedir"]).joinpath(config["obsID"]+".tar.gz")
    if which == "odf_reduction_script_fn":
        #if obsid == None: raise Exception("Need an observation for which='odf_reduction_script_fn'")
        return Path(config["DATA"]["basedir"]).joinpath(config["FILENAMES"]["odf_reduction_fn"]+config["obsID"]+".sh")
    if which == "sas_init_script":
        return path4(config, which="datadir").joinpath(config["XMM"]["SAS_init_script"])

    if which == "conf-file":
        return path4(config, which='datadir').joinpath(config["obsID"]).joinpath("xmmpy.conf")
    
# SCRIPTS

    if which == "pn_spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script_pn"])
    if which == "m1_spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script_m1"])
    if which == "m2_spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script_m2"])
    if which == "spec_script":
        return path4(config, which="datadir").joinpath(config["SPECTRA"]["script"])

    
    if which == "pn_event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script_pn"])
    if which == "m1_event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script_m1"])
    if which == "m2_event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script_m2"])
    if which == "event_script":
        return path4(config, which="datadir").joinpath(config["EVENTS"]["script"])
    
        
    if which == "pn_image_script":
        return path4(config, which="datadir").joinpath(config["IMAGES"]["script_pn"])
    if which == "m1_image_script":
        return path4(config, which="datadir").joinpath(config["IMAGES"]["script_m1"])
    if which == "m2_image_script":
        return path4(config, which="datadir").joinpath(config["IMAGES"]["script_m2"])
    if which == "image_script":
        return path4(config, which="datadir").joinpath(config["IMAGES"]["script"])
    
    
    if which == "pn_lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script_pn"])
    if which == "m1_lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script_m1"])
    if which == "m2_lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script_m2"])
    if which == "lc_script":
        return path4(config, which="datadir").joinpath(config["LIGHT CURVES"]["script"])
    
    if which == "rgs_script":
        return path4(config, which="datadir").joinpath(config["RGS"]["script"])

    if which == "ana_script":
        return path4(config, which="datadir").joinpath(config["FILENAMES"]["ana_script"])
    
    if which == "image_pdf":
        return path4(config, which="imgdir").joinpath("overview.pdf")
    
    #elif which == "src_reg":
        ##if obsid == None: raise Exception("Need an observation for which='src_reg'")
        #if src == None: raise Exception("Need a source name for which='src_reg'")
        #return path4(config, which="datadir", obsid=obsid).joinpath(config["Sources"][src][obsid]["src_reg"])
    
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
