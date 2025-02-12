from ..etc.io_helper import ofn_support
from ..etc.my_paths import path4
import os
import logging

script_content="""
echo "Running spectral extraction (EVTFILE, YYYYYYYY)"

evselect table=EVTFILE spectrumset=FLTSETSRC withspectrumset=yes energycolumn=PI withspecranges=yes specchannelmin=0 specchannelmax=AAAAAAAA spectralbinsize=BBBBBBBB expression="XXXXXXXX && region(SRCREGIONFITS, X, Y)"

evselect table=EVTFILE spectrumset=FLTSETBKG withspectrumset=yes energycolumn=PI withspecranges=yes specchannelmin=0 specchannelmax=AAAAAAAA spectralbinsize=BBBBBBBB expression="XXXXXXXX && region(BKGREGIONFITS, X, Y)"

backscale spectrumset=FLTSETSRC badpixlocation=EVTFILE
backscale spectrumset=FLTSETBKG badpixlocation=EVTFILE

rmfgen spectrumset=FLTSETSRC rmfset=SRCRMF
arfgen spectrumset=FLTSETSRC arfset=SRCARF withrmfset=yes rmfset=SRCRMF badpixlocation=EVTFILE detmaptype=psf


if [ -f BINSPEC ]; then
   echo "File 'BINSPEC' exists."
   mv BINSPEC BINSPEC_old
fi
grppha FLTSETSRC BINSPEC comm="chkey respfile SRCRMF & chkey backfile FLTSETBKG & chkey ancrfile SRCARF & YYYYYYYY & exit"

######################################

"""


@ofn_support
def spec_script(exp, postfix="", t0=None, t1=None):
    """
    Parameters
    -----------
    exp : xmmpy.Exposure instance
    t0, t1 : float
      begin and end of spectrum, in XMM seconds
    """
    #print(exp.config)
    filt = ""
    if exp.config["SPECTRA"]["use_filtered"] == True:
        filt="_filt"
    d = exp.det
    
    spec_dir = os.path.expanduser(str(path4(exp.config, "specdir")))
    if not os.path.exists(spec_dir): os.mkdir(spec_dir)
    
    evt_file = str(path4(exp.config, d+"_evt"+filt))
    src_spec_file = str(path4(exp.config, d+"_src_spec_file", postfix=postfix))
    bkg_spec_file = str(path4(exp.config, d+"_bkg_spec_file", postfix=postfix))

    src_reg = str(path4(exp.config, "src_reg"))
    bkg_reg = str(path4(exp.config, "bkg_"+d+"_reg"))
    rmf = str(path4(exp.config, d+"_rmf", postfix=postfix))
    arf = str(path4(exp.config, d+"_arf", postfix=postfix))
    binspec = str(path4(exp.config, d+"_bin", postfix=postfix))
    
    binning = exp.config["SPECTRA"]["binning_expression"]
    
    ll = logging.getLogger("xmmpy")
    ll.debug("Generating spec script for "+str(exp)+":")
    ll.debug("    spec_dir: "+str(spec_dir))
    ll.debug("    evt_file: "+str(evt_file))
    ll.debug("    src_reg: "+str(src_reg))
    ll.debug("    bkg_reg: "+str(bkg_reg))
    ll.debug("    binning: "+str(binning))      
          
    n_script = script_content.replace("EVTFILE",evt_file)
    n_script = n_script.replace("FLTSETSRC",src_spec_file) 
    n_script = n_script.replace("FLTSETBKG",bkg_spec_file) 

    n_script = n_script.replace("SRCREGIONFITS",src_reg) 
    n_script = n_script.replace("BKGREGIONFITS",bkg_reg) 

    n_script = n_script.replace("SRCARF",arf) 
    n_script = n_script.replace("SRCRMF",rmf) 

    n_script = n_script.replace("YYYYYYYY",binning) 

    n_script = n_script.replace("BINSPEC",binspec) 


    if d == "m1" or d == "m2":
        n_script = n_script.replace("XXXXXXXX", "#XMMEA_EM && (PATTERN<=12)XXXXXXXX") # expresion
        n_script = n_script.replace("BBBBBBBB", "15") # spectralbinsize
        n_script = n_script.replace("AAAAAAAA", "11999") # specchannelmax

    elif d == "pn":
        n_script = n_script.replace("XXXXXXXX", "(FLAG==0) && (PATTERN<=4)XXXXXXXX")
        n_script = n_script.replace("BBBBBBBB", "5") # spectralbinsize
        n_script = n_script.replace("AAAAAAAA", "20479") # specchannelmax
    else:
        raise Exception("spec_script - No suitable detector ("+str(d)+"); must be in [pn, m1, m2]")

    if t0 is None and t1 is None:
        n_script = n_script.replace("XXXXXXXX", "")
    elif t1 is None:
        n_script = n_script.replace("XXXXXXXX", str(" && (TIME>%10.3f)" % (t0)))
    elif t0 is None:
        n_script = n_script.replace("XXXXXXXX", str(" && (TIME<%10.3f)" % (t1)))
    else:
        n_script = n_script.replace("XXXXXXXX", str(" && (TIME>%10.3f) && (TIME<%10.3f)" % (t0, t1)))
    return n_script
    
    
#@ofn_support
#def spec_script(exp, verbose=10):
    #"""
    #"""
    ##print("XXX", exp.obs_start_time, exp.obs_stop_time)
    #if exp.config["SPECTRA"]["spec_time_bins"] and str(exp.config["SPECTRA"]["spec_time_bins"]).strip().lower() != "none":
        #tb = exp.config["SPECTRA"]["spec_time_bins"]
        #if verbose>1: print("Setting spec time bins", tb, type(tb))
        #if type(tb)==int:
            #if verbose>0: print("Producing spectra for ",tb," bins.")
        #return ""
    #else:
        #tb = None
        #return one_spec_script(exp)
   
  ##xmm_source_products.py 0892000501/ 'TOI-776'
 
 
