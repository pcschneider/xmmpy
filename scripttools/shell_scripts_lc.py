#from .io_helper import ofn_support
from ..etc.io_helper import ofn_support
from ..etc.my_paths import path4
import os

import logging
#from .my_logger import *

script_content="""
######################################
echo "Running light curve extraction (EVTFILE, binning: BINNING, energy: ENERGLO:ENERGHI)"

evselect table=EVTFILE withrateset=yes rateset="RATESETSRC" makeratecolumn=yes maketimecolumn=yes timecolumn=TIME timebinsize=BINNING expression="XXXXXXXX && (PI in [ENERGLO:ENERGHI]) && region(SRCREGIONFITS, X, Y)"
evselect table=EVTFILE withrateset=yes rateset=RATESETBKG makeratecolumn=yes maketimecolumn=yes timecolumn=TIME timebinsize=BINNING expression="XXXXXXXX && (PI in [ENERGLO:ENERGHI]) && region(BKGREGIONFITS, X, Y)"
epiclccorr srctslist=RATESETSRC eventlist=EVTFILE outset=RATESETCRR bkgtslist=RATESETBKG withbkgset=yes applyabsolutecorrections=no
"""

script_content="""
######################################
echo "Running light curve extraction (EVTFILE, binning: BINNING, energy: ENERGLO:ENERGHI)"

evselect table=EVTFILE withrateset=yes rateset="RATESETSRC" makeratecolumn=yes maketimecolumn=yes timecolumn=TIME timebinsize=BINNING expression="XXXXXXXX && (PI in [ENERGLO:ENERGHI]) && region(SRCREGIONFITS, X, Y)" ZZZZZZZZ
evselect table=EVTFILE withrateset=yes rateset=RATESETBKG makeratecolumn=yes maketimecolumn=yes timecolumn=TIME timebinsize=BINNING expression="XXXXXXXX && (PI in [ENERGLO:ENERGHI]) && region(BKGREGIONFITS, X, Y)" ZZZZZZZZ
epiclccorr srctslist=RATESETSRC eventlist=EVTFILE outset=RATESETCRR bkgtslist=RATESETBKG withbkgset=yes applyabsolutecorrections=no
"""

@ofn_support
def lc_script(exp):
    #print(exp.config)
    filt = ""
    if exp.config["SPECTRA"]["use_filtered"] == True:
        filt="_filt"
    d = exp.det
    
    lc_dir = path4(exp.config, "lcdir")
    if not os.path.exists(lc_dir): os.mkdir(lc_dir)
    
    out = ""
    
    ll = logging.getLogger("xmmpy")
    
    for band in exp.config["LIGHT CURVES"]["energies"]:
        lo, hi = band.split(":")
        e_expression = str(lo)+"-"+str(hi)+"eV"
        evt_file = str(path4(exp.config, d+"_evt"+filt))
        src_lc_file = str(path4(exp.config, d+"_src_lc", energy_range=e_expression))
        bkg_lc_file = str(path4(exp.config, d+"_bkg_lc", energy_range=e_expression))
        crr_lc_file = str(path4(exp.config, d+"_crr_lc", energy_range=e_expression))
        
        
        src_reg = str(path4(exp.config, "src_reg"))
        bkg_reg = str(path4(exp.config, "bkg_"+d+"_reg"))
        
        binning = exp.config["LIGHT CURVES"]["binning"]
        try:
            t0 = exp.config["LIGHT CURVES"]["t0"]
        except Exception as EE:
            #print(EE)
            #raise EE
            t0 = None
        try:
            t1 = exp.config["LIGHT CURVES"]["t1"]
        except:
            t1 = None
            
        
        ll.debug("Generating light curve script for "+str(exp)+":")
        ll.debug("    lc_dir: "+str(lc_dir))
        ll.debug("    evt_file: "+str(evt_file))
        ll.debug("    src_reg: "+str(src_reg))
        ll.debug("    bkg_reg: "+str(bkg_reg))
        ll.debug("    binning: "+str(binning))
        ll.debug("    to, t1:  "+str(t0)+" "+str(t1))
        ll.debug("    energies: "+e_expression)  
        
        
        n_script = script_content.replace("EVTFILE",evt_file)
        n_script = n_script.replace("RATESETSRC",src_lc_file) 
        n_script = n_script.replace("RATESETCRR",crr_lc_file) 

        n_script = n_script.replace("RATESETBKG",bkg_lc_file) 
        n_script = n_script.replace("BINNING",str(binning) )
        n_script = n_script.replace("ENERGLO",str(lo) )
        n_script = n_script.replace("ENERGHI",str(hi) )
        n_script = n_script.replace("SRCREGIONFITS",src_reg) 
        n_script = n_script.replace("BKGREGIONFITS",bkg_reg) 


        if d == "m1" or d == "m2":
            n_script = n_script.replace("XXXXXXXX", "#XMMEA_EM && (PATTERN<=12)")
        elif d == "pn":
            n_script = n_script.replace("XXXXXXXX", "(FLAG==0) && (PATTERN<=4)")        
        
        appendix = ""
        if t0:
            appendix+=" timemin="+str(t0)
        if t1:
            appendix+=" timemax="+str(t1)
        n_script = n_script.replace("ZZZZZZZZ",appendix)

        out+=n_script
    return out
    
    
