from ..etc.io_helper import ofn_support
from ..etc.my_paths import path4
import os
import logging

script_content="""
######################################
echo "Running event extraction (EVTFILE, YYYYYYYY)"

evselect table=EVTFILE filteredset=OFILESRC expression="region(SRCREGIONFITS, X, Y)"

evselect table=EVTFILE filteredset=OFILEBKG expression="region(BKGREGIONFITS, X, Y)"

"""


@ofn_support
def evt_script(exp):
    d = exp.det
    
    evt_dir = path4(exp.config, "evtdir")
    if not os.path.exists(evt_dir): os.mkdir(evt_dir)
    
    evt_file = str(path4(exp.config, d+"_evt"))
    src_evt_file = str(path4(exp.config, d+"_src_evt_file"))
    bkg_evt_file = str(path4(exp.config, d+"_bkg_evt_file"))

    src_reg = str(path4(exp.config, "src_reg"))
    bkg_reg = str(path4(exp.config, "bkg_"+d+"_reg"))
    
    ll = logging.getLogger("xmmpy")
    ll.debug("Generating spec script for "+str(exp)+":")
    ll.debug("    spec_dir: "+str(evt_dir))
    ll.debug("    evt_file: "+str(evt_file))
    ll.debug("    src_reg: "+str(src_reg))
    ll.debug("    bkg_reg: "+str(bkg_reg))
          
    n_script = script_content.replace("EVTFILE",evt_file)
    n_script = n_script.replace("OFILESRC",src_evt_file) 
    n_script = n_script.replace("OFILEBKG",bkg_evt_file) 

    n_script = n_script.replace("SRCREGIONFITS",src_reg) 
    n_script = n_script.replace("BKGREGIONFITS",bkg_reg) 

    return n_script
    
    
