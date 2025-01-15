from ..etc.io_helper import ofn_support
from ..etc.my_paths import path4
import os
import logging

script_content="""
######################################
echo "Running image creation (EVTFILE, YYYYYYYY)"

evselect table=EVTFILE  withimageset=yes imageset=OFILE xcolumn=X ycolumn=Y imagebinning=IMGBIN ximagesize=IMXSIZE yimagesize=IMYSIZE expression="(PI in [ENERGLO:ENERGHI])"

"""


@ofn_support
def img_script(exp):
    d = exp.det
    out = ""


    img_dir = path4(exp.config, "imgdir")
    if not os.path.exists(img_dir): os.mkdir(img_dir)
    
    evt_file = str(path4(exp.config, d+"_evt"))
    src_evt_file = str(path4(exp.config, d+"_src_evt_file"))
    #bkg_evt_file = str(path4(exp.config, d+"_bkg_evt_file"))
    
    ll = logging.getLogger("xmmpy")

    for band in exp.config["IMAGES"]["energies"]:
        lo, hi = band.split(":")
        e_expression = str(lo)+"-"+str(hi)+"eV"
            
        img_file = str(path4(exp.config, d+"_image", postfix=e_expression))

        xim, yim = exp.config["IMAGES"]["x_size"], exp.config["IMAGES"]["y_size"]
        bt = exp.config["IMAGES"]["binning_type"]

        #src_reg = str(path4(exp.config, "src_reg"))
        #bkg_reg = str(path4(exp.config, "bkg_"+d+"_reg"))
        
        ll.debug("Generating image script for "+str(exp)+":")
        ll.debug("    expression:"+str(band))
        ll.debug("    img_dir: "+str(img_dir))
        ll.debug("    evt_file: "+str(evt_file))
        ll.debug("    ofn: "+str(img_file))
        #ll.debug("    src_reg: "+str(src_reg))
        #ll.debug("    bkg_reg: "+str(bkg_reg))
            
        n_script = script_content.replace("EVTFILE",evt_file)
        n_script = n_script.replace("OFILE",img_file) 
        #n_script = n_script.replace("OFILEBKG",bkg_evt_file) 

        n_script = n_script.replace("IMGBIN",bt) 
        n_script = n_script.replace("IMXSIZE",str(xim))
        n_script = n_script.replace("IMYSIZE",str(yim))
        
        n_script = n_script.replace("ENERGLO", str(lo))
        n_script = n_script.replace("ENERGHI", str(hi))
    
        out+=n_script
    return out

    
