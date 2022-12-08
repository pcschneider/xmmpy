#from ..scripttools import spec_script as xs
#import xmm_lightcurve as xl
#import xmm_plots as xp
import glob
from .io_helper import ofn_support
from .shell_scripts_spec import *
from .shell_scripts_lc import *
import logging
import warnings
warnings.filterwarnings("ignore", message="The read_ds9 function is deprecated and may be removed in a future version.")

#logger = logging.basicConfig(filename='example.log', level=logging.DEBUG, format='%(asctime)s %(message)s', filemode='w')    
#logging.getLogger().addHandler(logging.StreamHandler())

#gen_scripts("GJ341/0892000201/odata", "test")



#@ofn_support
#def spec_script(exp):
    ##print(exp.config)
    #filt = ""
    #if exp.config["SPECTRA"]["use_filtered"] == True:
        #filt="_filt"
    #d = exp.det
    
    #evt_file = exp.config["FILENAMES"][d+"_evt_file"+filt]
    #print(evt_file)
    #src_reg = exp.config["REGIONS"]["src"]
    #bkg_reg = exp.config["REGIONS"]["bkg_"+d]
    #print("Regions: ",src_reg, bkg_reg)
    #binning = exp.config["SPECTRA"]["binning"]
    
        
def ana_scripts_pn_spec():
    r = ""
    return r

@ofn_support
def source_products_script(obs):
    # read/copy xmm EPIC ana script analog    
    for e in obs.exposures.values():
        print(e, e.det)
        # generate call signature
    return "x"


if __name__ == "__main__":
    logger = logging.basicConfig(filename='example.log', level=logging.DEBUG, format='%(asctime)s %(message)s', filemode='w')    
    logging.getLogger().addHandler(logging.StreamHandler())

    gen_scripts("GJ341/0892000201/odata", "test")
