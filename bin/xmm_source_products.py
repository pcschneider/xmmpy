#!/usr/bin/python3
import argparse
from argparse import RawDescriptionHelpFormatter, RawTextHelpFormatter
import glob
import os 
from difflib import SequenceMatcher
import tempfile

possible_data_products = ["rgs","spec","lc", "evt", "img"]

def scripts_from_config(config, tmp=None, products=None):
    r = "import logging\nfrom xmmpy.etc import default_config, read_config, path4\nfrom xmmpy import Obs\n\n"
    #r+= "conf = read_config(\""+config+"\")\n"
    r+= "o = Obs(conf_file = \""+config+"\")\n"
    r+= "cfn = str(path4(o.config, which=\"ana_script\"))\n"
    if products is None:
        r+= "o.shell_scripts(ofn = cfn)\n"
    else:
        tp = ", ".join([str(x)+"="+str(products[x]) for x in products])
        r+= "o.shell_scripts(ofn = cfn, "+tp+")\n"
    if tmp:
        r+= "with open(\""+tmp.name+"\", \"w\") as oo:\n"
        r+= "    oo.write(str(cfn))\n"
    r+= "ll = logging.getLogger(\"xmmpy\")\n"
    r+="ll.info(\"Source data product script:\"+cfn)\n"
    
    return r

def scripts_from_directory(directory, source=None, tmp=None, products=None):
    def get_best_matching_config(name):
        br, bc = 0, ""
        for pc in possible_configs:
            r = SequenceMatcher(None, name, pc).ratio()
            #print(pc, " -> ", r)
            if r > br: 
                bc=pc
                br=r
        if br>0:
            return bc
        else:
            raise Exception("Did not find a matching config for "+str(name)+" in "+directory)
        
    r = "import logging\nfrom xmmpy.etc import default_config, read_config, path4\nfrom xmmpy import Obs\n\n"
    #print("directory",directory)
    #r+= "conf = read_config(\""+config+"\")\n"
    possible_configs = glob.glob(os.path.join(directory,"*.conf"))
    #print("possible_configs: ", possible_configs)
    
    if source is not None:
        src_name = source.replace(" ","_")        
        try:
            config = get_best_matching_config(src_name)
        except:
            config = None
    else:
        try:
            config = get_best_matching_config("xmmpy.config")
        except:
            config = None
            
    if config is None:
        #Try to generate default config
        raise Exception("No config found in directory  \'"+str(directory)+"\', and auto-config not implemented yet...")
    r+= "o = Obs(conf_file = \""+config+"\")\n"
    r+= "cfn = str(path4(o.config, which=\"ana_script\"))\n"
    if products is None:
        r+= "o.shell_scripts(ofn = cfn)\n"  
    else:
        print("YYY",products)
        r+=""
    if tmp:
        r+= "with open(\""+tmp.name+"\", \"w\") as oo:\n"
        r+= "    oo.write(str(cfn))\n"
    r+= "ll = logging.getLogger(\"xmmpy\")\n"
    r+="ll.info(\"Source data product script:\"+cfn)\n"
    #if isinstance(filename , str):
        #r+= "o.shell_scripts(ofn = \""+str(filename)+"\")\n"
    #else:
        #r+= "o.shell_scripts(ofn = "+str(filename)+")\n"
    return r
    

def scripts(directory=None, config=None, source=None, pythoncall="/home/majestix/hdd/python/bin/python3.8", filename=None, products=None, dry=False):
    """
      Generate script that downloads and processes a particular obsID
    """
    def check_products_arg(args):
        if args is None: return None
        tp = []
        for x in args:
            if x.lower() in possible_data_products: tp.append(x.lower())
            else: print("Don't know what to do with 'product=", x,"'")
        dct = {p: True if p in tp else False for p in possible_data_products}
        if len(dct) == 0: return None
        return dct
    
    tmp = tempfile.NamedTemporaryFile(delete=True, dir='.')
    ret = None
    checkedp = check_products_arg(products)
    if config is not None:
       ret = pythoncall+ " <<< '\n"+scripts_from_config(config, tmp=tmp, products=checkedp)+"'\n"  
    elif directory is not None:
        ret = pythoncall+ " <<< '\n"+scripts_from_directory(directory, source=source, tmp=tmp, products=checkedp)+"'\n"  
    else:
        print("You need to provide either a directory, a directory and source name, or a config-file.")
    ret+="\n"
    ret+="fn=$(cat "+tmp.name+")\n"
    ret+="rm "+tmp.name+"\n"
    ret+="log_fn=${fn%.sh}.log\n"
    ret+="echo \"Running ${fn} with log-file ${log_fn}\"\n"
    if dry: ret+="#"
    ret+=". $fn 2>&1 | tee  $log_fn \n"
    if filename is not None:
        with open(filename,"w") as oo:
            oo.write(ret)
     
    return ret

if __name__ == "__main__":
    # requires:
    #
    # export xmmpy=/home/majestix/hdd/tools/xmmpy
    # export PATH=$PATH:${xmmpy}/bin
    # export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools

    parser = argparse.ArgumentParser(
                    prog = 'xmm_source_products',
                    description = 'Generate script that produces source products for a given directory.\n Only some parameters are exposed, fine-tuning by adjusting the xmmpy-config file.\n Default xmm-config-file will be gerenated upon first run if not existent.',
                    epilog = 'There are three options:\n 1) provide an xmmpy-config file\n 2) provide a directory, the script will generate a standard xmmpy-config file assuming \n       that src.reg\', \'bkg_pn.reg\', \'bkg_m1.reg\', and \'bkg_m2.reg\' exist in the odata directory.\n 3) provide a directory and source name in which case the same as in 2) applies, but \n       with \'_sourcename\' included in the regions files (e.g., \'src_sourcename.reg\', etc.).\n\nExample: \n\n Use at your own discretion...', formatter_class = RawTextHelpFormatter)
    
    #parser.add_argument('obsID')    
    #, nargs='?', default=os.getcwd()
    parser.add_argument("directory",  default=os.getcwd(),  help="Refers to \'basedir/obsID\'/odata, i.e., the directory containing the \'odata\'-folder.",nargs='?')
    parser.add_argument('source', default=None, help="Source-name; the script will try to use the config-file named \'directory\'/$ource_obsID.config",nargs='?') 
    
    parser.add_argument('--config', default=None,  help="An xmmpy config-file, path is relative to \'directory\', i.e., \'$directory/$config\'.")
    
    parser.add_argument('--script', default=None, help="Generate script-file, filename will be the value provided by the \'script\'-argument.")

    # spec=None, lc=None, evt=None, rgs=None
    parser.add_argument('-p', '--products', default=None, nargs='*', help="Optionally select specific data product(s); select 1+ from: "+str(possible_data_products)+".\nOverrides source product processing steps from config-file, i.e., only selected data products will be generated.")#, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-d', '--dry', default=False, action='store_true', help="You can call the script to generate the scripts, but they will  _not_ be run.")
    args = parser.parse_args()
    
    print(scripts(directory=args.directory,source=args.source, config=args.config, filename=args.script, products=args.products, dry=args.dry))
    
