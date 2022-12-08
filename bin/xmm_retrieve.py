#!/usr/bin/python3
import argparse
import tempfile
from pathlib import Path

def retriever(obsID, directory, pythoncall="/home/majestix/hdd/python/bin/python3.7", filename=None):
    """
      Generate script that downloads and processes a particular obsID
    """
    tmp1 = tempfile.NamedTemporaryFile(delete=True, dir='.')
    tmp2 = tempfile.NamedTemporaryFile(delete=True, dir='.')
    rr = pythoncall+" <<< \'"
    rr+="""
from xmmpy import Obs 
from xmmpy.etc import path4

o = Obs(obsID=\""""
    rr+=str(obsID)
    rr+='\", populate=False)\n'
    rr+='o.config["DATA"]["basedir"] = \"'
    rr+=str(directory)
    rr+="\""
    rr+="""
o.download_odf()
r = o.odf_reduction_script()

o.write_config(\""""+tmp1.name+"""\")
rs = path4(o.config, which="odf_reduction_script_fn")
with open(\""""+tmp2.name+"""\", "w") as oo:
  oo.write(str(rs))
"""
    rr+="\'\n"
    rr+="fn=$(cat "
    rr+=tmp2.name
    rr+=")\n"
    rr+="rm "+tmp2.name+"\n"
    rr+="log_fn=${fn%.sh}.log\n"
    rr+="echo \"Running ${fn} with log-file ${log_fn}\"\n"
    rr+=". $fn &> $log_fn \n"
    #rr+="mv "+tmp1.name+" "+str(Path(directory).joinpath(obsID))+"/"+obsID+".xmmpy"
    rr+="mv "+tmp1.name+" "+str(Path(directory).joinpath(obsID))+"/xmmpy"+obsID+".conf"
    
    if filename is not None:
        with open(filename,'w') as oo:
            oo.write(rr)
    
    return rr

if __name__ == "__main__":
    # requires:
    #
    # export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools
    # export PATH=$PATH:/home/majestix/hdd/tools/xmmpy/scripttools
    parser = argparse.ArgumentParser(
                    prog = 'xmm_retrieve',
                    description = 'Generate script that downloads and processes a particular obsID',
                    epilog = 'Use at your own discretion...')
    
    parser.add_argument('obsID')    
    parser.add_argument('directory', default='.', nargs='?')
    parser.add_argument('--script', default=None)
    
    args = parser.parse_args()
    
    print(retriever(args.obsID,str(Path(args.directory).resolve()),filename=args.script))
    
