from ..etc.io_helper import ofn_support
from ..etc.my_paths import path4
from ..etc.region_tools import source_center
import os
import logging

script_content="""
ZZZZZZZZ

source SASINIT

######################################
echo "Running RGS extraction (YYYYYYYY)"

obs=AAAAAAAA
#   observation type (S for scheduled observation, U for unscheduled)
otype=S
expno1=BBBBBBBB
expno2=CCCCCCCC
withsrc=y
srcra=EEEEEEEE
srcdec=FFFFFFFF
srcid=GGGGGGGG
srcnr=SCRNR
# Name extension for output
prfx="HHHHHHHH"

# Enable/disable background-correction 
# Note that when a separate background file is produced and
# used, this correction should be disabled (=n).
bkgcor=n

# apply GTI file if necessary
gti=n
gtifile=gti.fits

#-----------------------------------------------------
#     from here on no changes should be necessary    |
#-----------------------------------------------------


#   set binning for response matrices
#    set to one if one need for wavelength scale
#    other wise set to 4000
ebins=1

# ------ change directory --------------
pwd=$(pwd -P)
cd RGSDIR

# #--------rgsproc----------------------
echo '  running rgsproc, please be patient...'
  if [ "${gti}" == "n" ]; then
    rgsproc srcra=${srcra} srcdec=${srcdec} withsrc=${withsrc} srclabel=${srcid} bkgcorrect=${bkgcor} withprefix=y prefix=${prfx} withmlambdacolumn=yes RGSPROC_EXTRA_ARGS -V 2 # >& my_rgsproc_logfile;
  else
    rgsproc srcra=${srcra} srcdec=${srcdec} withsrc=${withsrc} srclabel=${srcid} bkgcorrect=${bkgcor} withprefix=y prefix=${prfx} auxgtitables=${gtifile} withmlambdacolumn=yes RGSPROC_EXTRA_ARGS -V 2 # >& my_rgsproc_logfile
  fi
  [ $? = 0 ] || die "rgsproc failed"
echo '  ...done'

#--------RGS background lightcurve and gti_file (reprocess data if necessary)  --------------------

# evselect table="${did}R1${otype}${expno1}EVENLI0000.FIT:EVENTS" makeratecolumn=yes maketimecolumn=yes timecolumn=TIME timebinsize=100 expression="(CCDNR == 9) && ((M_LAMBDA,XDSP_CORR) in REGION(${did}R1${otype}${expno1}SRCLI_0000.FIT:RGS1_BACKGROUND))" rateset="rgs1_bglc$nex.fits"
# tabgtigen table="rgs1_bglc$nex.fits" gtiset="gti_rgs1$nex.fits" expression="(RATE < 1.0)"

# evselect table="${did}R2${otype}${expno2}EVENLI0000.FIT:EVENTS" makeratecolumn=yes maketimecolumn=yes timecolumn=TIME timebinsize=100 expression="(CCDNR == 9) && ((M_LAMBDA,XDSP_CORR) in REGION(${did}R2${otype}${expno2}SRCLI_0000.FIT:RGS2_BACKGROUND))" rateset="rgs2_bglc$nex.fits"
# tabgtigen table="rgs2_bglc$nex.fits" gtiset="gti_rgs2$nex.fits" expression="(RATE < 1.0)"

#--------region and banana plot----------------------
echo
###   R1  ###
R1_table=$(ls ${prfx}R1*EVENLI*.FIT)
if [ ${#R1_table} -ge 1 ];
  then
    R1_srclst=$(ls ${prfx}R1*SRCLI*.FIT)
    echo "Using R1 table='${R1_table}' and source list='${R1_srclst}' and srcid='${srcnr}'"

    evselect table="${R1_table}:EVENTS" withimageset=yes imageset="${prfx}_spatial_R1.fits" xcolumn='M_LAMBDA' ycolumn='XDSP_CORR'
    evselect table="${R1_table}:EVENTS" withimageset=yes imageset="${prfx}_banana_R1.fits" xcolumn='M_LAMBDA' ycolumn='PI' withyranges=yes yimagemin=0 yimagemax=3000 expression="region(${R1_srclst}:RGS1_SRC${srcnr}_SPATIAL,M_LAMBDA,XDSP_CORR)"
    rgsimplot endispset="${prfx}_banana_R1.fits" spatialset="${prfx}_spatial_R1.fits" srcidlist="${srcnr}" srclistset="${R1_srclst}" withendispregionsets=yes withendispset=yes withspatialregionsets=yes withspatialset=yes device=/cps plotfile="${prfx}_region_R1.ps"
  else 
    echo "No R1 event table."
fi

###   R2  ###

R2_table=$(ls ${prfx}R2*EVENLI*.FIT)
if [ ${#R2_table} -ge 1 ];
  then
    R2_srclst=$(ls ${prfx}R2*SRCLI*.FIT)
    echo
    echo "Using R2 table='${R2_table}' and source list='${R2_srclst}' and srcid='${srcnr}'"

    evselect table="${R2_table}:EVENTS" withimageset=yes imageset="${prfx}_spatial_R2.fits" xcolumn='M_LAMBDA' ycolumn='XDSP_CORR'
    evselect table="${R2_table}:EVENTS" withimageset=yes imageset="${prfx}_banana_R2.fits" xcolumn='M_LAMBDA' ycolumn='PI' withyranges=yes yimagemin=0 yimagemax=3000 expression="region(${R1_srclst}:RGS1_SRC${srcnr}_SPATIAL,M_LAMBDA,XDSP_CORR)"
    rgsimplot endispset="${prfx}_banana_R2.fits" spatialset="${prfx}_spatial_R2.fits" srcidlist="${srcnr}" srclistset="${R2_srclst}" withendispregionsets=yes withendispset=yes withspatialregionsets=yes withspatialset=yes device=/cps plotfile="${prfx}_region_R2.ps"
  else 
    echo "No R2 event table."
fi


#-------- done ------------------
cd ${pwd}

"""


@ofn_support
def rgs_script(conf, rgsproc_extra_args=None):
    #print(exp.config)
    gti = conf["RGS"]["gti"]
        
    spec_dir = os.path.expanduser(str(path4(conf, "rgsdir")))
    if not os.path.exists(spec_dir): os.mkdir(spec_dir)
    ll = logging.getLogger("xmmpy")
    ll.debug("Generating rgs script for "+str(conf["DATA"]["source_name"])+":")
    
    wcs_det = conf["DATA"]["detectors"][0]
    
    wcs_ref_file = str(path4(conf, wcs_det+"_evt"))
    ll.debug("Using wcs-Detektor="+wcs_det+" and \'"+ wcs_ref_file+"\'")
    src_reg_file = str(path4(conf, "src_reg"))
    x, y = source_center(src_reg_file, wcs_ref_file=wcs_ref_file)
    prefix=str(conf["DATA"]["source_name"]).replace(" ","_")+"_"+str(conf["obsID"])
    # src_spec_file = str(path4(exp.config, d+"_src_spec_file"))
    # bkg_spec_file = str(path4(exp.config, d+"_bkg_spec_file"))

    # binning = exp.config["SPECTRA"]["binning_expression"]
    
    ll.debug("Source position from \'"+str(src_reg_file)+"\', which translates to src RA="+str(x)+", Dec="+str(y))
    ll.debug("    gti: "+str(gti))
    ll.debug("    rgsdir: "+str(spec_dir))
    ll.debug("    prefix: "+prefix)
    # ll.debug("    src_reg: "+str(src_reg))
    # ll.debug("    bkg_reg: "+str(bkg_reg))
    # ll.debug("    binning: "+str(binning))      
          
    n_script = script_content.replace("AAAAAAAA",conf["obsID"])
    n_script = n_script.replace("BBBBBBBB","4") 
    n_script = n_script.replace("CCCCCCCC","5") 
    # print("XXX", path4(conf,"SAS_init_script"), conf  )
    n_script = n_script.replace("SASINIT", str(os.path.abspath(path4(conf,"SAS_init_script"))))
    n_script = n_script.replace("EEEEEEEE",str(x)) 
    n_script = n_script.replace("FFFFFFFF",str(y)) 
    n_script = n_script.replace("HHHHHHHH",str(prefix))
    n_script = n_script.replace("SCRNR",str(3))
    n_script = n_script.replace("GGGGGGGG",str(conf["DATA"]["source_name"]).replace(" ","_")) 
    n_script = n_script.replace("RGSDIR", str(os.path.abspath(spec_dir)))
    n_script = n_script.replace("YYYYYYYY", conf["DATA"]["source_name"]+", obsID "+conf["obsID"]+" -> "+str(os.path.abspath(spec_dir)))
    n_script = n_script.replace("ZZZZZZZZ",str(conf["XMM"]["tools"])) 
    if rgsproc_extra_args is not None:
       n_script = n_script.replace("RGSPROC_EXTRA_ARGS", str(rgsproc_extra_args))
    else:
       n_script = n_script.replace("RGSPROC_EXTRA_ARGS","")       
    #print(n_script)
    return n_script
    
    