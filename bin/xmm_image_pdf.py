#!/usr/bin/python3
import sys
#from os.path import dirname
sys.path.append("/home/majestix/hdd/python/lib/python3.8/site-packages/")
import glob
from xmmpy.qc import ax4image_and_sourceregion
from xmmpy.etc import path4
from xmmpy import Obs
from photutils.aperture import aperture_photometry, SkyCircularAperture
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits as pyfits
from astropy.wcs import WCS
import astropy.units as u
import os.path
from astropy.coordinates import SkyCoord
#from my_filenames import image_file, target_conf, coverage_fn
from collections import namedtuple

def one_page(obs, band=None, conf_file=None, verbose=10):
    def one_panel(i, e, band=None):
        """
        Plot one detector (or rather exposure)
        """
        ax0 = plt.subplot(2,2,i+1)
        x0, y0, wx, wy = ax0.get_position().bounds
            
        if verbose>2: print("bbox: ", x0, y0, wx, wy)
            
        #image_fn = "tst.pdf"#image_file(oi, e.det)
        #band = obs.config["IMAGES"]["energies"][0]
        lo, hi = band.split(":")
        e_expression = str(lo)+"-"+str(hi)+"eV"
        image_fn = path4(obs.config, e.det+"_image", postfix=e_expression)

        src_reg_fn = path4(obs.config, "src_reg")
        bkg_reg_fn = path4(obs.config, "bkg_"+e.det+"_reg")

        # if str(bkg_reg_fn).lower().strip()=="none":
        #     bkg_reg_fn=None

        


        if verbose>2: 
            print("fns ",image_fn, src_reg_fn, bkg_reg_fn, type(bkg_reg_fn))
            print("call arguments: ",e.det,  image_fn, src_reg_fn, bkg_reg_fn)
        # tmp = ax4image_and_sourceregion(image_fn, region_fn = src_reg_fn, bkg_region_fn = [bkg_reg_fn,  str(bkg_reg_fn)[0:-5]+"_new.fits"], idx=i+1)
        if bkg_reg_fn:
            tmp = ax4image_and_sourceregion(image_fn, region_fn = src_reg_fn, bkg_region_fn = bkg_reg_fn, idx=i+1)
        else:
            tmp = ax4image_and_sourceregion(image_fn, region_fn = src_reg_fn, idx=i+1)
        axs.append(tmp)
            
        with pyfits.open(src_reg_fn) as rff:
            sx, sy, sr = rff[1].data["X"][0], rff[1].data["Y"][0], rff[1].data["R"][0]
            if verbose>2: print("src-reg: ",sx, sy, sr, end="")
        if bkg_reg_fn:    
            with pyfits.open(bkg_reg_fn) as rff:
                bx, by, br = rff[1].data["X"][0], rff[1].data["Y"][0], rff[1].data["R"][0]
            if verbose>2: print(" bkg-reg:", bx,by, br)
        else:
            if verbose>2: print()   
            
        axs[i].set_title(str(e.exp_id)+" ("+e.det+")")
        line_temp = str(" ; %15s" % str(e.exp_id))
        dets_temp = "                  %s        " % e.det
        with pyfits.open(image_fn) as ff:
            wcs = WCS(ff[0].header)
            wcsL = WCS(ff[0].header, key="L")

            tmp = wcsL.wcs_world2pix(sx, sy, 1)
            sky = wcs.wcs_pix2world(tmp[0], tmp[1], 1)
            if verbose>3: print(tmp, sky)
            coord = SkyCoord(sky[0], sky[1], unit=(u.degree, u.degree))
            ap_src = SkyCircularAperture(coord, r=sr/20*u.arcsec)  
            src_phot = aperture_photometry(ff[0].data, ap_src, wcs=wcs)
            if verbose>2: print("src-phot", src_phot)

            if bkg_reg_fn and str(bkg_reg_fn).lower().strip()!="none":
                #print("XXX", bkg_reg_fn, type(bkg_reg_fn))  
                tmp = wcsL.wcs_world2pix(bx, by, 1)
                sky = wcs.wcs_pix2world(tmp[0], tmp[1], 1)
                if verbose>3: print(tmp, sky)
                coord = SkyCoord(sky[0], sky[1], unit=(u.degree, u.degree))
                ap_src = SkyCircularAperture(coord, r=br/20*u.arcsec)  
                bkg_phot = aperture_photometry(ff[0].data, ap_src, wcs=wcs)
                if verbose>2: print("bkg-phot", bkg_phot)
                net_cts = src_phot["aperture_sum"] - (sr/br)**2 * bkg_phot["aperture_sum"]
            else:
                net_cts = 0

            if verbose>2: print(" net", net_cts)
            plt.annotate("%2s - exposure: %4.1f ks, filter: %s" % (e.det, ff[0].header["EXPOSURE"]/1000., ff[0].header["Filter"]), xy=(plt_txt_x, 0.2 + i*0.03), xycoords="figure fraction")
            bb = band.split(":")
            plt.annotate("epoch: %s, energies: %6s - %6s (eV)" % (ff[0].header["DATE-OBS"].split("T")[0], bb[0], bb[1]),  xy=(plt_txt_x, 0.396), xycoords="figure fraction")
            plt.annotate("%s: ds9 %s -region load %s -region load %s" % (e.det, image_fn, src_reg_fn, bkg_reg_fn), xy=(0.02, 0.02+i*0.03), xycoords="figure fraction", size=3, annotation_clip=False)
            if bkg_reg_fn:                  
                plt.annotate("%s - src: %7i    bkg: %7i    net: %8.1f  rate (cts/ks): %8.2f" % (e.det, src_phot["aperture_sum"], bkg_phot["aperture_sum"], net_cts, net_cts/ff[0].header["EXPOSURE"]*1000), xy=(plt_txt_x, 0.3+i*0.03), xycoords="figure fraction")
            else:
                plt.annotate("%s - src: %5i" % (e.det, src_phot["aperture_sum"]), xy=(plt_txt_x, 0.3+i*0.03), xycoords="figure fraction")

        return line_temp, dets_temp
        
    name = obs.config["DATA"]["source_name"]
    if conf_file is None:
        global cnf_file
        conf_file = cnf_file
    plt_txt_x = 0.52

    print("Using source name=",name)
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.2, top=0.95, wspace=0.3, hspace=0.3)
    axs = []
        
    line = name + ""
    dets = len(name)*" "

    plt.rcParams.update({'font.size': 6})    
    exposures = obs.exposures.values()
    
    if len(exposures) > 0:
        mp = {obs.exposures[x].det:x for x in obs.exposures}
        # print("mp",mp)
        i,j = 0, 0
        for d in ["pn", "m1", "m2"]:
            
            if d in mp:
                e = obs.exposures[mp[d]]
            else:
                e = namedtuple('exposure', ['det', 'exp_id'])
                e.det, e.exp_id = d, None

            lt, dt = one_panel(i,e, band) # expID , detector-type
            i+=1
            line+=lt
            dets+=dt
    else:
        plt.gca().axis('off')
        plt.annotate("no usable exposure", xy=(0.1,0.9), xycoords="figure fraction", size=14)
        for i, det in enumerate(["pn", "m1", "m2"]):
            fn = path4(o.config, det+"_evt")
            with pyfits.open(fn) as ff:
                plt.annotate(det+": "+ff[0].header["Filter"], xy=(0.12, 0.8-i*0.05), xycoords="figure fraction")
    oi = obs.config["obsID"]
    plt.annotate(name+" ("+oi+")", xy=(plt_txt_x, 0.48), xycoords="figure fraction", size=14)
    bb = band.split(":")
    #print("bb",bb)
    #plt.annotate("energies (eV): %5s - %5s" % (bb[0],bb[1]), xy=(plt_txt_x, 0.45), xycoords="figure fraction", size=8)
    if conf_file: plt.annotate(conf_file, xy=(plt_txt_x, 0.425), xycoords="figure fraction", size=6, color='0.5')
    if len(exposures)>0: plt.annotate(line, xy=(plt_txt_x, 0.15), xycoords="figure fraction")
    plt.annotate(dets, xy=(plt_txt_x, 0.13), xycoords="figure fraction", color='0.5')


if __name__ == "__main__":
    import argparse    
    parser = argparse.ArgumentParser(
                    prog='images_pdf',
                    description='create pdf for specific source',
                    epilog='')
    # parser.add_argument('target', help  ='Name of target or \'all\'')   
    parser.add_argument('conf',nargs=1, help  ='xmmpy-config file')   
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help  ='verbosity')   
    args = parser.parse_args()
    #print(    args.conf[0])
    #print(" verbose: ",args.verbose)
    cname = glob.glob(args.conf[0])
    #print(cname)
    if args.verbose:
        verbosity=10
    else:
        verbosity = 1
    assert len(cname) == 1
    cnt=0
    max_cnt = 1e6
    
    global conf_file
    cnf_file = cname[0]
    
    obs = Obs(conf_file=cname[0])
    ofn = path4(obs.config, "image_pdf")

    with PdfPages(ofn) as pdf:
        for b in obs.config["IMAGES"]["energies"]:
            one_page(obs, band=b, verbose=verbosity)
        
        pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    print("\npdf: ",ofn)
    print(" or: ",os.path.relpath(ofn))


