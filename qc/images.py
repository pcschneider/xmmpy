from regions import Regions, CirclePixelRegion, Region
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm
import matplotlib
import numpy as np
import copy
from pathlib import Path
from os.path import basename, abspath
import pathlib
import glob
from matplotlib.pyplot import cm
from ..etc import read_config, path4, Surrounding


def ax4image_and_sourceregion(image_fn, region_fn, bkg_region_fn=None, r_scaling=5, idx=1, single_ax=False, fig=None):
    """
    Parameters
    ----------
    bkg_region_fn : str, Path, or list thereof
    r_scaling : float
        Describes how much larger the image will be compared to the radius (does NOT affect how the regions are plotted)
    single_ax : boolean
        If just one axis shall be on the figure
    """
    def add_patch(fn):
        """
        NOTE: Plots only the first region entry in the file
        """
        with pyfits.open(fn) as rff:
            if rff[1].data[0][0].lower() == "circle":
                x, y, r = rff[1].data[0][1], rff[1].data[0][2], rff[1].data[0][3]
                pth = Circle((x,y), r, facecolor='none', edgecolor='r', transform=ax.get_transform('world')) 
                ax.add_patch(pth)
            else:            
                raise Exception("Only circles allowed at the moment.")
        return x,y,r, pth

    ff = pyfits.open(image_fn)
    wcs = WCS(ff[0].header, key="L")
    if fig is not None:
        #print("Using fig")
        if single_ax:
            idx = 1
            ax = fig.add_subplot(1,1, idx, projection=wcs)
        else:
            ax = fig.add_subplot(2,2, idx, projection=wcs)
        
    else:    
        if single_ax:
            idx = 1
            ax = plt.subplot(1,1, idx, projection=wcs)
        else:
            ax = plt.subplot(2,2, idx, projection=wcs)
            
    mcmap = copy.copy(matplotlib.cm.get_cmap('jet')) 
    mcmap.set_bad((0,0,0))
    im = ax.imshow(ff[0].data, norm=LogNorm(vmax=10), cmap = mcmap)
    cbar = plt.colorbar(im, ticks=[1,5,10])
    cbar.set_label("Counts")
    cbar.set_ticks([1,2,3,4,6,10])
    cbar.set_ticklabels(["1","2","3","4","6","10"])
    
    sx, sy, sr, _ = add_patch(region_fn)
    if bkg_region_fn:
        if isinstance(bkg_region_fn, str) or isinstance(bkg_region_fn, Path):
            bkg_region_fn = [bkg_region_fn]
        
        color = iter(cm.rainbow(np.linspace(0, 2, len(bkg_region_fn))))    
        for c, fn in zip(color, bkg_region_fn):
            bx, by, br, p = add_patch(fn)
            #print("bkg", bx, by, br)
            p.set(label=basename(str(fn)))
            p.set(edgecolor=c)
            #ax.annotate(fn, xy=(bx, by+br), transform=ax.get_transform('world'))
        dx, dy, dr = 0.5*(sx+bx), 0.5*(sy+by), br
        tlim1 = (dx-5*dr, dy-5*dr)
        tlim2 = (dx+5*dr, dy+5*dr)        
    else:
        tlim1 = (sx-r_scaling*sr, sy-r_scaling*sr)
        tlim2 = (sx+r_scaling*sr, sy+r_scaling*sr)        

    ll = ax.get_xticks()[::3]
    # print(ll)
    # nllx = [str("%i" % (x-ll[len(ll)//2])) for x in ll[::3]]
    # print("x ll",ll, ll[len(ll)//2], nllx)
    ax.set_xticks(ll, [str("%i" % x) for x in ll])#, nllx)

    ttx = ax.get_transform("world")
    cc0 = ttx.transform(tlim1)
    cc0 = ax.transData.inverted().transform(cc0)
    cc1 = ttx.transform(tlim2)
    cc1 = ax.transData.inverted().transform(cc1)
    ax.set_xlim(cc0[0], cc1[0])
    ax.set_ylim(cc0[1], cc1[1])
    # ax.annotate(ff[0].header["Filter"], xy=(0.2, 0.2), color='white', xycoords="axes fraction")
    ax.set_xlabel("x (physical)")
    ax.set_ylabel("y (physical)")
    ax.legend()
    return ax



def ax4image_and_regions(image_fn, region_fns, r_scaling=2, ax=None):
    """
    Parameters
    ----------
    bkg_region_fn : str, Path, or list thereof
    r_scaling : float
        Describes how much larger the image will be compared to the radius (does NOT affect how the regions are plotted)
    single_ax : boolean
        If just one axis shall be on the figure
    """
    def add_patches(fn):
        # print("fn", fn)
        rr = []
        if isinstance(fn, str) or isinstance(fn, pathlib.PosixPath):
            regs = Regions.read(fn)
        elif isinstance(fn, Region):            
            # print("is region", type(fn))
            regs = [fn]
        else:
            print(rr)
            raise Exception("XXX")
        for r in regs:
            art = r.as_artist()                
            if isinstance(r, CirclePixelRegion):
                # print("meta", r.meta, r.visual)
                ls = r.visual["linestyle"] if "linestyle" in r.visual else '-'
                clr = r.visual["color"] if "color" in r.visual else 'r'
                txt = r.meta["text"] if "text" in r.meta else None
                pth = Circle((r.center.x, r.center.y), r.radius,  linestyle=ls, facecolor='none', edgecolor=clr, transform=ax.get_transform('world'))
                # print(pth.get_transform())
                ax.add_patch(pth)
                rr.append({"x":r.center.x, "y":r.center.y, "r":r.radius, "text":txt})
            else:
                raise Exception("Don't know how to handle region of type ", r)
            if hasattr(r, "origin") and r.origin!="Simbad" and r.origin !="sourcelist": pth.set(label=txt)
            elif hasattr(r, "filename"): pth.set(label=r.filename)
            if "text" in r.meta:
                # if hasattr(r, "origin"):
                #     print("origin:",r.origin)
                if hasattr(r, "origin") and r.origin=="Simbad":
                    an = ax.annotate(r.meta["text"], xy=(r.center.x, r.center.y-r.radius), color=clr, xycoords=ax.get_transform('world'), ha='center', va='top')
                elif hasattr(r, "origin") and r.origin=="sourcelist":
                    an = ax.annotate(r.meta["text"], xy=(r.center.x, r.center.y+r.radius), color=clr, xycoords=ax.get_transform('world'), ha='center', va='bottom')
                else:
                    an = ax.annotate(r.meta["text"], xy=(r.center.x, r.center.y+r.radius), color=clr, xycoords=ax.get_transform('world'), ha='center', va='bottom')
                # print(an, dir(an), an.get_text(), an.get_transform() )
        return rr

    ff = pyfits.open(image_fn)

    # plt.imshow(ff[0].data, norm=LogNorm(vmax=30))
    # plt.show()
    wcs = WCS(ff[0].header, key="L")
    wcs0 = WCS(ff[0].header)
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(1,1, 1, projection=wcs)
    else:
        raise Exception("Providing an axis is currently not supported; and likely will never be due to axis-transformtion.")
    mcmap = copy.copy(matplotlib.cm.get_cmap('jet')) 
    mcmap.set_bad((0,0,0))
    # print(np.shape(ff[0].data))
    cutout = ff[0].data
    im = ax.imshow(cutout, norm=LogNorm(vmax=10), cmap = mcmap)
    cbar = plt.colorbar(im, ticks=[1,5,10])
    cbar.set_label("Counts")
    cbar.set_ticks([1,2,3,4,6,10])
    cbar.set_ticklabels(["1","2","3","4","6","10"])
    
    radii = []
    cx, cy = [], []
    for fn in region_fns:    
        p = add_patches(fn)
        for r in p:
            cx.append(r["x"]), cy.append(r["y"])
            radii.append(r["r"])
    
    tlim1 = (min(cx)-r_scaling*max(radii), min(cy)-r_scaling*max(radii))
    tlim2 = (max(cx)+r_scaling*max(radii), max(cy)+r_scaling*max(radii))
    
    ll = ax.get_xticks()[::3]
    ax.set_xticks(ll, [str("%i" % x) for x in ll])#, nllx)

    ttx = ax.get_transform("world")
    cc0 = ttx.transform(tlim1)
    cc1 = ttx.transform(tlim2)
    cc0i = ax.transData.inverted().transform(cc0)
    cc1i = ax.transData.inverted().transform(cc1)
    if False:
        print("cc0", cc0)
        print("cc1", cc1)
        print("tlim: ",tlim1, tlim2)
        print("wcs0", wcs0.all_world2pix([tlim1[0]], [tlim1[1]],1))
        print("wcs", wcs.all_world2pix([tlim1[0]], [tlim1[1]],1))
        print("wcs0", wcs0.all_pix2world([tlim1[0]], [tlim1[1]],1))
        print("wcs", wcs.all_pix2world([tlim1[0]], [tlim1[1]],1))
        print("ccXi:", cc0i, cc1i)
        print("im limits:", cc0i[0], cc1i[0], cc0i[1], cc1i[1])
    im_lims = np.array([cc0i[0],cc1i[0], cc0i[1],cc1i[1]]).astype(int)
    # print(r'{im_lims}')
    cutout = ff[0].data[im_lims[2]: im_lims[3], im_lims[0]:im_lims[1]]
    mx = np.max(cutout)
    ax.set_xlim(cc0i[0], cc1i[0])
    ax.set_ylim(cc0i[1], cc1i[1])
    ax.annotate("%6.2f ks" % (ff[0].header["EXPOSURE"]/1e3), xy=(0.1, 0.14), color='white', xycoords="axes fraction")
    ax.annotate(ff[0].header["Filter"], xy=(0.1, 0.1), color='white', xycoords="axes fraction")
    ax.set_title(basename(image_fn))
    ax.set_xlabel("x (physical)")
    ax.set_ylabel("y (physical)")
    ax.legend()
    im.set(norm=LogNorm(vmin=1, vmax=1.2*mx))
    if False:
        plt.title("Check image")
        plt.figure()
        plt.imshow(cutout, norm=LogNorm(vmax=0.8*mx), cmap = mcmap, origin='lower')
        plt.xlabel("pix")
        plt.ylabel("pix")
    return ax


def check_one_config(fn, band="0.5-2.0keV"):
    cnf = read_config(fn)
    sur = Surrounding(fn)
    dr = str(path4(cnf, which='odata'))
    print(cnf["DATA"]["detectors"])
    image_fn_mapper = {"pn":"_EPN_", "m1":"_EMOS1_", "m2":"_EMOS2_"}
    for d in cnf["DATA"]["detectors"]:
        glob_str = dr+"/"+"*"+cnf["obsID"]+image_fn_mapper[d]+"*_"+band+".fits"
        image_fn = glob.glob(glob_str)
        if len(image_fn)!=1:
            raise Exception("Did not find exactly one single matching image for "+str(glob_str))
        src_reg_fn = path4(cnf, "src_reg")
        bkg_reg_fn = path4(cnf, "bkg_"+d+"_reg")
        print(fn, glob_str)
        print(image_fn[0], src_reg_fn, bkg_reg_fn)
        regs = Regions.read(src_reg_fn)
        for r in regs:
            r.filename=src_reg_fn            
        for item in Regions.read(bkg_reg_fn):
            item.filename=bkg_reg_fn
            regs.append(item)
        for item in sur.regions(coord_system="physical"):
            regs.append(item)
        # regs.append 
        ax4image_and_regions(image_fn[0], regs)
        plt.show()
        print()

if __name__ == "__main__":    

    if False:    
        ifn = "/home/majestix/hdd/test/xmmpy/data2/0865070101/odata/3862_0865070101_EPN_S003_ImagingEvt_image_0.2-1.0keV.fits"
        rfn = "data2/0784240701/odata/HD_120690_reg_src.fits"
        bfn = "data2/0784240701/odata/HD_120690_reg_bkg_pn.fits"
        ax = ax4image_and_sourceregion(ifn, rfn, bkg_region_fn=bfn)
        plt.show()

    x = check_one_config("/home/majestix/hdd/Xdata/trex/data/0865040801/HL_Tau_0865040801.conf")

    #regions = Regions.read("/home/majestix/hdd/test/xmmpy/data2/0784240701/odata/HD_120690_reg_src.fits")
