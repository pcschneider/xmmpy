from regions import Regions
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from matplotlib.patches import Circle
from matplotlib.colors import LogNorm
import matplotlib
import numpy as np
import copy
from pathlib import Path
from os.path import basename
from matplotlib.pyplot import cm

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
        print("Using fig")
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
            #print("p",p)
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

if __name__ == "__main__":    
        
    ifn = "/home/majestix/hdd/test/xmmpy/data2/0865070101/odata/3862_0865070101_EPN_S003_ImagingEvt_image_0.2-1.0keV.fits"
    rfn = "data2/0784240701/odata/HD_120690_reg_src.fits"
    bfn = "data2/0784240701/odata/HD_120690_reg_bkg_pn.fits"
    ax = ax4image_and_sourceregion(ifn, rfn, bkg_region_fn=bfn)

    plt.show()


    #regions = Regions.read("/home/majestix/hdd/test/xmmpy/data2/0784240701/odata/HD_120690_reg_src.fits")
