Regions
=========

The XMM-SAS can use fits-style region files. 

.. automethod:: xmmpy.etc.fits_region

|   
 
.. automethod:: xmmpy.etc.aperture_from_fits


Some useful information about XMM-Newton
==========================================

RGS - Dispersion 
-----------------
Approx 10 mA per pixel, depending on wavelength


Times used by XMM-Newton
---------------------------

The MJDREF value is 50814.0 and pertains to XMM-Newton seconds = 0. Or 63.184 via::

    xmmtimeconv time='50814.0' format=MJD
    
Description here: 
  https://xmm-tools.cosmos.esa.int/external/xmm_user_support/documentation/uhb/reftime.html

Some checks can be performed by::

    from astropy.io import fits as pyfits
    from astropy.time import Time

    ff = pyfits.open("~/hdd/TTs/ULYSSES/XMM/RU_Lup/new/new_0882060501/odata/pn.fits")
    do = Time(ff[0].header["DATE-OBS"])
    dt = do - Time(50814.0, format='mjd') # 50814 equals 1998-01-01T00:00:00.000, and equlas MJDREF
    print("Seconds elapsed since 1998-01-01T00:00:00: %10.2f" % dt.sec)
    print("                   elapsed Seconds + 63.2: %10.2f" % (dt.sec + 63.184))
    print(" Time stamp first photon:                  %10.2f" % min(ff[1].data["TIME"]))
    print(" CXC seconds for DATE-OBS:                 %10.2f" % (do.cxcsec))
