Some config-file entries
=========================

* ``[SPECTRA][time_bins]`` -- must be in ["none", int, [[a0, a1], [b0, b1], ...]
      If "none", the full observation will be used (equal to a parameter value of 1)
      
      If int, the observation will be split in ``time_bins`` chunks/intervals. Postfix will be ``_Xks_binX``
      
      If list, the specific time intervals will be used. The unit is ks after start of exposure. Postfix will be ``t0_t1ks``.

.. * as      
         

Description of high-level scripts
======================================

The following scripts can be used to process data:

.. _xmm_retrieve:

xmm_retrieve
------------

Download data from the XMM-Newton archive and performs the first reduction steps. If run as::

  xmm_retrieve.py obsID directory --script=xt.sh

it generates a script (``xt.sh``), which can be sourced to download the data and to run the script::

  directory/reduce_odf_{obsID}.sh
  
which itself calls::

  directory/make_xmm.sh    # <- copy of the general reduction script
  directory/{obsID}/sas_{obsID}.sh   # <- usual file to set SAS_ODF etc.

to generate::

  directory/{obsID}/odata # with the event-files and images

and the conf-file::

  directory/{obsID}/xmmpy0109060301.conf # standard conf-file
  
as well as the log-files::

  directory/reduce_odf_{obsID}.log
  
and log-files for some SAS-tasks::

  directory/{obsID}/cifbuild.log
  directory/{obsID}/odfingest.log
  directory/{obsID}/emproc.log
  directory/{obsID}/epproc.log
  
Essentially, the script downloads the data and runs ``make_xmm.sh``.

(Requires SAS and partially runs in bash.)

.. _xmm_source_regions:

xmm_source_regions
------------------

Generates source and background regions. Run as::

  xmm_source_regions.py directory/{obsID} "source name"  --script=xr.sh
  
The source and background regions will be::

  directory/{obsID}/odata/*reg*.fits

The script looks up the source in Simbad, propagates the sky position to the observation date and generates *one* source region (``*_reg_src.fits``) as well as individual background regions for pn, m1, and m2 (``*_reg_bkg_{det}.fits``). Source regions are circles. 

Log-file: ``xmmpy.log``

.. note::

  This script requires a directory as input, although it really uses the xmmpy{obsID}.conf-file to find files.
  
.. _xmm_source_products:  
  
xmm_source_products
-------------------

Run::

  xmm_source_products.py directory/{obsID} "source name"  --script=sp.sh

For example::

  xmm_source_products.py --config=./RU_Lup_0882060501.conf --script="new_specs.sh"

You can also adjust the desired source products, e.g., light curves only as in the following example via::

  xmm_source_products.py --config=./RU_Lup_0882060501.conf --p=lc --script="new_specs.sh"
  
The optional source products are "rgs","spec","lc", "evt".
  

xmm_full_process
-----------------

This is a bash-script, which combines :ref:`xmm_retrieve`, :ref:`xmm_source_regions`, and :ref:`xmm_source_products`.

