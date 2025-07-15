Documentation of the class and scripts
--------------------------------------

.. _Time:

Time
-----
The 'time since start of obs' refers to the pn-detector.

An overview of the different times can be shown by running::

    p12 $xmmpy/bin/xmm_time_tools.py CONF-FILE

Internally, all times are expressed in Telescope seconds (=cxcsec).

Add new data products
------------------------
Steps needed to add new data products:
  1) Add suitable part to config:
      :meth:`xmmpy.etc.default_config`
  2) Make sure that config is adjusted for specific source:
      :meth:`xmmpy.etc.update_source_in_config`
  3) Add data product to shell scripts  (may require to add new methods to :class:`xmmpy.Obs` or its subordinate class/methods)
      :meth:`xmmpy.Obs.shell_scripts`
      
(May need to add pathes to :meth:`xmmpy.etc.path4`)



Functions
---------

.. automethod:: xmmpy.etc.path4

|

.. automethod:: xmmpy.etc.default_config

|

.. automethod:: xmmpy.etc.update_source_in_config

|

.. automethod:: xmmpy.etc.rewrite_config

|

.. automethod:: xmmpy.etc.update_value_in_config


.. _scripts:

Scripts
-------

.. automethod:: xmmpy.bin.page4obs


Classes
-------

.. autoclass:: xmmpy.Obs
   :members:
   :special-members:
   :exclude-members: __dict__,__weakref__
   
|
|
   
.. autoclass:: xmmpy.Exposure
   :members:
   :special-members:
   :exclude-members: __dict__,__weakref__   
