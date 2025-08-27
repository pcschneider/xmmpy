# from xmmpy import Obs
import xmmpy.obstools.xmm_obs as xo
from xmmpy.etc import path4
import os

class XMM_Observation_for_xspec:
    """
    Provides two options to generate a string for loading data in xspec, from an xmmpy-config file or using  glob

    Example
    -------
    cnf_file = 'xmmpy-config.conf
    xo = XMM_Observation_for_xspec(cnf_file)
    xs = xo.data_string()
    print(xs)
    """
    def __init__(self, config_fn=None):
        self.epochs_spec_fnames = None
        if config_fn is not None:
            self.obs = xo.Obs(conf_file = config_fn)
            self.data_string_from_xmmpy_conf(config_fn)

    def data_string_from_glob(self, dir, gstr, det_map={"m1":"m1", "m2":"m2","pn":"pn"}, verbose=1):
        """
        Generates xspec data-string using 'dir+"/"+gstr' and 'det_map'
        """
        spec_fnames = glob.glob(dir+"/"+gstr)
        det = det_map.keys()
        epochs = {}
        for fn in spec_fnames:
            for d in det:
                if det_map[d] in fn:
                    if verbose>2: print(fn, " is ", d)
                    bn = os.path.basename(fn)
                    prfx_idx0 = bn.rfind(det_map[d]) + len(det_map[d])
                    prfx_idx1 = bn.rfind(".fits")
                    epoch_ident = bn[prfx_idx0:prfx_idx1]
                    if verbose>3: print("bn",bn, epoch_ident)
                    if epoch_ident in epochs:
                        epochs[epoch_ident][d] = fn
                    else:
                        epochs[epoch_ident] = {d:fn}
        if verbose>2: print()
        if verbose>1:
            for e in epochs:
                print(e, epochs[e])
                print()
        print("epochs: ",epochs.keys())
        self.epochs_spec_fnames = epochs
        return epochs


    def data_string_from_xmmpy_conf(self, conf_fn=None):
        """
        String for loading data in xspec from the spec files as defined in the xmmpy config file
        """
        if conf_fn is None and self.epochs_spec_fnames is None: raise Exception("No config filename provided.")
        elif conf_fn is not None:
            self.epoch_dict(config_fn = conf_fn)    
        xstr = ""
        if self.epochs_spec_fnames is None: raise Exception("No epoch data, try to run \'epoch_dict(fn)\' to retrieve filenames.")
        for i, e in enumerate(self.epochs_spec_fnames):
            for j, d in enumerate(self.epochs_spec_fnames[e]):
                tmp = str("%i:%i %s " % (i+1, j+1, self.epochs_spec_fnames[e][d]))
                xstr+=tmp    
        # print(xstr)
        return xstr


    def data_string(self, keys=None):
        """
        String for loading data in xspec
        """
        if self.epochs_spec_fnames is None: raise Exception("No filenames provided. Use xmmpy-config file or glob to discover spectra.")
        # elif conf_fn is not None:
        #     self.epoch_dict(config_fn = conf_fn)
        xstr = ""
        if self.epochs_spec_fnames is None: raise Exception("No epoch data, try to run \'epoch_dict(fn)\' to retrieve filenames.")
        for i, e in enumerate(self.epochs_spec_fnames):
            if keys is not None and e not in keys: continue
            for j, d in enumerate(self.epochs_spec_fnames[e]):
                tmp = str("%i:%i %s " % (i+1, j+1, self.epochs_spec_fnames[e][d]))
                xstr+=tmp
        return xstr

    def epoch_dict(self, config_fn=None, verbose=1):
        """
        Construct dictionary of spec_fn for epochs
        """
        if config_fn is not None:
            o = xo.Obs(conf_file=config_fn)
            self.obs = o
        else: 
            o = self.obs
        specs = o.spec_files()
        # print(specs)
        
        epochs = {}
        for d in specs:
            # print(d, specs[d])
            spec_fn_prefix =str( path4(o.config, which=d+"_bin"))
            idx = spec_fn_prefix.rfind(".fits")
            if verbose>2: print(spec_fn_prefix, spec_fn_prefix[0:idx])
            for fn in specs[d]:
                idx_fits = str(fn).rfind(".fits")
                postfix = str(fn)[idx:idx_fits]
                if verbose>2: print(fn, os.path.exists(fn), postfix)
                if postfix not in epochs: epochs[postfix] = {d:fn}
                else: epochs[postfix][d] = fn
                
            # print()
        if verbose>2:
            for e in epochs:
                print(e,epochs[e].keys())
        self.epochs_spec_fnames = epochs
        return epochs