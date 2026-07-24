# from xmmpy import Obs
import xmmpy.obstools.xmm_obs as xo
from xmmpy.etc import path4
import numpy as np
import pandas as pd
import os
try:
    import xspec # requires 'heainit'
except: 
    print("xmmpy - No xspec; run '> heainit' first")

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
        """Initialize; if config_fn is given, load epoch/spectrum filenames from the xmmpy config."""
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


    def data_string(self, keys=None, dets=None):
        """
        String for loading data in xspec

        Parameters
        ----------
        keys : list of epoch identifiers, "" if only one epoch
        dets : list of detectors, must be in ["pn", "m1", "m2"]
        """
        if self.epochs_spec_fnames is None: raise Exception("No filenames provided. Use xmmpy-config file or glob to discover spectra.")
        # elif conf_fn is not None:
        #     self.epoch_dict(config_fn = conf_fn)
        xstr = ""
        if self.epochs_spec_fnames is None: raise Exception("No epoch data, try to run \'epoch_dict(fn)\' to retrieve filenames.")
        for i, e in enumerate(self.epochs_spec_fnames):
            if keys is not None and e not in keys: continue
            for j, d in enumerate(self.epochs_spec_fnames[e]):
                if not os.path.exists(self.epochs_spec_fnames[e][d]): 
                    print("File does not exist: ",self.epochs_spec_fnames[e][d])
                    continue
                if dets is None:
                    tmp = str("%i:%i %s " % (i+1, j+1, self.epochs_spec_fnames[e][d]))
                    xstr+=tmp
                elif d in dets:
                    tmp = str("%i:%i %s " % (i+1, j+1, self.epochs_spec_fnames[e][d]))
                    xstr+=tmp
        if xstr =="": return None
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
    


def untie_epochs(ref_comp=1, ignore_names=[]):
    """
    Unfreeze all parameters across all data groups.

    Parameters
    ----------
    ref_comp : int
        Untie links to this component
    ignore_names : list of str
        Ignore parameters with these names

    Returns
    -------
    list of int
        Parameter indices that were unfrozen.
    """

    def _comp4link(link):
        for i0 in range(xspec.AllData.nGroups):
            i = i0+1
            idx0, idx1 = xspec.AllModels(i).startParIndex, xspec.AllModels(i).startParIndex+xspec.AllModels(i).nParameters-1
            if idx0 <= int(link) <= idx1:
                return i
        return None

    i = 0
    idx0, idx1 = xspec.AllModels(i+1).startParIndex, xspec.AllModels(i+1).startParIndex+xspec.AllModels(i+1).nParameters-1
    print("nGroups",xspec.AllData.nGroups, "base index", idx0, idx1)
    for i_tmp in range(xspec.AllData.nGroups-1):
        i = i_tmp + 2
        print(i)
        p0, p1 = xspec.AllModels(i).startParIndex, xspec.AllModels(i).startParIndex+xspec.AllModels(i).nParameters-1
        print(i, " parameters: ",p0, p1)
        for p in range(p1-p0+1):
            print(i, p)
            link = xspec.AllModels(i)(p+1).link
            l0 = pIdx4link(link)
            comp = _comp4link(l0)
            print(p, xspec.AllModels(i)(p+1).name, "link: ", link, l0, " component", comp)
            if comp == ref_comp: 
                if xspec.AllModels(i)(p+1).name in ignore_names:
                    print("Ignoring", xspec.AllModels(i)(p+1).name, " for group", i, "parameter", p+1)
                    continue
                xspec.AllModels(i)(p+1).link = ""
                ref_is_free = xspec.AllModels(i)(p+1).frozen == 0
                print("ref_is_free", ref_is_free, "(comp: ",comp, "param: ",l0,")")
                if ref_is_free:
                    if xspec.AllModels(i)(p+1).frozen:
                        xspec.AllModels(i)(p+1).frozen = 0

                    

            

def untie_based_on_name(xm, xd, name="norm", ignore_links=False,  model_comps = None, first_only=False):
    """
    Unfreeze all parameters matching `name` across all data groups.

    Parameters
    ----------
    xm : xspec.AllModels
    xd : xspec.AllData
    name : str
        Parameter name to match (e.g. "norm").
    ignore_links : bool
        If True, clear the link string before unfreezing.
    model_comps : list or None
        Restrict to these 1-based component indices; None means all components.
    first_only : bool
        Stop after the first matching parameter is found.

    Returns
    -------
    list of int
        Parameter indices that were unfrozen.
    """
    found = False
    indices = []
    for i in range(xd.nGroups):
        idx0, idx1 = xm(i+1).startParIndex, xm(i+1).startParIndex+xm(i+1).nParameters-1
        idx = idx0
        comp_idx = 1
        for j, comp_name in enumerate(xm(i+1).componentNames):
            print("  Data group", i, " model component:",j+1, "(",comp_name,")")#idx0, comp_idx, comp_name)#, xspec.AllModels(i+1))#(j+1().nParameters)
            comp = getattr(xm(i+1), comp_name)
            # print(comp)
            for k, pn in enumerate(comp.parameterNames):
                # print(k, pn)
                param = xm(i+1)(comp_idx) 
                
                print("     %12s (%3i), link: %7s,  frozen: %5s," % (pn,idx,param.link, param.frozen), end="")
                print(" values:", param.values)
                if model_comps is not None and j+1 not in model_comps: 
                    pass
                else:
                    if pn == name: 
                        if param.link != "" and ignore_links:
                            param.link = ""
                        print("\'%s\'" % param.link)                        
                        if param.link == "": 
                            param.frozen = False
                            found = True
                            indices.append(idx)
                if first_only and found: return indices
                idx+=1
                comp_idx+=1
                # print("param:", param)
    return indices  

def link_based_on_names(xm, xd, names=["He"], model_comps = None):
    """Call link_based_on_name for each name in `names`."""
    for n in names:
        link_based_on_name(xm, xd, name=n, model_comps=model_comps)

def link_based_on_name(xm, xd, name="He", model_comps = None):
    """
    Link all occurrences of parameter `name` to its first occurrence across data groups.

    Parameters
    ----------
    xm : xspec.AllModels
    xd : xspec.AllData
    name : str
        Parameter name to link (e.g. abundance element symbol).
    model_comps : list or None
        Restrict to these 1-based component indices; None means all components.
    """
    for i in range(xd.nGroups):
        idx0, idx1 = xm(i+1).startParIndex, xm(i+1).startParIndex+xm(i+1).nParameters-1
        idx = idx0
        comp_idx = 1
        first = None
        for j, comp_name in enumerate(xm(i+1).componentNames):
            print("  Data group", i, " model component:",j+1, "(",comp_name,")")#idx0, comp_idx, comp_name)#, xspec.AllModels(i+1))#(j+1().nParameters)
            comp = getattr(xm(i+1), comp_name)
            # print(comp)
            for k, pn in enumerate(comp.parameterNames):
                # print(k, pn)
                param = xm(i+1)(comp_idx)

                print("     %12s (%3i), link: %7s,  frozen: %5s," % (pn,idx,param.link, param.frozen), end="")
                print(" values:", param.values)
                if model_comps is not None and j+1 not in model_comps:
                    pass
                else:
                    if pn == name:
                        if first is not None:
                            link_str = "%s" % first
                            print("idx", idx, " link_str", link_str, "first: ",first)
                            param.link = link_str
                            # param.link = ""
                        else:
                            first = idx
                        print("\'%s\'" % param.link)
                        # param.frozen = False
                idx+=1
                comp_idx+=1

def link_pattern(xm, xd, ref="", dependend=[], model_comps = None):
    """
    Link `ref` parameter across data groups to its first occurrence, and link all
    parameters in `dependend` to that same reference index.

    Parameters
    ----------
    xm : xspec.AllModels
    xd : xspec.AllData
    ref : str
        Name of the reference parameter (linked to itself across groups).
    dependend : list of str
        Parameter names to link to the reference parameter's first index.
    model_comps : list or None
        Restrict to these 1-based component indices; None means all components.
    """
    def get_first_idx(tmp_name, grp):
        idx0, idx1 = xm(grp+1).startParIndex, xm(grp+1).startParIndex+xm(grp+1).nParameters-1
        idx = idx0
        comp_idx = 1
        first = None
        for j, comp_name in enumerate(xm(i+1).componentNames):
            print("  Data group", i, " model component:",j+1, "(",comp_name,")")#idx0, comp_idx, comp_name)#, xspec.AllModels(i+1))#(j+1().nParameters)
            comp = getattr(xm(i+1), comp_name)
            # print(comp)
            for k, pn in enumerate(comp.parameterNames):
                # print(k, pn)
                param = xm(i+1)(comp_idx)

                print("     %12s (%3i), link: %7s,  frozen: %5s," % (pn,idx,param.link, param.frozen), end="")
                print(" values:", param.values)
                if model_comps is not None and j+1 not in model_comps:
                    pass
                else:
                    if pn == tmp_name:
                        if first is None:
                            return idx
                idx+=1
                comp_idx+=1

    for i in range(xd.nGroups):
        idx0, idx1 = xm(i+1).startParIndex, xm(i+1).startParIndex+xm(i+1).nParameters-1
        idx = idx0
        comp_idx = 1
        first = get_first_idx(ref, i)
        for j, comp_name in enumerate(xm(i+1).componentNames):
            print("  Data group", i, " model component:",j+1, "(",comp_name,")")#idx0, comp_idx, comp_name)#, xspec.AllModels(i+1))#(j+1().nParameters)
            comp = getattr(xm(i+1), comp_name)
            # print(comp)
            for k, pn in enumerate(comp.parameterNames):
                # print(k, pn)
                param = xm(i+1)(comp_idx)

                print("     %12s (%3i), link: %7s,  frozen: %5s," % (pn,idx,param.link, param.frozen), end="")
                print(" values:", param.values)
                if model_comps is not None and j+1 not in model_comps:
                    pass
                else:
                    if pn == ref and idx != first:
                        link_str = "%s" % first
                        print("idx", idx, " link_str", link_str, "first: ",first)
                        param.link = link_str
                            # param.link = ""
                    if pn in dependend:
                        link_str = "%s" % first
                        print("idx", idx, " link_str", link_str, "first: ",first)
                        param.link = link_str

                    print("\'%s\'" % param.link)

                        # param.frozen = False
                idx+=1
                comp_idx+=1

def FIP_pattern(xm, xd, model_comps=None):
    """
    Apply FIP-bias linking: group elements by First Ionization Potential and link
    abundances within each group to a single reference element.

    Low-FIP (ref Fe): Mg, Si, Al, Ca, Ni — mid-FIP (ref C): N, O, Si, S — high-FIP (ref Ne): Ar.

    Returns list of reference element names ['C', 'Fe', 'Ne'].
    """
    link_pattern(xm, xd, ref = "C", dependend = ["N", "O", "Si", "S"])

    link_pattern(xm, xd, ref = "Fe", dependend = ["Mg", "Si", "Al", "Ca", "Ni"])
    link_pattern(xm, xd, ref = "Ne", dependend = ["Ar"])
    return ["C", "Fe", "Ne"]    

def value_tuple():
    """Placeholder — not yet implemented."""
    pass


def value_str(value, error, default_fmt="%5.2f"):
    """
    Format a value with asymmetric error bars as a LaTeX string.

    Returns plain ``"V.VV"`` if `error` is NaN or both bounds are zero,
    otherwise ``"V.VV$^{+hi}_{-lo}$"``.
    """
    out = str(default_fmt % value)
    if type(error) == float and np.isnan(error): return out
    error_tuple = error[1:-1].split(",")
    lo, hi = float(error_tuple[0]), float(error_tuple[1])
    tmp = "$^{+"+default_fmt+"}_{-"+default_fmt+"}$"
    out+=str(tmp % (hi-value, value-lo))
    if lo == 0.0 and hi == 0.0:
        out = str("%5.2f" % value)
    return out.replace(" ","")

def pIdx4link(link):
    """Extract the integer parameter index from an xspec link string (e.g. ``'= p42'`` → 42)."""
    if type(link) == float and np.isnan(link): return None
    idx = link.rfind("p")
    nr = link[idx+1:]
    # print(nr)
    return int(nr)

def range4param(df, param_name):
    """Return the min/max values of a parameter."""
    return df[df["parameter_name"] == param_name]["value"].min(), df[df["parameter_name"] == param_name]["value"].max()

def range4param_idx(df, param_name):
    """Return the min/max values of a parameter."""
    return df[df["parameter_index"] == param_name]["value"].min(), df[df["parameter_index"] == param_name]["value"].max()


def tex4params(ifn, params=None, dg_mapper=None, verbose=1):
    """
    Generate LaTeX table rows from a write_model CSV file.

    One row per data group is produced: the first row is a header with parameter
    names, followed by one data row per group.

    Parameters
    -----------
    ifn : str
        Path to a CSV file written by FitProperties.write() / write_model.
    params : list of str or None
        Parameter names (strings) to include as columns, taken from the
        ``parameter_name`` column of the CSV.  The same name may appear
        multiple times if the model has several components with that name
        (e.g. ``"kT"`` in a 3T APEC model) — all matching indices are
        included in order.  Integers (parameter indices) are not supported;
        always pass names.  If None, all unique parameter indices present in
        the first data group are included.
    dg_mapper : dict or None
        Maps data-group integers to the label printed in the first column,
        e.g. ``{1: "0764100501"}``.  Unmapped groups fall back to
        ``"Data group N"``.  If None, all groups use the fallback label.
    verbose : int
        Verbosity level (default 1).

    Returns
    -------
    list of str
        LaTeX row strings.  Element 0 is the header row (column names);
        subsequent elements are one data row per data group.

    Examples
    --------
    # All free parameters, data-group number as row label:
    rows = tex4params("0764100501_spec_params_3T_frozen.csv")

    # Select specific parameters, map group 1 to an ObsID string:
    rows = tex4params(
        "0764100501_spec_params_3T_frozen.csv",
        params=["lg10Flux", "kT", "kT_4", "kT_5"],
        dg_mapper={1: "0764100501"},
    )
    print("\\n".join(rows))
    """
    def _select_fmt(lo, hi, margin=1):
        log_lo_floor = int(np.floor(np.log10(abs(lo)))) + margin
        if log_lo_floor < 0:
            log_lo_floor*=-1
            return "%" + str(log_lo_floor+3) + "."+str(log_lo_floor+1)+"f"
        return "%" + str(log_lo_floor+2) + "."+str(1+margin)+"f"


    if verbose>1: print("Reading file", ifn)
    df = pd.read_csv(ifn, delimiter=',', comment='#')
    df.drop_duplicates(inplace=True )
    rows = []
    bold = "\\textbf{"
    nobold = "{"
    fett = nobold
    dgrps = np.unique(df['data_group'])
    top_row, pop_top = " & ", True
    param_names = {}
    for dg in dgrps:
        gi = np.where(df["data_group"] == dg)[0]
        if len(gi) == 0: continue
        if dg_mapper is not None:
            if dg in dg_mapper.keys():
                row = fett +"%s} & " % dg_mapper[dg]
            else:
                row = fett+"Data group %i} & " % dg
        else:
            row = fett+"Data group %i} & " % dg
        # row = " & "
        if params is None:
            param_indices = np.unique(df["parameter_index"].iloc[gi])    
            params = [df["parameter_name"].iloc[gi[np.where(df["parameter_index"] == p)[0]]].values for p in param_indices]
        elif type(params) == list:
            param_indices = []
            if all(type(p)==str for p in params):
                for p in params:
                    pi = np.where(df["parameter_name"].iloc[gi] == p)[0]
                    tmp = df["parameter_index"].iloc[gi[pi]]
                    if verbose>1:print("index for ", p, " in data group", dg, ":", pi, " tmp:", tmp)
                    for t in tmp: param_indices.append(t)

        ranges = {} # key: param_index, value: (lo, hi)
        for p in param_indices:
                name = df["parameter_name"].iloc[np.where(df["parameter_index"] == p)[0]].values[0]
                param_names[p] =name
        for p in param_indices:
            ranges[p] = range4param_idx(df, p)
            print("range for param", param_names[p], "/", p," in data group", dg, ":", ranges[p])

        for p in param_indices:
            pi = np.where(df["parameter_index"] == p)[0]
            if len(pi) == 0: 
                print("Parameter index %i not found in data group %i" % (p, dg))
                continue
            if pop_top:
                if verbose>1: print("pop_top, for param_index p", p)
                top_row += df["parameter_name"].iloc[pi[0]] + " & "

            # row += str(p) + " & "
            # row += " & "
            fmt= _select_fmt(*ranges[p])
            print(p, param_names[p], fmt)
            v = df["value"].iloc[pi[0]]
            e = df["error"].iloc[pi[0]]
            row += value_str(v, e, default_fmt=fmt) + " & "
        if pop_top:
            pop_top = False
            top_row = top_row[:-3] + "\\\\\n\\hline"
            rows.append(top_row)
        row = row[:-3] + "\\\\"
        rows.append(row)
    return rows

class FitProperties:
    """
    Snapshot of the current xspec model state.

    Per-parameter attributes (one entry per model parameter across all data groups):
        paridx, datagrp, compidx, compname, parname, link, frozen, value, error
        ``error`` holds the raw xspec tuple ``(lo, hi, flags_str)`` per parameter.

    Per-data-group attributes:
        grp_ids, exposure, rate_net, rate_err, rate_gross, rate_pred

    Fit-statistic scalars:
        statistic, dof, statmethod, expression
    """

    def __init__(self):
        # per-parameter
        self.paridx   = []
        self.datagrp  = []
        self.compidx  = []
        self.compname = []
        self.parname  = []
        self.link     = []
        self.frozen   = []
        self.values    = []
        self.error    = []
        # per-data-group
        self.grp_ids    = []
        self.exposure   = []
        self.rate_net   = []
        self.rate_err   = []
        self.rate_gross = []
        self.rate_pred  = []
        # fit statistics
        self.statistic  = None
        self.dof        = None
        self.statmethod = None
        self.expression = None
        # per-spectrum (one entry per loaded spectrum across all data groups)
        self.specnums   = []   # spectrum number within its data group (1-based)
        self.filenames  = []
        self.ignored    = []
        # flux command output
        self.flux_emin  = None   # lower energy bound used for flux calculation (keV)
        self.flux_emax  = None   # upper energy bound
        self.flux       = {}     # {grp_id: (phot_flux, phot_err, ener_flux, ener_err)}

    def __getitem__(self, key):
        """
        Retrieve parameter records by name and optional data group.

        Key forms
        ---------
        fp["kT"]          – all records whose parameter_name is "kT"
        fp["kT", 1]       – records named "kT" in data group 1
        fp[("kT", 1)]     – same

        Returns
        -------
        list of dict, each with keys:
            paridx, datagrp, compidx, compname, parname,
            link, frozen, value, error
        Raises KeyError when no match is found.
        """
        if isinstance(key, tuple):
            name, grp = key
        else:
            name, grp = key, None
        
        results = [
            dict(paridx=pi, datagrp=dg, compidx=ci, compname=cn,
                 parname=pn, link=lk, frozen=fr, value=va, error=er)
            for pi, dg, ci, cn, pn, lk, fr, va, er in zip(
                self.paridx, self.datagrp, self.compidx, self.compname,
                self.parname, self.link, self.frozen, self.values, self.error)
            if pn == name and (grp is None or dg == grp)
        ]

        if not results:
            raise KeyError(
                "No parameter %r found%s" % (name, "" if grp is None else " in data group %i" % grp))
        return results

    def value(self, key, grp=None):
        """
        Return parameter value(s) by name and optional data group.

        Mirrors the ``__getitem__`` call signature but returns only the
        best-fit value float(s) instead of the full record dicts.

        Parameters
        ----------
        key : str or tuple (name, grp)
            Parameter name, or a ``(name, grp)`` tuple.
        grp : int or None
            Data group (1-based).  Ignored when *key* is already a tuple.

        Returns
        -------
        list of float
            One entry per matching parameter record (e.g. three entries for
            ``kT`` in a 3T model).  Raises ``KeyError`` when nothing matches.

        Examples
        --------
        fp.value("lg10Flux")        # → [float]
        fp.value("kT", 1)           # → [0.2, 0.6, 1.0]  (three apec components)
        fp.value(("kT", 1))         # → same
        """
        if isinstance(key, tuple):
            name, grp = key
        else:
            name = key
        return [r["value"] for r in self[name, grp]]

    @staticmethod
    def calc_errors(errors=True, error_delta=2.706):
        """
        Calculate confidence errors for free, unlinked model parameters in xspec.

        Iterates over all data groups, identifies free and unlinked parameters
        within the requested index range, and calls xspec.Fit.error() for each.
        After this method returns, param.error in xspec memory is up-to-date.

        Parameters
        ----------
        errors : bool or list
            If True, calculate errors for all free parameters.
            If a list of (1-based) parameter indices, only calculate for those.
        error_delta : float
            Delta fit statistic defining the confidence interval
            (default 2.706 corresponds to 90% single-parameter confidence).
        """
        if not errors:
            return

        calc_indices = np.arange(10_000)
        if isinstance(errors, list):
            calc_indices = list(errors)
            errors = True
        elif not isinstance(errors, bool):
            errors = False

        error_calc_strings = {}
        for i in range(xspec.AllData.nGroups):
            idx = xspec.AllModels(i + 1).startParIndex
            comp_idx = 1
            for comp_name in xspec.AllModels(i + 1).componentNames:
                comp = getattr(xspec.AllModels(i + 1), comp_name)
                for pn in comp.parameterNames:
                    param = xspec.AllModels(i + 1)(comp_idx)
                    if len(str(param.link)) == 0 and not param.frozen and idx in calc_indices:
                        error_calc_strings[idx] = f"{error_delta} {idx}"
                    idx += 1
                    comp_idx += 1

        for idx, es in error_calc_strings.items():
            print(f"Calculating error for {es}")
            xspec.Fit.error(es)

    @classmethod
    def collect(cls, errors=False, error_delta=2.706, flux_range=(0.2,2.4)):
        """
        Calculate errors (if requested) and snapshot the current xspec model state.

        Parameters
        ----------
        errors : bool or list
            Passed to calc_errors.
        error_delta : float
            Passed to calc_errors.
        flux_range : tuple
            Energy limits for flux calculation (keV).

        Returns
        -------
        FitProperties
        """
        cls.calc_errors(errors, error_delta)
        fp = cls()

        for i in range(xspec.AllData.nGroups):
            idx0 = xspec.AllModels(i + 1).startParIndex
            idx1 = xspec.AllModels(i + 1).startParIndex + xspec.AllModels(i + 1).nParameters - 1
            print("Data group", i, " parameter indices", idx0, idx1)
            idx = idx0
            comp_idx = 1
            for j, comp_name in enumerate(xspec.AllModels(i + 1).componentNames):
                print("  Data group", i, " model component:", j + 1, "(", comp_name, ")")
                comp = getattr(xspec.AllModels(i + 1), comp_name)
                for pn in comp.parameterNames:
                    param = xspec.AllModels(i + 1)(comp_idx)
                    print("     %12s (%3i), link: %7s,  frozen: %5s," % (pn, idx, param.link, param.frozen), end="")
                    print(" values:", param.values)
                    fp.paridx.append(idx)
                    fp.datagrp.append(i + 1)
                    fp.compidx.append(comp_idx)
                    fp.compname.append(comp_name)
                    fp.parname.append(pn)
                    fp.link.append(str(param.link))
                    fp.frozen.append(bool(param.frozen))
                    fp.values.append(param.values[0])
                    fp.error.append(param.error)
                    idx += 1
                    comp_idx += 1

        spec_counts   = {}
        first_spec_idx = {}   # {grp: first AllData 1-based index for that group}
        for i in range(1, xspec.AllData.nSpectra + 1):
            spec     = xspec.AllData(i)
            grp      = spec.dataGroup
            spec_counts[grp] = spec_counts.get(grp, 0) + 1
            if grp not in first_spec_idx:
                first_spec_idx[grp] = i
            r        = spec.rate
            exposure = spec.exposure
            print("spectrum %i (group %i/%i) rates %s  exposure %s" % (
                i, grp, spec_counts[grp], str(r), str(exposure)))
            fp.grp_ids.append(grp)
            fp.specnums.append(spec_counts[grp])
            fp.filenames.append(str(spec.fileName))
            fp.ignored.append(str(spec.ignoredString()))
            fp.exposure.append(float(exposure))
            fp.rate_net.append(float(r[0]))
            fp.rate_err.append(float(r[1]))
            fp.rate_gross.append(float(r[2]))
            fp.rate_pred.append(float(r[3]))

        fp.statistic  = xspec.Fit.statistic
        fp.dof        = xspec.Fit.dof
        fp.statmethod = xspec.Fit.statMethod
        fp.expression = xspec.AllModels(1).expression

        # flux: stored in Spectrum objects (not Model objects) when spectra are loaded.
        # Tuple layout: (erg_val, erg_lo, erg_hi, phot_val, phot_lo, phot_hi) per model.
        fp.flux_emin = flux_range[0]
        fp.flux_emax = flux_range[1]
        xspec.AllModels.calcFlux("%g %g" % flux_range)
        for grp_id, spec_i in sorted(first_spec_idx.items()):
            fv = xspec.AllData(spec_i).flux   # 6-tuple per model
            fp.flux[grp_id] = (fv[3], fv[4], fv[5], fv[0], fv[1], fv[2])
            # stored as (phot_val, phot_lo, phot_hi, erg_val, erg_lo, erg_hi)

        return fp

    def adapt(self, ifn=None):
        """
        Apply the stored parameter values and frozen flags to the current xspec model.

        Iterates over every stored parameter.  Linked parameters are skipped
        because their values are derived from the link target.  All other
        parameters have their current value and frozen flag updated to match
        the snapshot held in this instance.

        The model topology (number of groups, components, parameters) must
        match the stored snapshot.

        Parameters
        ----------
        ifn : str, optional
            Input filename (CSV or FITS). Adapt parameters from this file.            
        """
        if ifn is not None:
            if os.path.exists(ifn):
                self.read(ifn)

        # --- topology check ---
        n_xspec  = xspec.AllData.nGroups
        stored_groups = sorted(set(self.datagrp))
        if n_xspec != len(stored_groups):
            raise ValueError(
                "adapt: data group count mismatch — "
                "xspec has %i group(s), snapshot has %i" % (n_xspec, len(stored_groups)))
        for grp in stored_groups:
            xspec_names = []
            for comp_name in xspec.AllModels(grp).componentNames:
                xspec_names.extend(getattr(xspec.AllModels(grp), comp_name).parameterNames)
            stored_names = [self.parname[i] for i, g in enumerate(self.datagrp) if g == grp]
            if xspec_names != stored_names:
                raise ValueError(
                    "adapt: parameter name mismatch in data group %i —\n"
                    "  xspec:    %s\n"
                    "  snapshot: %s" % (grp, xspec_names, stored_names))
        
        # Adapt parameters
        for grp, par_idx, val, lnk, frz in zip(
                self.datagrp, self.paridx, self.values, self.link, self.frozen):
            start = xspec.AllModels(grp).startParIndex
            param = xspec.AllModels(grp)(par_idx - start + 1)
            try:
                lnk_str = "" if (lnk is None or (isinstance(lnk, float) and np.isnan(lnk))) \
                             else str(lnk).strip()
            except (TypeError, ValueError):
                lnk_str = ""
            if lnk_str:
                param.link = lnk_str
            else:
                param.link   = ""
                param.values = float(val)
                param.frozen = bool(frz)

    def restore(self, ifn=None):
        """
        Fully restore the xspec session from a snapshot.

        Reads from ``ifn`` if given, then:
        - loads ``self.filenames`` via ``xspec.AllData`` using explicit
          ``"{grp}:{specnum} filename"`` tokens where specnum counts from 1
          within each data group
        - ignores bad channels (``xspec.AllData.ignore("bad")``)
        - re-applies the per-spectrum ignored string stored in ``self.ignored``
        - sets ``xspec.Fit.statMethod = "cstat"``
        - creates the model from ``self.expression``
        - calls ``self.adapt()`` to restore all parameter values, frozen flags,
          and link expressions

        Parameters
        ----------
        ifn : str, optional
            CSV or FITS file to read before restoring.  If None the instance
            must already be populated (e.g. via ``FitProperties.read``).
        """
        if ifn is not None:
            fp = FitProperties.read(ifn)
            self.__dict__.update(fp.__dict__)

        if not self.filenames:
            raise ValueError("restore: no filenames stored — call read() first or pass ifn=")

        tokens = ["%i:%i %s" % (grp, snum, fn)
                  for grp, snum, fn in zip(self.grp_ids, self.specnums, self.filenames)]
        xspec.AllData(" ".join(tokens))

        xspec.AllData.ignore("bad")

        for seq_idx, ign in enumerate(self.ignored, start=1):
            if ign and ign.strip():
                xspec.AllData(seq_idx).ignore(ign.strip())

        xspec.Fit.statMethod = "cstat"

        if self.expression:
            xspec.Model(self.expression)
            self.adapt()

    @staticmethod
    def _read_csv(ifn):
        """
        Read a CSV file written by _write_csv.

        Returns a dict with keys ``params`` (DataFrame), ``fit`` (dict or None),
        and ``rates`` (list of dicts or None).
        """
        import re
        df = pd.read_csv(ifn, comment='#')
        df.drop_duplicates(inplace=True)

        _param_defaults = {
            "parameter_index": 0,     "data_group":     0,
            "component_index": 0,     "component_name": "",
            "parameter_name":  "",    "link":           "",
            "frozen":          False, "value":          float("nan"),
            "error":           float("nan"),
        }
        for col, default in _param_defaults.items():
            if col not in df.columns:
                df[col] = default
        # pandas reads empty CSV cells as NaN floats; coerce to the expected types
        df["link"]   = df["link"].fillna("").astype(str).str.strip()
        df["frozen"] = df["frozen"].map(
            lambda x: str(x).strip().lower() in ("true", "1", "yes")).astype(bool)

        fit = {}
        rates = {}
        with open(ifn) as fh:
            for line in fh:
                if not line.startswith('#'):
                    continue
                m = re.match(r'# model expression: (.+)', line)
                if m:
                    fit["expression"] = m.group(1).rstrip()
                    continue
                m = re.match(
                    r'# fit statistic: ([\d.eE+-]+), degrees of freedom: ([\d.eE+-]+), method: (\S+)',
                    line)
                if m:
                    fit.update(statistic=float(m.group(1)), dof=float(m.group(2)),
                               statmethod=m.group(3).rstrip())
                    continue
                # new format: # data group G, spectrum S, <key>: <val>
                m = re.match(r'# data group (\d+), spectrum (\d+), file: (.+)', line)
                if m:
                    key = (int(m.group(1)), int(m.group(2)))
                    rates.setdefault(key, {"datagrp": key[0], "specnum": key[1]})["filename"] = m.group(3).rstrip()
                    continue
                m = re.match(r'# data group (\d+), spectrum (\d+), ignored: (.*)', line)
                if m:
                    key = (int(m.group(1)), int(m.group(2)))
                    rates.setdefault(key, {"datagrp": key[0], "specnum": key[1]})["ignored"] = m.group(3).rstrip()
                    continue
                m = re.match(r'# data group (\d+), spectrum (\d+), rates.*?: \(([^)]+)\)', line)
                if m:
                    key = (int(m.group(1)), int(m.group(2)))
                    vals = [float(v.strip()) for v in m.group(3).split(',')]
                    rates.setdefault(key, {"datagrp": key[0], "specnum": key[1]}).update(
                        rate_net=vals[0], rate_err=vals[1], rate_gross=vals[2], rate_pred=vals[3])
                    continue
                m = re.match(r'# data group (\d+), spectrum (\d+), exposure: ([\d.eE+-]+)', line)
                if m:
                    key = (int(m.group(1)), int(m.group(2)))
                    rates.setdefault(key, {"datagrp": key[0], "specnum": key[1]})["exposure"] = float(m.group(3))
                    continue
                # old format (one spectrum per group, no specnum): backward compat
                m = re.match(r'# data group (\d+), file: (.+)', line)
                if m:
                    key = (int(m.group(1)), 1)
                    rates.setdefault(key, {"datagrp": key[0], "specnum": 1})["filename"] = m.group(2).rstrip()
                    continue
                m = re.match(r'# data group (\d+), ignored: (.*)', line)
                if m:
                    key = (int(m.group(1)), 1)
                    rates.setdefault(key, {"datagrp": key[0], "specnum": 1})["ignored"] = m.group(2).rstrip()
                    continue
                m = re.match(r'# data group (\d+), rates.*?: \(([^)]+)\)', line)
                if m:
                    key = (int(m.group(1)), 1)
                    vals = [float(v.strip()) for v in m.group(2).split(',')]
                    rates.setdefault(key, {"datagrp": key[0], "specnum": 1}).update(
                        rate_net=vals[0], rate_err=vals[1], rate_gross=vals[2], rate_pred=vals[3])
                    continue
                m = re.match(r'# data group (\d+), exposure: ([\d.eE+-]+)', line)
                if m:
                    key = (int(m.group(1)), 1)
                    rates.setdefault(key, {"datagrp": key[0], "specnum": 1})["exposure"] = float(m.group(2))
                    continue
                # flux energy range header
                m = re.match(r'# flux energy range: ([\d.eE+-]+) ([\d.eE+-]+) keV', line)
                if m:
                    fit["flux_emin"] = float(m.group(1))
                    fit["flux_emax"] = float(m.group(2))
                    continue
                # per-group flux values
                m = re.match(r'# data group (\d+), flux \([^)]+\): \(([^)]+)\)', line)
                if m:
                    grp = int(m.group(1))
                    vals = [float(v.strip()) for v in m.group(2).split(',')]
                    fit.setdefault("flux", {})[grp] = tuple(vals)

        rates_list = [rates[k] for k in sorted(rates)] if rates else None
        return dict(params=df, fit=fit or None, rates=rates_list)

    @staticmethod
    def _read_fits(ifn):
        """
        Read a FITS file written by _write_fits.

        Returns a dict with keys ``params`` (DataFrame), ``fit`` (dict),
        and ``rates`` (list of dicts).  The ``error`` column in ``params``
        uses the same ``(lo, hi, flags)`` string convention as the CSV format
        so the same display logic applies to both.
        """
        from astropy.io import fits as pyfits

        def _col(tbl, name, default):
            """Return column values, or a list of ``default`` if the column is absent."""
            if name not in tbl.names:
                return [default] * len(tbl)
            return tbl[name]

        with pyfits.open(ifn) as hdul:
            phdr = hdul[0].header
            fit = dict(
                statistic =phdr.get("FITSTAT"),
                dof       =phdr.get("FITDOF"),
                statmethod=phdr.get("STATMETH"),
                expression=phdr.get("MODELEXP"),
                flux_emin =phdr.get("FLUXEMIN"),
                flux_emax =phdr.get("FLUXEMAX"),
            )

            try:
                pt = hdul["PARAMS"].data
                nan = float("nan")
                err_col = []
                for lo, hi, flags in zip(_col(pt, "ERR_LO",    nan),
                                         _col(pt, "ERR_HI",    nan),
                                         _col(pt, "ERR_FLAGS", "")):
                    if np.isnan(float(lo)) or np.isnan(float(hi)):
                        err_col.append(nan)
                    else:
                        err_col.append(f"({lo}, {hi}, {flags})")
                df = pd.DataFrame({
                    "parameter_index": [int(v)   for v in _col(pt, "PARIDX",  0)],
                    "data_group":      [int(v)   for v in _col(pt, "DATAGRP", 0)],
                    "component_index": [int(v)   for v in _col(pt, "COMPIDX", 0)],
                    "component_name":  [str(v).strip() for v in _col(pt, "COMPNAME", "")],
                    "parameter_name":  [str(v).strip() for v in _col(pt, "PARNAME",  "")],
                    "link":            [str(v).strip() for v in _col(pt, "LINK",     "")],
                    "frozen":          [bool(v)  for v in _col(pt, "FROZEN", False)],
                    "value":           [float(v) for v in _col(pt, "VALUE",  nan)],
                    "error":           err_col,
                })
            except KeyError:
                df = pd.DataFrame()

            try:
                rt = hdul["RATES"].data
                rates = [dict(
                    datagrp      =int(  _col(rt, "DATAGRP",       0)   [i]),
                    specnum      =int(  _col(rt, "SPECNUM",       1)   [i]),
                    filename     =str(  _col(rt, "FILENAME",      "")  [i]).strip(),
                    ignored      =str(  _col(rt, "IGNORED",       "")  [i]).strip(),
                    exposure     =float(_col(rt, "EXPOSURE",      nan) [i]),
                    rate_net     =float(_col(rt, "RATE_NET",      nan) [i]),
                    rate_err     =float(_col(rt, "RATE_ERR",      nan) [i]),
                    rate_gross   =float(_col(rt, "RATE_GROSS",    nan) [i]),
                    rate_pred    =float(_col(rt, "RATE_PRED",     nan) [i]),
                    flux_phot   =float(_col(rt, "FLUX_PHOT",    nan) [i]),
                    flux_phot_lo=float(_col(rt, "FLUX_PHOT_LO", nan) [i]),
                    flux_phot_hi=float(_col(rt, "FLUX_PHOT_HI", nan) [i]),
                    flux_ener   =float(_col(rt, "FLUX_ENER",    nan) [i]),
                    flux_ener_lo=float(_col(rt, "FLUX_ENER_LO", nan) [i]),
                    flux_ener_hi=float(_col(rt, "FLUX_ENER_HI", nan) [i]),
                ) for i in range(len(rt))]
            except KeyError:
                rates = []

        return dict(params=df, fit=fit, rates=rates or None)

    @classmethod
    def read(cls, ifn):
        """
        Read a model file and return a populated FitProperties instance.

        The format is selected by file extension: ``.fits`` / ``.fit`` uses
        _read_fits, anything else uses _read_csv.

        Parameters
        ----------
        ifn : str
            Input filename (CSV or FITS).

        Returns
        -------
        FitProperties
        """
        if os.path.exists(ifn) is False:
            raise IOError("FitProperties::read - File not found: %s" % ifn)
        
        if ifn.lower().endswith(('.fits', '.fit')):
            data = cls._read_fits(ifn)
        else:
            data = cls._read_csv(ifn)

        fp = cls()
        df = data["params"]

        fp.paridx   = df["parameter_index"].tolist()
        fp.datagrp  = df["data_group"].tolist()
        fp.compidx  = df["component_index"].tolist()
        fp.compname = df["component_name"].tolist()
        fp.parname  = df["parameter_name"].tolist()
        fp.link     = df["link"].tolist()
        fp.frozen   = df["frozen"].tolist()
        fp.values    = df["value"].tolist()
        fp.error    = df["error"].tolist()

        fit = data["fit"]
        if fit:
            fp.statistic  = fit.get("statistic")
            fp.dof        = fit.get("dof")
            fp.statmethod = fit.get("statmethod")
            fp.expression = fit.get("expression")
            fp.flux_emin  = fit.get("flux_emin")
            fp.flux_emax  = fit.get("flux_emax")
            fp.flux       = fit.get("flux", {})

        seen_grps = set()
        for r in (data["rates"] or []):
            fp.grp_ids.append(r["datagrp"])
            fp.specnums.append(r.get("specnum", 1))
            fp.filenames.append(r.get("filename", ""))
            fp.ignored.append(r.get("ignored", ""))
            fp.exposure.append(r["exposure"])
            fp.rate_net.append(r["rate_net"])
            fp.rate_err.append(r["rate_err"])
            fp.rate_gross.append(r["rate_gross"])
            fp.rate_pred.append(r["rate_pred"])
            # reconstruct per-group flux from FITS RATES columns (first spectrum per group wins)
            grp = r["datagrp"]
            if grp not in seen_grps and "flux_phot" in r:
                seen_grps.add(grp)
                fp.flux[grp] = (r["flux_phot"], r["flux_phot_lo"], r["flux_phot_hi"],
                                r["flux_ener"], r["flux_ener_lo"], r["flux_ener_hi"])

        return fp

    @staticmethod
    def display(ifn):
        """
        Print a human-readable summary of a model file written by write().

        Accepts both CSV (.csv) and FITS (.fits/.fit) files.  Shows fit
        statistics, parameter values grouped by data group and component,
        and per-group count rates where available.
        """
        if ifn.lower().endswith(('.fits', '.fit')):
            data = FitProperties._read_fits(ifn)
        else:
            data = FitProperties._read_csv(ifn)

        df   = data["params"]
        fit  = data["fit"]
        rates = data["rates"]

        if fit:
            if fit.get("expression"):
                print("Model: %s" % fit["expression"])
            if fit.get("statistic") is not None:
                print("Fit statistic: %.4f   dof: %g   method: %s" % (
                    fit["statistic"], fit["dof"], fit["statmethod"]))
            print()

        for i in np.unique(df['data_group']):
            gi = np.where(df["data_group"] == i)[0]
            print("data group", i)
            for comp in np.unique(df["component_name"].iloc[gi]):
                ci = np.where((df["data_group"] == i) & (df["component_name"] == comp))[0]
                print("  component", comp)
                for j in ci:
                    print("  ", df["parameter_index"].iloc[j], end="")
                    print(" ", df["parameter_name"].iloc[j], end=" ")
                    print(" (frozen", df["frozen"].iloc[j], end=") ")
                    print("", value_str(df["value"].iloc[j], df["error"].iloc[j]), end=" ")
                    lnk = str(df["link"].iloc[j]).strip()
                    if lnk:
                        lidx = pIdx4link(lnk)
                        if lidx is None:
                            print()
                            continue
                        ri = np.where(df["parameter_index"] == lidx)[0][0]
                        print("    link: %i (group %i: %s)" % (
                            lidx, df["data_group"].iloc[ri], df["parameter_name"].iloc[ri]), end=" ")
                    print()
            print()

        if rates:
            print("Count rates:")
            for r in rates:
                ign = r.get("ignored", "")
                ign_str = ("  ignored: %s" % ign) if ign else ""
                print("  group %i/%i  %s  exp=%.1fs  net=%.4f±%.4f  gross=%.4f  pred=%.4f  ct/s%s" % (
                    r["datagrp"], r.get("specnum", 1), r.get("filename", ""),
                    r["exposure"], r["rate_net"], r["rate_err"], r["rate_gross"], r["rate_pred"],
                    ign_str))
            print()

        flux_dict = data.get("fit", {}) or {}
        stored_flux = flux_dict.get("flux", {})
        if stored_flux:
            emin = flux_dict.get("flux_emin", "?")
            emax = flux_dict.get("flux_emax", "?")
            print("Flux (%s-%s keV):" % (emin, emax))
            for grp, fv in sorted(stored_flux.items()):
                # fv = (phot_val, phot_lo, phot_hi, erg_val, erg_lo, erg_hi)
                print("  group %i: phot=%.3e (%.3e-%.3e) ph/cm2/s  erg=%.3e (%.3e-%.3e) erg/cm2/s" % (
                    grp, fv[0], fv[1], fv[2], fv[3], fv[4], fv[5]))
            print()

    @classmethod
    def write(cls, ofn, errors=False, error_delta=2.706, flux_range=(0.2, 2.4)):
        """
        Collect the current xspec model state and write it to ``ofn``.

        The output format is selected by file extension: ``.fits`` / ``.fit``
        produce a FITS file; any other extension produces a CSV.

        Parameters
        ----------
        ofn : str
            Output filename.
        errors : bool or list
            Passed to collect / calc_errors.
        error_delta : float
            Passed to collect / calc_errors.
        flux_range : tuple of (float, float)
            Energy range in keV for the xspec flux command (default 0.2–2.4 keV).
        """
        fp = cls.collect(errors, error_delta, flux_range=flux_range)
        if ofn.lower().endswith(('.fits', '.fit')):
            fp._write_fits(ofn)
        elif ofn.lower().endswith('.csv'):
            fp._write_csv(ofn)
        else:
            raise ValueError("FitProperties::write - Unrecognized file extension: %s" % ofn)

    def _write_csv(self, ofn):
        print("Writing to", ofn)
        df = pd.DataFrame({
            "parameter_index": self.paridx,
            "data_group":      self.datagrp,
            "component_index": self.compidx,
            "component_name":  self.compname,
            "parameter_name":  self.parname,
            "link":            self.link,
            "frozen":          self.frozen,
            "value":           self.values,
            "error":           self.error,
        })
        df.to_csv(ofn, index=False)
        with open(ofn, 'a') as oo:
            oo.write("# model expression: %s\n" % (self.expression or ""))
            oo.write("# fit statistic: %f, degrees of freedom: %f, method: %s\n" % (
                self.statistic, self.dof, self.statmethod))
            for i, (grp, snum) in enumerate(zip(self.grp_ids, self.specnums)):
                r = (self.rate_net[i], self.rate_err[i], self.rate_gross[i], self.rate_pred[i])
                oo.write("# data group %i, spectrum %i, file: %s\n" % (grp, snum, self.filenames[i]))
                oo.write("# data group %i, spectrum %i, ignored: %s\n" % (grp, snum, self.ignored[i] if i < len(self.ignored) else ""))
                oo.write("# data group %i, spectrum %i, rates (bkg subtracted, net uncertainty, gross rate, predicted): %s\n" % (grp, snum, str(r)))
                oo.write("# data group %i, spectrum %i, exposure: %s, counts (bkg subtracted, net uncertainty, gross rate, predicted): %s \n" % (
                    grp, snum, str(self.exposure[i]), str([ri * self.exposure[i] for ri in r])))
            if self.flux:
                oo.write("# flux energy range: %g %g keV\n" % (self.flux_emin, self.flux_emax))
                for grp, fv in sorted(self.flux.items()):
                    # fv = (phot_val, phot_lo, phot_hi, erg_val, erg_lo, erg_hi)
                    oo.write("# data group %i, flux (%g-%g keV, phot (value, lo, hi), erg (value, lo, hi)): (%g, %g, %g, %g, %g, %g)\n" % (
                        grp, self.flux_emin, self.flux_emax,
                        fv[0], fv[1], fv[2], fv[3], fv[4], fv[5]))

    def _write_fits(self, ofn):
        """
        Output extensions
        -----------------
        PRIMARY   header keywords: FITSTAT, FITDOF, STATMETH, MODELEXP,
                                   FLUXEMIN, FLUXEMAX
        PARAMS    BinTable – one row per model parameter:
                      PARIDX     (I)  global parameter index
                      DATAGRP    (I)  data group (1-based)
                      COMPIDX    (I)  component index within data group
                      COMPNAME   (A)  component name (e.g. 'phabs', 'apec')
                      PARNAME    (A)  parameter name (e.g. 'nH', 'kT')
                      LINK       (A)  xspec link string, empty if unlinked
                      FROZEN     (L)  True if parameter is frozen
                      VALUE      (D)  best-fit value
                      ERR_LO     (D)  lower confidence bound (absolute)
                      ERR_HI     (D)  upper confidence bound (absolute)
                      ERR_FLAGS  (A)  xspec error flags string
        RATES     BinTable – one row per loaded spectrum:
                      DATAGRP       (I)  source data group (1-based)
                      SPECNUM       (I)  spectrum number within data group (1-based)
                      FILENAME      (A)  spectrum filename
                      IGNORED       (A)  xspec ignoredString for this spectrum
                      EXPOSURE      (D)  exposure time (s)
                      RATE_NET      (D)  background-subtracted net count rate
                      RATE_ERR      (D)  net rate uncertainty
                      RATE_GROSS    (D)  gross count rate
                      RATE_PRED     (D)  model-predicted count rate
                      FLUX_PHOT     (D)  photon flux (ph/cm²/s) over FLUXEMIN–FLUXEMAX
                      FLUX_PHOT_LO  (D)  lower 1-sigma bound on photon flux
                      FLUX_PHOT_HI  (D)  upper 1-sigma bound on photon flux
                      FLUX_ENER     (D)  energy flux (erg/cm²/s) over FLUXEMIN–FLUXEMAX
                      FLUX_ENER_LO  (D)  lower 1-sigma bound on energy flux
                      FLUX_ENER_HI  (D)  upper 1-sigma bound on energy flux
        """
        from astropy.io import fits as pyfits

        err_lo    = [float(e[0]) if e else float("nan") for e in self.error]
        err_hi    = [float(e[1]) if e else float("nan") for e in self.error]
        err_flags = [str(e[2])   if len(e) > 2          else ""  for e in self.error]

        max_cname = max((len(s) for s in self.compname), default=1)
        max_pname = max((len(s) for s in self.parname),  default=1)
        max_link  = max((len(s) for s in self.link),     default=1)
        max_flags = max((len(s) for s in err_flags),     default=1)

        params_hdu = pyfits.BinTableHDU.from_columns([
            pyfits.Column(name="PARIDX",    format="I",             array=np.array(self.paridx,  dtype=np.int16)),
            pyfits.Column(name="DATAGRP",   format="I",             array=np.array(self.datagrp, dtype=np.int16)),
            pyfits.Column(name="COMPIDX",   format="I",             array=np.array(self.compidx, dtype=np.int16)),
            pyfits.Column(name="COMPNAME",  format=f"{max_cname}A", array=np.array(self.compname)),
            pyfits.Column(name="PARNAME",   format=f"{max_pname}A", array=np.array(self.parname)),
            pyfits.Column(name="LINK",      format=f"{max_link}A",  array=np.array(self.link)),
            pyfits.Column(name="FROZEN",    format="L",             array=np.array(self.frozen)),
            pyfits.Column(name="VALUE",     format="D",             array=np.array(self.values)),
            pyfits.Column(name="ERR_LO",    format="D",             array=np.array(err_lo)),
            pyfits.Column(name="ERR_HI",    format="D",             array=np.array(err_hi)),
            pyfits.Column(name="ERR_FLAGS", format=f"{max_flags}A", array=np.array(err_flags)),
        ])
        params_hdu.name = "PARAMS"

        nan = float("nan")
        _fnan = (nan,) * 6
        # fv layout: (phot_val, phot_lo, phot_hi, erg_val, erg_lo, erg_hi)
        flux_phot      = [self.flux.get(grp, _fnan)[0] for grp in self.grp_ids]
        flux_phot_lo   = [self.flux.get(grp, _fnan)[1] for grp in self.grp_ids]
        flux_phot_hi   = [self.flux.get(grp, _fnan)[2] for grp in self.grp_ids]
        flux_ener      = [self.flux.get(grp, _fnan)[3] for grp in self.grp_ids]
        flux_ener_lo   = [self.flux.get(grp, _fnan)[4] for grp in self.grp_ids]
        flux_ener_hi   = [self.flux.get(grp, _fnan)[5] for grp in self.grp_ids]

        max_fname = max((len(s) for s in self.filenames), default=1)
        max_ign   = max((len(s) for s in self.ignored),   default=1)
        rates_hdu = pyfits.BinTableHDU.from_columns([
            pyfits.Column(name="DATAGRP",      format="I",             array=np.array(self.grp_ids,   dtype=np.int16)),
            pyfits.Column(name="SPECNUM",      format="I",             array=np.array(self.specnums,  dtype=np.int16)),
            pyfits.Column(name="FILENAME",     format=f"{max_fname}A", array=np.array(self.filenames)),
            pyfits.Column(name="IGNORED",      format=f"{max_ign}A",   array=np.array(self.ignored)),
            pyfits.Column(name="EXPOSURE",     format="D",             array=np.array(self.exposure)),
            pyfits.Column(name="RATE_NET",     format="D",             array=np.array(self.rate_net)),
            pyfits.Column(name="RATE_ERR",     format="D",             array=np.array(self.rate_err)),
            pyfits.Column(name="RATE_GROSS",   format="D",             array=np.array(self.rate_gross)),
            pyfits.Column(name="RATE_PRED",    format="D",             array=np.array(self.rate_pred)),
            pyfits.Column(name="FLUX_PHOT",    format="D",             array=np.array(flux_phot)),
            pyfits.Column(name="FLUX_PHOT_LO", format="D",             array=np.array(flux_phot_lo)),
            pyfits.Column(name="FLUX_PHOT_HI", format="D",             array=np.array(flux_phot_hi)),
            pyfits.Column(name="FLUX_ENER",    format="D",             array=np.array(flux_ener)),
            pyfits.Column(name="FLUX_ENER_LO", format="D",             array=np.array(flux_ener_lo)),
            pyfits.Column(name="FLUX_ENER_HI", format="D",             array=np.array(flux_ener_hi)),
        ])
        rates_hdu.name = "RATES"

        primary_hdu = pyfits.PrimaryHDU()
        primary_hdu.header["FITSTAT"]  = (self.statistic,        "fit statistic value")
        primary_hdu.header["FITDOF"]   = (self.dof,              "degrees of freedom")
        primary_hdu.header["STATMETH"] = (self.statmethod,       "fit statistic method")
        primary_hdu.header["MODELEXP"] = (self.expression or "", "xspec model expression")
        primary_hdu.header["FLUXEMIN"] = (self.flux_emin if self.flux_emin is not None else nan,
                                          "flux lower energy bound (keV)")
        primary_hdu.header["FLUXEMAX"] = (self.flux_emax if self.flux_emax is not None else nan,
                                          "flux upper energy bound (keV)")

        hdul = pyfits.HDUList([primary_hdu, params_hdu, rates_hdu])
        hdul.writeto(ofn, overwrite=True)
        print(f"Written: {ofn}")


def photon_percentage_range(percentage, index=1, group=1):
    """
    Calculate the energy range containing a given percentage of the total source photons based on the current xspec model.

    Parameters
    ----------
    percentage : float
        Fraction of total photons to enclose, e.g. 0.9 for 90%.
    index : int
        xspec data group index (1-based) used to retrieve the energy grid.

    Returns
    -------
    tuple of (float, float)
        (e_low, e_high) in keV enclosing `percentage` of the total photon flux.
    """
    margin = (1.0 - percentage) / 2.0

    last_rebin = xspec._pyXspec.getplotSettings("lastrebin")
    min_sig = last_rebin[0]
    max_bins = last_rebin[1]
    groups = last_rebin[2]
    err = last_rebin[3]
    xspec.Plot.setRebin(minSig=0, maxBins=1, groupNum=-1, errType='quad')

    device = xspec.Plot.device
    xspec.Plot.device='/null'

    xspec.Plot("data")
    energies = np.array(xspec.Plot.x(index))
    e_lo = energies[:-1]
    e_hi = energies[1:]
    xspec.Plot("data model")
    folded_values = np.array(xspec.Plot.model(index))

    total = folded_values.sum()
    cumulative = np.cumsum(folded_values) / total
    if total == 0:
        raise ValueError("Model photon flux is zero — cannot compute percentage range.")

    idx_lo = np.searchsorted(cumulative, margin)
    idx_hi = np.searchsorted(cumulative, 1.0 - margin)
    # print("commands: ",xspec.Plot.commands)
    xspec.Plot.device=device
    xspec.Plot.setRebin(minSig=min_sig, maxBins=max_bins, groupNum=groups, errType=err)
    return float(e_lo[idx_lo]), float(e_hi[min(idx_hi, len(e_hi) - 1)])


def fit_and_photon_percentage(data_str, ignore_ranges, model_fn, stat_method="cstat"):
    """
    Load data into xspec, apply energy ignores, fit a model, and return spectral results.

    Parameters
    ----------
    data_str : str
        xspec data string (filename or multi-spectrum "1:1 file1 1:2 file2" format).
    ignore_ranges : list of str
        xspec ignore expressions applied in order, e.g. ["**-0.2", "5.-**", "bad"].
    model_fn : callable
        Function that builds and sets the xspec model.
    stat_method : str
        xspec fit statistic (default: "cstat").

    Returns
    -------
    list of float
        [e_lo, e_hi, cflux] — energy range enclosing 80% of model photons (keV)
        and log10 of the absorbed flux from the cflux component.
    """
    xspec.AllData(data_str)
    for r in ignore_ranges:
        xspec.AllData.ignore(r)
    xspec.Plot.device = '/null'
    xspec.Plot.xAxis = "keV"
    xspec.Fit.statMethod = stat_method
    model_fn()
    xspec.Fit.nIterations = 100
    xspec.Fit.perform()
    a, b = photon_percentage_range(0.8)
    cflux = xspec.AllModels(1).cflux.lg10Flux.values[0]
    print("e_lo, e_hi", a, b, "cflux", cflux)
    return [a, b, cflux]
