from xmmpy import Obs
import argparse
import astropy.units as aunits

def ks2time(conf, ks=None, margin_sec=1, format='cxcsec', verbose=1):
    obs = Obs(conf_file=conf, verbose=verbose)
    t0 = obs.observation_start - margin_sec*aunits.second
    ret = t0 + ks*aunits.second * 1000.
    if format == 'cxcsec':
        ret = ret.cxcsec
    try:
        ret.format = format
    except:
        raise Exception("Don't understand format: "+format)
    return ret

def time2ks(conf, t=None, margin_sec=1, verbose=1):
    obs = Obs(conf_file=conf, verbose=verbose)
    t0 = obs.obsservation_start - margin_sec*aunits.second
    ks = (t-t0.cxcsec)/aunits.second
    ks = ks/1000.
    return ks.value
    ret = t0 + ks*aunits.second * 1000.
    if format == 'cxcsec':
        ret = ret.cxcsec
    return ret

def time_overview(conf, verbose=1):
    obs = Obs(conf_file=conf, verbose=verbose)
    print("Obs %s has start-time: %8.1f (%s)" % (obs.config["obsID"], obs.observation_start.cxcsec, obs.observation_start.iso))
    for e in obs.exposures.values():
        print("Exposure: ", e)
        t0, duration = e['start'], e["exposure"]
        print("      starting at: %s; mjd: %.4f; cxcsec: %.2f; ontime: %.2f ks" % (t0, t0.mjd, t0.cxcsec, duration/1e3))
    time_bins, tb_iso = obs._time_bins(), obs._time_bins(ref_system='iso')
    for b, c in zip(time_bins, tb_iso):
        print("Time bin: ", b, "(Telescope time) or ",c, "Duration: ",(b[1]-b[0])/1000, "ks")

        
if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument(metavar="conf-file", type=str, dest='conf', help='Config-file')
    args.add_argument('-ks', type=float, help='Time in ks from start of obs')
    args.add_argument('-t', type=float, help='Time telescope time to ks from start of obs')
    args.add_argument('-f', type=str, default='cxcsec', help='Output format')
    args.add_argument('-o', '--overview', action='store_true', default=False, help='Print overview')
    args.add_argument('-v', type=int, default=1, help='Verbosity')
    args = args.parse_args()
    if args.overview:
        time_overview(args.conf, verbose=args.v)
        exit()
    if args.ks is not None and args.t is None:
        print(ks2time(args.conf, ks=args.ks, format=args.f, verbose=args.v))
    elif args.t is not None and args.ks is None:
        print(time2ks(args.conf, t=args.t, verbose=args.v))
    else:
        print("Need to specify either -ks or -t")