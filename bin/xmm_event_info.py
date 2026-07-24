#!/usr/bin/python3
import argparse
from pathlib import Path

KEYS = ["INSTRUME", "EXPIDSTR", "OBJECT", "FILTER", "OBS_MODE", "SUBMODE", "DATE-OBS", "ONTIME"]


def find_event_files(path, pattern="*ImagingEvts.ds", recursive=True):
    p = Path(path)
    if p.is_file():
        return [p]
    globber = p.rglob if recursive else p.glob
    return sorted(globber(pattern))


def read_event_info(filename):
    from astropy.io import fits as pyfits
    info = {"file": filename}
    with pyfits.open(filename) as ff:
        headers = [h.header for h in ff]
        for k in KEYS:
            info[k] = next((h[k] for h in headers if k in h), None)
    return info


def report(path, pattern="*ImagingEvts.ds", recursive=True):
    files = find_event_files(path, pattern=pattern, recursive=recursive)
    if not files:
        print("No event files found matching %r under %s" % (pattern, path))
        return []
    rows = []
    for fn in files:
        try:
            rows.append(read_event_info(fn))
        except Exception as e:
            print("  Could not read %s: %s" % (fn, e))
    return rows


def print_table(rows, relative_to=None):
    if not rows:
        return
    cols = ["file"] + KEYS
    display = []
    for r in rows:
        d = dict(r)
        if relative_to is not None:
            try:
                d["file"] = str(Path(r["file"]).relative_to(relative_to))
            except ValueError:
                d["file"] = str(r["file"])
        else:
            d["file"] = str(r["file"])
        display.append(d)

    widths = {c: max(len(c), *(len(str(d[c])) for d in display)) for c in cols}
    header = "  ".join(c.ljust(widths[c]) for c in cols)
    print(header)
    print("-" * len(header))
    for d in display:
        print("  ".join(str(d[c]).ljust(widths[c]) for c in cols))


if __name__ == "__main__":
    # requires:
    #
    # export PYTHONPATH=$PYTHONPATH:/home/majestix/hdd/tools
    parser = argparse.ArgumentParser(
                    prog='xmm_event_info.py',
                    description='Read XMM event-files and report basic header information (target, filter, ontime, etc).',
                    epilog='Use at your own discretion...')

    parser.add_argument('path', help="Directory to search (recursively, by default), or a single event file.")
    parser.add_argument('--pattern', default="*ImagingEvts.ds", help="Glob pattern for event files when 'path' is a directory (default: *ImagingEvts.ds)")
    parser.add_argument('--no-recursive', dest='recursive', action='store_false', help="Only search the given directory itself, not its subdirectories.")

    args = parser.parse_args()

    rows = report(args.path, pattern=args.pattern, recursive=args.recursive)
    print_table(rows, relative_to=args.path if Path(args.path).is_dir() else None)
