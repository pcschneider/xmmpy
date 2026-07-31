#!/home/majestix/hdd/python/bin/python3.8
"""
Plain verification script for Obs._add_exposure() (obstools/xmm_obs.py):

- two exposures sharing the same raw exp_id (different detectors) still get
  disambiguated via the "+50" offset, as before
- three (or more) exposures sharing the same raw exp_id are no longer an
  error: the collision-avoidance loop keeps walking (+50, +100, ...) until
  it finds a free slot, instead of raising after a single fallback attempt
  (this was the actual observed bug: a background with pn/m1/m2 all sharing
  one EXP_ID would drop the m2 exposure entirely)
- re-adding an exposure with the same evt_filename (e.g. a repeated
  populate_exposures() call) reuses the existing exp_id instead of minting
  a new synthetic one, so the exposures dict doesn't grow with duplicates
- the "already in list, replacing" warning goes through the named "xmmpy"
  logger, not the root logger

Run directly:
    /home/majestix/hdd/python/bin/python3.8 tests/test_add_exposure.py
"""
import os
import sys
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from xmmpy.obstools import Obs

failures = []


def check(name, condition):
    status = "PASS" if condition else "FAIL"
    print("[%s] %s" % (status, name))
    if not condition:
        failures.append(name)


class FakeExposure:
    """
    Minimal stand-in for xmmpy.obstools.xmm_exp.Exposure, exposing only the
    attributes _add_exposure() actually touches (exp_id, det, evt_filename),
    so these checks don't need real event files on disk.
    """
    def __init__(self, exp_id, det, evt_filename):
        self.exp_id = exp_id
        self.det = det
        self.evt_filename = evt_filename

    def __str__(self):
        return "Exp: "+str(self.exp_id)+" ("+str(self.det)+" -> "+str(self.evt_filename)+")"


class _RecordCollector(logging.Handler):
    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append(record)


def new_obs():
    return Obs(obsID="0000000000", populate=False, verbose=0)


def test_two_way_collision():
    o = new_obs()
    o._add_exposure(FakeExposure("0000000000002", "pn", "pn.fits"))
    o._add_exposure(FakeExposure("0000000000002", "m1", "m1.fits"))
    check("two-way collision: both exposures present", len(o.exposures) == 2)
    check("two-way collision: second exposure shifted by +50",
          "0000000000052" in o.exposures and o.exposures["0000000000052"].det == "m1")


def test_three_way_collision():
    # Reproduces the real-world bug: pn, m1, m2 event files that all share
    # the same raw EXP_ID header value.
    o = new_obs()
    o._add_exposure(FakeExposure("0000000000002", "pn", "pn.fits"))
    o._add_exposure(FakeExposure("0000000000002", "m1", "m1.fits"))
    o._add_exposure(FakeExposure("0000000000002", "m2", "m2.fits"))
    check("three-way collision: all three exposures present (no dropped m2)", len(o.exposures) == 3)
    check("three-way collision: dets are pn/m1/m2",
          set(e.det for e in o.exposures.values()) == {"pn", "m1", "m2"})
    check("three-way collision: m2 walked past the occupied +50 slot to +100",
          "0000000000102" in o.exposures and o.exposures["0000000000102"].det == "m2")


def test_repeated_add_same_filename_no_duplicate():
    o = new_obs()
    o._add_exposure(FakeExposure("0000000000002", "pn", "pn.fits"))
    o._add_exposure(FakeExposure("0000000000002", "m1", "m1.fits"))
    o._add_exposure(FakeExposure("0000000000002", "m2", "m2.fits"))

    # Simulate a second populate_exposures() pass: fresh Exposure objects,
    # same evt_filename, exp_id reset back to the raw header value.
    o._add_exposure(FakeExposure("0000000000002", "pn", "pn.fits"))
    o._add_exposure(FakeExposure("0000000000002", "m1", "m1.fits"))
    o._add_exposure(FakeExposure("0000000000002", "m2", "m2.fits"))

    check("repeated add: exposure count stays at 3 (no duplicates)", len(o.exposures) == 3)
    check("repeated add: original exp_ids are reused",
          set(o.exposures.keys()) == {"0000000000002", "0000000000052", "0000000000102"})


def test_warning_uses_named_logger():
    collector = _RecordCollector()
    logger = logging.getLogger("xmmpy")
    prev_level = logger.level
    logger.addHandler(collector)
    logger.setLevel(logging.WARNING)
    try:
        o = new_obs()
        o._add_exposure(FakeExposure("0000000000002", "pn", "pn.fits"))
        o._add_exposure(FakeExposure("0000000000002", "pn", "pn.fits"))  # same filename -> replace + warn
    finally:
        logger.removeHandler(collector)
        logger.setLevel(prev_level)

    warn_records = [r for r in collector.records if r.levelno == logging.WARNING]
    check("replacement warning was emitted exactly once", len(warn_records) == 1)
    check("replacement warning uses the named 'xmmpy' logger, not root",
          bool(warn_records) and warn_records[0].name == "xmmpy")


def main():
    test_two_way_collision()
    test_three_way_collision()
    test_repeated_add_same_filename_no_duplicate()
    test_warning_uses_named_logger()

    print()
    if failures:
        print("%d check(s) FAILED: %s" % (len(failures), ", ".join(failures)))
        sys.exit(1)
    else:
        print("All checks passed.")


if __name__ == "__main__":
    main()
