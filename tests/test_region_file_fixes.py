#!/home/majestix/hdd/python/bin/python3.8
"""
Plain verification script for the region-file fixes in
bin/xmm_correct_region_file_format.py and etc/io_helper.py:

- multi-region ds9 text files no longer silently drop all but the last
  region (fixed in xmm_correct_region_file_format.py, which used to
  reassign rather than append to its region list)
- annulus/ellipse regions are no longer rejected by fits_multi_region_writer
  due to case-sensitive shape matching ("Annulus"/"Ellipse" vs. the
  lowercase shape names ds9 actually writes)
- arcsec (") / arcmin (') radius unit suffixes are still converted to
  decimal degrees
- single-region files still convert cleanly (regression check)

Run directly:
    /home/majestix/hdd/python/bin/python3.8 tests/test_region_file_fixes.py
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from astropy.io import fits as pyfits
from xmmpy.etc import fits_multi_region_writer
from xmmpy.bin.xmm_correct_region_file_format import parse_region_lines

failures = []


def check(name, condition):
    status = "PASS" if condition else "FAIL"
    print("[%s] %s" % (status, name))
    if not condition:
        failures.append(name)


def test_multi_circle_writer(tmpdir):
    ofn = os.path.join(tmpdir, "two_circles.fits")
    lines = ["circle(24845.502,24864.19,200.00009)", "circle(22811.648,25228.416,700.00024)"]
    fits_multi_region_writer(lines, ofn)
    with pyfits.open(ofn) as ff:
        data = ff[1].data
        check("multi-circle writer produces 2 rows", len(data) == 2)
        check("multi-circle writer keeps first circle's values",
              data["SHAPE"][0].strip().upper() == "CIRCLE" and abs(data["X"][0]-24845.502) < 1e-2)
        check("multi-circle writer keeps second circle's values",
              abs(data["X"][1]-22811.648) < 1e-2 and abs(data["R"][1]-700.00024) < 1e-2)


def test_annulus_writer(tmpdir):
    ofn = os.path.join(tmpdir, "annulus.fits")
    fits_multi_region_writer(["annulus(24845.502,24864.19,200.0,300.0)"], ofn)
    with pyfits.open(ofn) as ff:
        data = ff[1].data
        check("annulus (lowercase, as ds9 writes it) no longer rejected", len(data) == 1)
        check("annulus shape recorded correctly", data["SHAPE"][0].strip().lower() == "annulus")


def test_ellipse_writer(tmpdir):
    ofn = os.path.join(tmpdir, "ellipse.fits")
    fits_multi_region_writer(["ellipse(24845.502,24864.19,200.0,100.0,45.0)"], ofn)
    with pyfits.open(ofn) as ff:
        data = ff[1].data
        check("ellipse (lowercase, as ds9 writes it) no longer rejected", len(data) == 1)


def test_unit_suffix_conversion(tmpdir):
    reg_fn = os.path.join(tmpdir, "sky_units.reg")
    with open(reg_fn, "w") as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=green\n")
        f.write("fk5\n")
        f.write('circle(266.99915,-60.087312,5.0")\n')
        f.write("circle(266.99915,-60.087312,2.0')\n")
    lines = parse_region_lines(reg_fn)
    check("unit-suffix parsing yields both lines", len(lines) == 2)
    check("arcsec suffix converted to degrees (5.0\"/3600)",
          abs(float(lines[0].split(",")[-1].rstrip(")")) - 5.0/3600) < 1e-9)
    check("arcmin suffix converted to degrees (2.0'/60)",
          abs(float(lines[1].split(",")[-1].rstrip(")")) - 2.0/60) < 1e-9)


def test_correct_regions_multi_and_single(tmpdir):
    # Reproduces the exact real-world failure: a background region file with
    # two circles (see LHS_1140_reg_bkg_pn_2circs.fits, which ended up with
    # only 1 row due to this bug) should now yield a genuine 2-row FITS file.
    two_fn = os.path.join(tmpdir, "two_circles.reg")
    with open(two_fn, "w") as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("physical\n")
        f.write("circle(24845.502,24864.19,200.00009)\n")
        f.write("circle(22811.648,25228.416,700.00024)\n")

    one_fn = os.path.join(tmpdir, "one_circle.reg")
    with open(one_fn, "w") as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("physical\n")
        f.write("circle(24845.502,24864.19,200.00009)\n")

    two_lines = parse_region_lines(two_fn)
    one_lines = parse_region_lines(one_fn)

    two_out = os.path.join(tmpdir, "two_circles_out.fits")
    one_out = os.path.join(tmpdir, "one_circle_out.fits")
    fits_multi_region_writer(two_lines, two_out)
    fits_multi_region_writer(one_lines, one_out)

    with pyfits.open(two_out) as ff:
        check("two-region .reg file produces 2-row FITS (the original bug)", len(ff[1].data) == 2)
    with pyfits.open(one_out) as ff:
        check("single-region .reg file still produces 1-row FITS (regression check)", len(ff[1].data) == 1)


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        test_multi_circle_writer(tmpdir)
        test_annulus_writer(tmpdir)
        test_ellipse_writer(tmpdir)
        test_unit_suffix_conversion(tmpdir)
        test_correct_regions_multi_and_single(tmpdir)

    print()
    if failures:
        print("%d check(s) FAILED: %s" % (len(failures), ", ".join(failures)))
        sys.exit(1)
    else:
        print("All checks passed.")


if __name__ == "__main__":
    main()
