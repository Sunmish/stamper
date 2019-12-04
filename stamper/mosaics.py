# A utility that uses montage (http://montage.ipac.caltech.edu/) to create
# mosaics.
# ----------------
# Initial functionality is only for DSS2.


import numpy as np
import os

from subprocess import Popen
from astropy.io import fits

from stamper.coord import hms_dms_dd
from stamper.postage_stamps import DSS2


def angular_distance(coords1, coords2):
    """Get the angular distance between a set of RA, DEC coordinates in [dd]."""

    cos_angle = math.sin(math.radians(coords1[1])) * \
                math.sin(math.radians(coords2[1])) + \
                math.cos(math.radians(coords1[1])) * \
                math.cos(math.radians(coords2[1])) * \
                math.cos(math.radians(coords1[0]) - math.radians(coords2[0]))

    try:
        gamma = math.degrees(math.acos(cos_angle))
    except ValueError:
        gamma = math.degrees(math.acos(min(1, max(cos_angle, -1))))

    return gamma



def coord_gen(coords, fov, max_size=5.0):
    """Generate coords necessary for covering fov.

    Coordinages (`coords`) in decimal degrees. `fov` in degrees.
    """


    





def GLEAM_mosaic(coords, fov, bands=["170-231"], indir="./", outname="mosaic"):
    """
    """









def DSS_mosaic(coords, fov, bands=["DSS2B", "DSS2R", "DSS2IR"], indir="./", \
			   outname="mosaic"):
	"""
	"""

	if not indir.endswith("/"): indir += "/"
	if not os.path.exists(indir): os.mkdir(indir)

	try: coords = hms_dms_dd(coords[0], coords[1])
	except Exception: pass

	for band in bands:

		if not os.path.exists(indir+"raw/"): os.mkdir(indir+"raw/")
		if not os.path.exists(indir+"proj/"): os.mkdir(indir+"proj/")

		Popen("cd {0}raw/".format(indir), shell=True).wait()
		Popen("mArchiveList dss {0} {1} {2} {3} {3} remote.tbl".format(\
			  band, coords[0], coords[1], fov), shell=True).wait()

		archive_list = np.genfromtxt(indir+"raw/remote.tbl", dtype=object, \
			                         skip_header=3)

		ra_list, dec_list = archive_list[:, 12], archive_list[:, 13]
		fnames = archive_list[:, -1]

		for i in range(len(fnames)):
			DSS2(float(ra_list[i]), float(dec_list[i]), 0.5, indir+"raw/",
				 fnames[i])

		Popen("cd ..", shell=True).wait()

		Popen("mImgtbl {0}raw rimages.tbl".format(indir), shell=True).wait()
		Popen("mMakeHdr rimages.tbl {0}hdr.hdr".format(indir), shell=True).wait()

		Popen("mProjExec -p {0}raw rimages.tbl {0}hdr.hdr {0}proj stats.tbl".format(\
			  indir), shell=True).wait()
		Popen("mImgtbl {0}proj pimage.tbl".format(indir), shell=True).wait()

		Popen("mAdd -p {0}proj pimages.tbl {0}hdr.hdr {1}_{2}.fits".format(\
			  indir, outname, band))
