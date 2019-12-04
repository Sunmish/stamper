#! /usr/bin/env python

from stamper import postage_stamps, coord
from argparse import ArgumentParser


def main():
    """
    """

    ps = ArgumentParser(description="Download DSS2 stamps.")
    ps.add_argument("-r", "--ra", default=None, required=True, 
                    )
    ps.add_argument("-d", "--dec", default=None, required=True,
                    )
    ps.add_argument("-s", "--size", "--fov", dest="fov", default=1.,
                    type=float, help="Image size in degrees. Max 1 degree.")
    ps.add_argument("-n", "--name", default=None, type=str)
    ps.add_argument("-D", "--directory", default="./", type=str)
    ps.add_argument("-b", "--band", default="blue", type=str)
    ps.add_argument("-o", "--overwrite", action="store_true")

    args = ps.parse_args()

    try:
        coords = (float(args.ra), float(args.dec))
    except ValueError:
        coords = coord.hms_dms_dd(args.ra, args.dec)


    if args.name is None:
        args.name = "J{:.2f}{:+.2f}_dss_{}".format(coords[0], coords[1], args.band)
        print(">>> Setting name to {}".format(args.name))

    postage_stamps.DSS2(coords[0], coords[1], args.fov, args.directory, args.name,
                        band=args.band, overwrite=args.overwrite)


if __name__ == "__main__":
    main()

