#! /usr/bin/env python

from stamper import postage_stamps, coord
from argparse import ArgumentParser


def main():
    """
    """

    ps = ArgumentParser(description="Download GLEAM stamps.")
    ps.add_argument("-r", "--ra", default=None, required=True, 
                    )
    ps.add_argument("-d", "--dec", default=None, required=True,
                    )
    ps.add_argument("-s", "--size", "--fov", dest="fov", default=5.,
                    type=float, help="Image size in degrees. Max 5 degrees.")
    ps.add_argument("-n", "--name", default=None, type=str)
    ps.add_argument("-D", "--directory", default="./", type=str)
    ps.add_argument("-b", "--band", default=["170-231"], type=str, nargs="*")
    ps.add_argument("-o", "--overwrite", action="store_true")

    args = ps.parse_args()

    try:
        coords = (float(args.ra), float(args.dec))
    except ValueError:
        coords = coord.hms_dms_dd(args.ra, args.dec)

    for band in args.band:
        if args.name is None:
            args.name = "J{:.2f}{:+.2f}".format(coords[0], coords[1])
            print("Setting name to {}".format(args.name))


        postage_stamps.GLEAM(coords[0], coords[1], args.fov, args.directory, args.name,
                             band=band, overwrite=args.overwrite)


if __name__ == "__main__":
    main()