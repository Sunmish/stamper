'''
# Postage stamp downloader. ---------------------------------------------------
#
# So far this includes:
#     NVSS
#     SUMSS
#     TGSS (ADR1)
#     GLEAM
#     DSS2
#     PS1 (DR1)
#
# TODO:
# `ALL` function to allow downloading all relevant stamps for a source.
# Incorporate Montage to facillitate downloading of fits files that are larger
# than the limits of the postage stamp servers for e.g. PS1 (NVSS?).
# 
# Add "query name instead of coorinates" option
# FIX NVSS coordinate input - at the moment it doesn't accept ':' or any
# delimiter other than ' '.
# Change coord to SkyCoord from astropy.coordinates. Possibly. Nah this works.
# Move `PanSTARRS_mosaic` to a different module? 
# PS1_mosaic 
#
# POSTAGE STAMPS to add:
#      WISE
#      RASS
#      VLSSr ?
#      SDSS (when there is a proper postage stamp server...)
#      IRAS  ?
#      
#
# -----------------------------------------------------------------------------
'''


from __future__ import print_function, division

# For scraping the postage stamp servers and cleaning up.
from mechanize import Browser
import shutil
import os

from subprocess import Popen

import math
import numpy

# For editing FITS files to make them more usable.
from astropy.io import fits 

# For coordinate conversion (especially for NVSS, SUMSS, and TGSS):
from stamper.coord import hms_dms_dd, dd_hms_dms
from stamper.serrors import DeclinationError

#import matplotlib.pyplot as plt
#import aplpy

import logging
logging.basicConfig(format='%(levelname)s (%(module)s): %(message)s', \
    level=logging.DEBUG)



def fix_eq(image):
    '''Add `EQUINOX` key to FITS images.

    Images such as those from PS1 and XMM-Newton archive do not have built-in
    `EQUINOX` keys and some software is confused by this.'''

    if isinstance(image, basestring):
        hdulist = fits.open(image, mode='update')
        updatefile = True
    elif isinstance(image, fits.HDUList):
        hdulist = image
        updatefile = False
    else:
        raise TypeError('Image must be a filename or an astropy.io.fits.HDUList object.')

    hdulist[0].header['EQUINOX'] = 2000

    if updatefile:
        hdulist.flush()
        hdulist.close()
    else:
        return hdulist


def fix_cd(image):
    '''Changes CD keys to CDELT for GLEAM.'''

    hdulist = fits.open(image, mode='update')

    hdulist[0].header['CDELT1'] = hdulist[0].header['CD1_1']
    hdulist[0].header['CDELT2'] = hdulist[0].header['CD2_2']
    hdulist[0].header['CROTA1'] = hdulist[0].header['CD1_2']
    hdulist[0].header['CROTA2'] = hdulist[0].header['CD2_1']

    for key in ['CD1_1', 'CD2_2', 'CD1_2', 'CD2_1']:
        del hdulist[0].header[key]

    hdulist.flush()
    hdulist.close()   


def add_cdelt3(image, stokes="I"):
    """Add CDELT3 keys for FREQ parameters."""

    hdulist = fits.open(image, mode="update")

    hdulist[0].header["CRVAL3"] = stokes
    hdulist[0].header["CDELT3"] = 1.0
    hdulist[0].header["CTYPE3"] = "STOKES"
    hdulist[0].header["CUNIT3"] = "STOKES"
    hdulist[0].header["CRPIX3"] = 1.0
    hdulist[0].header["CROTA3"] = 0.0

    hdulist.flush()
    hdulist.close()


def add_cdelt4(image):
    """Add CDELT4 keys for STOKES parameters."""

    hdulist = fits.open(image, mode="update")

    hdulist[0].header["CRVAL4"] = hdulist[0].header["FREQ"]
    hdulist[0].header["CDELT4"] = 1.0
    hdulist[0].header["CTYPE4"] = "FREQ"
    hdulist[0].header["CUNIT4"] = "Hz"
    hdulist[0].header["CRPIX4"] = 1.0
    hdulist[0].header["CROTA4"] = 0.0

    hdulist.flush()
    hdulist.close()


def fix_ps(image):
    '''Convenience function to change FITS keywords of the PS1 data.

    From the documentation:
    "One quirk of these images is that they use the obsolete WCS 
    keywords PC001001, PC001002, PC002001, PC002002 instead of 
    the FITS standard keywords CD1_1, CD1_2, etc.  Many software 
    packages automatically handle the old PC keywords, but some 
    require special processing (e.g., in IDL you should call 
    the fits_cd_fix procedure to modify the header).  If you 
    do not use these keywords, you will find that the RA pixel 
    spacing has the wrong sign (it is positive instead of negative)."
    '''

    hdulist = fits.open(image, mode='update')
    hdulist[0].header['CD1_1'] = hdulist[0].header['CDELT1'] * \
                                 hdulist[0].header['PC001001']
    hdulist[0].header['CD1_2'] = hdulist[0].header['PC001002']
    hdulist[0].header['CD2_2'] = hdulist[0].header['CDELT2'] * \
                                 hdulist[0].header['PC002002']
    hdulist[0].header['CD2_1'] = hdulist[0].header['PC002001']

    del hdulist[0].header['PC001001']
    del hdulist[0].header['PC001002']
    del hdulist[0].header['PC002002']
    del hdulist[0].header['PC002001']
    del hdulist[0].header['CDELT1']
    del hdulist[0].header['CDELT2']

    hdulist.flush()
    hdulist.close()


def fix_minmax(image):
    '''Recalculate the datamin/datamax.'''

    hdulist = fits.open(image, mode='update')
    datamin, datamax = numpy.nanmin(hdulist[0].data), numpy.nanmax(hdulist[0].data)
    hdulist[0].header['DATAMIN'] = datamin
    hdulist[0].header['DATAMAX'] = datamax

    hdulist.flush()
    hdulist.close()


def cleanup(indir='./'):
    '''Cleans temporary files created by stamper. 

    This should only be needed if a function fails to complete.
    WARNING: if you have other files on your system in the working 
    directory with the name of one of the temporary files THEY WILL
    BE DELETED. 
    '''

    if not indir.endswith('/'): indir += '/'

    files = ['gleam.html', 'FIT', 'sumss.html', 'tgss.html', 'PanSTARRS.html']

    for f in files:
        try: os.remove(indir+f)
        except OSError: pass

    logging.info('Temporary files deleted.')


def DSS2_mosaic(RA, DEC, fov, outdir, outname, band="b", overwrite="c", \
                montage=True, cleanup=True, subtraction=False):
    """Make a mosaic of multiple DSS2 images.

    Use only if image is to be greater than 1 degree. Otherwise use 
    normal `DSS2` function which handles images up to 1 degree.

    """

    try: ra, dec = float(RA), float(DEC)
    except ValueError: ra, dec = hms_dms_dd(RA, DEC) 

    extra_fov = fov + 0.5  # Allow padding so that we can subtract a nice square image.
    size = extra_fov/2.0

    coords = {dec: [ra]}

    n_overlaps = int(math.ceil((size - 0.7)/0.6))

    for i in range(n_overlaps):
        coords[dec+0.6*(i+1)] = [ra]
        coords[dec-0.6*(i+1)] = [ra]

    for d in coords.keys():

        for i in range(0, n_overlaps):

            if (ra + (0.6*(1 + i) / numpy.cos(numpy.radians(d)))) >= 360.0:

                coords[d].append((ra + (0.6*(1 + i) / \
                                  numpy.cos(numpy.radians(d)))) - 360.0)

            else:

                coords[d].append(ra + (0.6*(1 + i) / \
                                 numpy.cos(numpy.radians(d))))


            if (ra - (0.3*(1 + i) / numpy.cos(numpy.radians(d)))) < 0.0:

                coords[d].append(ra - (0.6*(1 + i) / \
                                 numpy.cos(numpy.radians(d))) + 360.0)

            else:

                coords[d].append(ra - (0.6*(1 + i) / \
                                 numpy.cos(numpy.radians(d))))

    logging.info('>>> Downloading %i stamps to form mosaic.' % \
                 int(len(coords.keys())*n_overlaps))

    count = 0

    for d in coords.keys():

        for i in range(0, n_overlaps):

            if overwrite == "c":
                try:
                    DSS2(coords[d][i], d, 0.8, outdir, '%s_%i_%.2f_%s' \
                          % (outname, i, round(d, 1), band), band, overwrite='r')
                except IOError:
                    pass

            else:

                DSS2(coords[d][i], d, 0.8, outdir, '%s_%i_%.2f_%s' \
                          % (outname, i, round(d, 1), band), band, overwrite)
            count += 1

    logging.info('>>> Downloaded %i stamps.' % count)

    if montage:

        for b in band:
            work_dir = os.path.abspath(outdir)
            work_dir = '%s/%s' % (work_dir, b)
            raw_dir  = '%s/raw' % work_dir
            prj_dir  = '%s/projected' % work_dir

            if not os.path.exists(work_dir):
                os.mkdir(work_dir)
            if not os.path.exists(raw_dir): 
                os.mkdir(raw_dir)
            else: 
                logging.warning('>>> Raw directory already exists: working in there.') 
            if not os.path.exists(prj_dir): 
                os.mkdir(prj_dir)
            else: 
                logging.warning('>>> Projected directory already exists: '\
                                'working in there.')  

            for spec in os.listdir(outdir):
                if spec.endswith('_%s.fits' % b):
                    shutil.copy2('%s/%s' % (outdir, spec), \
                                 '%s/%s' % (raw_dir, spec))    


            logging.info('>>> DSS2: making optimal header...')
            Popen('mImgtbl %s %s/rimages.tbl' % \
                 (raw_dir, work_dir), shell=True).wait()
            Popen('mMakeHdr %s/rimages.tbl %s/%s.hdr' % \
                 (work_dir, work_dir, outname), shell=True).wait()

            logging.info('>>> DSS2: re-projecting to optimal header...')
            Popen('mProjExec -p %s %s/rimages.tbl %s/%s.hdr %s %s/stats.tbl' % \
                 (raw_dir, work_dir, work_dir, outname, prj_dir, work_dir), \
                 shell=True).wait()
            Popen('mImgtbl %s %s/pimages.tbl' % (prj_dir, work_dir), \
                  shell=True).wait()

            if subtraction:
                # Add montage background subtraction routine.
                # raise ValueError('Subtraction not yet working.')
                final_dir = prj_dir
                final_tbl = '%s/pimages.tbl' % work_dir            

            else:

                final_dir = prj_dir
                final_tbl = '%s/pimages.tbl' % work_dir

            logging.info('>>> DSS2: adding images to form mosaic...')
            Popen('mAdd -e -p %s %s %s/%s.hdr %s/%s.fits' % \
                 (final_dir, final_tbl, work_dir, outname, work_dir, outname), \
                 shell=True).wait()

        
            if cleanup:

                logging.info('>>> DSS2: deleting extra FITS directories.')
                shutil.rmtree(raw_dir)
                shutil.rmtree(prj_dir)

            logging.info('>>> DSS2: trimming image...')
            Popen('mSubimage %s/%s.fits %s/%s_%s.fits %f %f %f %f' % \
                  (work_dir, outname, work_dir, outname, b, ra, dec, \
                  fov, fov), shell=True).wait()



def PS1_mosaic(RA, DEC, fov, outdir, outname, band=['g'], overwrite='r', \
    montage=True, cleanup=True, subtraction=False):
    '''Make a small mosaic of PanSTARRS data to make sure the region of
    interest is covered by the downloaded stamps.'''


    try: ra, dec = float(RA), float(DEC)
    except ValueError: ra, dec = hms_dms_dd(RA, DEC) 


    # if fov < 1.0:
    #     raise ValueError('>>> FOV for mosaics should be at least 1 degree.')
    if dec < -40:
        raise ValueError('>>> Pan-STARRS has no data below -40 declination.')
    if -40 < dec < -30:                  
        logging.warning('>>> Pan-STARRS data -40 < dec < -30 is low-quality.')

    fov = fov + 0.3  # Allow padding to centre images.
    size = (fov/2.0)

    coords = {dec: [ra]}
 
    n_overlaps = int(math.ceil((size - 0.1) / 0.15))
    print(n_overlaps)

    for i in range(n_overlaps):
    
        coords[dec+0.15*(i+1)] = [ra]
        coords[dec-0.15*(i+1)] = [ra]


    for d in coords.keys():

        for i in range(0, n_overlaps):

            if (ra + (0.15*(1 + i) / numpy.cos(numpy.radians(d)))) >= 360.0:

                coords[d].append((ra + (0.15*(1 + i) / \
                                  numpy.cos(numpy.radians(d)))) - 360.0)

            else:

                coords[d].append(ra + (0.15*(1 + i) / \
                                 numpy.cos(numpy.radians(d))))


            if (ra - (0.3*(1 + i) / numpy.cos(numpy.radians(d)))) < 0.0:

                coords[d].append(ra - (0.15*(1 + i) / \
                                 numpy.cos(numpy.radians(d))) + 360.0)

            else:

                coords[d].append(ra - (0.15*(1 + i) / \
                                 numpy.cos(numpy.radians(d))))

    logging.info('>>> Downloading %i stamps to form mosaic.' % \
                 int(len(coords.keys())*n_overlaps))

    count = 0

    for d in coords.keys():

        for i in range(0, n_overlaps):

            if overwrite == 'c':
                try:
                    PS1(coords[d][i], d, 0.3, outdir, '%s_%i_%.2f' \
                        % (outname, i, round(d, 1)), band, overwrite='r')
                except IOError:
                    pass

            else:

                PS1(coords[d][i], d, 0.3, outdir, '%s_%i_%.2f' \
                          % (outname, i, round(d, 1)), band, overwrite)
            count += 1

    logging.info('>>> Downloaded %i stamps.' % count)


    if montage:

        for b in band:
            work_dir = os.path.abspath(outdir)
            work_dir = '%s/%s' % (work_dir, b)
            raw_dir  = '%s/raw' % work_dir
            prj_dir  = '%s/projected' % work_dir

            if not os.path.exists(work_dir):
                os.mkdir(work_dir)
            if not os.path.exists(raw_dir): 
                os.mkdir(raw_dir)
            else: 
                logging.warning('>>> Raw directory already exists: working in there.') 
            if not os.path.exists(prj_dir): 
                os.mkdir(prj_dir)
            else: 
                logging.warning('>>> Projected directory already exists: '\
                                'working in there.')  

            for spec in os.listdir(outdir):
                if spec.endswith('_%s.fits' % b):
                    shutil.copy2('%s/%s' % (outdir, spec), \
                                 '%s/%s' % (raw_dir, spec))    


            logging.info('>>> Pan-STARRS: making optimal header...')
            Popen('mImgtbl %s %s/rimages.tbl' % \
                 (raw_dir, work_dir), shell=True).wait()
            Popen('mMakeHdr %s/rimages.tbl %s/%s.hdr' % \
                 (work_dir, work_dir, outname), shell=True).wait()

            logging.info('>>> Pan-STARRS: re-projecting to optimal header...')
            Popen('mProjExec -p %s %s/rimages.tbl %s/%s.hdr %s %s/stats.tbl' % \
                 (raw_dir, work_dir, work_dir, outname, prj_dir, work_dir), \
                 shell=True).wait()
            Popen('mImgtbl %s %s/pimages.tbl' % (prj_dir, work_dir), \
                  shell=True).wait()

            if subtraction:
                # Add montage background subtraction routine.
                # raise ValueError('Subtraction not yet working.')
                final_dir = prj_dir
                final_tbl = '%s/pimages.tbl' % work_dir            

            else:

                final_dir = prj_dir
                final_tbl = '%s/pimages.tbl' % work_dir

            logging.info('>>> Pan-STARRS: adding images to form mosaic...')
            Popen('mAdd -e -p %s %s %s/%s.hdr %s/%s.fits' % \
                 (final_dir, final_tbl, work_dir, outname, work_dir, outname), \
                 shell=True).wait()

        
            if cleanup:

                logging.info('>>> Pan-STARRS: deleting extra FITS directories.')
                shutil.rmtree(raw_dir)
                shutil.rmtree(prj_dir)

            logging.info('>>> Pan-STARRS: trimming image...')
            Popen('mSubimage %s/%s.fits %s/%s_%s.fits %f %f %f %f' % \
                  (work_dir, outname, work_dir, outname, b, ra, dec, \
                  (fov-0.4), (fov-0.4)), shell=True).wait()





def PS1(RA, DEC, fov, outdir, outname, band=['g'], overwrite='r', \
        fix_fits=True, cutout=True, fix_ps_keys=True):
    '''Download Pan-STARRS data from the postage stamp server.

    Parameters
    ----------
    RA        : string or float
                Right ascension - hms:dms or dd:dd.
    DEC       : string or float
                Declination - hms:dms or dd:dd.
    fov       : float
                Field of view - max size is ~ 0.42 degrees.
    outdir    : string
                Output directory. Use "./" for working directory.
    outname   : string
                Output basename. Band name is appended.
    band      : string, optional
                Pan-STARRS band. List or one of ["g", "r", "i", "y", "z"].
    overwrite : {'r', 'w'}, optional
                `r` does not overwriting. `w` does overwriting.
    fix_ps_key: bool, optional
                Change PS1 pizel spacing keys to be consistent with current
                conventions.

    '''


    def downloadlink(l):

        with open(l.text, 'w') as f:
            br.follow_link(l)
            f.write(br.response().read())


    try: ra, dec = float(RA), float(DEC)
    except ValueError: ra, dec = hms_dms_dd(RA, DEC) 


    # Perform some checks:
    # Declination must be > -40:
    if dec < -40:
        raise ValueError('>>> Pan-STARRS has no data below -40 declination.')
    if -40 < dec < -30:                  
        logging.warning('>>> Pan-STARRS data -40 < dec < -30 is low-quality and ' \
                        'not currently available. Stamp may not download.')
        # raise ValueError('>>> Pan-STARRS data -40 < dec < -30 is low-quality and ' \
        #                  'not current available.')


    # Check fov - if greater than 6000 pixels default to 6000 pixels:
    pov = int(fov * 3600.0 * 4.0)
    if (pov > 6000.0) and cutout:
        logging.warning('>>> Maximum size for Pan-STARRS images is 6000 pixels.')
        logging.warning('>>> Setting FOV to 6000 pixels (~0.4167 degrees).')
        pov = 6000
    

    # Check that the band is valid:
    if not isinstance(band, list): band = [band]
    for b in band:
        if b not in ['g', 'r', 'i', 'y', 'z']:
            raise ValueError('>>> `%s` is not a valid band for Pan-STARRS. \n'\
                             '>>> Valid bands are: [`g`, `r`, `i`, `y`, `z`].'\
                             % b )


    if not outdir.endswith('/'): outdir += '/'
    if not os.path.exists(outdir): os.mkdir(outdir)


    # Download each band separately.
    for b in band:

        flag = False

        if os.path.exists('%s%s_%s.fits' % (outdir, outname, b)):
            if overwrite == 'r':
                logging.warning('>>> %s_%s.fits already exists.' % (outname, b))
                flag = True
            elif overwrite == 'w':
                logging.warn('>>> %s_%s.fits already exists. '\
                             'Overwriting...' % (outname, b))
                os.remove('%s%s_%s.fits' % (outdir, outname, b))

        try:

            if not flag:

                logging.info('>>> Querying Pan-STARRS for %f, %f...' % (ra, dec))

                br = Browser()

                br.set_handle_robots(False)
                br.addheaders = [('User-agent', 'Firefox')]

                br.open('http://ps1images.stsci.edu/cgi-bin/ps1cutouts')

                br.select_form(nr=0)
                
                br['pos'] = '%f %f' % (ra, dec)

                control = br.find_control('filter')
                for item in control.items:
                    if item.name == '%s' % b:
                        item.selected = True
                    else:
                        item.selected = False

                br['size'] = '%i' % pov

                result = br.submit()
                content = result.read()

                with open('PanSTARRS.html', 'w') as f:
                    f.write(content)

                if cutout: lname = 'FITS-cutout'
                else: lname = 'FITS' 

                for link in br.links():
                    if link.text == lname:
                        stamp_name = link.text
                        downloadlink(l=link)
                        break

                shutil.copy2('./%s' % stamp_name, '%s%s_%s.fits' % (outdir, outname, b))

                os.remove(stamp_name)
                os.remove('PanSTARRS.html')

                # PS1 uses obsolete keywords for pixel spacing. We should replace those.
                if fix_ps_keys:

                    fix_ps('%s%s_%s.fits' % (outdir, outname, b))

                    logging.info('>>> Pan-STARRS stamp acquired and fixed.')

                else:

                    logging.info('>>> Pan-STARRS stamp acquired.')


        except Exception:
            raise



def GLEAM(RA, DEC, fov, outdir, outname, projection='ZEA_regrid', \
        band=['170-231'], overwrite='r', mask=True, robust=-1, \
        robust0_uname=None, robust0_pword=None, fix=True):
    '''Function to download GLEAM postage stamp(s) from the GLEAM VO server.

    Parameters
    ----------
    RA         : string or float
               Right ascension. Can be in many formats.
    DEC        : string or float
               Declination. Can be in many formats. Must be < +25 degrees.
    fov        : float
               Field size. Must be < 5.0 degrees.
    outdir     : string
               Output directory. `outdir` is created if it does not exist.
    outname    : string
               Output filename. If `outname` exists in `outdir` function stops.
    projection : {'ZEA', 'ZEA_regrid', 'SIN_regrid'}, optional
               Projection for fits file. Default is 'ZEA_regrid'.
    band       : {'white', '072-080', '080-088', '088-095', '095-103',
                '103-111', '103-134', '111-118', '118-126', '126-134',
                '139-147', '139-170', '147-154', '154-162', '162-170',
                '170-177', '170-231', '177-185', '185-193', '193-200',
                '200-208', '208-216', '216-223', '223-231', all}, optional
               Freq. band. Default is '170-231' - the stacked white image.
    overwrite  : {'w', 'r'}, optional
               Determine whether to overwrite existing files of same name.

    Raises
    ------
    UnboundLocalError
        If this occurs it is likely there is no stamp for the queried region.

    '''


    gleam_bands = ['072-080', '080-088', '088-095', '095-103',
                   '103-111', '103-134', '111-118', '118-126', '126-134',
                   '139-147', '139-170', '147-154', '154-162', '162-170',
                   '170-177', '170-231', '177-185', '185-193', '193-200',
                   '200-208', '208-216', '216-223', '223-231', '072-103']

    narrow_bands = ['072-080', '080-088', '088-095', '095-103',
                   '103-111', '111-118', '118-126', '126-134',
                   '139-147', '147-154', '154-162', '162-170',
                   '170-177', '177-185', '185-193', '193-200',
                   '200-208', '208-216', '216-223', '223-231']


    # The actual web-scraping function:
    def gleam_stamp(coords, fov, outdir, outname, b, mask, \
        projection='SIN_regrid', overwrite='r', robust=-1, \
        robust0_uname=None, robust0_pword=None, fix=True):
        '''Accessing the GLEAM server. 

        This is a separate function simply for ease of looping over all 
        bands if needed.
        '''


        # Download link must be a function or the FITS files become corrupted.
        # I don't know why.
        def downloadlink(l):

            with open(l.text, 'wb') as f:
                br.follow_link(l)
                f.write(br.response().read())

        str_coords = "%f, %f" % (coords[0], coords[1])

        logging.info('Querying GLEAM for {:.2f}d, {:.2f}d at {}...'.format(
                     coords[0], coords[1], b))

        br = Browser()
        br.set_handle_robots(False)
        br.addheaders = [('User-agent', 'Firefox')]
        if robust == -1:
            br.open('http://gleam-vo.icrar.org/gleam_postage/q/form')
        elif robust == 0:
            # no longer available?
            raise ValueError("no robust 0 image available.")
            # Requires username and password.
            br.add_password('http://gleam-vo.icrar.org/gleam_postage/q/form', \
                            robust0_uname, robust0_pword)
            br.open('http://gleam-vo.icrar.org/gleam_postage/q/form')

        br.select_form(nr=0)
        br['pos'] = str_coords
        br['size'] = str(fov)
        control1 = br.find_control('freq')
        for item in control1.items:
            if item.name == b:
                item.selected = True
        control2 = br.find_control('proj_opt')
        for item in control2.items:
            if item.name == projection:
                item.selected = True

        try:

            result = br.submit()
            content = result.read()

            with open('gleamer.html', 'wb') as f:
                f.write(content)

            for link in br.links():
                if link.text == 'FITS':
                    stamp_name = link.text
                    if not mask:
                        link.url = link.url+'&gleamer=1'
                        downloadlink(l=link)
                    else:
                        downloadlink(l=link)

            shutil.copy2(stamp_name, outdir+outname+'.fits')

            if fix:

                fix_minmax(outdir+outname+'.fits')
                fix_cd(outdir+outname+'.fits')

            os.remove(stamp_name)
            os.remove('gleamer.html')

            logging.info('GLEAM stamp acquired.')

        except Exception:
            logging.warn("GLEAM stamp for band {0} not acquired. HTTP Error 504 most likely.".format(b))
           
        
        br.close()




    if band == 'all':
        band = gleam_bands
    elif band == "narrow":
        band = narrow_bands

    # A number of checks are performed - if any fail the function will
    # terminate and return None.

    if not outname.endswith('.fits'): outname += '.fits'
    if not os.path.exists(outdir): os.makedirs(outdir)
    if not outdir.endswith('/'): outdir += '/'

    # Convert coordinates to [dd] - not necessary but it makes things easier.
    try: ra, dec = float(RA), float(DEC)
    except ValueError: ra, dec = hms_dms_dd(RA, DEC) 


    # Make sure projection is valid:
    if projection not in ['ZEA', 'ZEA_regrid', 'SIN_regrid']:
        raise ValueError('GLEAM: projection must be one of ' \
                         ' [`ZEA`, `ZEA_regrid`, `SIN_regrid`].')
    # Check to make sure the image won't be larger than 5.0 degrees:
    if float(fov) > 5.0:
        logging.warn('GLEAM: maximum size for GLEAM stamps is 5.0 degrees. \n'\
                     'GLEAM: setting size to 5.0 degrees.')
    # Check to see if the requested image is actually in the GLEAM survey area.
    if dec > 30.0:
        raise ValueError('GLEAM: no good GLEAM data above +30 degrees.')

    coords = (ra, dec)

    if isinstance(band, list):
        for b in band:

            # Make sure the band is valid:
            if b == 'white':
                b = '170-231'
            if b not in gleam_bands:
                raise ValueError('GLEAM: %s is not a valid band.' % b)

                # Check to see if the stamp already exists:
            if os.path.exists(outdir+outname[:-5]+'_'+b+'.fits'):
                if overwrite == 'r' or overwrite == False:
                    logging.warn('GLEAM: %s already exists. Passing.' % outname)
                elif overwrite == 'w' or overwrite == True:
                    logging.info('GLEAM: %s already exists - overwriting.' \
                                 % outname)
                    os.remove(outdir+outname[:-5]+'_'+b+'.fits')    

                    gleam_stamp(coords, fov, outdir, outname=outname[:-5]+'_'+b, \
                        b=b, projection=projection, overwrite=overwrite, mask=mask, \
                        fix=fix)
            else:
                 gleam_stamp(coords, fov, outdir, outname=outname[:-5]+'_'+b, \
                        b=b, projection=projection, overwrite=overwrite, mask=mask, \
                        fix=fix)

    else:
        if band == 'white':
            band = '170-231'
        gleam_stamp(coords, fov, outdir, outname[:-5]+'_'+band, \
            projection=projection, b=band, overwrite=overwrite, mask=mask, \
            fix=fix)





def TGSS(RA, DEC, fov, outdir, outname, overwrite='r'):
    '''Function to download a TGSS postage stamp(s) from the TGSS ADR server.

    Parameters
    ----------
    RA      : string or float
            Right ascension. Is converted to decimal degrees if not already.
    DEC     : string or float
            Declination. As above with RA, but must be > -55.0 degrees.
    fov     : float
            Field size. Must be < 1.0 degrees.
    outdir  : string
            Output directory. `outdir` is created if it does not exist.
    outname : string
            Output filename. If `outname` exists in `outdir` function stops.
    overwrite  : {'w', 'r'}, optional
               Determine whether to overwrite existing files of same name.

    '''


    def downloadlink(l):

        f = open(l.text, 'wb')
        br.follow_link(l)
        f.write(br.response().read())


    try: ra, dec = float(RA), float(DEC)
    except ValueError: ra, dec = hms_dms_dd(RA, DEC)


    # A number of checks are performed - if any fail the function will
    # terminate and return None.

    # Check to see if the stamp already exists:
    if os.path.exists(outdir+'/'+outname+'.fits'):
        if overwrite == 'r' or overwrite == False:
            print('"'+outname+'" already exists.')
            return None
        elif overwrite == 'w' or overwrite == True:
            logging.info('Existing "'+outname+'".fits will be overwritten.')
            os.remove(outdir+'/'+outname+'.fits')

    if dec < -55.0:
        raise ValueError('TGSS: %f, %f is out of the TGSS survey area.' % \
                         (ra, dec))

    # Check to make sure the image won't be larger than 1.0 degrees:
    if float(fov) > 1.0:
        raise ValueError('TGSS: Stamps size too big. Maximum size is 1 degree.')

    coords = '%f, %f' % (ra, dec)

    try:

        logging.info('Querying TGSS for {:.2f}d, {:.2f}d...'.format(ra, dec))

        br = Browser()

        br.set_handle_robots(False)
        br.addheaders = [('User-agent', 'Firefox')]

        br.open('http://vo.astron.nl/tgssadr/q_fits/cutout/form')

        br.select_form(nr=0)
        br['hPOS'] = coords
        br['hSIZE'] = str(fov)
        control = br.find_control('hINTERSECT')
        for item in control.items:
            if item.name == 'COVERS':
                item.selected = True

        result = br.submit()
        content = result.read()

        with open('tgss.html', 'wb') as f:
            f.write(content)

        for link in br.links():
            if '.FITS' in link.text:
                stamp_name = link.text
                downloadlink(l=link)
                break

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        shutil.copy2('./'+stamp_name, \
                    outdir+'/'+outname+'.fits')

        os.remove(stamp_name)
        os.remove('tgss.html')

        logging.info('TGSS stamp acquired.')

    except Exception:

        if os.path.exists('tgss.html'):
            os.remove('tgss.html')
        logging.warn('TGSS: something went wrong.')
        with open('problem_stamps.txt', 'a+') as f:
            f.write('No stamp in TGSS for recovered for: \n')
            f.write('{0}\n'.format(outname))


def NVSS_mosaic(RA, DEC, fov, outdir, outname, projection='TAN', overwrite='r', \
                montage=True):
    '''Make NVSS using montage.

    Parameters
    ----------
    RA         : string or float
               Right ascension. Must be hh:mm:ss, but is converted
               automatically with coord.py.
    DEC        : string or float
               Declination. Must be dd:mm:ss, but is converted automatically
               with coord.py. Must be > -40 degrees.
    fov        : float
               Field size. Must be < 2.0 degrees.
    outdir     : string
               Output directory. `outdir` is created if it does not exist.
    outname    : string
               Output filename. If `outname` exists in `outdir` function stops.
    projection : {'SIN', 'TAN', 'ARC', 'NCP', 'GLS', 'MER', 'AIT', 'STG'},
                optional
               Map projection. Default is 'TAN'.
    overwrite  : {'w', 'r'}, optional
               Determine whether to overwrite existing files of same name.

    '''


    try: ra, dec = float(RA), float(DEC)
    except ValueError: ra, dec = hms_dms_dd(RA, DEC) 


    if dec < -40:
        raise ValueError('>>> NVSS has no data below -40 declination.')

    fov = fov + 2.0
    size = (fov / 2.0)

    coords = {dec: [ra]}

    n_overlaps = int(math.ceil((size - 0.5) / 0.5))

    for i in range(n_overlaps):

        coords[dec+0.5*(i+1)] = [ra]
        coords[dec-0.5*(i+1)] = [ra]

    for d in coords.keys():

        for i in range(n_overlaps):

            if (ra + (0.5*(i+1) / numpy.cos(numpy.radians(d)))) >= 360.0:

                 coords[d].append((ra + (0.5*(1 + i) / \
                                  numpy.cos(numpy.radians(d)))) - 360.0)

            else:

                coords[d].append(ra + (0.5*(1 + i) / \
                                 numpy.cos(numpy.radians(d))))


            if (ra - (0.3*(1 + i) / numpy.cos(numpy.radians(d)))) < 0.0:

                coords[d].append(ra - (0.5*(1 + i) / \
                                 numpy.cos(numpy.radians(d))) + 360.0)

            else:

                coords[d].append(ra - (0.5*(1 + i) / \
                                 numpy.cos(numpy.radians(d))))     


    logging.info('>>> Downloading %i stamps to form an NVSS mosaic.' % \
                  int(n_overlaps*n_overlaps + 1))

    image_count = 0

    for d in coords.keys():

        for i in range(n_overlaps):

            NVSS(coords[d][i], d, 2.0, outdir, '%s_%i_%.2f' \
                 % (outname, i, round(d, 1)), projection, overwrite=overwrite, fix=True)

            image_count += 1

    print(image_count)

    if montage:

        work_dir = os.path.abspath(outdir)
        raw_dir  = '%s/raw' % work_dir
        prj_dir  = '%s/projected' % work_dir

        if not os.path.exists(work_dir):
            os.mkdir(work_dir)
        if not os.path.exists(raw_dir): 
            os.mkdir(raw_dir)
        else: 
            logging.warning('>>> Raw directory already exists: working in there.') 
        if not os.path.exists(prj_dir): 
            os.mkdir(prj_dir)
        else: 
            logging.warning('>>> Projected directory already exists: '\
                            'working in there.')  

        for spec in os.listdir(outdir):
            if spec.endswith('.fits'):
                shutil.copy2('%s/%s' % (outdir, spec), \
                             '%s/%s' % (raw_dir, spec))

        logging.info('>>> NVSS: making optimal header...')
        Popen('mImgtbl %s %s/rimages.tbl' % \
             (raw_dir, work_dir), shell=True).wait()
        Popen('mMakeHdr %s/rimages.tbl %s/%s.hdr' % \
             (work_dir, work_dir, outname), shell=True).wait()
        logging.info('>>> NVSS: re-projecting to optimal header...')
        Popen('mProjExec -p %s %s/rimages.tbl %s/%s.hdr %s %s/stats.tbl' % \
             (raw_dir, work_dir, work_dir, outname, prj_dir, work_dir), \
             shell=True).wait()
        Popen('mImgtbl %s %s/pimages.tbl' % (prj_dir, work_dir), \
              shell=True).wait()

        final_dir = prj_dir
        final_tbl = '%s/pimages.tbl' % work_dir

        logging.info('>>> NVSS: adding images to form mosaic...')
        Popen('mAdd -e -p %s %s %s/%s.hdr %s/%s.fits' % \
             (final_dir, final_tbl, work_dir, outname, work_dir, outname), \
             shell=True).wait()

        if cleanup:

            logging.info('>>> NVSS: deleting extra FITS directories.')
            shutil.rmtree(raw_dir)
            shutil.rmtree(prj_dir)

        logging.info('>>> NVSS: trimming image...')
        Popen('mSubimage %s/%s.fits %s/%s.fits %f %f %f %f' % \
              (work_dir, outname, work_dir, outname, ra, dec, \
              (fov-0.4), (fov-0.4)), shell=True).wait()

        logging.info('>>> NVSS mosaic complete.')




def NVSS(RA, DEC, fov, outdir, outname, projection='TAN', overwrite='r', fix=True):
    '''Function to download NVSS postage stamp(s).

    Projection cannot be 'NCP' or 'SIN' otherwise annotations in kvis will
    be all over the place if not in the image.

    Will use SkyView if size is > 2 degrees. (One day...) 

    Parameters
    ----------
    RA         : string or float
               Right ascension. Must be hh:mm:ss, but is converted
               automatically with coord.py.
    DEC        : string or float
               Declination. Must be dd:mm:ss, but is converted automatically
               with coord.py. Must be > -40 degrees.
    fov        : float
               Field size. Must be < 2.0 degrees.
    outdir     : string
               Output directory. `outdir` is created if it does not exist.
    outname    : string
               Output filename. If `outname` exists in `outdir` function stops.
    projection : {'SIN', 'TAN', 'ARC', 'NCP', 'GLS', 'MER', 'AIT', 'STG'},
                optional
               Map projection. Default is 'TAN'.
    overwrite  : {'w', 'r'}, optional
               Determine whether to overwrite existing files of same name.

    '''


    def fix_NVSS(image):
        '''Add clean BMAJ and BMIN to header.'''

        hdulist = fits.open(image, mode='update')
        hdulist[0].header['BMAJ'] = 1.25e-2
        hdulist[0].header['BMIN'] = 1.25e-2
        hdulist[0].header['BPA']  = 0.0

        hdulist.flush()
        hdulist.close()

    def downloadlink(l):

        f = open(l.text, 'wb')
        br.follow_link(l)
        f.write(br.response().read())


    # A number of checks are performed - if any fail the function will
    # terminate and return None.

    # Check to see if the stamp already exists:
    if os.path.exists(outdir+'/'+outname+'.fits'):
        if overwrite == 'r' or overwrite == False:
            raise IOError('NVSS: %s already exists.' % outname)
        elif overwrite == 'w' or overwrite == True:
            logging.info('NVSS: Existing %s already exists - overwriting.' \
                % outname)
            os.remove(outdir+'/'+outname+'.fits')

    # Make sure projection is valid:
    if (fov <= 2.0) and (projection not in ['SIN', 'TAN', 'ARC', 'NCP', 'GLS', \
        'MER', 'AIT', 'STG']):
        raise ValueError('NVSS: %s is an invalid projection for NVSS at this size.' \
                         % projection)
    elif (fov > 2.0) and (projection not in ["SIN", "TAN", "CAR", "AIT", "ZEA"]):
        raise ValueError("NVSS: %s is an invalid projection for NVSS stamps" \
                         " downloaded from SkyView." % projection)

    # Make sure the coordinates are in HMS, DMS format and in survey region:
    try:
        RA_test = float(RA)
        ra, dec = dd_hms_dms(RA, DEC)
    except ValueError:
        ra, dec = RA, DEC
        RA, DEC = hms_dms_dd(RA, DEC)

    if '-' or '+' in dec:
        DEC_limit = float(dec[0:3])
    else:
        DEC_limit = float(dec[0:2])

    if DEC_limit < -40.0:
        raise ValueError('NVSS: %s, %s is out of the survey area.' \
                         % (ra, dec))

    # Check to make sure the image won't be larger than 2.0 degrees:
    if 2.0 < float(fov) <= 5.0:
        logging.warn("NVSS: Using SkyView for NVSS stamp.")
        skyview = True
    elif float(fov) > 5.0:
        raise ValueError('NVSS: maximum size for NVSS stamps is 5.0 degrees.')
    else:
        skyview = False

    if skyview:
        raise ValueError("NVSS: SkyView not currently supported.")

    if not outdir.endswith('/'): outdir += '/'
    if not os.path.exists(outdir): os.makedirs(outdir)

    try:

        logging.info('NVSS: querying postage stamp server for {:.2f}d, {:.2f}d...'.format(RA, DEC))

        radius = '%f, %f' % (fov, fov)

        br = Browser()

        # NVSS doesn't like robots:
        br.set_handle_robots(False)
        br.addheaders = [('User-agent', 'Firefox')]

        if not skyview:
            # Access the NVSS postage stamp server:
            br.open('http://www.cv.nrao.edu/nvss/postage.shtml')
            br.select_form(nr=0)
            br['RA'] = ra
            br['Dec'] = dec
            br['Size'] = radius
            control1 = br.find_control('Type')
            for item in control1.items:
                if item.name == 'application/octet-stream':
                    item.selected = True
            control2 = br.find_control('MAPROJ')
            for item in control2.items:
                if item.name == projection:
                    item.selected = True

            result = br.submit()
            content = result.read()

            with open(outdir+outname+'.fits', 'wb') as f:
                f.write(content)

        else:
            # Access SkyView to get large NVSS image:
            br.open("https://skyview.gsfc.nasa.gov/current/cgi/query.pl")
            br.select_form(nr=0)
            coord_string = "{0}, {1}".format(ra, dec)
            br["Position"] = coord_string
            control1 = br.find_control("skv_survey", id="Radio:GHz")
            for item in control1.items:
                if item.name == "NVSS":
                    item.selected = True
            control2 = br.find_control("projection")
            for item in control2.items:
                if item.name == projection.title():
                    item.selected = True
            imsize = 3.0*fov/0.0125
            br["pixels"] = "{0}".format(imsize)
            br["imscale"] = "{0}".format(fov)

            br.submit(nr=0, type="submit")
            # br.click(type="submit")
            # return result
            for link in br.links():
                print(link.text)
            return None
            # content = result.read()

            with open("skyview.html", "wb") as f:
                f.write(content)

            for link in br.links():
                if "FITS" in link.text:
                    stamp_name = link.text
                    downloadlink(l=link)

            shutil.copy2(stamp_name, outdir+outname+".fits")
            os.remove(stamp_name)
            os.remove("skyview.html")


        if fix:

            fix_NVSS(outdir+outname+'.fits')

        logging.info('>>> NVSS stamp acquired.')

    except Exception:
        raise
        # logging.warn('>>> NVSS: something went wrong.')


def SUMSS(RA, DEC, fov, outdir, outname, projection='SIN', overwrite='r', \
    samedir=False):
    '''Function to download SUMSS postage stamp(s).

    Parameters
    ----------
    RA         : string or float
               Right ascension. Must be hh:mm:ss, but is converted
               automatically with coord.py.
    DEC        : string or float
               Declination. Must be dd:mm:ss, but is converted automatically
               with coord.py. Must be < -30 degrees.
    fov        : float
               Field size. Must be < 2.0 degrees.
    outdir     : string
               Output directory. `outdir` is created if it does not exist.
    outname    : string
               Output filename. If `outname` exists in `outdir` function stops.
    projection : {'SIN', 'TAN', 'ARC', 'NCP', 'GLS', 'AIT', 'STG'}, optional
               Map projection. Default is 'TAN'.
    overwrite  : {'w', 'r'}, optional
               Determine whether to overwrite existing files of same name.
    samedir    : bool, optional
               If the SUMSS stamp is to go to the working directory, this 
               should be set to True. (This is needed for the shutil functions.)
    '''


    def downloadlink(l):

        f = open(l.text, 'wb')
        br.follow_link(l)
        f.write(br.response().read())


    # A number of checks are performed - if any fail the function will
    # terminate and return None.

    if not samedir:
        if not outdir.endswith('/'):
            outdir += '/'

    # Check to see if the stamp already exists:
    if os.path.exists(outdir+outname+'.fits'):
        if overwrite == 'r' or overwrite == False:
            raise IOError(outname+' already exists.')
        elif overwrite == 'w' or overwrite == True:
            logging.info('SUMSS: Existing '+outname+' will be overwritten.')
            os.remove(outdir+outname+'.fits')

    # Make sure projection is valid:
    if projection not in ['NCP', 'SIN', 'TAN', 'TAN', 'ARC', 'GLS', 'AIT', \
        'STG']:
        raise ValueError(\
            '{0} is an invalid projection for SUMSS. Type help(SUMSS)'
            ' for options.'.format(projection))

    # Make sure the coordinates are in HMS, DMS format and in survey region:
    try:
        RA_test = float(RA)
        RA, DEC = dd_hms_dms(RA, DEC)
        ra, dec = float(RA), float(DEC)
    except ValueError:
        ra, dec = hms_dms_dd(RA, DEC)

    if '-' or '+' in DEC:
        DEC_limit = float(DEC[0:3])
    else:
        DEC_limit = float(DEC[0:2])

    if DEC_limit > -30.0:
        raise DeclinationError(\
            '{0}, {1} is out of SUMSS\'s survey area: "{2}" \n' \
            'No SUMSS data above DEC -30 degrees.'.format(RA, DEC, outname), DEC_limit)

    # Check to make sure the image won't be larger than 5.0 degrees:
    if float(fov) > 2.0:
        raise ValueError(\
            '{0} is too large for SUMSS: {1} \n' \
            'SUMSS stamps cannot be larger than 2.0 (ish) '\
            'degrees.'.format(fov, outname))

    try:

        logging.info('Querying SUMSS for {:.2f}d, {:.2f}d ...'.format(
            ra, dec))
        radius = str('{0} {0}'.format(fov, fov))

        br = Browser()

        # SUMSS does not like robots:
        br.set_handle_robots(False)
        br.addheaders = [('User-agent', 'Firefox')]

        br.open('http://www.astrop.physics.usyd.edu.au/cgi-bin/postage.pl')
        br.select_form(nr=0)
        br['RA'] = RA
        br['DEC'] = DEC
        br['fieldsize'] = radius
        control = br.find_control('projection')
        for item in control.items:
            if item.name == projection:
                item.selected = True

        result = br.submit()
        content = result.read()

        with open('sumss.html', 'wb') as f:
            f.write(content)

        for link in br.links():
            if 'fits' in link.text:
                stamp_name = link.text
                downloadlink(l=link)

        if not samedir:

            if not os.path.exists(outdir):
                os.makedirs(outdir)

            shutil.copy2(stamp_name, outdir+outname+'.fits')

        else:

            os.rename(stamp_name, outname+'.fits')

        os.remove(stamp_name)
        os.remove('sumss.html')

        logging.info('SUMSS stamp acquired.')

    except Exception:
        raise

        if os.path.exists('sumss.html'):
            os.remove('sumss.html')
        print('Something went wrong: '+outname+' is ' \
              'probably out of the survey area.')


def DSS2(RA, DEC, fov, outdir, outname, band='b', overwrite='r'):
    '''Function to download DSS2 postage stamp(s).

    Parameters
    ----------
    RA         : string or float
               Right ascension. Can be literally any format.
    DEC        : string or float
               Declination. Can be literally any format.
    fov        : float
               Field size. Must be < 1.0 degrees.
    outdir     : string
               Output directory. `outdir` is created if it does not exist.
    outname    : string
               Output filename. If `outname` exists in `outdir` function stops.
    band       : string, optional
               Color band for DSS2. Options are blue, red, or IR.
    overwrite  : {'w', 'r'}, optional
               Determine whether to overwrite existing files of same name.

    '''

    # A number of checks are performed - if any fail the function will
    # terminate and return None.

    # Check to see if the stamp already exists:
    if os.path.exists(outdir+'/'+outname+'.fits'):
        if overwrite == 'r' or overwrite == False:
            logging.info('"'+outname+'" already exists.')
            return None
        elif overwrite == 'w' or overwrite == True:
            logging.info('Existing "'+outname+'".fits will be overwritten.')
            os.remove(outdir+'/'+outname+'.fits')

    # Check to make sure the image won't be larger than 1.0 degrees:
    if float(fov) > 1.0:
        raise ValueError('{0} is too large for DSS2 ({1}). DSS2 stamps cannot ' \
                         'be larger than 1.0 degrees.'.format(fov, outname))
        return None

    # Determine which band to use (blue or red):
    if band in ['b', 'B', 'blue', 'Blue', 'BLUE', "DSS2B"]:
        band = 'poss2ukstu_blue'
        mess = 'blue'
    elif band in ['r', 'R', 'red', 'Red', 'RED', "DSS2R"]:
        band = 'poss2ukstu_red'
        mess = 'red'
    elif band in ['ir', 'iR', 'IR', 'Ir', 'infrared', 'Infrared', 'INFRARED', "DSS2IR"]:
        band = 'poss2ukstu_ir'
        mess = 'infrared'
    else:
        raise ValueError('DSS2 band must be "blue", "red", or "infrared".')

    try:

        logging.info('Querying DSS2...')

        br = Browser()

        br.set_handle_robots(False)
        br.addheaders = [('Users-agent', 'Firefox')]

        br.open('http://archive.stsci.edu/cgi-bin/dss_form')
        br.select_form(nr=1)
        control = br.find_control('v')
        for item in control.items:
            if item.name == band:
                item.selected = True
        br['r'] = str(RA)
        br['d'] = str(DEC)
        br['h'] = str(fov*60.0)
        br['w'] = str(fov*60.0)
        br.find_control('s').items[0].selected = True

        result = br.submit()
        content = result.read()

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        with open(outdir+'/'+outname+'.fits', 'wb') as f:
            f.write(content)

        logging.info('DSS2 ' + mess + ' stamp acquired.')

    except Exception as e:

        print(e)
        logging.info("Something went wrong. DSS2 covers most of the sky, so "
                     "this could be a connectivity problem or you are "
                     "somewhere near the galactic centre.")
        pass


    return None
