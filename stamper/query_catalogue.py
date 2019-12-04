import numpy as np

from mechanize import Browser
import shutil
import os

from astroquery.irsa import Irsa
from astroquery.vizier import Vizier
from astropy import units as u
from astropy.constants import c

from stamper.coord import dd_hms_dms, hms_dms_dd

from curve_fitting import loglog

import logging
logging.basicConfig(format='%(levelname)s (%(module)s): %(message)s', \
    level=logging.INFO)



def angular_distance(coords1, coords2):
    '''Get the angular distance between a set of RA, DEC coordinates in [dd].'''

    cos_angle = math.sin(math.radians(coords1[1])) * \
            math.sin(math.radians(coords2[1])) + \
            math.cos(math.radians(coords1[1])) * \
            math.cos(math.radians(coords2[1])) * \
            math.cos(math.radians(coords1[0] - coords2[0]))

    try:
        gamma = math.degrees(math.acos(cos_angle))
    except ValueError:
        gamma = math.degrees(math.acos(min(1, max(cos_angle, -1))))

    return gamma


def get_galex_results(coords, search_limit=2.0, verbose=True, flux=True):
    '''Query GALEX through VizieR and return FUV and NUV results.

    Optionally get measurements in Jy.

    Parameters:
    -----------
    coords  : tuple
    search_limit : float, optional
                   Search radius in arcmin.
    verbose : bool, optional
              If True magnitudes (and fluxes) are printed to the terminal.
    flux    : bool, optional
              If True, calculate flux densities in Jy.

    '''

    try: ra, dec = float(coords[0]), float(coords[1])
    except ValueError: ra, dec = hms_dms_dd(coords[0], coords[1])

    logging.info('>>> Querying VizieR for GALEX results...')

    try:
        table = Vizier.query_region('%f %f' % (ra, dec), \
                                    radius=search_limit*u.arcmin, \
                                    catalog='II/312/ais')[0]

    
        flag = False
        for i in range(0, len(table['b'])):

            if table['b'][i] == 3:

                mNUV, mFUV, emNUV, emFUV = table['NUV'][i], table['FUV'][i], \
                                           table['e_NUV'][i], table['e_FUV'][i]
                fNUV, fFUV, efNUV, efFUV = table['Nflux'][i], table['Fflux'][i], \
                                           table['e_Nflux'][i], table['e_Fflux'][i]

                flag = True
                break


        if not flag:

            mNUV, mFUV, emNUV, emFUV = table['NUV'][0], table['FUV'][0], \
                                       table['e_NUV'][0], table['e_FUV'][0]
            fNUV, fFUV, efNUV, efFUV = table['Nflux'][0], table['Fflux'][0], \
                                       table['e_Nflux'][0], table['e_Fflux'][0]





        if verbose:
            print('MAGNITUDES: \n' \
                  '       NUV = %f +/- %f mag \n' \
                  '       FUV = %f +/- %f mag \n' \
                  '    FLUXES: \n' \
                  '       NUV = %f +/- %f mag \n' \
                  '       FUV = %f +/- %f mag \n' % \
                  (mNUV, emNUV, mFUV, emFUV, fNUV, efNUV, fFUV, efFUV))

        return mNUV, mFUV, emNUV, emFUV, fNUV, fFUV, efNUV, efFUV

    except IndexError:
        raise


def get_wise_results(coords=None, search_limit=2.0, flux=False, correct=True, \
                    verbose=False, W=None, eW=None):
    """Get all 4 WISE band measurements in MAG (and flux densitiy).

    """


    if W is None:


        try: ra, dec = float(coords[0]), float(coords[1])
        except ValueError: ra, dec = coords[0], coords[1]

        logging.info(">>> Querying VizieR for WISE results...")

        table = Vizier.query_region('%f %f' % (ra, dec), \
                                    radius=search_limit*u.arcmin, \
                                    catalog='II/328/allwise')[0]

        flag = False
        for i in range(0, len(table['AllWISE'])):
            if 'U' not in table['qph'][i]:

                W1, E1, W2, E2, W3, E3, W4, E4 = table['W1mag'][i], table['e_W1mag'][i], \
                                                 table['W2mag'][i], table['e_W2mag'][i], \
                                                 table['W3mag'][i], table['e_W3mag'][i], \
                                                 table['W4mag'][i], table['e_W4mag'][i]

                flag = True
                break

        if not flag:
            raise ValueError(">>> No results found with all bands.")

        W = [W1, W2, W3, W4]
        eW = [E1, E2, E3, E4]    

    if verbose:
        for i in range(len(W)):
            print("W%i = %.2f +/- %.2f" % (i+1, W[i], eW[i]))

    if flux:

        f0 = [306.682, 170.663, 29.0448, 8.2839]
        ef0 = [4.600, 2.600, 0.436, 0.124]

        if correct:

            f_temp = []
            for i in range(0, 4):
                f_temp.append(f0[i]*10**(-0.4*W[i]))

            m, m_unc, cc, c_unc = loglog(x=[c.value/3.4e-6, c.value/4.6e-6, \
                c.value/12.0e-6, c.value/22.0e-6], y=f_temp, \
                unc=[0.05*f_temp[0], 0.05*f_temp[1], 0.05*f_temp[2], \
                0.05*f_temp[3]], \
                absolute_unc=True, show=False)

            if -0.5 < m < 0.5:
                fc1, fc2, fc3, fc4 = 0.9907, 0.9935, 0.9169, 0.9905
            elif -1.5 < m <= -0.5:
                fc1, fc2, fc3, fc4 = 0.9921, 0.9943, 0.9373, 0.9926
            elif -2.5 < m <= -1.5:
                fc1, fc2, fc3, fc4 = 1.0, 1.0, 1.0, 1.0
            elif -3.5 < m <= -2.5:
                fc1, fc2, fc3, fc4 = 1.0142, 1.0107, 1.1081, 1.0130
            elif m <= -3.5:
                fc1, fc2, fc3, fc4 = 1.0347, 1.0265, 1.2687, 1.0319
            elif 0.5 <= m < 1.5:
                fc1, fc2, fc3, fc4 = 0.9961, 0.9976, 0.9393, 0.9934
            elif 1.5 <= m < 2.5:
                fc1, fc2, fc3, fc4 = 1.0084, 1.0066, 1.0088, 1.0013
            elif m >= 2.5:
                fc1, fc2, fc3, fc4 = 1.0283, 1.0206, 1.1344, 1.0142

        else:
            fc1, fc2, fc3, fc4 = 1.0, 1.0, 1.0, 1.0

        fc = [fc1, fc2, fc3, fc4]
        F = []
        eF = []

        for i in range(0, 4):

            F.append(fc[i] * f0[i] * 10**(-0.4*W[i]))
            eF.append(\
                fc[i] * 10**(-0.4*W[i]) * np.sqrt(\
                (ef0[i]**2) + (np.log(10**(-0.4)) * f0[i] * eW[i])**2))

        if verbose:
            for i in range(len(F)):
                print("F%i = %.2f +/- %.2f [mJy]" % (i+1, F[i]*1000.0, eF[i]*1000.0))

    else:
        F, eF = [], []


    return W, eW, F, eF



def get_galex_mag(coords, search_limit=2.0, verbose=True, flux=False, \
                  path_to_galex='/Users/duchesst/GOOGLE/scripts/galex_ais_full.tsv'):
    '''Get FUV and NUV GALEX AIS data for a source.

    Optionally get measurements in Jy.

    Parameters:
    -----------
    coords  : tuple
    search_limit : float, optional
                   Search radius in arcmin.
    verbose : bool, optional
              If True magnitudes (and fluxes) are printed to the terminal.
    flux    : bool, optional
              If True, calculate flux densities in Jy.
    path_to_galex : str
              Location of the GALEX AIS table from VizieR.

    '''

    try: ra, dec = float(coords[0]), float(coords[1])
    except ValueError: ra, dec = coords[0], coords[1]

    galex = np.genfromtxt(path_to_galex, \
                          names='_RAJ2000,_DEJ2000,RAJ2000,DEJ2000,FUV,e_FUV,' \
                          'NUV,e_NUV,Fflux,e_Fflux,Nflux,e_Nflux,Fr,Nr')

    search = True
    while search:

        c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        galex_catalogue = SkyCoord((galex['RAJ2000'], galex['DEJ2000']), \
                                   unit=(u.deg, u.deg))
        i = c.match_to_catalog_sky(galex_catalogue)[0]

        dist = angular_distance((ra, dec), (galex['RAJ2000'][i], \
                                galex['DEJ2000'][i])) * 60.0

        if (galex['b'][i] == 3) and (dist <= search_radius):

            # 3 -> both FUV and NUV measurements exits:
            search = False

        elif dist <= search_radius:

            # We assume that if we don't have have both NUV and FUV that our 
            # source is actually further down the list.
            galex = np.delete(galex, (i), axis=0)

        else:

            raise ValueError('>>> No source found with both NUV and FUV within' \
                 ' %f arcmin of search coordinates.' % search_radius)


    return galax[:, i]







def get_wise_mag(coords, verbose=True, flux=False, correct=True):
    '''Get all 4 WISE band measurements in MAG.

    Optionally get measurements in Janskies, corrected based on the SED
    between the 4 band measurements.

    Parameters:
    coords  : tuple
    verbose : bool, optional
              If True magnitudes (and fluxes) are printed to the terminal.
    flux    : bool, optional
              If True, calculate flux densities in Jy.
    correct : bool, optional
              If True, attempt to apply corrections to fluxes based on SED.

    '''

    try: ra, dec = dd_hms_dms(coords[0], coords[1])   # Have to be in hms:dms.
    except TypeError: ra, dec = coords[0], coords[1] 

    if verbose:
        print('Looking around {0}, {1}...'.format(ra, dec))

    try:

        table = Irsa.query_region(ra+' '+dec, catalog='allwise_p3as_psd', \
            spatial='Box', width=1*u.arcsec)

        W1, W2, W3, W4 = table['w1mpro'][0], table['w2mpro'][0], table['w3mpro'][0], \
            table['w4mpro'][0]
        eW1, eW2, eW3, eW4 = table['w1sigmpro'][0], table['w2sigmpro'][0], \
            table['w3sigmpro'][0], table['w4sigmpro'][0]

    except IndexError:

        logging.warning('No sources within 1 arcsec of input coordinates.')
        logging.warning('Looking within a 2 arcsec box...')

        table = Irsa.query_region(ra+' '+dec, catalog='allwise_p3as_psd', \
            spatial='Box', width=2*u.arcsec)

        W1, W2, W3, W4 = table['w1mpro'][0], table['w2mpro'][0], table['w3mpro'][0], \
            table['w4mpro'][0]
        eW1, eW2, eW3, eW4 = table['w1sigmpro'][0], table['w2sigmpro'][0], \
            table['w3sigmpro'][0], table['w4sigmpro'][0]


    W = [W1, W2, W3, W4]
    eW = [eW1, eW2, eW3, eW4]

    if verbose:
        for i in range(0, 4):
            print('W{0} = {1} +/- {2}'.format(i+1, W[i], eW[i]))

    if flux:

        f0 = [306.682, 170.663, 29.0448, 8.2839]
        ef0 = [4.600, 2.600, 0.436, 0.124]

        if correct:

            f_temp = []
            for i in range(0, 4):
                f_temp.append(f0[i]*10**(-0.4*W[i]))

            m, m_unc, cc, c_unc = loglog(x=[c.value/3.4e-6, c.value/4.6e-6, \
                c.value/12.0e-6, c.value/22.0e-6], y=f_temp, \
                unc=[0.05*f_temp[0], 0.05*f_temp[1], 0.05*f_temp[2], \
                0.05*f_temp[3]], \
                absolute_unc=True, show=False)

            # p = np.polyfit(x=[c.value/3.4e-6, c.value/4.6e-6, \
            #   c.value/12.0e-6, c.value/22.0e-6], y=f_temp, deg=1)
            # m = p[0]

            if -0.5 < m < 0.5:
                fc1, fc2, fc3, fc4 = 0.9907, 0.9935, 0.9169, 0.9905
            elif -1.5 < m <= -0.5:
                fc1, fc2, fc3, fc4 = 0.9921, 0.9943, 0.9373, 0.9926
            elif -2.5 < m <= -1.5:
                fc1, fc2, fc3, fc4 = 1.0, 1.0, 1.0, 1.0
            elif -3.5 < m <= -2.5:
                fc1, fc2, fc3, fc4 = 1.0142, 1.0107, 1.1081, 1.0130
            elif m <= -3.5:
                fc1, fc2, fc3, fc4 = 1.0347, 1.0265, 1.2687, 1.0319
            elif 0.5 <= m < 1.5:
                fc1, fc2, fc3, fc4 = 0.9961, 0.9976, 0.9393, 0.9934
            elif 1.5 <= m < 2.5:
                fc1, fc2, fc3, fc4 = 1.0084, 1.0066, 1.0088, 1.0013
            elif m >= 2.5:
                fc1, fc2, fc3, fc4 = 1.0283, 1.0206, 1.1344, 1.0142

            if verbose: print('SED is calculated to be F ~ v^{0}'.format(m))

        else:

            fc1, fc2, fc3, fc4 = 1.0, 1.0, 1.0, 1.0

        fc = [fc1, fc2, fc3, fc4]
        F = []
        eF = []

        for i in range(0, 4):

            F.append(fc[i] * f0[i] * 10**(-0.4*W[i]))
            eF.append(\
                fc[i] * 10**(-0.4*W[i]) * np.sqrt(\
                (ef0[i]**2) + (np.log(10**(-0.4)) * f0[i] * eW[i])**2))

        if verbose:
            for i in range(0, 4):
                print('F{0} = {1} +/- {2} [Jy]'.format(i+1, F[i], eF[i]))

    else:

        F, eF = [], []


    return W, eW, F, eF
