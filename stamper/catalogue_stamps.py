# Use a VizieR-queried catalogue as a reference list to download postage stamps
# around each entry in the catalogue.
#
# Requires:
#    - astroquery
#    - astropy
#    - numpy
# ------------------------------------------------------------------------------


from stamper.postage_stamps import GLEAM, NVSS, SUMSS, TGSS, DSS2
from stamper.coord import hms_dms_dd, dd_hms_dms

import numpy

import logging
logging.basicConfig(format='%(levelname)s (%(module)s): %(message)s', \
	level=logging.INFO)

from astroquery.vizier import Vizier


def get_catalogue(catalogue, row_limit=-1):
	'''Queries VizieR for the catalogue given and returns the Table object.

	Parameters
	----------
	catalogue : string
				The catalogue code, e.g. "I/322A/out". 

	Returns
	-------
	cat       : astropy.table.Table
	'''

	Vizier.ROW_LIMIT = -1  # To give all results of a table. Default is 50.

	cats = Vizier.get_catalogs(catalogue)
	try:
		cat = cats[catalogue]
	except KeyError:
		raise KeyError(' (catalogue_stamps): '\
			'The catalogue name/index is not correct. \n' \
			'    {0} is not a valid catalogue name.'.format(catalogue))

	return cat



def get_stamps(catalogue, name, outdir='./', gleam=False, nvss=False, \
	sumss=False, tgss=False, dss2=False, row_limit=-1):
	'''Download postage stamps relevant to each entry in a VizieR-queried table.

	Parameters
	----------
	catalogue : string
				The catalogue code, e.g. "I/322A/out". 
	name      : string
				Name scheme for output files. This string will be used to search
				for what column in the VizieR table gives the name.
	outdir    : string, optional
			    Output directory. Note that each set of stamps is output into 
			    their own directory structures within this output directory.
	gleam     : {True, False, string}, optional
			    If string then the string should specify the GLEAM band to use.
	nvss      : bool, optional
	sumss     : bool, optional
	tgss      : bool, optional
	dss2      : {True, False, string}, optional
				If string then the string should specify the DSS2 band to use.
	row_limit : int, optional
				Number of rows for VizieR to return. -1 for infinite.
	'''


	if not outdir.endswith('/'): outdir += '/'


	cat = get_catalogue(catalogue, row_limit)

	# Try to determine the naming scheme to use for output files:
	for n in cat.colnames:
		if name == n: flag = True
	if not flag:
		for n in cat.colnames:
			if name in n: 
				flag = True
				name = n


	RA, DEC = cat['_RAJ2000'], cat['_DEJ2000']  # These are VizieR computed.

	for i in range(0, len(RA)):

		filename = name + cat[name][i]
		logging.info('{0}: Starting...'.format(filename))
		# print 'INFO: {0}: Starting...'.format(filename)

		if nvss:
			try: NVSS(RA[i], DEC[i], 1.0, outdir+'nvss', filename+'_nvss')
			except Exception: pass

		if sumss:
			try: SUMSS(RA[i], DEC[i], 1.0, outdir+'sumss', filename+'sumss')
			except Exception: pass

		if tgss:
			try: TGSS(RA[i], DEC[i], 1.0, outdir+'tgss', filename+'_tgss')
			except Exception: pass

		if gleam != False:
			if gleam == 'all':
				out = outdir+'gleam/'+filename
				band = 'all'
			elif gleam == True:
			 	out = outdir+'gleam'
			 	band = '170-231'
			else:
			 	out = outdir+'gleam'
			 	band = gleam
			try: GLEAM(RA[i], DEC[i], 3.0, out, filename, band=band, \
				overwrite='r')
			except Exception: pass

		if dss2 != False:
			if dss2 == True:
				band = 'b'
			else:
				band = dss2
			DSS2(RA[i], DEC[i], 1.0, outdir+'dss', filename+'_'+banb, band=band)

		logging.info('{0}: Finished.'.format(filename))








