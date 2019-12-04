'''stamper is for downloading postage stamps from various web servers.

The following postage stamp servers are currently supported:
    GLEAM      -- GaLactic and Extragalactic All-sky MWA survey
    TGSS ADR   -- TFIR GMRT Sky Survey Alternate Data Release
    NVSS       -- NRAO VLA Sky Survey
    SUMSS      -- Sydney Univeristy Molonglo Sky Survey
    DSS2       -- Digitized Sky Survey 2 (blue, red, or IR)
Servers that should be added (but probably won't!):
    RASS       -- ROSAT X-Ray All-Sky Survey
    VLSSr      -- VLA Low-frequency Sky Survey redux.

Also included are the functions "hms_dms_dd" and "dd_hms_dms" which convert
RA and DEC between HMS, DMS and decimal degree formats. These are needed as
some postage stamp servers will only take decimal degrees (e.g., TGSS ADR) and
others will only take HMS, DMS (e.g., NVSS and SUMSS). 

PACKAGE REQUIREMENTS
--------------------
mechanize.Browser
numpy
math
shutil
os

TODO:
 -- include basic NED query to make annotation files for kvis and DS9
 -- include WISE postage stamps

'''

__all__ = ['hms_dms_dd', 'dd_hms_dms', 'GLEAM', 'TGSS', 'NVSS', 'SUMSS', \
    'DSS2', 'PS1', 'PS1_mosaic', 'cleanup']

# The coord module needs to be accessible to the postage_stamps module.
# from coord import hms_dms_dd, dd_hms_dms

# For convenience.
# from postage_stamps import GLEAM, TGSS, NVSS, SUMSS, DSS2, PS1, PS1_mosaic, cleanup


