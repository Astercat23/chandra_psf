
import numpy as np
from astropy.io import fits
from astropy.io import ascii
import errors as err

## Updated 2018.02.14 for general use
#----------------------------------------------
## The Profile object

class RadProfile(object):
    def __init__(self, from_fits_file=None, key='SUR_BRI'):
        self.rleft  = 0.0
        self.rright = 0.0
        self.value  = 0.0
        self.value_err  = 0.0
        self.r_unit = None
        self.value_unit = None
        if from_fits_file is not None:
            assert isinstance(from_fits_file, str)
            assert isinstance(key, str)
            self._read_from_fits(from_fits_file, key)
        return

    def _read_from_fits(filename, key):
        hdu_list = fits.open(filename)
        data     = hdu_list[1].data
        self.rleft = data['R'][:,0]
        self.rright = data['R'][:,1]
        self.surbri = data[key]
        self.surbri_err = data[key + '_ERR']
        self.value_unit = data.columns[key].unit
        return

    @property
    def rmid(self):
        return 0.5 * (self.rleft + self.rright)

    @property
    def area(self):
        return np.pi * (self.rright**2 - self.rleft**2) # pix^2

    def __getslice__(self, i,j):
        result = RadProfile()
        result.rleft  = self.rleft[i:j]
        result.rright = self.rright[i:j]
        result.surbri = self.surbri[i:j]
        result.surbri_err = self.surbri_err[i:j]
        return result

    def __getitem__(self, ivals):
        result = RadProfile()
        result.rleft  = self.rleft[ivals]
        result.rright = self.rright[ivals]
        result.surbri = self.surbri[ivals]
        result.surbri_err = self.surbri_err[ivals]
        return result

    def minus(self, value, value_err=0):
        oldsb     = self.surbri
        oldsb_err = self.surbri_err
        self.surbri = oldsb - value
        self.surbri_err = err.prop_add( oldsb_err, value_err )
        return

    def plus(self, value, value_err=0):
        oldsb     = self.surbri
        oldsb_err = self.surbri_err
        self.surbri = oldsb + value
        self.surbri_err = err.prop_add( oldsb_err, value_err )
        return

    def divide(self, value, value_err=0):
        oldsb     = self.surbri
        oldsb_err = self.surbri_err
        self.surbri = oldsb / value
        self.surbri_err = err.prop_div( oldsb, value, oldsb_err, value_err )
        return

    def multiply(self, value, value_err=0):
        oldsb     = self.surbri
        oldsb_err = self.surbri_err
        self.surbri = oldsb * value
        self.surbri_err = err.prop_mult( oldsb, value, oldsb_err, value_err )
        return

    def write(self, filename, indices='all', sci_note=False):
        if indices == 'all':
            indices = range(len(self.rmid)-1)

        FORMAT = "%f \t%f \t%f \t%f\n"
        if sci_note:
            FORMAT = "%e \t%e \t%e \t%e\n"

        f = open(filename, 'w')
        f.write( "# Bin_left\tBin_right\tSurbri\tSurbri_err\n" )
        for i in indices:
            f.write( FORMAT % \
            (self.rleft[i], self.rright[i], self.surbri[i], self.surbri_err[i]) )
        f.close()
        return

#----------------------------------------------
## Useful functions

def copy_profile(profile):
	result = RadProfile()
	result.rleft  = np.array( profile.rleft )
	result.rright = np.array( profile.rright )
	result.surbri = np.array( profile.surbri )
	result.surbri_err = np.array( profile.surbri_err )
	return result

def get_profile(filename):
    result = RadProfile()
    data   = ascii.read( filename )
    keys   = data.keys()
    result.rleft  = data[keys[0]]
    result.rright = data[keys[1]]
    result.surbri = data[keys[2]]
    result.surbri_err = data[keys[3]]
    return result


def add_profile(profile1, profile2=RadProfile(), weight1=1.0, weight2=1.0):
    result = RadProfile()
#    if profile1.rleft != profile2.rleft or profile1.rright != profile2.rright:
#        print 'Error: Profile bins need to match up'
#        return
    result.surbri = profile1.surbri * weight1 + profile2.surbri * weight2
    result.surbri_err = np.sqrt( profile1.surbri_err**2 * weight1**2 + profile2.surbri_err**2 * weight2**2 )
    result.rleft  = profile1.rleft
    result.rright = profile1.rright
    return result

def make_bkg_profile(template, bkg_value, bkg_err=0.0):
    result = RadProfile()
    result.rleft  = template.rleft
    result.rright = template.rright
    result.surbri = np.zeros( len(template.rleft) ) + bkg_value
    result.surbri_err = np.zeros( len(template.rleft) ) + bkg_err
    return result

#----------------------------------------------
## Added Feb 5, 2013 : More useful functions

def add_bkg(profile, bkg_counts, bkg_area):
    ## To subtract, put a - sign before bkg_counts
    bkg_surbri = bkg_counts / bkg_area
    bkg_err    = np.sqrt( bkg_counts ) / bkg_area
    sbnew = profile.surbri + bkg_surbri
    sbnew_err = np.sqrt( profile.surbri_err**2 + bkg_err**2 )
    profile.surbri = sbnew
    profile.surbri_err = sbnew_err
    return

def residual_profile( profile, model ):
    result = RadProfile()
    result.rleft  = np.array(profile.rleft)
    result.rright = np.array(profile.rright)
    result.surbri = profile.surbri - model.surbri
    result.surbri_err = np.sqrt( profile.surbri_err**2 + model.surbri_err**2 )
    return result
