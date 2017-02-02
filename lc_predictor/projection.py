import numpy as np
from empca import empca
import pickle
from scipy.stats import nanmean
from sklearn.pls import PLSRegression
import scipy.interpolate
import savgol
import pkgutil
import os

def calculate_empca_projections(derivative, weight,empca_dict_file=None):
    if empca_dict_file == None:
        empca_dict = pickle.loads(pkgutil.get_data(
            'lc_predictor', 'trained_data/trained_empca.dict'))
    else:
        filehandler = open(empca_dict_file, 'r')
        empca_dict = pickle.load(filehandler)
        filehandler.close()
    mean_derivative = empca_dict['mean_derivative']
    eigvec = empca_dict['eigvec']
    empca_model = empca(np.zeros([1,np.size(mean_derivative)]), niter=0)
    empca_model.eigvec = eigvec
    centered_derivative = np.nan_to_num(derivative - mean_derivative)
    empca_model.set_data(np.array([centered_derivative]), np.array([weight]))
    return empca_model.coeff[0]
    
def calculate_pls_curve(coeff, pls_band='Bmag', pls_dict_file=None): 

    def dict2plsca(dict):
        plsca = PLSRegression(n_components=np.shape(dict['coefs'])[0])
        plsca.x_mean_ = dict['x_mean']
        plsca.y_mean_ = dict['y_mean']
        plsca.coefs = dict['coefs'] 
        return plsca
   
    def load_pls(bands='Bmag', pls_dict_file=None):
        if pls_dict_file == None:
            dict = pickle.loads(pkgutil.get_data(
                'lc_predictor', 'trained_data/pls_%s.dict' % bands))
        else:
            filehandler3 = open(pls_dict_file, 'r')
            dict = pickle.load(filehandler3)
            filehandler3.close()
        plsca = dict2plsca(dict)
        epochs = [None]
        if 'epochs' in dict:
            epochs = dict['epochs'] 
        return dict, epochs
    
    def dict2mean(X, dict):
        plsca = PLSRegression(n_components=np.shape(dict['coefs'])[0])
        plsca.x_mean_ = dict['x_mean']
        plsca.y_mean_ = dict['y_mean']
        plsca.coefs = dict['coefs'] 
        return plsca.predict(X)

    dict_pls, epochs = load_pls(bands=pls_band, pls_dict_file=pls_dict_file) 
    phot = dict2mean(np.array([coeff]), dict_pls)
    return epochs, phot[0]

def read_input(input='input_data/spectra_sn2008Z/',
               SN_METADATA='sn_metadata.dat', SPECTRA_METADATA='metadata.dat'):
    if not input[0] in ['/', '~', '.']:
        input = pkgutil.get_loader('lc_predictor').filename + '/' + input
    dirin = False
    if os.path.isdir(input):
        dirin = True
        in_dir = input
        input = input + SN_METADATA
    label_attr_map = {
    'Helio_Redshift' : ['zhelio', float],
    'Date_Bmax' : ['day_Bmax', float],
    'Deredshift_spectra' : ['DEREDSHIFT', int],
    'Spectra_Metadata' : ['spectra_meta', str],
    'Spectra_Dir' : ['spectra_dir', str],
    'SN_Name' : ['sn_name', str]
    }
    in_param = {}

    if not input[0] in ['/', '~', '.']:
        in_file = (pkgutil.get_data('lc_predictor', input)).splitlines()
    else:
        in_file = open(input, 'rb').read().splitlines()
    for line in in_file:
        row = line.split()
        if row != []:
            if row[0] != '#':
                label = row[0]
                data = row[1:]
                attr = label_attr_map[label][0]
                datatypes = label_attr_map[label][1:]
                values = [
                    (datatypes[i](data[i])) for i in range(len(data)) ]
                in_param[attr] = values
    if dirin:
        in_param['spectra_dir'] = [in_dir]
    if not in_param['spectra_dir'][0][0] in ['/', '~', '.']:
        in_param['spectra_dir'][0] = (
            pkgutil.get_loader('lc_predictor').filename
            + '/' + in_param['spectra_dir'][0])
    if not 'spectra_meta' in in_param.keys():
        in_param['spectra_meta'] = [
            in_param['spectra_dir'][0] + SPECTRA_METADATA ]    
    if not in_param['spectra_meta'][0][0] in ['/', '~', '.']:
        in_param['spectra_meta'][0] = (
            pkgutil.get_loader('lc_predictor').filename
            + '/' + in_param['spectra_meta'][0])
    in_param['spectra_list'] = np.atleast_1d(np.genfromtxt(
        in_param['spectra_meta'][0], dtype=None))
    return in_param

def calculate_derivative(input_data, 
                         empca_dict_file=None,
                         redd_corr=None):
    if empca_dict_file == None:
        empca_dict = pickle.loads(pkgutil.get_data('lc_predictor',
                                  'trained_data/trained_empca.dict'))
    else:
        filehandler = open(empca_dict_file, 'r')
        empca_dict = pickle.load(filehandler)
        filehandler.close()
    spectra_dir = input_data['spectra_dir'][0]
    day_Bmax = input_data['day_Bmax'][0]
    DEREDSHIFT = input_data['DEREDSHIFT'][0]
    zhelio = input_data['zhelio'][0]
    spectra_list = input_data['spectra_list']
    bins = empca_dict['bins']
    wave_range = empca_dict['wave_range']
    vel_savgol = empca_dict['sg_param'] 
    savgol_par = [int(0.00408*vel_savgol[0]*np.sqrt(vel_savgol[2]))*2-1,
                  vel_savgol[1]]

    def mean_sp(wave, log, weight):
        wave_out = wave[0]
        if np.any(np.mean(wave,0) != wave[0]):
            print 'error incoherent wavelengths'
        spectrum = np.sum(
            np.nan_to_num(np.array(log))*np.nan_to_num(np.array(weight)),
            0)/np.nan_to_num(np.sum(weight,0)) 
        w_out = (np.sum(np.nan_to_num(np.array(weight)), 0)) 
        return wave_out, spectrum, w_out

    def resample_sp(in_list):
        wave, sp, w = in_list
        wave_in = wave[0]
        new_index_0 = (wave_range[0] - wave_in)/2.
        new_index_1 = (wave_range[1] - wave_in)/2.
        new_wave = range(wave_range[0],wave_range[1], wave_range[2])
        new_sp = np.zeros((wave_range[1]-wave_range[0])/wave_range[2])
        new_w = np.zeros((wave_range[1]-wave_range[0])/wave_range[2])
        for i in range(-wave_range[2]/4+1, wave_range[2]/4+1):
            # The mean can be done better. 
            new_sp += sp[new_index_0 +i:new_index_1 +i:wave_range[2]/2]
            new_w += 1./ w[new_index_0 +i:new_index_1 +i:wave_range[2]/2]
        new_sp = new_sp/((wave_range[2]/2))
        # propagated weight on the mean 
        new_w = 1./(new_w/((wave_range[2]/2)**2)) 
        return new_wave, new_sp, new_w

    def calculate_res(spectrum_wave, spectrum_values):
        def rebin_(x, y, wave_out):
            f = scipy.interpolate.interp1d(x, y)
            min_i = min([
                i for i in range(np.size(wave_out)) if wave_out[i] > x[0] ])
            max_i = max([
                i for i in range(np.size(wave_out)) if wave_out[i] < x[-1] ])
            ynew = np.zeros(np.size(wave_out))
            if max_i > min_i:
                ynew[min_i:max_i] = f(wave_out[min_i:max_i])
            return np.array(ynew)

        w_spectrum_values = np.zeros(np.size(spectrum_values))
        w_spectrum_values[:] = np.nan
        i_3800 = 0
        i_8000 = -1
        for i_ in range(np.size(spectrum_wave)):
            if spectrum_wave[i_] < 3800.:
                i_3800 = i_
            if spectrum_wave[i_] < 8000.:
                i_8000 = i_
        if i_3800 != i_8000:
            max_value = max(np.nan_to_num(spectrum_values[i_3800:i_8000]))
        else:
            max_value = max(np.nan_to_num(spectrum_values))
            print 'do not use spectra out of the wavelength range 3000-8000'
        spectrum_values = (spectrum_values)/max_value
        # This needs to be simplified. 
        b = 1000
        c = 15000
        n_of_points = 40
        factor_resid = 1.e-4
        wavelength_array = np.array(range(b,c,2))
        spectrum_values_tmp = rebin_(spectrum_wave, spectrum_values,
                                    wavelength_array)
        print 'check equal ', np.size(
            spectrum_values_tmp), np.size(wavelength_array)
        w_spectrum_values_tmp = np.zeros(np.size(spectrum_values_tmp))
        resid_vec = [] 
        for j in range(np.size(wavelength_array)/n_of_points):
            wave_ = wavelength_array[j*n_of_points:(j+1)*n_of_points]
            flux_ = spectrum_values_tmp[j*n_of_points:(j+1)*n_of_points]
            pol_ = np.polyfit(wave_,flux_,4,full=True)
            resid = pol_[1][0]
            pol = pol_[0]
            resid_vec.append(resid)
        mean_resid = nanmean(resid_vec)
        found = 0
        for j in range(np.size(wavelength_array)/n_of_points):
            if resid_vec[j] > factor_resid*mean_resid:
                w_spectrum_values_tmp[
                    j*n_of_points:(j+1)*n_of_points] = 1./np.sqrt(resid_vec[j])
        return wavelength_array, spectrum_values_tmp, w_spectrum_values_tmp

    def load_spectra(range_, spectra_list=None, Bmax=None, der_s=False, 
                     zhelio=None, spectra_dir=None):
        list_spectra = [ spectrum_name for spectrum_name,
                        dayspectrum in spectra_list if range_[0] <=
                        (float(dayspectrum)- Bmax)/(1.+zhelio) and 
                        (float(dayspectrum)- Bmax)/(1.+zhelio) < range_[1] ]
        print ''
        print 'Computing epochs: ', range_
        print 'Spectra included: ', list_spectra
        wave = []
        value = []
        weight = []
        for ii, spectrum in enumerate(list_spectra):
            sp = np.genfromtxt(spectra_dir+spectrum).T
            _value_ = sp[1]
            if der_s:
                _wave_ = sp[0]/(1.+zhelio)
            else:
                _wave_ = sp[0]
            if redd_corr != 0.:
                _value_ = (redd_law(_wave_,_value_, redd=redd_corr))
            _wave_,_value_,_weight_ = calculate_res(_wave_, _value_)
            wave.append(_wave_)
            value.append(_value_)
            weight.append(_weight_)
        return wave, value, weight

    def redd_law(wave, spectrum, redd=0.):
        out_spectrum = np.array(spectrum)*10**(
            0.4*(redd*(3.1+2.002*((10000./np.array(wave))-(1.0/0.55)))))
        return out_spectrum

    def clean_neg_val(vec, vec2):
        for i in range(np.size(vec)):
            if vec[i] < 0.:
                vec[i] = 0.
                vec2[i] = 0.
        return vec, vec2

    if redd_corr == None:
        redd_corr = 0
    size_time_bin = ((wave_range[1]-wave_range[0])/wave_range[2])
    row = np.zeros(size_time_bin*np.size(bins[:-1]))
    w_row = np.zeros(size_time_bin*np.size(bins[:-1]))
    wave_row = np.zeros(size_time_bin*np.size(bins[:-1]))
    row[:] = np.NaN
    wave_row[:] = np.NaN
    for j, day in enumerate(bins[:-1]):
        binning = bins[j+1]-day
        spectrum_wave, spectrum_values, w_spectrum_values = load_spectra(
            [day,day+binning], spectra_list=spectra_list, Bmax=day_Bmax,
            zhelio=zhelio, der_s=DEREDSHIFT, spectra_dir=spectra_dir)
        if not spectrum_wave == []:
            spectrum_der = []
            spectrum_values_smooth = []
            for k, sp in enumerate(spectrum_values):
                sp1, w_spectrum_values[k] = clean_neg_val(sp,
                                                          w_spectrum_values[k])
                spectrum_der.append(savgol.savitzky_golay(sp1,savgol_par[0], 
                                    savgol_par[1],deriv=1))
                spectrum_values_smooth.append(savgol.savitzky_golay(
                    sp,savgol_par[0],savgol_par[1],deriv=0))
            wave_res, der_res, w_res = resample_sp(mean_sp(spectrum_wave, 
                spectrum_der, w_spectrum_values))
            wave_res, values_res, w_res = resample_sp(mean_sp(spectrum_wave, 
                spectrum_values_smooth, w_spectrum_values))
            w_res = np.nan_to_num((w_res)*(((2.3026)*values_res)**2))
            values_res = np.log10(np.e)*der_res/values_res
            wave_row[j*size_time_bin:(j+1)*size_time_bin] = wave_res
            row[j*size_time_bin:(j+1)*size_time_bin] = values_res 
            w_row[j*size_time_bin:(j+1)*size_time_bin] = w_res
    return row, w_row, wave_row


