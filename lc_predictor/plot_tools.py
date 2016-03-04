import matplotlib
matplotlib.use('Agg')
import numpy as np
import pickle
from scipy.stats import nanmean
import matplotlib.pyplot as plt
import os
import pkgutil

def create_dir(path):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    return

def spectra(
        sn_coeff, input_data,
        empca_dict_file=None, out_dir=None,
        one_plot=False):
    if out_dir!= None:
        create_dir(out_dir)
    if empca_dict_file == None:
        empca_dict = pickle.loads(pkgutil.get_data(
            'lc_predictor', 'trained_data/trained_empca.dict'))
    else:
        filehandler = open(empca_dict_file, 'r')
        empca_dict = pickle.load(filehandler)
        filehandler.close()
    mean_derivative = empca_dict['mean_derivative']
    eigvec = empca_dict['eigvec']
    sn_name = input_data['sn_name'][0]
    reconstruct_deriv = np.sum(eigvec.T*sn_coeff,1)+mean_derivative
    spectra_dir = input_data['spectra_dir'][0]
    spectra_list = input_data['spectra_list']
    Bmax = input_data['day_Bmax'][0]
    der_s = input_data['DEREDSHIFT'][0]
    zhelio = input_data['zhelio'][0]
    bins = empca_dict['bins']
    wave_range = empca_dict['wave_range']

    def integrate_spectrum(vec):
        out_vec = [0]
        for val in vec:
            out_vec.append(out_vec[-1]-val*wave_range[2]/2.)
        return np.array(out_vec)

    if one_plot:
        plt.figure(figsize=(8,1.2*np.size(spectra_list)))
    sort_index = np.argsort([ spectrum[1] for spectrum in spectra_list])
    for n_spectrum, spectrum in enumerate(spectra_list[sort_index]):
        print spectrum
        dayspectrum = float(spectrum[1])
        spectrum = spectrum[0]
        for j in range(np.size(bins[:-1])):
            if bins[j] <= (dayspectrum- Bmax)/(1.+zhelio) and (
                    dayspectrum- Bmax)/(1.+zhelio) < bins[j+1]:
                sp = np.genfromtxt(spectra_dir+spectrum).T
                _value_ = sp[1]
                if der_s:
                    _wave_ = sp[0]/(1.+zhelio)
                else:
                    _wave_ = sp[0]
                range_1 = np.arange(
                    wave_range[0]-wave_range[2]/2., wave_range[1] +
                    wave_range[2]/2.,wave_range[2])
                min_common = max(min(_wave_), min(range_1))
                max_common = min(max(_wave_), max(range_1))
                if min_common < max_common:
                    log_const = nanmean([
                        np.log10(_value_)[i] for i in range(np.size(_value_))
                        if min_common < _wave_[i] and _wave_[i] < max_common ])
                    if one_plot:
                        i_0 = min([ i_w for i_w in range(np.size(_wave_))
                                   if _wave_[i_w] > 3000.])
                        i_1 = max([ i_w for i_w in range(np.size(_wave_)) 
                                   if _wave_[i_w] < 10000.])
                        plt.plot(_wave_[i_0:i_1], np.log10(
                                _value_[i_0:i_1])-log_const-0.6*n_spectrum,'g')
                    else:
                        plt.figure()
                        plt.plot(_wave_, _value_, label='%s' % (spectrum))
                    range_0 = np.size(reconstruct_deriv)/(np.size(bins)-1)
                    log_reconstr = integrate_spectrum(
                        reconstruct_deriv[range_0*j:range_0*(j+1)])
                    log_reconstr -= np.mean([
                        log_reconstr[i] for i in range(np.size(range_1))
                        if min_common < range_1[i] and range_1[i] < max_common
                        ])
                    log_reconstr += log_const 
                    if one_plot:
                        plt.plot(range_1,
                                 log_reconstr-log_const - 0.6*n_spectrum, 'k')
                        plt.text(
                            range_1[-1],
                            log_reconstr[-1]-log_const + 0.1 - 0.6*n_spectrum,
                            '%.2f' % ((dayspectrum- Bmax)/(1.+zhelio)))
                    else:
                        plt.plot(range_1, 10**log_reconstr,
                                 label='PCA reconstr.')
                        plt.legend()
                        plt.title(
                            'epoch = %.2f' % ((dayspectrum-Bmax)/(1.+zhelio)))
                        
                        plt.xlabel('wavelenght (\AA)')
                        plt.ylabel('flux')
                        plt.xlim((3000,8000))
                    #plt.ylim((0,1.2*max(10**log_reconstr)))
                    if out_dir != None and not one_plot:
                        plt.savefig(out_dir+'epoch_%.2fdays.pdf' %
                                    ((dayspectrum-Bmax)/(1.+zhelio)))
    if one_plot:
        plt.plot([],[],'g',label='%s' % sn_name)
        plt.plot([],[],'k',label='PCA reconstr.')
        plt.legend(loc=3)
        plt.xlabel('wavelenght (\AA)')
        plt.ylabel('log(flux) + const')
        plt.xlim((3000,8000))
        if out_dir != None:
            plt.savefig(out_dir+'spectra.pdf')
    return 

def LC(epoch, B_lc, out_dir=None):
    if out_dir != None:
        create_dir(out_dir)
    # plot LC
    plt.figure()
    plt.plot(epoch, B_lc)
    plt.xlabel('epoch since $B$ max (day)')
    plt.ylabel('predicted $B$mag (H_0=70)')
    plt.gca().invert_yaxis()
    if out_dir != None:
        plt.savefig(out_dir+'predicted_Bmag.pdf')
    return

def color(epoch, color_curve, out_dir=None):
    if out_dir != None:
        create_dir(out_dir)
    plt.figure()
    plt.plot(epoch, color_curve)
    plt.xlabel('epoch since $B$ max (day)')
    plt.ylabel('predicted $B-V$ color')
    if out_dir != None:
        plt.savefig(out_dir+'predicted_BmV.pdf')
    return

def load_training_coeff(pkl_file='trained_data/training_coeff.pkl'):
    if not pkl_file[0] in ['/', '~', '.']:
        pkl_file = pkgutil.get_loader(
            'lc_predictor').filename + '/' + pkl_file
    # load training sample coeff
    filehandler = open(pkl_file, 'r')
    training_sne, training_coeff = pickle.load(filehandler)
    filehandler.close()
    return training_sne, training_coeff

def scatter_plot(names, training, sn_name=None, sn_coeff=None, plot_text=False,
                 PC_x=1, PC_y=2, out_dir=None):
    if out_dir != None:
        create_dir(out_dir)
    plt.figure()
    plt.plot(training.T[PC_x-1], training.T[PC_y-1],'o', label='training set')
    if plot_text:
        [plt.text(training.T[PC_x-1][i], training.T[PC_y-1][i])
         for i in range(np.size(names))]
    plt.plot(sn_coeff[PC_x-1], sn_coeff[PC_y-1], 'Dr', label=sn_name)
    plt.legend()
    plt.xlabel('PC %d' % (PC_x))
    plt.ylabel('PC %d' % (PC_y))
    if out_dir != None:
        plt.savefig(out_dir+'PC%d_PC%d_scatter_plot.pdf' % (PC_x, PC_y))
    return
