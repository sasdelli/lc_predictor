#!/usr/bin/python

from lc_predictor import projection
from lc_predictor import plot_tools

# This uses as input an example directory in the module.
input_directory = 'input_data/spectra_sn2008Z/'
#input_directory = '../input_data/spectra_sn2011fe/'
#input_directory = 'input_data/2016zd/'
# This is the output directory for the pdf images.
output_directory = './out_dir/'

code_dir = '/home/wfh/PCA_codes/lc_predictor/lc_predictor/'

def main():
    """
    This example code shows the capability of these libraries on few
    example SNe. 
    With the input_file option is possible to load a SN input file
    different from the default 'input_data/sn_metadata.dat'. E.g:

    input_param = projection.read_input(
        input='/PATH/TO/INPUT/DATA/')

    """
    # It loads the input from the input_directory variable. 
    input_param = projection.read_input(input=input_directory)
    # It creates an input vector.
    derivative, weight, wave = projection.calculate_derivative(
        input_data=input_param)
    # It calculates the PCA projections.
    sn_coeff = projection.calculate_empca_projections(derivative, weight, empca_dict_file=code_dir+'trained_data/trained_empca_150919.dict')
    # It calculates the PLS prediction of the B mag LC.
    epoch, B_lc = projection.calculate_pls_curve(
        sn_coeff, pls_band='Bmag', pls_dict_file=code_dir+'trained_data/pls_Bmag_150919.dict')
    # This is the PLS prediction of the B-V color curve
    epoch, color = projection.calculate_pls_curve(
        sn_coeff, pls_band='BmV', pls_dict_file=code_dir+'trained_data/pls_BmV_150919.dict')
    # It plots the scatter plot of the SN and of the training set
    plot_tools.spectra(sn_coeff, input_data=input_param,
                       out_dir=output_directory, one_plot=True)
    # # This plots the light curve and color curve predicted from the
    # # spectra.
    # plot_tools.LC(epoch, B_lc, out_dir=output_directory)
    # plot_tools.color(epoch, color, out_dir=output_directory)
    import numpy as np
    print np.shape(epoch)
    print np.shape(B_lc)

    # plot the PCA coefficients of the training set
    training_sne, training_coeff = plot_tools.load_training_coeff(pkl_file=code_dir+'trained_data/training_coeff_150919.pkl')
    plot_tools.scatter_plot(training_sne, training_coeff,
                            sn_name=input_param['sn_name'][0],
                            sn_coeff=sn_coeff, out_dir=output_directory)
if __name__ == '__main__':
    main()
