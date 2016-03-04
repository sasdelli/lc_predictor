#!/usr/bin/python

from lc_predictor import projection
from lc_predictor import plot_tools

# This uses as input an example directory in the module.
input_directory = 'input_data/spectra_sn2008Z/'
#input_directory = '../input_data/spectra_sn2011fe/'
#input_directory = 'input_data/2016zd/'
# This is the output directory for the pdf images.
output_directory = './out_dir/'

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
    sn_coeff = projection.calculate_empca_projections(derivative, weight)
    # It calculates the PLS prediction of the B mag LC.
    epoch, B_lc = projection.calculate_pls_curve(
        sn_coeff, pls_band='Bmag')
    # This is the PLS prediction of the B-V color curve
    epoch, color = projection.calculate_pls_curve(
        sn_coeff, pls_band='BmV')
    # It plots the scatter plot of the SN and of the training set
    plot_tools.spectra(sn_coeff, input_data=input_param,
                       out_dir=output_directory, one_plot=True)
    # This plots the light curve and color curve predicted from the
    # spectra.
    plot_tools.LC(epoch, B_lc, out_dir=output_directory)
    plot_tools.color(epoch, color, out_dir=output_directory)

    # plot the PCA coefficients of the training set
    training_sne, training_coeff = plot_tools.load_training_coeff()
    plot_tools.scatter_plot(training_sne, training_coeff,
                            sn_name=input_param['sn_name'][0],
                            sn_coeff=sn_coeff, out_dir=output_directory)
if __name__ == '__main__':
    main()
