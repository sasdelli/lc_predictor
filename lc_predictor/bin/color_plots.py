#!/usr/bin/env python

from lc_predictor import projection
#from lc_predictor import plot_tools
#from lc_predictor.plot_tools import load_training_coeff

import numpy as np
import re
from numpy import isnan
import matplotlib.pyplot as plt
import pickle

from mpl_toolkits.mplot3d import Axes3D

input_directory = 'input_data/spectra_sn2008Z/'


code_dir = '/home/wfh/PCA_codes/lc_predictor/lc_predictor/'
PLS_code_dir = '/home/wfh/PCA_codes/PCA/'
#training_sne, training_coeff = load_training_coeff(pkl_file=code_dir+'trained_data/training_coeff_150919.pkl')

pkl_file=code_dir+'trained_data/training_coeff_150919.pkl'
filehandler = open(pkl_file, 'r')
training_sne, training_coeff = pickle.load(filehandler)
filehandler.close()

# It loads the input from the input_directory variable. 
input_param = projection.read_input(input=input_directory)
# It creates an input vector.
derivative, weight, wave = projection.calculate_derivative(
    input_data=input_param)
# It calculates the PCA projections.
sn_coeff = projection.calculate_empca_projections(derivative, weight, empca_dict_file=code_dir+'trained_data/trained_empca_150919.dict')


def create_class_vec(new_name):
    cfa_dir=PLS_code_dir+"data/cfaspec_snIa/"
    SNe_data=np.loadtxt(cfa_dir+'/cfasnIa_param_mod.dat', dtype={'names': ('SN_name', 'zhel', 'tMaxB', 'err_tMaxB', 'ref', 'Dm15', 'err_Dm15', 'ref2', 'M_B', 'err_M_B', "BmV", "err_BmV", "BmmVm", "err_BmmVm", "Phot_ref"),'formats': ('S15', "f8", "f8","f8", "S15", "f8", "f8","S15","f8" , "f8","f8", "f8","f8", "f8","S15")})
    spectra_data=np.loadtxt(cfa_dir+'/cfasnIa_mjdspec.dat', dtype={'names': ('spectrum_name', 'time'),'formats': ('S40', "f8")})
    SNe_BranchWang_class=np.loadtxt(cfa_dir+'/branchwangclass_mod.dat', dtype={'names': ('SN_name', 'pEW5972', 'pEW6355', 'vabs6355', 'phase', 'Branch', 'Wang'),'formats': ('S15', "f8", "f8","f8",  "f8","S15","S15")})
    name_regex = re.compile('(.+)\-\d+\.\d+')
    name_vector=[]
    for spectrum_name in enumerate(spectra_data["spectrum_name"]):
        name_vector.append(name_regex.search(spectrum_name[1]).group(1))
    #It creates the vectors of the classification of Branch and Wang
    #SN_name_vec=[]
    pEW5972_vec=[]
    pEW6355_vec=[]
    vabs6355_vec=[]
    Branch_vec=[]
    Wang_vec=[]
    for i, supernova in enumerate(new_name):
        pEW5972_tmp=np.nan
        pEW6355_tmp=np.nan
        vabs6355_tmp=np.nan
        Branch_tmp=np.nan
        Wang_tmp=np.nan
        for name_sn in enumerate(SNe_BranchWang_class["SN_name"]):
            if name_sn[1] ==  supernova:
                SN_name_tmp, pEW5972_tmp, pEW6355_tmp, vabs6355_tmp, phase_tmp, Branch_tmp, Wang_tmp= SNe_BranchWang_class[name_sn[0]]
        #SN_name_vec.append(SN_name_tmp)
        pEW5972_vec.append(pEW5972_tmp)
        pEW6355_vec.append(pEW6355_tmp)
        vabs6355_vec.append(vabs6355_tmp)
        Branch_vec.append(Branch_tmp)
        Wang_vec.append(Wang_tmp)
    #color plot for Branch 
    color_plot_Branch=[]
    for i in  range(0,np.size(new_name)):
        if Branch_vec[i]=="CN":
            color_plot_Branch.append('r')
        elif  Branch_vec[i]=="SS":
            color_plot_Branch.append('g')
        elif  Branch_vec[i]=="BL":
            color_plot_Branch.append('b')
        elif  Branch_vec[i]=="CL":
            color_plot_Branch.append('y')
        else:
            color_plot_Branch.append('w')
    #color plot for Wang 
    color_plot_Wang=[]
    for i in  range(0,np.size(new_name)):
        if Wang_vec[i]=="91T":
            color_plot_Wang.append('r')
        elif  Wang_vec[i]=="N" :
            color_plot_Wang.append('g')
        elif  Wang_vec[i]=="pec":
            color_plot_Wang.append('b')
        elif  Wang_vec[i]=="HV":
            color_plot_Wang.append('y')
        elif  Wang_vec[i]=="91bg":
            color_plot_Wang.append('c')
        else:
            color_plot_Wang.append('w')
    return color_plot_Wang, color_plot_Branch
color_plot_Wang, color_plot_Branch = create_class_vec(training_sne)

Wang_keys = [
    ["ro", "$91T$"],
    ["go","$Norm$"],
    ["bo","$Pec$"],
    ["yo","$HV$"],
    ["co","$91bg$"],
             ] 
Branch_keys = [
    ["ro", "$CN$"],
    ["go","$SS$"],
    ["bo","$BL$"],
    ["yo","$CL$"],
             ] 

def plot_3D(m_coeff, color_plot, keys=[], components=[0,1,2], SN=None):
    """
    Make a 3D plot of components and an additional SN. The default is the first three PCs.
    """

    #plot 3d components
    #fig = plt.figure(figsize=(16, 10))
    fig= plt.figure()
    ax = Axes3D(fig)
    if SN != None:
        ax.scatter(SN[components[0]], SN[components[1]], SN[components[2]], c='k', marker='D', s=30)

    ax.scatter(m_coeff[:,components[0]], m_coeff[:,components[1]], m_coeff[:,components[2]], c=color_plot, marker='o')
    #ax.scatter([0],[0],[0],marker='x')
    ax.set_xlabel('component %d' % (components[0]+1))
    ax.set_ylabel('component %d' % (components[1]+1))
    ax.set_zlabel('component %d' % (components[2]+1))
    for i in keys:
        ax.plot([],[],i[0],label=i[1])
    ax.legend(loc= 'lower right', bbox_to_anchor = (1, 0))


# This is an example of how to use the plotting script.
plot_3D(training_coeff, color_plot_Wang, keys = Wang_keys, SN=sn_coeff)
plot_3D(training_coeff, color_plot_Wang, keys = Wang_keys, SN=sn_coeff, components=[0,1,3])
plot_3D(training_coeff, color_plot_Wang, keys = Wang_keys, SN=sn_coeff, components=[0,2,3])
plot_3D(training_coeff, color_plot_Wang, keys = Wang_keys, SN=sn_coeff, components=[1,2,3])

plt.show()
    
    
    
