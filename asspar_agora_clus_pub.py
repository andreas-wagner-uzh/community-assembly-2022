# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 12:50:22 2021

@author: aw

To hold parameters and to print those parameters 
by writing to an output file
"""
#to be able to write to stdout.
import sys

from datetime import datetime
dateFORMAT = '%d-%m-%Y'
foutname="test_"+ datetime.now().strftime(dateFORMAT) + "_out.txt"

#for the random number generator
seed=65377

#path to a file where a list of models to be used is stored, together with some model properties, such as biomass growth rate
modellist_inpath = "metmodels\\Agora1_03_without_mucins\\"
modellist_infname = "agora_modelsample_spec_100_08-09-2021_out.txt"

#the directory where all agora models will be stored
model_inpath="metmodels\\Agora1_03_without_mucins\\mat\\"

med_default=0  #default initial concentration of all medium components
#the diet that will define the environment in which the community exists
diet='western'
diet_infile="metmodels\\Agora1_03_without_mucins\\Magnusdottir_2017_TableS12_diets.csv"

#all nutrient concentrations below this value will be set to zero
minconc=1e-09

initbiomass=0.01 #initial biomass in grams for all invading species

#the fraction of total biomass below which a strain is considered extinct
ext_thresh=1e-03

D=0.017     #dilution rate, dimensionless
deltat=0.1  #time step of simulations, in hours
tmax=5000   #maximal simulation time for each assembly sequence, in hours
n_ass=1     #number of assembly sequences to try 
tint_eq=24  #test a community for equilibrium every tint_eq hours
twin_eq=4   #the time window into the past that is considered to test for equilibrium
#maximal deviation of max-min from mean tolerated to consider a community in equilibrium
devtol_eq=0.01


#the following parameters are to test communities for sensitivity to medium perturbations
#unused because perturbation analysis too costly
#expected number of species perturbations per hour
m_spec_pert=0
#expected number of environmental perturbations per hour
m_env_pert=0
env_pert_strength=1.5 #factor by which to perturb media and fresh medium nutrient concentrations 
spec_pert_strength=1.5 #factor by which to perturb species biomasses 
n_pert=3 #number of trial perturbation to examine stability of each community
tmax_eq=256 #maximal time allowed until a community must have reached equilibrium after a perturbation
tint_eq_pert=24  #perturbed community is examined every tint_eq_pert hours for equilibrium, used for equilibration after perturbation 

#standard (reference) Michaelis Menten constants for nutrient uptake
Vmaxref=20 #mmol / g / hr
Kmref=0.05 #mmol / l (concentration at half-maximal velocity)

#minimal biomass growth flux difference that cross-feeding must cause in the recipient
#to be considered cross-feeding
minfluxdiff=1e-09

#to print parameters to stdout and an output file
def parprint(outfilename):
    fout = open(outfilename, 'a')
    #to write to both the dedicated outfile and the console
    outarr=[fout, sys.stdout]
    for outfile in outarr:
        
        print('\nlist of models to be loaded ', modellist_infname, end ='', file = outfile)        
        print('\noutput file', foutname, end ='', file = outfile)
        print('\ndiet input file', diet_infile, end ='', file = outfile)
        print('\ndiet ', diet, end ='', file = outfile)

        print('\nseed=', seed, end ='', file = outfile)
        print('\nmed_default={0:.4f}'.format(med_default), end ='', file = outfile)
        print('\nminconc={0:1.2e}'.format(minconc), end ='', file = outfile, flush = True)

    
        print('\ninitbiomass={0:.4f}'.format(initbiomass), end ='', file = outfile)
        
        print('\next_thresh={0:.4e}'.format(ext_thresh), end ='', file = outfile)
        
        print('\nD={0:.4f}'.format(D), end ='', file = outfile)
        print('\ndeltat={0:.4f}'.format(deltat), end ='', file = outfile)
        print('\ntmax={0:.4f}'.format(tmax), end ='', file = outfile)
        print('\nn_ass={0:.4f}'.format(n_ass), end ='', file = outfile)
        print('\ntint_eq={0:.4f}'.format(tint_eq), end ='', file = outfile)
        print('\ntwin_eq={0:.4f}'.format(twin_eq), end ='', file = outfile)
        print('\ndevtol_eq={0:.4e}'.format(devtol_eq), end ='', file = outfile)
        
        print('\nm_spec_pert={0:.4f}'.format(m_spec_pert), end ='', file = outfile)
        print('\nm_env_pert={0:.4f}'.format(m_env_pert), end ='', file = outfile)

        
        print('\nenv_pert_strength={0:.4f}'.format(env_pert_strength), end ='', file = outfile)
        print('\nspec_pert_strength={0:.4f}'.format(spec_pert_strength), end ='', file = outfile)

        print('\nn_pert={0:.4f}'.format(n_pert), end ='', file = outfile)
        print('\ntmax_eq={0:.4f}'.format(tmax_eq), end ='', file = outfile)
        print('\ntint_eq_pert={0:.4f}'.format(tint_eq_pert), end ='', file = outfile)
    
        print('\nVmaxref={0:.4f}'.format(Vmaxref), end ='', file = outfile)
        print('\nKmref={0:.4f}'.format(Kmref), end ='', file = outfile, flush = True)
        
        print('\nminfluxdiff={0:1.2e}'.format(minfluxdiff), end ='', file = outfile, flush = True)

    
    #do not close outfile here or we'll get into trouble with the file handle in the 
    #calling program
    



