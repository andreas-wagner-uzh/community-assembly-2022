# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 08:49:02 2021

@author: aw
"""

"""
loads multiple agora metabolic models and assembles them into communities
through serial invasion of single species, using dFBA, records the sequence 
of invasions, extinctions, and equilibrium species communities reached, as well
as multiple community statistics.

Needs to be called as 'assembly_cstat_agora_cluscode_pub.py [jobid]', where
jobid is a positive integer that will determine the random number seed
to be used, intended to run on a cluster where multiple jobs can be 
run simultaneously.

"""


"""
SET PARAMETERS AND WRITE TO OUTPUT FILE
"""

#read the command line argument jobid
import sys
if len(sys.argv)<2:
    print('error aw: called program without a job index')
    sys.exit()
jobid=sys.argv[1]

#the parameter file and a function to output the parameters

import asspar_agora_clus_pub as par

import random
#set the seed to a number that is a function of the job id
currseed=par.seed + 100*int(jobid)
random.seed(currseed)

import numpy as np
#to write to stdout at the end of this file
import sys
import copy
import cobra

#a file with my own metabolic functions
import metfuncs_aw_pub as mfunc

#give the output file a unique job id 
outfilename='job_'+jobid+'_'+par.foutname
fout = open(outfilename, 'w')
#a bit circuituous, but otherwise we overwrite for some reason the parameter list
fout.close()
fout = open(outfilename, 'a')

#print parameter file to the same file with prefix par_,  
paroutfilename='par_'+par.foutname
par.parprint(paroutfilename)

#write output both to a file and to stdout for monitoring
outarr=[fout, sys.stdout]
#neader for global output
for outfile in outarr:
    print("jobid\tcurrseed\t", end='', file=outfile)
    print("n_ass\tn_comm\tt_curr\tn_spec\tEShannon\tESimpson\tcomm_id\t", end='', file=outfile)
    print("doms\tbiom_dom\tmu_doms_west\tmu_doms_highfib\tnexcr_doms_west\tnexcr_doms_highfib\t", end='', file=outfile)       
    print("t_since_last\tn_inv\tn_ext\tbiom\ttot_n_nut\tnutconc\tc_conc\t", end='', file=outfile)
    print("spec_pert_sens\tenv_pert_sens\tfvia_fresh\tfcomp_all\tfneut_all\tffaci_all\tfcomp_pw\tfneut_pw\tffaci_pw\tf_cfeed\t", end='', file=outfile)

                        
"""
LOAD MODELS
"""


infile=par.modellist_inpath + par.modellist_infname

#useful utils to load entire data tables
import pandas as pd


df =  pd.read_csv(infile, sep='\s+', index_col=False)

#now associate with each species the file name from which it should be loaded,
#as well as various statistics, through hashes keyed by a species index
speciesfile={}
sp_file={}
sp_reactions={}
sp_exch_r={}
sp_blocked_r={}
sp_metabolites={}
sp_mu_western={}
sp_mu_highfib={}
sp_excr_western={}
sp_excr_highfib={}

for i in df.index:
    #numerical species identifier
    sp_id=df.at[i, 'i']
    sp_file[sp_id]=df.at[i, 'species']
    sp_reactions[sp_id]=df.at[i, 'reactions']
    sp_exch_r[sp_id]=df.at[i, 'exch_r']
    sp_blocked_r[sp_id]=df.at[i, 'blocked_r']
    sp_metabolites[sp_id]=df.at[i, 'metabolites']
    sp_mu_western[sp_id]=df.at[i, 'mu_western']
    sp_mu_highfib[sp_id]=df.at[i, 'mu_highfib']
    sp_excr_western[sp_id]=df.at[i, 'excr_western']
    sp_excr_highfib[sp_id]=df.at[i, 'excr_highfib']




print('\nagora data set contains ', len(sp_file), ' distinct species from file ', infile )
print('\n loading mat files from directory ', par.model_inpath)
#a dictionary of species keyed by an integer index for the species
species={}
#an integer id that will be a shorthand for the species
species_index=0
#a dictionary (map) from the index to the file name holding the species
index_to_name={}

for i in sorted(sp_file.keys()):
    
    species[i]=cobra.io.load_matlab_model(par.model_inpath+sp_file[i])

    index_to_name[i]=sp_file[i]
    print('loading', par.model_inpath+sp_file[i])

    #close all exchange fluxes for import, to be on the safe side
    #the function for ranvia should also work for agora models
    mfunc.set_all_exchange_fluxes_ranvia_SR2018(species[i], 0, 1000)

print('loaded ', len(species), 'agora species')



"""
CREATE MEDIUM FROM THE UNION OF EXTERNAL METABOLITES
create a dictionary of metabolite concentrations keyed by all the 
metabolites that are external to at least one loaded model
"""
#load a diet as a dictionary keyed by metabolites, and giving lower bounds for imports
#of the metabolite
diet_lb=mfunc.load_diet_agora(par.diet_infile, par.diet)
print("diet:", par.diet)
print("testing dietary requirements for loaded species ")

spec_for_del=[]
for i in species.keys():
    mfunc.set_uptake_from_diet_agora(species[i], diet_lb)
    

#all calculations for the medium below implicitly assume a constant volume of 1 liter
print("\ninitializing medium for all species")

#first translate the loaded diet dict, which gives lower bounds for uptake rates
#into a concentration file. I  assume that conc=-uptake_lower_bound. 
#Rationale: if there is 1g of biomass
#per l medium, and only species i exists in the community, then the maxcimal 
#uptake rate of each metabolite for this species is given 
#by the value in the loaded diet file. Thus, the concentration of this metabolite must be at least
#as high as the minus the uptake rate

#to begin with, map exchange reactions in the diet file to external metabolites
#note that not every species may have the needed exchange reaction or metabolite
#the returned dictionary is keyed by exchange reactions and valued by external 
#metabolites, and contains key-value a pair only if at least one model has the pair
media_metabolites=mfunc.map_exchange_reactions_onto_external_metabolites(species, list(diet_lb.keys()))

diet_conc={}
for exr in media_metabolites.keys():
    diet_conc[media_metabolites[exr]]= - diet_lb[exr]

medium=mfunc.init_medium_multiple_models_agora(species, par.med_default, diet_conc)
#define what an unspent medium looks like
medfresh = mfunc.init_medium_multiple_models_agora(species, par.med_default, diet_conc)

#a backup copy of medfresh, for simulations wehere medfresh may change dover time
medfresh_init=copy.deepcopy(medfresh)


#a 2D dictionary keyed by carbon source and modelname that holds the maximum uptake rate Vmax and Km for each nutrient
#for each species
Vmax={}
Km={}
for nutrient in medium.keys():
    Vmax[nutrient]={}
    Km[nutrient]={}
    for name in species.keys():           
        Vmax[nutrient][name]=par.Vmaxref
        Km[nutrient][name]=par.Kmref

       

"""
CONNECT MEDIA COMPONENTS TO EXCHANGE REACTIONS
"""
#dictionary of dictionary of exchange reactions for each model
print("\nmapping media components onto exchange reactions")
ex_r={}
#identify exchange reactions for each species and medium component
for name in species.keys():
    ex_r[name]=mfunc.map_external_metabolites_onto_exchange_reactions(species[name], medium) 

"""
MASTER ARRAYS FOR DATA STORAGE 
"""

#array of arrays indexed by assembly number, array[i] will hold the time  
#that elapsed between the i-th and i+1st stable assembly
time_bet_ass=[]
#analogous for extinctions and invasions between assemblies
ext_bet_ass=[]
inv_bet_ass=[]

#analogus for species per assembly
spec_per_ass=[]

#array of arrays indexed by assembly number, array[i] will hold
#the perturbation resistance of the i-th assembly
env_pert_sens_per_ass=[]
spec_pert_sens_per_ass=[]

#an array indexed by assembly sequence that will hold a dictionary.
#the dictionary will be keyed by time and each entry will hold an 
#entire community object 
comm=[]
for i in range(0, par.n_ass):
    #initialize with an empty dictionary
    comm.append({})

"""
START INDIVIDUAL ASSEMBLY SEQUENCES
"""
#master loop over multiple assembly sequences
for ass_ctr in range(0, par.n_ass):
    print("\nASSEMBLY SEQUENCE ", ass_ctr)
    
    #to monitor how the medium changes over time
    #a dictionary keyed by media component whose values are arrays of nutrient concentrations
    medium_evolve={}
    #one key points to an array of time steps
    medium_evolve['t']=[]
    medium_evolve['t'].append(0)
    #the others hold an array of nutrient concentrations
    for m in medium.keys():
        medium_evolve[m]=[]
        medium_evolve[m].append(medium[m])
    
    
    
    #set up a dictionary keyed with species id (an integer) to biomass for each species, 
    #this is easy to work with but not memory efficient, because most species will not be present at most times
    biomass={} #biomass for each model
    for name in species.keys():
        #initial biomass in grams / l
        biomass[name]=0
    #chose a random species to be present from the beginning
    firstspecies=random.choice(list(species.keys()))
    
    #to store state of random number generator to ensure reproducibility,a  remnant
    #of a debugging stage of the code, but left in here and below for code stability
    rndstate=random.getstate()
    print('initializing with species ', firstspecies,  index_to_name[firstspecies])
    biomass[firstspecies]=par.initbiomass
    
    #a dictionary keyed by time that indicates when changes in a community have taken place
    #these can be species extinctions, additions of species, and the attainment of an 
    #equilibrium without extinction.   
    ass_seq={}
    
    #keyed by time, will hold an array of species that were introduced at that time 
    inv_seq={}
    inv_seq[0]=[firstspecies]
    #keyed by time, will hold an array of species that went extinct at that time 
    ext_seq={}
    #dictionary keyed by time that will hold the environmental perturbation sensitivity of the assembly at certain times
    env_pert_sens_seq={}
    #dictionary keyed by time that will hold the species (biomass) perturbation sensitivity of the assembly at certain times
    spec_pert_sens_seq={}
    
    
    biomassarr={}   #a dictionary keyed with species name that points to an array
                    #with model-specific biomass as a function of time
    for name in species.keys():           
        biomassarr[name]=[]
        biomassarr[name].append(biomass[name])
    
    
    #initial sum of biomass values for all species
    totbiomass=0
    for name in species.keys():           
        totbiomass = totbiomass+biomass[name]
    totbiomassarr=[] #to store how biomass changes over time
    totbiomassarr.append(totbiomass)
             
    
   
    """
    LOOP OVER TIME STEPS OF THE CURRENT ASSEMBLY SEQUENCE
    """
    #the number of persisting equilibrium communities that have become established in this assembly sequence
    comm_ctr=0
    for timestep in range(1, int(par.tmax/par.deltat)+1): 
        tcurr=timestep*par.deltat
        
        #apply species perturbations at random time intervals
        if random.random() < par.m_spec_pert*par.deltat:
            biomass = mfunc.ran_pert_spec(biomass, par.spec_pert_strength)
            
        #apply environmental perturbations, also to media feed, at random time intervals
        if random.random() < par.m_env_pert*par.deltat:
            medfresh=mfunc.ran_pert_env_medfresh(medfresh, medfresh_init, par.env_pert_strength)
            
                
        #now, for each nutrient and species, determine the maximal amount that the species can uptake
        for nutrient in medium.keys():
            
            #notice that the partition nutrient function will only return an uptake rate for species
            #with biomass greater than zero that have an exchange reaction for the nutrient
            uptake = mfunc.partition_nutrients(nutrient, medium[nutrient], ex_r, species, biomass, Vmax[nutrient], Km[nutrient],  par.deltat)
            #do the following only for those species names in uptake, which are those with 
            #biomass greater than zero that have an exchange reaction for the nutrient
            for name in uptake.keys(): 
                    #identify the exchange reaction for the current model and nutrient
                    exchange=species[name].reactions.get_by_id(ex_r[name][nutrient])
                    #note that uptake fluxes are negative, hence the minus sign
                    #also note that FBA works with fluxes of mmoles/g/h, hence the division
                    exchange.lower_bound= - uptake[name]/(par.deltat*biomass[name])
      
        solution={} #dictionary for FBA solution of each species
        biomass_added={} #currently added biomass of each species
        mu_est={} #will hold estimated growth rate mu for each species
        t2_est={} #estimated doubling time for each species
        totbiomass=0
        for (name, model) in species.items():
            #only do FBA if a model actually has positive biomass
            if biomass[name]>0:
                solution[name] = model.optimize()
    
        
        #now update the concentration of all media components
        #this must come BEFORE the new biomass is computed, because the 
        #fluxes are based on the old biomass values
        #note: a negative exchange flux means that a molecule leaves the medium
        for m in medium.keys():
            for (name, model) in species.items():
                #only update media components for species that are actually present            
                if biomass[name]>0:
                    #in case a specific model does not have metabolite m
                    try:
                        model.metabolites.get_by_id(m)
                        if solution[name].status == 'optimal':
                            medium[m]=medium[m] + solution[name].fluxes[ex_r[name][m]]*biomass[name]*par.deltat
                        else:
                            print('warning from assembly_agora_ranvia_v5.py: solver status is', \
                                  solution[name].status, ' for species ', name, ' leaving medium unchanged ')
                            
                    except KeyError:
                        pass #do nothing
                        
        #set all medium components that have fallen below a threshold value to zero, to avoid
        #solver problems with very small values
        mfunc.reset_medium(medium, par.minconc)
    
        #now add to the arrays showing time change in media components
        medium_evolve['t'].append(timestep*par.deltat)
        for m in medium.keys():
            medium_evolve[m].append(medium[m])
                  
        
            
        for (name, model) in species.items():
            # only update media components for species that are actually present 
            if biomass[name]==0:
                biomass_added[name]=0
            else:
                if solution[name].status == 'optimal':
                    biomass_added[name] = (solution[name].objective_value)*biomass[name]*par.deltat
                    mu_est[name]=solution[name].objective_value   
                else:
                    print('warning from assembly_cstat_agora_cluscode_pub.py: solver status is', \
                          solution[name].status, ' for species ', name, ' adding zero biomass ')
                   
                    biomass_added[name]=0
                    mu_est[name]=0
            
            
            #update all arrays regardless of whether a species is present or not
            biomass[name] = biomass[name] + biomass_added[name]
            biomassarr[name].append(biomass[name])
            totbiomass=totbiomass+biomass[name]
            #now add to the biomass array
        totbiomassarr.append(totbiomass)
        
        #check for equilibrium and extinctions every tint_eq hours
        if(tcurr>=par.tint_eq and int(tcurr/par.tint_eq)==tcurr/par.tint_eq):
            print("checking for equilibrium, current time: ", tcurr)
            
            #check whether we have reached an equilibrium
            equil=mfunc.equilibrium_cstat(species, biomassarr, par.deltat, par.twin_eq, par.devtol_eq)
            
                         
            #check whether a species has gone extinct
            extspecies=mfunc.extinction1(species, biomassarr, totbiomassarr, par.ext_thresh)
            
            if len(extspecies)>0:
                print(tcurr, ' species ', extspecies, ' went extinct ')
                ext_seq[tcurr]=extspecies
            for extname in extspecies:                
                #adjust total biomass and then set the species biomass to zero
                totbiomass=totbiomass-biomass[name]
                biomass[extname]=0
                biomassarr[extname][-1]=0
                
            #now for the special case where all species have gone extinct
            #in this case immediately introduce a new species, and revert to
            #non-equilibrium
            sctr=0
            for s in biomass.keys():
                  if biomass[s]>0: 
                      sctr+=1
            if sctr==0:
                print('warning from assembly_cstat_agora_cluscode_pub.py: all species extinct')
                random.setstate(rndstate)
                name=random.choice(sorted(species.keys()))
                rndstate=random.getstate()
                print("\n", tcurr, 'introducing species after all went extinct: species id ', name, index_to_name[name], flush = True)
                biomass[name]=par.initbiomass
                inv_seq[tcurr]=[name] 
                equil=0
              
            
            #update the assembly sequence only if an equilibrium has been reached
            #and if it has changed since last time
            if(equil==1):
                
                #in case of a change to the assembly, output various statistics for the
                #community, and introduce a new species, 
                #the case where this is the first equilibrium community needs to be treated as a special case
                currspecies=mfunc.species_present(species, biomass)
                #if this is the first equilibrium community 
                if len(ass_seq)==0:
                    lasttime=0
                #if this is NOT the first equilibrium community
                else:
                    lasttime=np.max(list(ass_seq.keys()))
                #the first condition below applies if this is the first equilbrium community
                if(len(ass_seq)==0 or ass_seq[lasttime] != currspecies):
                    comm_ctr+=1
                    ass_seq[tcurr] = currspecies
                                      
                    #store the current equilibrium community for future analysis
                    comm[ass_ctr][tcurr]=mfunc.CCMcomm(species, biomass, medium)
                    
                            
                    
                    
                    #if a stable community that is different from the previous stable community 
                    #has been established, analyze the community and output the results
                    for outfile in outarr:
                        
                        currcomm=comm[ass_ctr][tcurr]
                        
                        #specific for CLUSTER jobs
                        print("\n", jobid, "\t", end ='', file = outfile)
                        print(currseed, "\t", end ='', file = outfile)

                        print('{0:2d}\t'.format(ass_ctr), end ='', file = outfile)
                        print('{0:2d}\t'.format(comm_ctr), end ='', file = outfile)
                        print('{0:3.1f}\t'.format(tcurr), end ='', file = outfile)
                        #number of species in the current community
                        print('{0:2d}\t'.format(len(currspecies)), end ='', file = outfile)
                        
                        
                        #evenness measures, first normalized Shannon entropy EShannon (Shannon evenness)
                        n_currspec=len(currcomm.species.keys())   
                        if n_currspec>=2:
                            expH = mfunc.diversity_Hill(species, biomass, 1)
                            EShannon = np.log(expH)/np.log(n_currspec)
                            print('{0:4.3f}\t'.format(EShannon), end ='', file = outfile, flush = True)
                        else:
                            print('na\t', end ='', file = outfile, flush = True)
                        
                        #evenness measures, first normalized Simpson index
                        n_currspec=len(currcomm.species.keys())   
                        if n_currspec>=2:
                            Simpson = mfunc.diversity_Hill(species, biomass, 2)
                            ESimpson = Simpson/n_currspec
                            print('{0:4.3f}\t'.format(ESimpson), end ='', file = outfile, flush = True)
                        else:
                            print('na\t', end ='', file = outfile, flush = True)
                        
                        
                        #a community id consisting of the concatenated species ids 
                        commid='S'
                        for name in currspecies:
                            commid = commid+'.'+str(name)
                        print('{0:20s}\t'.format(commid), end ='', file = outfile)
                         
                        #id and biomass growth of the dominant species
                        if n_currspec>0:
                            biom_dom=0
                            for name in currspecies:
                                if biomass[name]>biom_dom:
                                    biom_dom=biomass[name]
                                    domspec=name
                            print('{0:5s}\t'.format(str(domspec)), end ='', file = outfile)
                            print('{0:3.2f}\t'.format(biom_dom), end ='', file = outfile)
                            print('{0:3.2f}\t'.format(sp_mu_western[domspec]), end ='', file = outfile)
                            print('{0:3.2f}\t'.format(sp_mu_highfib[domspec]), end ='', file = outfile)
                            print('{0:3.2f}\t'.format(sp_excr_western[domspec]), end ='', file = outfile)
                            print('{0:3.2f}\t'.format(sp_excr_highfib[domspec]), end ='', file = outfile)

                        else:
                            print ('na\tna\tna\tna\t', end ='', file = outfile)

                        #time since last stable community
                        print('{0:3.1f}\t'.format(tcurr-lasttime), end ='', file = outfile)
                        
                        #number of invasion attempts since last stable community
                        inv_since_last=0
                        for t in sorted(inv_seq.keys()):
                            #count the attempted invasion that happened at the last community
                            #(but not the one we are about to attemopt)
                            if t>= lasttime and t<tcurr:
                               inv_since_last+=1
                        print('{0:3d}\t'.format(inv_since_last), end ='', file = outfile)
                        
                        #number of extinctions since last stable community
                        ext_since_last=0
                        for t in sorted(ext_seq.keys()):
                            #count only extinctions that happened since the last stable 
                            #community, but including the present extinction, if any
                            if t>lasttime and t<=tcurr:
                               ext_since_last+=1
                        print('{0:3d}\t'.format(ext_since_last), end ='', file = outfile)
                        
                        #total community biomass
                        print('{0:3.2f}\t'.format(totbiomass), end ='', file = outfile)
                        
                        #total number of nutrients
                        tot_n_nut=0
                        for nut in medium.values():
                            if (nut>0):
                                tot_n_nut+=1
                        print('{0:3d}\t'.format(tot_n_nut), end ='', file = outfile)
                        
                        #total nutrient concentration
                        totnut=0
                        for nut in medium.values():
                            totnut+=nut
                        print('{0:3.2f}\t'.format(totnut), end ='', file = outfile)
                        
                        #total concentration of carbon atoms for those molecules where a formula exists
                        # the third argument is a list of metabolites to be excluded from the calculation
                        totCconc = mfunc.tot_C_conc_ranvia(species, medium, ['co2[e]'])
                        print('{0:6.2f}\t'.format(totCconc), end ='', file = outfile)
                        
                        
                        #perturbation sensitivity, calculate only for communities of at least two species
                        #computationally much too expensive to do routinely
                        if len(currcomm.species.keys())>=2:
                            env_pert_sens=mfunc.env_pert_sensitivity_medfresh_cstat(species, biomass, medium, medfresh, medfresh_init, par.D, ex_r, Vmax, Km, \
                                           par.deltat, par.tint_eq_pert, par.tmax_eq, par.devtol_eq, \
                                           par.ext_thresh, par.env_pert_strength, par.n_pert, par.minconc)               
                            env_pert_sens_seq[tcurr]=env_pert_sens
                            
                            spec_pert_sens=mfunc.spec_pert_sensitivity_cstat(species, biomass, medium, medfresh, \
                                                               par.D, ex_r, Vmax, Km, par.deltat, \
                                                               par.tint_eq_pert, par.tmax_eq, par.devtol_eq, \
                                                               par.ext_thresh, par.spec_pert_strength, par.n_pert, par.minconc)                
                            spec_pert_sens_seq[tcurr]=spec_pert_sens
                            
                            print('{0:3.2f}\t'.format(spec_pert_sens), end ='', file = outfile)
                            print('{0:3.2f}\t'.format(env_pert_sens), end ='', file = outfile, flush = True)
                        else:
                            print('na\tna\t', end ='', file = outfile, flush = True)      
                    
                        
                        #the fraction of species that are viable by themselves on the fresh medium feed
                        fviable_medfresh=mfunc.frac_viable_on_medfresh(currcomm, medfresh, par.D, ex_r, Vmax, Km, par.deltat)
                        print('{0:3.2f}\t'.format(fviable_medfresh), end ='', file = outfile, flush = True)
                        
                        #the fraction of species that experience the rest of the community as competitive, neutral, or facilitative
                        if len(currcomm.species.keys())>=2:
                            [f_comp, f_neut, f_faci]=mfunc.frac_comp_neut_fac_all_postdil(currcomm, medfresh, par.D, ex_r, Vmax, Km, par.deltat)
                            print('{0:3.2f}\t'.format(f_comp), end ='', file = outfile, flush = True)
                            print('{0:3.2f}\t'.format(f_neut), end ='', file = outfile, flush = True)
                            print('{0:3.2f}\t'.format(f_faci), end ='', file = outfile, flush = True)
                        else:
                            print('na\tna\tna\t', end ='', file = outfile, flush = True)

                        #the fraction of species that experience the rest of the community as competitive, neutral, or facilitative
                        if len(currcomm.species.keys())>=2:
                            [f_comp, f_neut, f_faci]=mfunc.frac_comp_neut_fac_pairwise_postdil(currcomm, medfresh, par.D, ex_r, Vmax, Km, par.deltat)
                            print('{0:3.2f}\t'.format(f_comp), end ='', file = outfile, flush = True)
                            print('{0:3.2f}\t'.format(f_neut), end ='', file = outfile, flush = True)
                            print('{0:3.2f}\t'.format(f_faci), end ='', file = outfile, flush = True)
                        else:
                            print('na\tna\tna\t', end ='', file = outfile, flush = True)
                            
                        #the fraction of species pairs that show cross-feeding on at least one metabolite
                        if len(currcomm.species.keys())>=2:
                            #the minimal difference in biomass growth flux that cross-feeding must cause
                            f_cf=mfunc.frac_crossfeed(currcomm, par.D, ex_r, Vmax, Km, par.deltat, par.minfluxdiff)
                           
                            print('{0:3.2f}\t'.format(f_cf), end ='', file = outfile, flush = True)
                        else:
                            print('na\t', end ='', file = outfile, flush = True)
                            
                         
                         
    
                            

                #if an equilibrium has been reached, regardless of whether
                #the community that has been established is new,
                #introduce a new species that is not already present
                #if all species have previously gone extinct, this will introduce a new species
                #setstate is to ensure reproducibility of random number generation 
                random.setstate(rndstate)
                name=random.choice(sorted(species.keys()))
                rndstate=random.getstate()

                while(biomass[name]>0):
                     random.setstate(rndstate)
                     name=random.choice(sorted(species.keys()))
                     rndstate=random.getstate()

                print("\n", tcurr, 'introducing species', name, index_to_name[name], flush = True)
                biomass[name]=par.initbiomass
                inv_seq[tcurr]=[name] 
                
        
        # now apply a dilution flux, i.e., replace a fraction D of the current 
        # culture with a fresh medium without biomass    
        # replace a fraction D*deltat of the culture, which means that a fraction
        # 1-D*deltat is preserved
        [biomass, post_medium] = mfunc.dilute_multiple_models(biomass, medium, medfresh, 1-(par.D*par.deltat))
        for m in medium.keys():
            medium[m]=post_medium[m]
       
   

 