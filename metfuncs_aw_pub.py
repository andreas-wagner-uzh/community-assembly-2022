# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 16:36:05 2021

@author: aw
"""

from cobra import Model, Reaction, Metabolite

import random

import numpy as np

#for regular expression searches to map metabolite KEGG ids onto metabolite names
import re

import copy

import pandas as pd


#compute hamming distance of two vectors
def hamming(a, b):
    return len([i for i in filter(lambda x: x[0] != x[1], zip(a, b))])
            
# takes a list of file names containing agora species names
# returns two lists of file names, one that holds only a single random representative
# of each species 
# another that holds only a single random representative of each genus
def random_species_agora(filelist):
    one_per_species=[]
    #list of species files such that each agora genus occurs only once
    one_per_genus=[]
    #shuffle a copy of the species file list (leave the passed original unchanged)
    ranlist=[]
    for x in filelist:
        ranlist.append(x)
    #make sure that we chose a random genus or species
    random.shuffle(ranlist)
    for sname in ranlist:
        tmpmatch=re.match('([a-zA-Z]+)_([a-zA-Z]+)', sname)
        genus=tmpmatch.group(1)
        species=tmpmatch.group(2)
        tmpspecname = genus + '_' + species
        flag = 0
        #if the species is already in our array, don't add it
        for x in one_per_species:
            if tmpspecname in x:
                flag=1
                break
        if flag == 0:
            one_per_species.append(sname)
        #now do the same for the genera
        flag = 0
        for x in one_per_genus:
            if genus in x:
                flag=1
                break
        if flag==0:
            one_per_genus.append(sname)
    
    return [one_per_genus, one_per_species]       



#sets all exchange fluxes in a model to values between a specific 
#lower and upper bound
#assumes that an exchange reaction is a reaction that has only
#an educt (no product) and where the educt is an external metabolite
#assumes that the compartment variable of the external metablite is 
#set to 'e', tested for the reaction universe from San Roman 2018 
#and should thus also be ok for random viable networks derived from it.
def open_exchange_fluxes(metmodel, lower, upper):
    #set exchange reactions for each external 
    #metabolite with generous influx-outflux boundaries
    for reaction in metmodel.reactions:
        flag=0
        for met in reaction.metabolites:
            #is this metabolite internal? If so, this is not an exchange reaction
            if met.compartment != 'e':
                flag=1
                break
            #is this metabolite a product? If so, this cannot be an exchange reaction either
            if reaction.get_coefficient(met.id)>0:
                flag=1
                break
        #if flag is still zero, the reaction has only external metabolites
        #and no products
        if flag==0:
            reaction.lower_bound = lower  
            reaction.upper_bound = upper
            
#a function to set all exchange fluxes of random viable networks derived from the
#San Roman 2018 universe to specific values            
def set_all_exchange_fluxes_ranvia_SR2018(metmodel, lower, upper):
    tmpctr=0
    for reaction in metmodel.reactions:
        if('EX_' in reaction.id):
            reaction.lower_bound = lower  
            reaction.upper_bound = upper
            tmpctr+=1
            
    
#a function to set the exchange fluxes of metabolites that are keys of 
#dict_lower_bound to the value specified by the value of dict_lower_bound            
def set_some_lower_bounds_ranvia_SR2018(metmodel, dict_lower_bound):
    metlist=sorted(dict_lower_bound.keys())
    
    for reaction in metmodel.reactions:
        if 'EX_' in reaction.id:
            for met in reaction.metabolites:
                if met.id in metlist:
                    reaction.lower_bound = dict_lower_bound[met.id]  

#takes a dictionary keyed with exchange reactions and valued with lower bounds for
#these reactions. Sets those exchange reactions that the model contains to the 
#value in the dictionary
def set_uptake_from_diet_agora(metmodel, diet):
    for reaction in metmodel.reactions:
        if 'EX_' in reaction.id:
            if reaction.id in diet.keys():
                reaction.lower_bound = diet[reaction.id] 
                
# will cycle over all 'EX_' reactions in a model, store their 
# upper and lower bounds in a dictionary, and return an array of these 
# dictionaries
# will work as intended only if exchange fluxes are prefixed by EX_, 
# as in the agora models
def backup_exchange_fluxes(model):
    lower={}
    upper={}
    for r in model.reactions:
        if 'EX_' in r.id: 
            lower[r.id]=r.lower_bound
            upper[r.id]=r.upper_bound
    return [lower, upper]

#will cycle over all 'EX_' reactions in a model, and restore their 
#upper and lower bounds from the dictionaries lower and upper, keyed with r.id  
#will work as intended only if exchange fluxes are prefixed by EX_, as in the agora models
def restore_exchange_fluxes(model, lower, upper):
    for r in model.reactions:
        if 'EX_' in r.id:
            r.lower_bound=lower[r.id]
            r.upper_bound=upper[r.id]
            
#expects a comma-delimited infile whose format is that of Magnusdottir_2017_TableS12_diets.csv
#diet is a string that should be either highfib or western

#Note that in the original Magnusdottir diet list, the metabolite id is NOT exactly the same as in
#the agora models. For example, the compartment has been stripped from the id, and
#there also seem to be some minor name changes, e.g., from glc_D to glc
#this means that exchange reactions need to be mapped onto metabolite ids
#in a separate step/function

#will return a dictionary keyed by exchange reaction and the lower bound
#(negative value) of the corresponding exchange reaction
def load_diet_agora(infile, diet):
    
    df =  pd.read_csv(infile, sep=',')
    
    western_lb={}
    highfib_lb={}
    
    for i in df.index:
        western_lb[df.iloc[i]['Exchange_reaction']] = - df.iloc[i]['western']
        highfib_lb[df.iloc[i]['Exchange_reaction']] = - df.iloc[i]['high_fiber']
    if diet == 'western':
        return western_lb
    elif diet == 'highfib':
        return highfib_lb
    else:
        print('error_aw: invalid diet in load_diet_agora')

 
       

#initializes the environment (medium) by setting the concentrations of 
#all external metabolites of a model to initconc, creates a dictionary with 
#(key, value)=(metabolite cobra id (NOT name!), concentration) for all external metabolites
#requires that external metabolites in the model are in compartment 'e'
def init_medium(metmodel, initconc):
    medium={}
    for met in metmodel.metabolites:
        #don't include these as media components
        if met.id == "o2_e" or met.id == "co2_e"  or met.id == "h_e" or met.id == "h2o_e" or met.id == "nh4_e" or met.id == "pi_e": 
            continue
        if met.compartment == 'e':
            medium[met.id]=initconc
    return medium

#initializes the environment (medium) by setting the concentrations of 
#all external metabolites that occur in at least one of the passed 
#models to initconc, modeldict is a dictionary of models keyed by model name
#creates and returns a dictionary with 
#(key, value)=(metabolite cobra id (NOT name!), concentration) for all external metab
def init_medium_multiple_models(modeldict, initconc):
    medium={}
    for modelname, metmodel in modeldict.items():
        for met in metmodel.metabolites:
            #don't include these as media components
            if met.id == "o2_e" or met.id == "co2_e"  or met.id == "h_e" or met.id == "h2o_e" or met.id == "nh4_e" or met.id == "pi_e": 
                continue
            #if we found an external metabolite that is not already in the medium
            if met.compartment == 'e' and met not in medium:
                medium[met.id]=initconc
    return medium


#initializes the environment (medium) by setting the concentrations of 
#all external metabolites that occur in at least one of the passed 
#models to defaultconc or to an entry of the dictionary conc (keyed
#with metabolite name and valued with a concentration), 
#EXCEPT if the metabolite name is in the passed
#list non_media_metabolites, which should not be considered as part of the 
#medium
    
#modeldict is a dictionary of models keyed by model name
#creates and returns a dictionary with 
#(key, value)=(metabolite cobra id (NOT name!), concentration) for all external metabolites
def init_medium_multiple_models_ranvia(modeldict, non_media_metabolites, defaultconc, conc):
    medium={}
    for metmodel in modeldict.values():
        for met in metmodel.metabolites:
            #if the metabolite should not be added to the medium
            if met.id in non_media_metabolites:
                continue
            #if we found an external metabolite that is not already in the medium
            if met.compartment == 'e' and met not in medium:
                if met.id in conc.keys():
                    medium[met.id]=conc[met.id]
                else: 
                    medium[met.id]=defaultconc
    return medium


#defines a medium in which all external metabolites that exist in at least one 
#agora model have a concentration of defaultconc, except those in a dictionary conc, which
#are set to the concentration in that dictionary
def init_medium_multiple_models_agora(modeldict, defaultconc, conc):
    medium={}
    for metmodel in modeldict.values():
        for met in metmodel.metabolites:
            #if we found an external metabolite that is not already in the medium
            if met.compartment == 'e' and met not in medium:
                if met.id in conc.keys():
                    medium[met.id]=conc[met.id]
                else: 
                    medium[met.id]=defaultconc
    return medium


#set all those nutrient concentrations that lie below minconc to zero
def reset_medium(medium, minconc):
    for met in medium.keys():
        if medium[met]<=minconc:
            medium[met]=0
    


#takes a dictionary medium (m.id, conc), 
#then calculates the total concentration of 
#carbon molecules, based on those metabolites for which a 
#formula is available, does explicitly exclude from this calculation
#those metabolites whose KEGG IDs are in the list 'exclude'

#this should also work for agora models
def tot_C_conc_ranvia(speciesdict, medium, exclude):   
    tot_C_conc=0
    for m in medium.keys():
        if medium[m]>0 and m not in exclude:
            #now try to find the metabolite in at least one of the passed species
            found=0
            for i in speciesdict.keys():
                try:
                    mtmp=speciesdict[i].metabolites.get_by_id(m)
                    found=1
                    break
                #if the metabolite is not found, get by id will return an error
                except KeyError:
                    found=0
            if found==0:
                print('warning tot_C_conc_ranvia: did not find metabolite ', m, ' in any species, returning na')
                return 'na'
            if mtmp.formula=='':
                print('warning tot_c_conc_ranvia: metabolite with no formula: ', mtmp )
                continue
            else:
                #now search for a C followed by a number or by no number
                tmpmatch=re.match('C([0-9]{0,3})', mtmp.formula)
                #the metabolite has no carbon atoms
                if tmpmatch == None:
                    Catoms=0
                else:
                    Catomstmp=tmpmatch.group(1)
                    #if the C is not followed by a numeral, then the compound has 1 
                    #C atoms
                    if Catomstmp == '':
                        Catoms=1
                    else: 
                        Catoms=int(Catomstmp)
                tot_C_conc += Catoms*medium[m]
    return tot_C_conc

#performs FBA on a passed network (uptake flux bounds must already have been set) 
#determines the uptake fluxes for all the carbon sources in C sources (a list of metabolites), 
#calculated, based on the C1 units in those carbon sources
#then calculates the total excretion flux of all carbon sources
    
#in the calculations of both uptake and excretion fluxes, excludes those
#carbon sources in exclude

#flags species for excessive excretion of it excretes more carbon than it 
#imports

def excessive_excretion_ranvia(species, Csources, exclude): 
    
    sol=species.optimize()
    if sol.objective_value==0:
        print('warning_aw: excessive excretion ranvia: no biomass growth flux ')
    tot_C1_uptake=0
    for m in Csources:
        if m not in exclude:
            #now try to find the metabolite in the passed species           
            try:
                mtmp=species.metabolites.get_by_id(m)
                found=1
            #if the metabolite is not found, get by id will return an error
            except KeyError:               
                found=0
            #in this case move on to the next metabolite
            if found==0:
                print('warning excessive excretion ranvia: did not find C source ', m, ' in the species')
                continue
            
            if mtmp.formula=='':
                print('warning excessive excretion ranvia: metabolite with no formula: ', mtmp )
                #move on to next metabolite
                continue
            else:
                #now search for a C followed by a number or by no number
                #print(mtmp.formula)
                tmpmatch=re.match('C([0-9]{0,3})', mtmp.formula)
                #the metabolite has no carbon atoms
                if tmpmatch == None:
                    Catoms=0
                else:
                    Catomstmp=tmpmatch.group(1)
                    #if the C is not followed by a numeral, then the compound has 1 
                    #C atom
                    if Catomstmp == '':
                        Catoms=1
                    else: 
                        Catoms=int(Catomstmp)
                
            
                #next identify the exchange reaction for this metabolite
                for rid in sol.fluxes.keys():
                    rtmp=species.reactions.get_by_id(rid)
                    #if the reaction is an exchange reaction for carbon source m with a negative flux, which means the metabolite is imported
                    if 'EX_' in rid and len(rtmp.metabolites)==1 and m in rtmp.reaction and sol.fluxes[rid] < 0 :
                        C1_uptake = - sol.fluxes[rid]*Catoms #note that uptake fluxes are negative, so we need to invert them here
                        tot_C1_uptake += C1_uptake

        
    #now identify all exchange reactions and determine the carbon flux through those with a positive (excretion) flux
    tot_C1_excretion=0
    for rid in sol.fluxes.keys():
        rtmp=species.reactions.get_by_id(rid)
        #if the reaction is an exchange reaction with a nonzero flux
        if 'EX_' in rid and len(rtmp.metabolites)==1 and sol.fluxes[rid]>0 :
            #this is the metabolite in question
            mtmp=list(rtmp.metabolites.keys())[0]
            #do any calculations only if the metabolite id is not to be excluded
            if mtmp.id not in exclude:
                #if it has a formula, extract the number of C atoms
                if mtmp.formula != '':
                    tmpmatch=re.match('C([0-9]{0,3})', mtmp.formula)
                    #the metabolite has no carbon atoms
                    if tmpmatch == None:
                        Catoms=0
                    else:
                        Catomstmp=tmpmatch.group(1)
                        #if the C is not followed by a numeral, then the compound has 1 
                        #C atom
                        if Catomstmp == '':
                            Catoms=1
                        else: 
                            Catoms=int(Catomstmp)                           
                met_excretion = Catoms * sol.fluxes[rid]
                tot_C1_excretion += met_excretion
    print("excessive excretion ranvia: total C1 uptake: ", tot_C1_uptake , "total C1 excretion: ", tot_C1_excretion)
    if tot_C1_uptake < tot_C1_excretion:
        return 1
    else:
        return 0

#takes a random viable model in which all exchange flux bounds already need to be set
#performs FBA, and asks whether the total excretion of one carbon units (based on those
#metabolites that have a formula) is smaller or equal than the total uptake of one carbon units.
#this must be the case if the biomass growth flux is greater than zero, otherwise mass 
#conservation is violated. Prints a warning and returns 0 if that condition is violated, and
#returns 1 otherwise.
def C1flux_test_massconservation(species):   
   
    sol=species.optimize()
    if sol.objective_value==0:
        print('warning_aw: C1flux_test_massconservation: no biomass growth flux ')
     
    totC1flux=0
    for rid in sol.fluxes.keys():
        rtmp=species.reactions.get_by_id(rid)
        #if the reaction is an exchange reaction with a nonzero flux
        if 'EX_' in rid and len(rtmp.metabolites)==1 and sol.fluxes[rid] !=0 :
            #this is the metabolite to be exchanged
            mtmp=list(rtmp.metabolites.keys())[0]
            #if it has a formula, extract the number of C atoms
            if mtmp.formula != '':
                tmpmatch=re.match('C([0-9]{0,3})', mtmp.formula)
                #the metabolite has no carbon atoms
                if tmpmatch == None:
                    Catoms=0
                else:
                    Catomstmp=tmpmatch.group(1)
                    #if the C is not followed by a numeral, then the compound has 1 
                    #C atom
                    if Catomstmp == '':
                        Catoms=1
                    else: 
                        Catoms=int(Catomstmp)
            C1flux = Catoms * sol.fluxes[rid]            
            
            totC1flux += C1flux
    #imported metabolites make a negative contribution to this flux, and 
    #excretions a positive one, excretions must be at most equal to uptakes,
    #so this flux should never be positive
    if totC1flux <= 0: 
        return 1 
    else: 
        print('warning_aw: C1flux_test_massconservation: more C1 excreted than consumed ')
        return 0 


#map media components on the corresponding exchange reaction in a model
#assumes that media components are in compartment 'e' and will fail otherwise
#this could take some time for a large model and should thus be done sparingly
#medium is a hash whose keys are the metabolite ids (NOT names!) of media components

#returns a dictionary exchange_reaction with media components as keys and 
#exchange reaction ids (NOT names!) as values
    

def map_external_metabolites_onto_exchange_reactions(metmodel, medium):
    exchange_reaction={}
    for r in metmodel.reactions:
        #the reaction is the right one if it contains the 
        #metabolite as educt and no other product
        for met1 in medium.keys():
            flag=0
            for met2 in r.metabolites:
            #does the reaction involve a metabolite other than met1
            #if so, the reaction cannot be the exchange reaction for met1?
                if(met2.id != met1):
                    flag=1
                    break
            if flag==0:
                exchange_reaction[met1]=r.id
    return exchange_reaction

#maps exchange reactions for multiple models onto media components/metabolite ids for these exchange
#reactions, needed for agora models, where there may be non-standard metabolite ids
#in the diet files 

#passes a list of exchange reaction ids ex_r and identifies the corresponding media
#component

#returns a dictionary media_met with exchange reactions as keys and media
#component ids (not names!) as values
def map_exchange_reactions_onto_external_metabolites(modeldict, ex_r):
    media_met={}
    for exr in ex_r:
        flag=0
        for mod in modeldict.keys():
            try:
                r=modeldict[mod].reactions.get_by_id(exr)
                if len(r.metabolites)>1:
                    print('error_aw: echange reaction with more than one metabolite')
                #there should be only one metabolite in this reaction, if the reaction
                #is truly external
                for m in r.metabolites:
                    media_met[exr]=m.id
                flag=1
                #we have found the exchange reaction, so skip the other models -- inelegant but works
                break
            #if this model does not have the exchange reaction, try the next one
            except KeyError:
                flag=0
                pass #do nothing
        #if the external reaction exists in none of the loaded models
        if flag==0:
            print('warning aw: external reaction', exr, ' not found in any loaded model')
    return media_met



#computes a nutrient uptake rate whose dimensions depend on the dimension of Vmax
#usually mmol /g /hr. Assumes that the uptake rate is determined by Michelis
#Menten kinetics with Vmax and Km. Note that a Km=0 effectively amounts to a 
#step function where for any concentration conc, the uptake rate is at its maximum

#the function also ensures that no more nutrient is taken up than exists in the 
#medium, which requires information about the existing biomass and the time interval
#deltat over which the uptake is to take place.

#nutrient [mmoles] is the total amount of nutrient in the medium
#Km is the nutrient amount [mmoles] at half-maximal velocity
#Vmax is the maximal uptake rate [usually mmoles /g /hr]
#biomass is the total amount of biomass in the culture
#deltat is the time interval over which the uptake is to take place
    
#the returned value is the maximal uptake in mmoles (NOT per unit time and unit biomass!)
def max_absolute_nutrient_uptake_MM(nutrient, Vmax, Km, biomass, deltat):
    #MM transport limitation for the nutrient uptake flux in mmoles/g/hr
    maxuptake=Vmax*(nutrient/(Km + nutrient))
    #total mmoles of nutrient imported by the amount of present biomass in the time interval deltat [mmoles] 
    maxuptake=maxuptake*biomass*deltat
    #as defined above, maxuptake could be larger than the total amount of nutrient 
    #present, and the following prevents this problem  
    maxuptake=np.minimum(nutrient, maxuptake)
    return maxuptake

#will map the first occurrence of CDDDD (Kegg compound identifier) in a string
#onto a common compound name, KID_to_name needs to be a dictionary 
#with (KEGGID: common_name) key-value pairs 
def map_KID_to_name(KID_to_name, string): 
    KID=re.findall('C[0-9]{5}', string)
    if len(KID)>0 and KID[0] in KID_to_name.keys(): 
        metname=KID_to_name[KID[0]] 
    else: 
        metname = ''
    return metname


#returns an array of those species that have nonzero biomass
#species is a dict keyed by a species id whose value is a metabolic model
#biomass is a dict keyed by the same ids whose value is the current biomass
def species_present(species, biomass):
    present=[]
    for name in sorted(species.keys()):
        if biomass[name]>0:
            present.append(name)
    return present

#calculates a very general diversity measure due to Hill,
#based on Mittelbach and McGill, p 22
def diversity_Hill(species, biomass, q):
    #relative abundance of species in the community
    relbiom=[]
    for name in species.keys():
        if biomass[name]>0:
            relbiom.append(biomass[name])
    #don't do the calculation if we have fewer than two species
    if len(relbiom)<2: 
        return None
    #calculate relative abundances
    totbiom=0
    for x in relbiom: totbiom += x
    for i in range(0,len(relbiom)):
        relbiom[i] /= totbiom

    
    Dq=0
    #for q=1 the Hill diversity should be the exponential of the Shannon entropy
    if q == 1:
        for x in relbiom:
            Dq -= x*np.log(x)
        Dq = np.exp(Dq)
    else:
        for x in relbiom:
            Dq += x**q
        Dq=Dq**(1/(1-q))
    return Dq


    

"""
computes uptake flux bounds for a nutrient based on the biomass of each strain
that is present, procedure similar to Chio, PloS Comp Bio 2014, simpler than 
San Roman 2018 for multi-species communities

input includes a metabolite and its concentration, multiple metabolic models, their respective biomasses, 
model-specific Michaelis constants Vmax and Km for the nutrient's maximal uptake
 
metid is the id of the nutrient component in question

conc is a scalar that represents the nutrient concentration in mmoles / l (all
calculations are based on a reaction volume of 1 l)

ex_r is a 2d dictionary keyed with modelname and metabolite that identifies which
medium components actually have an exchange reaction in a given model, important
to apply the nutrient partitioning only to strains (models) that can actually 
import the nutrient

models,biomass, Vmax, Km, recent_consumption 
should all be dictionaries keyed with model names

returns  dictionary 'actual_import_bound' keyed by model (strain) that contains the absolute number of 
mmoles that the strain can uptake, considering its present biomass, in the time
interval deltat. This is NOT the uptake per gram biomass and hour but the total uptake.
Note that the dictionary will only have as many keys as there are models that can uptake the 
nutrient
"""
def partition_nutrients(metid, conc, ex_r, models, biomass, Vmax, Km,  deltat):
    #identify an indicator function indicating which models have an exchange function
    #for the nutrient and biomass different from zero, apply the rest of the 
    #function only to them. If there are none, will return an empty dictionary
    valid_model={}
    for modelname, metmodel in models.items():
        if metid in ex_r[modelname].keys() and biomass[modelname]>0:
            valid_model[modelname]=1
        else:
            valid_model[modelname]=0
    
    #compute the absolute Michaelis-Menten transport limit in mmoles of the nutrient for each
    #strain 'modelname', based on the strain's biomass, assumes that Vmax is in
    #mmoles / g / hr
    max_transport={}
    for modelname, metmodel in models.items():
        if valid_model[modelname]==1:
            max_transport[modelname]=Vmax[modelname]*(conc/(Km[modelname] + conc))
            # total mmoles of nutrient that can be imported by the amount of present 
            # biomass of the strain modelnam,e in the time interval deltat [mmoles] 
            max_transport[modelname]=max_transport[modelname]*biomass[modelname]*deltat
    
    #a dictionary to hold the maximal amount of nutrient available to each strain
    max_availability={}
    #will hold the weighted sum of biomasses for all models that consumed the nutrient
    sumbiomass=0
    for modelname, metmodel in models.items():
        if valid_model[modelname]==1:
            sumbiomass=sumbiomass+biomass[modelname]
    for modelname, metmodel in models.items():
        if valid_model[modelname]==1:
            #this is the maximal availability per gram of total biomass (all species) and unit time 
            max_availability[modelname]=conc/(sumbiomass*deltat)
            #calculate the maximal availability in terms of absolute mmoles of nutrient
            #to be imported by the strain modelname in the time interval deltat
            max_availability[modelname]=max_availability[modelname]*biomass[modelname]*deltat    
            #division and multiplication by deltat is redundant, which may make it easier
            #understand the logic of this computation
    
    
    #now set the actual upper bound for import for FBA
    actual_import_bound={}    
    #calculate the minimum of the transport limit and the availability limit, 
    for modelname, metmodel in models.items(): 
        if valid_model[modelname]==1:
            actual_import_bound[modelname]=np.minimum(max_transport[modelname], max_availability[modelname])


    #if the sum of the actual import bounds exceeds the total nutrient concentration, then
    #rescale it such that only at most the total amount of nutrient could be imported 
    #(note that the actual amount imported is determined by FBA)
    #notice also that the amounts calculated here are positive, whereas 
    #uptake fluxes in FBA would be negative. This needs to be set correctly outside this routine.
    max_sum_import=0
    for modelname, metmodel in models.items(): 
         if valid_model[modelname]==1:
             max_sum_import=max_sum_import+actual_import_bound[modelname]
           
    if max_sum_import>conc:
        for modelname, metmodel in models.items():
            if valid_model[modelname]==1:
                actual_import_bound[modelname]=actual_import_bound[modelname]*(conc/max_sum_import)
    #the newly normalized values for the actual import bound add up to at most conc
    return actual_import_bound
    

#simulates a single replacement event of a fraction frac_transfer of an existing 
#culture volume with spent (old) medium by an equal volume of fresh medium, where
#spent medium is represented by a dictionary containing the concentrations of media components
#biomass is a scalar, the amount of biomass in the spent culture
#returns an array with two elements containing (i) a dictionary for the 
#medium after transfer, (ii) the amount of biomass after transfer
def dilute(biomass_old, old, fresh, frac_transfer):
    biomass_new=frac_transfer*biomass_old
    new={}
    #loop over all metabolites in the spent medium
    for met in old.keys():
        new[met] = frac_transfer*old[met]+(1-frac_transfer)*fresh[met]
    return [biomass_new, new]

#same as dilute but for multiple models, where biomass now becomes a 
#dictionary keyed by model name, returns an array with two entries of the
#updated biomass and medium dictionaries
def dilute_multiple_models(biomass_old, medium_old, medium_fresh, frac_transfer):
    biomass_new={}
    for model in biomass_old.keys():
        biomass_new[model]=frac_transfer*biomass_old[model]
    medium_new={}
    #loop over all metabolites in the spent medium
    for met in medium_old.keys():
        medium_new[met] = frac_transfer*medium_old[met]+(1-frac_transfer)*medium_fresh[met]
    return [biomass_new, medium_new]

#environmental perturbation that changes every media chemical with a concentration c>0
#to a value that is c/f with probability 0.5, and c*f with probability 0.5 
def ran_pert_env(medium, frac):
    pert_medium={}
    for met in medium.keys():
        conc=medium[met]
        if conc>0:
            if random.random()>0.5:
                pert_medium[met]=conc*frac
            else:
                pert_medium[met]=conc/frac
        #it's important here that the new medium has the same number of molecules as the old
        #even if their concentrations are zero
        else:
            pert_medium[met]=0
    return pert_medium

#environmental perturbation that changes the concentration of every chemical in the 
#medium supply of the chemostat to a uniformly distributed 
#value in an interval (c/frac, c*frac) 
def ran_pert_env_medfresh(medfresh, medfresh_init, frac):
    pert_medium={}
    for met in medfresh.keys():
        conc=medfresh_init[met]
        if conc>0:
            pert_medium[met]= random.uniform(conc/frac, conc*frac)
        #it's important that the new medium has the same number of molecules as the old
        #even if their concentrations are zero
        else:
            pert_medium[met]=0
    return pert_medium


#perturbation that changes the biomass b of every species that is present
#to a value that is b/f with probability 0.5, and b*f with probability 0.5 
def ran_pert_spec(biomass, fac):
    #set up a new biomass dictionary so as not to touch the old one
    pert_biomass={}
    for spec_id in biomass.keys():
        biom=biomass[spec_id]
        if biom>0:
            if random.random()>0.5:
                pert_biomass[spec_id]=biom*fac
            else:
                pert_biomass[spec_id]=biom/fac
        #if the species is not currently present, do nothing
        else:
            pert_biomass[spec_id]=0
    return pert_biomass



#evaluates for each strain whether (i) its biomass as a fraction of total biomass
#has fallen below ext_thresh or if (ii) if it is the only strain present and 
#its total biomass has fallen below 1e-13. If so, appends it to a list of strain
#names slated for extinction. Does NOT set the biomass to zero, this should
#be done in the calling routine

#biomassarr is a dictionary keyed with strain ids that points to an array that holds a history of biomass values
#totbiomassarr is an array that holds the history of total biomass    

#returns a list of those strain ids that have gone extinct
def extinction1(strains, biomassarr, totbiomassarr, ext_threshold):
    extinctions=[]
    #first determine which strains are actually present
    present=[]    
    for name in strains.keys():
        if biomassarr[name][-1]>0:
            present.append(name)
    #if only one strain is present, then do not allow its biomass to fall below that
    #of one ecoli cell
    if len(present)==1:
        if biomassarr[present[0]][-1]<1e-13:
            extinctions.append(present[0])
    elif len(present)>1:
        for name in present:      
        #execute only for species where the most recent biomass value is >0
            if (biomassarr[name][-1]/totbiomassarr[-1])<ext_threshold:
                extinctions.append(name)
    else:
        print("error in extinctions 1: no strain present")
    return extinctions


#evaluates for each strain whether its biomass as a fraction of total biomass
#has fallen below ext_thresh and if so, appends it to a list of strain
#names slated for extinction. Does NOT set the biomass to zero, this should
#be done in the calling routine
    
#biomass is a dictionary keyed with strain ids that points to current biomass of the named strain
#totbiomass holds the current total biomass
    
#similar to extinction 1 except for the passed data structures

#returns a list of those strain ids that have gone extinct    
def extinction2(strains, biomass, ext_threshold):
    extinctions=[]
    #first calculate total biomass
    totbiomass=0
    for name in strains.keys():
        totbiomass=totbiomass+biomass[name]
    for name in strains.keys():
        #execute only for species where the most recent biomass value is >0
        if biomass[name]>0 and (biomass[name]/totbiomass)<ext_threshold:
            extinctions.append(name)
    return extinctions            
                    

#evaluates whether a collection of species (dict keyed by a species id)
#has reached a biomass equilibrium within the last timewin hours
#biomass is a dict of biomass values over time keyed by species id
#each value is an array showing the species biomass evolution over time
#totbiomass is an array with total biomass values of all species
#deltat is the time step width for the simulation

#tests for each species that accounts for more than 0.1 percent of total biomass
#whether the difference between the maximum and the minimum of 
#the recent biomass does not exceed a fraction devtol of the mean recent biomass
#the 0.1 percent threshold exists because otherwise species with tiny contributions 
#to biomass may make it seem as if the total community was not in equilibrium

#returns one if for each species an equilibrium has been reached, zero otherwise

def equilibrium_cstat(strains, biomass, deltat, timewin, devtol):
    #first build a new biomass array for each species that contains 
    #only the last deltat*timewin timesteps
    n_steps=int(timewin/deltat)

    #if the population history is too short, then we cannot evaluate equilibration   
    firstkey=list(biomass.keys())[0]
    if len(biomass[firstkey])<n_steps:
        return 0
    
    #the total biomass present at the current time
    totbiomass=0
    #a dictionary that will only hold recent biomass values for each species
    recent_biomass={}
    for name in strains.keys():
        #do the calculation only for species where the most recent biomass value is >0
        if biomass[name][-1]>0:
            #this may be a shallow copy, but it should not matter, because 
            #we are not editing the object in question
            recent_biomass[name]=biomass[name][-n_steps:]
            totbiomass += biomass[name][-1]
    #notice that the quantity below should be the same when calculated for absolute
    #biomass or for relative biomass
    for name in recent_biomass.keys():
        stat=(np.max(recent_biomass[name])- np.min(recent_biomass[name]))/np.mean(recent_biomass[name])
        if stat>devtol and biomass[name][-1]/totbiomass > 0.001:
            return 0
    return 1




        
#takes a set of species, with a given biomass distribution and a given medium
#and performs dFBA for a maximum time interval of maxtime in time steps of deltat
#with a dilution rate D by fresh medium
        
#checks every testint hours whether equilibrium has been reached
#based on the 'equilibrium_cstat' function with deviation devtol
#ex_r is a data structure to map media components onto
#exchange reactions
#Vmax, and Km hold Michaelis transport constants for each species and nutrient 
 
#resets the medium to render all values below threshold minconc to zero
       
#returns an array of three entries, a flag that indicates whether (1) or not (0) equilibrium
#has been reached, the new biomass of all species, and the new medium composition (at the time of return)
def equilibrate_cstat(species, biomass, medium, medfresh, D, ex_r, Vmax, Km, deltat, testint, maxtime, devtol, minconc):

    #first set up new arrays for biomass and medium so that we don't have to alter the
    #passed ones
    #also set up arrays keyed by species and nutrient to record biomass and 
    #medium over time       
    biom={}
    for s in biomass.keys():
        biom[s]=biomass[s]
    biom_evolve={}   #a dictionary keyed with species name that points to an array
                     #with species biomass as a function of time
    for s in biomass.keys():           
        biom_evolve[s]=[]
        biom_evolve[s].append(biom[s])

    med={}
    for m in medium.keys():
        med[m]=medium[m]
        
    #set all medium components that have fallen below a threshold value to zero, to avoid
    #solver problems with very small values
    reset_medium(med, minconc)
     
    
    med_evolve={}
    #one key points to an array of time steps
    med_evolve['t']=[]
    med_evolve['t'].append(0)
    #the others point to an array of nutrient concentrations
    for m in med.keys():
        med_evolve[m]=[]
        med_evolve[m].append(med[m])
        
    #master loop over time steps
    for timestep in range(1, int(maxtime/deltat)+1): 
        tcurr=timestep*deltat
        for m in med.keys():
            #calculate uptake bound of nutrient m by each species
            uptake = partition_nutrients(m, med[m], ex_r, species, biom, Vmax[m], Km[m],  deltat)
            #do the following only for those species names in uptake, which are those with 
            #biomass greater than zero that have an exchange reaction for the nutrient
            for s in uptake.keys(): 
                #identify the exchange reaction for the current model and nutrient
                exchange=species[s].reactions.get_by_id(ex_r[s][m])
                #note that uptake fluxes are negative, hence the minus sign
                #also note that FBA works with fluxes of mmoles/g/h, hence the division
                exchange.lower_bound= - uptake[s]/(deltat*biom[s])
                
        solution={} #dictionary for FBA solution of each species
        biom_added={} #currently added biomass of each species
        for (s, model) in species.items():
            #only do FBA if a model actually has positive biomass
            if biom[s]>0:
                solution[s] = model.optimize()
                #if we cannot get to an optimal solution for whatever reason
                #abandon equilibrium and return a non-equilibrium state
                if solution[s].status != 'optimal':
                    print('warning from equilibrate_cstat: no optimal solution found for species', s)
                    return [0, biom, med] 
                
        #now update the concentration of all media components
        #this should come BEFORE the new biomass is computed, because the 
        #fluxes are based on the old biomass values
        #note: a negative exchange flux means that a molecule leaves the medium
        for m in med.keys():
            for (s, model) in species.items():
                #only update media components for species that are actually present            
                if biom[s]>0:
                    #if a specific model does not have metabolite m, retrieving the metabolite 
                    #could give an error            
                    try:
                        model.metabolites.get_by_id(m)
                        #note that this will only be seen if no error occurred
                        med[m]=med[m] + solution[s].fluxes[ex_r[s][m]]*biom[s]*deltat
                    except KeyError:
                        pass #do nothing
                        
        #reset negligible medium components to zero
        reset_medium(med, minconc)

        #now add to the arrays showing time change in media components
        med_evolve['t'].append(timestep*deltat)
        for m in med.keys():
            med_evolve[m].append(med[m])
                                     
        #now update biomass   
        for (s, model) in species.items():
            # the below assumes that the biomass growth flux as a result of FBA
            # gives the actual grams of biomass added per g and hour, and needs no conversion
            # factor
            if biom[s]>0:
                biom_added[s] = (solution[s].objective_value)*biom[s]*deltat
                #leave biomass unchanged, but still update the arrays
            else:
                biom_added[s]=0
            #update all arrays regardless of whether a species is present or not
            biom[s] = biom[s] + biom_added[s]
            biom_evolve[s].append(biom[s])
                
        #now for equilibrium testing
        #begin testing after testint hours have elapsed and at every integer multiple of testint
        if(tcurr>=testint and int(tcurr/testint)==tcurr/testint):
                    
            #check whether we have reached an equilibrium
            #use data from the entire time interval since the last test
            equil=equilibrium_cstat(species, biom_evolve, deltat, testint, devtol)
            if equil==1:
                return [1, biom, med]
        #now apply a dilution flux, i.e., replace a fraction D of the current 
        [biom, post_med] = dilute_multiple_models(biom, med, medfresh, 1-(D*deltat))  
        for m in med.keys():
            med[m]=post_med[m]
       
    #if we get here, then we have never reached equilibrium
    return [0, biom, med]  



#applies n_pert random perturbations of a factor pert to all media components and to the medium feed
#equilibrates the community after each perturbation, and computes the fraction 
#of new equilibria where at least one species has gone to extinction
#returns -1 if no perturbation led to an equilibrium, or if the community is empty
def env_pert_sensitivity_cstat(species, biomass, medium, medfresh, D, ex_r, Vmax, Km, \
                           deltat, testint, maxtime, devtol, ext_thresh, pert, n_pert, minconc):                
    
    #determine all species that are present, if no species is present, then 
    #return -1
    presentctr=0
    for s in biomass.keys():
        if biomass[s]>0:
            presentctr+=1
    if presentctr==0:
        return -1
    
    #a successful perturbation is one where a population has reached equilibrium
    succ_ctr=0
    #number of perturbations that lead to an extinction
    ext_ctr=0
    for pert_ctr in range(0,n_pert):
        pert_medium = ran_pert_env(medium, pert)
        pert_medfresh =  ran_pert_env(medfresh, pert)
                
        [eqflag, eqbiom, eqmed]=equilibrate_cstat(species, biomass, pert_medium, pert_medfresh, \
                                                    D, ex_r, Vmax, Km, deltat, testint, maxtime, devtol, minconc)
        if eqflag==1: 
            succ_ctr += 1
            extspec_pert=extinction2(species, eqbiom, ext_thresh)
            if len(extspec_pert)>0:
                ext_ctr += 1
    if succ_ctr>0:
        frac_ext=ext_ctr/succ_ctr
    #if no perturbation lead to an equilibrium then the fraction of perturbations leading to 
    #extinctions extinctions cannot be calculated
    #so return an impossible value
    else:
        frac_ext=-1
    return frac_ext 



#applies n_pert random perturbations to the medium feed, where
#each perturbation creates a new uniformly distributed
#random concentration of each nonzero nutrient around the 
#loaded parameter value of the medium feed
#equilibrates the community after each perturbation, and computes the fraction 
#of new equilibria where at least one species has gone to extinction
#returns -1 if no perturbation led to an equilibrium, or if the community is empty
def env_pert_sensitivity_medfresh_cstat(species, biomass, medium, medfresh, medfresh_init, D, ex_r, Vmax, Km, \
                           deltat, testint, maxtime, devtol, ext_thresh, pert, n_pert, minconc):                
    
    #determine all species that are present, if no species is present, then 
    #return -1
    presentctr=0
    for s in biomass.keys():
        if biomass[s]>0:
            presentctr+=1
    if presentctr==0:
        return -1
    
    #a successful perturbation is one where a population has reached equilibrium
    succ_ctr=0
    #number of perturbations that lead to an extinction
    ext_ctr=0
    for pert_ctr in range(0,n_pert):
        #key difference to the other environmental perturbation routine
        pert_medfresh=ran_pert_env_medfresh(medfresh, medfresh_init, pert)
        #do not perturb the medium itself
        pert_medium=medium

        [eqflag, eqbiom, eqmed]=equilibrate_cstat(species, biomass, pert_medium, pert_medfresh, \
                                                    D, ex_r, Vmax, Km, deltat, testint, maxtime, devtol, minconc)
        if eqflag==1: 
            succ_ctr += 1
            #determine which species became extinct
            extspec_pert=extinction2(species, eqbiom, ext_thresh)
            if len(extspec_pert)>0:
                ext_ctr += 1
    if succ_ctr>0:
        frac_ext=ext_ctr/succ_ctr
    #if no perturbation lead to an equilibrium then the fraction of perturbations leading to 
    #extinctions extinctions cannot be calculated
    #so return a meaningless value
    else:
        frac_ext=-1
    return frac_ext   
  
                    
#applies n_pert random perturbations of a factor pert to all species biomasses
#equilibrates the community after each perturbation, and computes the fraction 
#of new equilibria where at least one species has gone to extinction
#returns -1 if no perturbation led to an equilibrium, or of the community is empty
def spec_pert_sensitivity_cstat(species, biomass, medium, medfresh, D, ex_r, Vmax, Km, \
                           deltat, testint, maxtime, devtol, ext_thresh, pert, n_pert, minconc):

    #determine all species that are present, if no species is present, then 
    #return -1
    presentctr=0
    for s in biomass.keys():
        if biomass[s]>0:
            presentctr+=1
    if presentctr==0:
        return -1

                
    #a successful perturbation is one where a population has reached equilibrium
    succ_ctr=0
    #number of perturbations that lead to an extinction
    ext_ctr=0
    for pert_ctr in range(0,n_pert):
       
        pert_biomass= ran_pert_spec(biomass, pert)

        [eqflag, eqbiom, eqmed]=equilibrate_cstat(species, pert_biomass, medium, medfresh, \
                                                    D, ex_r, Vmax, Km, deltat, testint, maxtime, devtol, minconc)
        if eqflag==1: 
            succ_ctr += 1
            #determine which species became extinct
            extspec_pert=extinction2(species, eqbiom, ext_thresh)
            
            if len(extspec_pert)>0:
                ext_ctr += 1
    if succ_ctr>0:
        frac_ext=ext_ctr/succ_ctr
    #if no perturbation lead to an equilibrium then the fraction of perturbations leading to 
    #extinctions extinctions cannot be calculated
    #so return a meaningless value
    else:
        frac_ext=-1
    return frac_ext   
                    

#given a community, what fraction of its species are viable on the fresh medium 
#when grown there alone, returns -1 if the community contains no species
def frac_viable_on_medfresh(comm, medfresh, D, ex_r, Vmax, Km, deltat):
    #create a new community, since we will be manipulating exchange fluxes 
    commtmp=copy.deepcopy(comm)
    biomtmp=copy.deepcopy(comm.biomass)
    
    
    #determine all species that are present
    present=[]
    for s in commtmp.biomass.keys():
        if commtmp.biomass[s]>0:
            present.append(s)
    #if no species are present, return -1
    if len(present)==0:
        return -1
        
        
    #number of species that are viable by themselves on fresh medium
    n_viable=0
    #cycle over all species that are present
    for s1 in present:
        #set the biomass of all other species to zero
        biomtmp[s1]=commtmp.biomass[s1]
        for s2 in present:
            if s2 != s1:
                biomtmp[s2]=0
        #now determine the uptake rate of the one remaining species s1
        #for the nutrients in the fresh medium, allocate all the nutrients
        #to this species
        for m in commtmp.medium.keys():
            uptake = partition_nutrients(m, medfresh[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
            for s3 in uptake.keys(): 
                #identify the exchange reaction for the current model and nutrient
                exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                #note that uptake fluxes are negative, hence the minus sign
                #also note that FBA works with fluxes of mmoles/g/h, hence the division
                exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])
        #now ask whether the species alone is viable
        sol = commtmp.species[s1].slim_optimize()
        if sol>0:  n_viable+=1
        
    return n_viable/len(present)
        
        


    
#given a community of at least two species, what fraction of its species grows more slowly, 
#equally fast, or faster when ALL other species are present as opposed to 
#when the species is alone
#returns an array of three elements [competition, neutral, facilitation] or [] if 
#called for a community with fewer than two pecies 
def frac_comp_neut_fac_all(comm, D, ex_r, Vmax, Km, deltat):
    if len(comm.species.keys())<2:
        print('error in frac_comp_neut_fac_all')
        return []
    
    #create a new community, since we will be manipulating exchange fluxes 
    commtmp=copy.deepcopy(comm)
    biomtmp=copy.deepcopy(comm.biomass)
    
    #determine all species that are present
    present=[]
    for s in commtmp.biomass.keys():
        if commtmp.biomass[s]>0:
            present.append(s)
            
    #first determine the growth rate of each species when all other species are
    #present
    jointgrowth={}
    for m in commtmp.medium.keys():
        uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, commtmp.biomass, Vmax[m], Km[m],  deltat)
        for s in uptake.keys(): 
            #identify the exchange reaction for the current model and nutrient
            exchange=commtmp.species[s].reactions.get_by_id(ex_r[s][m])
            #note that uptake fluxes are negative, hence the minus sign
            #also note that FBA works with fluxes of mmoles/g/h, hence the division
            exchange.lower_bound= - uptake[s]/(deltat*commtmp.biomass[s])
    for s in commtmp.biomass.keys():
        jointgrowth[s] = commtmp.species[s].slim_optimize()        
    
        #number of interactions that are competitive, neutral, or facilitative
    n_comp=0
    n_neut=0
    n_faci=0
    #cycle over all species that are present
    for s1 in present:
        #set the biomass of all other species to zero
        biomtmp[s1]=commtmp.biomass[s1]
        for s2 in present:
            if s2 != s1:
                biomtmp[s2]=0
        #now determine the uptake rate of the one remaining species s1
        #for the nutrients in the fresh medium, allocate all the nutrients
        #to this species
        for m in commtmp.medium.keys():
            uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
            for s3 in uptake.keys(): 
                #identify the exchange reaction for the current model and nutrient
                exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                #note that uptake fluxes are negative, hence the minus sign
                #also note that FBA works with fluxes of mmoles/g/h, hence the division
                exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])
        #now ask whether the species alone grows more slowly, equally fast, or faster 
        #than when all others are present
        sol = commtmp.species[s1].slim_optimize()
        if sol>jointgrowth[s1]:  
            n_comp+=1
        elif sol==jointgrowth[s1]: 
            n_neut+=1
        elif sol<jointgrowth[s1]: 
            n_faci+=1    
        
    return [n_comp/len(present), n_neut/len(present), n_faci/len(present)]

#like frac_comp_neut_fac_all, but applies one dilution step to the medium
#so that the community has some fresh nutrient, otherwise
#the glucose concentration is effectively zero in glucose minimal medium.

#given a community of at least two species, what fraction of its species grows more slowly, 
#equally fast, or faster when ALL other species are present as opposed to 
#when the species is alone?
    
#returns an array of three elements [competition, neutral, facilitation] or [] if the
#call was for a community with fewer than two pecies 
def frac_comp_neut_fac_all_postdil(comm, medfresh, D, ex_r, Vmax, Km, deltat):
    if len(comm.species.keys())<2:
        print('error in frac_comp_neut_fac_all')
        return []
    
    #create a new community, since we will be manipulating exchange fluxes 
    commtmp=copy.deepcopy(comm)
    
    #now apply a dilution step to add a squirt of nutrients
    [commtmp.biomass, commtmp.medium] = dilute_multiple_models(commtmp.biomass, commtmp.medium, medfresh, 1-(D*deltat))
    biomtmp=copy.deepcopy(commtmp.biomass)

    
    #determine all species that are present
    present=[]
    for s in commtmp.biomass.keys():
        if commtmp.biomass[s]>0:
            present.append(s)
            
    #first determine the growth rate of each species when all other species are
    #present
    jointgrowth={}
    for m in commtmp.medium.keys():
        uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, commtmp.biomass, Vmax[m], Km[m],  deltat)
        for s in uptake.keys(): 
            #identify the exchange reaction for the current model and nutrient
            exchange=commtmp.species[s].reactions.get_by_id(ex_r[s][m])
            #note that uptake fluxes are negative, hence the minus sign
            #also note that FBA works with fluxes of mmoles/g/h, hence the division
            exchange.lower_bound= - uptake[s]/(deltat*commtmp.biomass[s])
    for s in commtmp.biomass.keys():
        jointgrowth[s] = commtmp.species[s].slim_optimize()        
    
    #number of interactions that are competitive, neutral, or facilitative
    n_comp=0
    n_neut=0
    n_faci=0
    #cycle over all species that are present
    for s1 in present:
        #set the biomass of all other species to zero
        biomtmp[s1]=commtmp.biomass[s1]
        for s2 in present:
            if s2 != s1:
                biomtmp[s2]=0
        #now determine the uptake rate of the one remaining species s1
        #for the nutrients in the fresh medium, allocate all the nutrients
        #to this species
        for m in commtmp.medium.keys():
            uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
            for s3 in uptake.keys(): 
                #identify the exchange reaction for the current model and nutrient
                exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                #note that uptake fluxes are negative, hence the minus sign
                #also note that FBA works with fluxes of mmoles/g/h, hence the division
                exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])
        #now ask whether the species alone is viable
        sol = commtmp.species[s1].slim_optimize()
        if sol>jointgrowth[s1]:  
            n_comp+=1
        elif sol==jointgrowth[s1]: 
            n_neut+=1
        elif sol<jointgrowth[s1]: 
            n_faci+=1    
        
    return [n_comp/len(present), n_neut/len(present), n_faci/len(present)]


#given a community of at least two species, what fraction of its species grows more slowly, 
#equally fast, or faster when paired with each other species as opposed to  
#when the species is alone
    
#returns an array of three elements [competition, neutral, facilitation] or [] if the
#call was for a community with fewer than two species
def frac_comp_neut_fac_pairwise(comm, D, ex_r, Vmax, Km, deltat):
    if len(comm.species.keys())<2:
        print('error in frac_comp_neut_fac_pairwise')
        return []
    
    #create a new community, since we will be manipulating exchange fluxes 
    commtmp=copy.deepcopy(comm)
    biomtmp=copy.deepcopy(comm.biomass)
    
   
    #determine all species that are present
    present=[]
    for s in commtmp.biomass.keys():
        if commtmp.biomass[s]>0:
            present.append(s)
            
    
    growthalone={}
    #cycle over all species that are present and determine their growth alone
    for s1 in present:
        #set the biomass of all other species to zero
        biomtmp[s1]=commtmp.biomass[s1]
        for s2 in present:
            if s2 != s1:
                biomtmp[s2]=0
        #now determine the uptake rate of the one remaining species s1
        #for the nutrients in the fresh medium, allocate all the nutrients
        #to this species
        for m in commtmp.medium.keys():
            uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
            for s3 in uptake.keys(): 
                #identify the exchange reaction for the current model and nutrient
                exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                #note that uptake fluxes are negative, hence the minus sign
                #also note that FBA works with fluxes of mmoles/g/h, hence the division
                exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])
        #now ask whether the species alone is viable
        growthalone[s1] = commtmp.species[s1].slim_optimize()

    #now cycle over all species PAIRS, and set only the biomass corresponding to these
    #species to a value different from zero
    n_comp=0
    n_neut=0
    n_faci=0
    for s1 in present:
        for s2 in present:
             if s1 != s2:
                 biomtmp[s1]=commtmp.biomass[s1]
                 biomtmp[s2]=commtmp.biomass[s2]
                 #set the biomass of all other species to zero
                 for s3 in present:
                     if s3 != s1 and s3 !=s2:
                             biomtmp[s3]=0
                 #now determine the uptake rate of the species pair s1 and s2
                 #for the nutrients in the fresh medium, allocate all the nutrients
                 #to this species
                 for m in commtmp.medium.keys():
                    uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
                    for s3 in uptake.keys(): 
                        #identify the exchange reaction for the current model and nutrient
                        exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                        #note that uptake fluxes are negative, hence the minus sign
                        #also note that FBA works with fluxes of mmoles/g/h, hence the division
                        exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])
                 #now determine growth of s1 in the presence of s2
                 sol = commtmp.species[s1].slim_optimize()
                 if sol<growthalone[s1]:  
                     n_comp+=1
                 elif sol==growthalone[s1]: 
                     n_neut+=1
                 elif sol>growthalone[s1]: 
                     n_faci+=1    
    n_pairs=len(present)*(len(present)-1)
    return [n_comp/n_pairs, n_neut/n_pairs, n_faci/n_pairs]

#like frac_comp_neut_fac_pairwise, but applies one dilution step to the medium
#so that the community has some fresh nutrient, otherwise
#the glucose concentration is effectively zero in glucose minimal medium.
       
#given a community of at least two species, what fraction of its species grows more slowly, 
#equally fast, or faster when paired with each other species as opposed to  
#when the species is alone
    
#returns an array of three elements [competition, neutral, facilitation] or [] if the
#call was for a community with fewer than two species
    
def frac_comp_neut_fac_pairwise_postdil(comm, medfresh, D, ex_r, Vmax, Km, deltat):
    if len(comm.species.keys())<2:
        print('error in frac_comp_neut_fac_pairwise')
        return []
    
    
    #create a new community, since we will be manipulating exchange fluxes 
    commtmp=copy.deepcopy(comm)
    
    #now apply a dilution step to add a squirt of nutrients
    [commtmp.biomass, commtmp.medium] = dilute_multiple_models(commtmp.biomass, commtmp.medium, medfresh, 1-(D*deltat))
    biomtmp=copy.deepcopy(commtmp.biomass)

    #determine all species that are present
    present=[]
    for s in commtmp.biomass.keys():
        if commtmp.biomass[s]>0:
            present.append(s)
            
    
    growthalone={}
    #cycle over all species that are present and determine their growth alone
    for s1 in present:
        #set the biomass of all other species to zero
        biomtmp[s1]=commtmp.biomass[s1]
        for s2 in present:
            if s2 != s1:
                biomtmp[s2]=0
        #now determine the uptake rate of the one remaining species s1
        #for the nutrients in the fresh medium, allocate all the nutrients
        #to this species
        for m in commtmp.medium.keys():
            uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
            for s3 in uptake.keys(): 
                #identify the exchange reaction for the current model and nutrient
                exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                #note that uptake fluxes are negative, hence the minus sign
                #also note that FBA works with fluxes of mmoles/g/h, hence the division
                exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])
        #now ask whether the species alone is viable
        growthalone[s1] = commtmp.species[s1].slim_optimize()

    #now cycle over all species PAIRS, and set only the biomass corresponding to these
    #species to a value different from zero
    n_comp=0
    n_neut=0
    n_faci=0
    for s1 in present:
        for s2 in present:
             if s1 != s2:
                 biomtmp[s1]=commtmp.biomass[s1]
                 biomtmp[s2]=commtmp.biomass[s2]
                 #set the biomass of all other species to zero
                 for s3 in present:
                     if s3 != s1 and s3 !=s2:
                             biomtmp[s3]=0
                 #now determine the uptake rate of the species pair s1 and s2
                 #for the nutrients in the fresh medium, allocate all the nutrients
                 #to this species
                 for m in commtmp.medium.keys():
                    uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
                    for s3 in uptake.keys(): 
                        #identify the exchange reaction for the current model and nutrient
                        exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                        #note that uptake fluxes are negative, hence the minus sign
                        #also note that FBA works with fluxes of mmoles/g/h, hence the division
                        exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])
                 #now determine growth of s1 in the presence of s2
                 sol = commtmp.species[s1].slim_optimize()
                 if sol<growthalone[s1]:  
                     n_comp+=1
                 elif sol==growthalone[s1]: 
                     n_neut+=1
                 elif sol>growthalone[s1]: 
                     n_faci+=1    
    n_pairs=len(present)*(len(present)-1)
    return [n_comp/n_pairs, n_neut/n_pairs, n_faci/n_pairs]
        

                 
#for a community of species, determines the fraction of species pairs s1 s2 where s2
#crossfeeds on a waste product of s1. Criteria: s1 must excrete m, s2 must import it, 
#and if we set the uptake rate of m in s2 to zero, that must reduce the growth rate of s2
#by at least minfluxdiff
def frac_crossfeed(comm, D, ex_r, Vmax, Km, deltat, minfluxdiff):
    if len(comm.species.keys())<2:
        print('error in frac_crossfeed')
        return []
    
    #create a new community, since we will be manipulating exchange fluxes 
    commtmp=copy.deepcopy(comm)
    biomtmp=copy.deepcopy(comm.biomass)
    
    #determine all species that are present
    present=[]
    for s in commtmp.biomass.keys():
        if commtmp.biomass[s]>0:
            present.append(s)
    #number of cross-feeding interactions
    n_crossfeed=0
    #biomass growth when both species are present
    jointgrowth={}
    for s1 in present:
        for s2 in present:
            #cross-feeding is not necessarily symmetrical
            #so the following has two parts, for crossfeeding s1->s2, and s2->s1
             if s1 < s2:
                 #first determine the uptake rate of all present species
                 #for the nutrients in the present medium
                 for m in commtmp.medium.keys():
                     
                    uptake = partition_nutrients(m, commtmp.medium[m], ex_r, commtmp.species, biomtmp, Vmax[m], Km[m],  deltat)
                    for s3 in uptake.keys(): 
                        #identify the exchange reaction for the current model and nutrient
                        exchange=commtmp.species[s3].reactions.get_by_id(ex_r[s3][m])
                        #note that uptake fluxes are negative, hence the minus sign
                        #also note that FBA works with fluxes of mmoles/g/h, hence the division
                        exchange.lower_bound= - uptake[s3]/(deltat*biomtmp[s3])      
                 jointgrowth[s1]=commtmp.species[s1].optimize()
                 jointgrowth[s2]=commtmp.species[s2].optimize()
                 
                 #now cycle over all the metabolites excreted by s1
                 crossfeedflag=0
                 for m in commtmp.medium.keys():
                     #do the following ONLY if both species have an exchange reaction for the metabolite
                     if m in ex_r[s1].keys() and m in ex_r[s2].keys():
                         exchange1=commtmp.species[s1].reactions.get_by_id(ex_r[s1][m])
                         if jointgrowth[s1].fluxes[exchange1.id]>0:
                             #does s2 import m?
                             exchange2=commtmp.species[s2].reactions.get_by_id(ex_r[s2][m])
                             if jointgrowth[s2].fluxes[exchange2.id]<0:
                                 #set the import of m to zero
                                 exchange2_backup=exchange2.lower_bound
                                 exchange2.lower_bound=0
                                 sol=commtmp.species[s2].slim_optimize()
                                 exchange2.lower_bound=exchange2_backup
                                 #if the new biomass growth flux is smaller than the old flux by a factor that 
                                 #exceeds the minimal flux difference threshold, then cross-feeding occurs
                                 if jointgrowth[s2].objective_value - sol > minfluxdiff:
                                     #if at least one metabolite is involved in cross-feeding call
                                     #break out of the loop, because we know that cross feeding takes place
                                     crossfeedflag=1
                                     break
                 if crossfeedflag==1:
                     n_crossfeed+=1
                     
                 #now do the same but for s1 and s2 reversed
                 #cycle over all the metabolites excreted by s2
                 crossfeedflag=0
                 for m in commtmp.medium.keys():
                     #do the following ONLY if both species have an exchange reaction for the metabolite
                     if m in ex_r[s1].keys() and m in ex_r[s2].keys():
                         exchange2=commtmp.species[s2].reactions.get_by_id(ex_r[s2][m])
                         #does species 2 excrete m
                         if jointgrowth[s2].fluxes[exchange2.id]>0:
                             #does s1 import m?
                             exchange1=commtmp.species[s1].reactions.get_by_id(ex_r[s1][m])
                             if jointgrowth[s1].fluxes[exchange1.id]<0:
                                 #set the import of m to zero
                                 exchange1_backup=exchange1.lower_bound
                                 exchange1.lower_bound=0
                                 sol=commtmp.species[s1].slim_optimize()
                                 exchange1.lower_bound=exchange1_backup
                                 #if the new biomass growth flux is smaller than the old flux by a factor that 
                                 #exceeds the minimal flux difference threshold, then cross-feeding occurs
                                 if jointgrowth[s1].objective_value - sol > minfluxdiff:                                 
                                     #if at least one metabolite is involved in cross-feeding call
                                     #break out of the loop, because we know that cross feeding takes place
                                     crossfeedflag=1
                                     break
                 if crossfeedflag==1:
                     n_crossfeed+=1
                
     
    n_pairs=len(present)*(len(present)-1)
    return n_crossfeed/n_pairs
             

#prints some basic statistics from the current assembly sequence

#the passed arguments are dictionaries keyed by time that point to 
#species composition at time t, the species invading at time t, and the species going extinct
#at time t
def print_assembly_stats(assembly_seq, invasion, extinction, pert_sens):
    #this is an array of times at which the assembly changed
    assembly_change_t=sorted(assembly_seq.keys())
    #determine the number of species present in each assembly
    print("number of species in successive assemblies")
    n_species=[]
    for t in assembly_change_t:
        n_species.append(len(assembly_seq[t]))
    for s in range(0, len(n_species)):
        print("assembly", s, "species", n_species[s])
    
    #determine the time intervals between successive assembly changes
    print("time between successive assemblies")
    for t_index in range(1, len(assembly_change_t)):
        t_between = assembly_change_t[t_index] - assembly_change_t[t_index-1]
        print("assembly", t_index-1, "time to next stable assembly", t_between)
    
    inv_t=sorted(invasion.keys())
    ext_t=sorted(extinction.keys())
    #loop over the time index at which assembly changes 
    print("invasions between successive assemblies")
    for t_index in range(1, len(assembly_change_t)):
        #now loop over the actual time points at which invasions were attempted 
        invctr=0
        for t in inv_t:
            if t>=assembly_change_t[t_index-1] and t<assembly_change_t[t_index]:
                invctr=invctr+1
        print("assembly", t_index-1, "invasions", invctr)
        
    print("extinctions between successive assemblies")
    for t_index in range(1, len(assembly_change_t)):
        #now loop over the actual time points at which invasions were attempted 
        extctr=0
        for t in ext_t:
            if t>assembly_change_t[t_index-1] and t<=assembly_change_t[t_index]:
                extctr=extctr+1
        print("assembly", t_index-1, "extinctions", extctr) 
        
    print("perturbation sensitivity of successive assemblies")
    ass_ctr=0
    for sens in sorted(pert_sens.keys()):
        print('assembly ', ass_ctr, 'sensitivity', pert_sens[sens] )
        ass_ctr += 1
    
#compute some basic statistics from the current assembly sequence
#ass_seq is the community's species composition at times t where this composition changed
def store_assembly_stats(ass_seq, inv_seq, ext_seq, env_pert_sens_seq, spec_pert_sens_seq, spec_per_ass, \
                         time_bet_ass, ext_bet_ass, inv_bet_ass,  \
                         env_pert_sens_per_ass, spec_pert_sens_per_ass):

    
    #this is an array of times at which the assembly changed
    assembly_change_t=sorted(ass_seq.keys())
    #determine the number of species present in each assembly
    n_species=[]
    for t in assembly_change_t:
        n_species.append(len(ass_seq[t]))
    #for s in range(0, len(n_species)):   
    #s is the s-th assembly in the sequence
    for s in range(0, len(n_species)):
        #if the array spec_per_ass does not yet have an entry for the s-th
        #assembly (because no assembly sequence this long has been encountered before)
        #then create it
        if len(spec_per_ass)<=s:
            spec_per_ass.append([])         
            spec_per_ass[s].append(n_species[s])
        #if the array element exists already, append to it
        else:
            spec_per_ass[s].append(n_species[s])   
    
   
    
    #determine the time intervals between successive assembly changes
    for t_index in range(0, len(assembly_change_t)-1):
        t_between = assembly_change_t[t_index+1] - assembly_change_t[t_index]
   
        #if the array spec_per_ass does not yet have an entry for the s-th
        #assembly (because no assembly sequence this long has been encountered before)
        #then create it
        if len(time_bet_ass)<=t_index:
            time_bet_ass.append([])         
            time_bet_ass[t_index].append(t_between)
        #if the array element exists already, append to it
        else:
            time_bet_ass[t_index].append(t_between)
    
    #list of times at which invasions were attempted
    inv_t=sorted(inv_seq.keys())
    #loop over the time index at which assembly changes 
    for t_index in range(0, len(assembly_change_t)-1):
        #now loop over the actual time points at which invasions were attempted 
        invctr=0
        for t in inv_t:
            if t>=assembly_change_t[t_index] and t<assembly_change_t[t_index+1]:
                invctr=invctr+1
        #then create one
        if len(inv_bet_ass)<=t_index:
            inv_bet_ass.append([])         
            inv_bet_ass[t_index].append(invctr)
        #if the array element exists already, append to it
        else:
            inv_bet_ass[t_index].append(invctr)
        
    ext_t=sorted(ext_seq.keys())
    for t_index in range(0, len(assembly_change_t)-1):
        #now loop over the actual time points at which invasions were attempted 
        extctr=0
        for t in ext_t:
            if t>=assembly_change_t[t_index] and t<assembly_change_t[t_index+1]:
                extctr=extctr+1
        if len(ext_bet_ass)<=t_index:
            ext_bet_ass.append([])         
            ext_bet_ass[t_index].append(extctr)
        #if the array element exists already, append to it
        else:
            ext_bet_ass[t_index].append(extctr)
        
    ass_ctr=0
    #the i-th element will hold the sensitivity of the i-th assembly
    env_pert_sens_arr=[]
    for sens in sorted(env_pert_sens_seq.keys()):
        env_pert_sens_arr.append(env_pert_sens_seq[sens])
        ass_ctr += 1
    #s is the s-th assembly in the sequence
    for s in range(0, len(n_species)):
        #if the array env_pert_sens_per_ass does not yet have an entry for the s-th
        #assembly (because no assembly sequence this long has been encountered before)
        #then create it
        if len(env_pert_sens_per_ass)<=s:
            env_pert_sens_per_ass.append([])         
            env_pert_sens_per_ass[s].append(env_pert_sens_arr[s])
        #if the array element exists already, append to it
        else:
            env_pert_sens_per_ass[s].append(env_pert_sens_arr[s])   
   
    #now the same for biomass/species perturbations
    ass_ctr=0
    #the i-th element will hold the sensitivity of the i-th assembly
    spec_pert_sens_arr=[]
    for sens in sorted(spec_pert_sens_seq.keys()):
        spec_pert_sens_arr.append(spec_pert_sens_seq[sens])
        ass_ctr += 1
    #s is the s-th assembly in the sequence
    for s in range(0, len(n_species)):
        #if the array env_pert_sens_per_ass does not yet have an entry for the s-th
        #assembly (because no assembly sequence this long has been encountered before)
        #then create it
        if len(spec_pert_sens_per_ass)<=s:
            spec_pert_sens_per_ass.append([])         
            spec_pert_sens_per_ass[s].append(spec_pert_sens_arr[s])
        #if the array element exists already, append to it
        else:
            spec_pert_sens_per_ass[s].append(spec_pert_sens_arr[s])

#a class that will hold a community of models, where each model has a nonzero biomass
#and a medium in which the community exists (but may also apply to other objects)
#spec and biom are dictionaries keyed with model names whose values are the models and the biomass
#med is a dictionary keyed with (external) metabolite and holds its concentration
            
#designed for CCM models but works for others too
class CCMcomm:
    """
    def __init__(self, spec, biom, med):
        self.species={}
        self.biomass={}
        self.medium={}
    """
    
    def __init__(self, spec, biom, med):
        self.species={}
        self.biomass={}
        for specid in spec.keys():
            if biom[specid]>0:
                self.biomass[specid]=biom[specid]
                #uses the cobrapy deepcopy copy command to avoid problems with passing by reference
                self.species[specid]=spec[specid].copy()
        #now create a medium that contains only the metabolites that are present in the passed medium
        self.medium={}
        for m in med.keys():
            if med[m]>0:
                self.medium[m]=med[m]

    def commprint(self):
        print("species biomass")
        for specid in self.species.keys():
            print(specid, self.biomass[specid])
        print("metabolite concentration")
        for m in self.medium.keys():
            print(m, self.medium[m])

    
#a highly inefficient deepcopy routine for metabolic models, exists because
#both python's deep copy and cobrapy's deepcopy yield unpredictable results when 
#one chooses randomly reactions from the copy, i.e., even with the random number
#seed set, random.choice does not give predictable results. 
#I aim to alleviate this by building
#the copy from a sorted list of metabolite/reaction ids
#should be used only rarely because of the likely time overhead
def metmodel_deepcopy_aw(oldmod):
    
    #note that if we give the old model id as an argument here, the whole
    #old model will apparently be copied, i.e., it appears that the new
    #model will become a reference to an old model.
    newmod=Model()
    newmod.id=oldmod.id  
    newmod.solver=oldmod.solver
    
    midarr=[]
    for m in oldmod.metabolites:
        midarr.append(m.id)
    newmlist=[]
    #sort metabolite ids, in the hope that it will make their order and thus random choices
    #from the copied network more reproducible
    for mid in sorted(midarr):
        mtmp=oldmod.metabolites.get_by_id(mid)
        newmtmp=Metabolite(mtmp.id,formula=mtmp.formula,name=mtmp.name,compartment=mtmp.compartment)           
        newmlist.append(newmtmp) 
    newmod.add_metabolites(newmlist)
    #print("now done copying the metabolites")
    #now do the same thing with the reactions
    ridarr=[] 
    for r in oldmod.reactions:
        ridarr.append(r.id)
    for rid in sorted(ridarr):
        rtmp=oldmod.reactions.get_by_id(rid)
        rnewtmp=Reaction(rtmp.id)
        rnewtmp.name = rtmp.name
        rnewtmp.subsystem = rtmp.subsystem
        #note that the reaction has to be added to the model before one
        #can set the stoichiometry, or the metabolites of the model will not have been defined
        newmod.add_reactions([rnewtmp]) 
        #also, it is important to set the bounds last, or the fba solution
        #will differ between the old and the new model
        rnewtmp.reaction = rtmp.reaction  
        #also, it seems rather important to set the bounds last, or the fba solution
        #will differ between the old and the new model
        rnewtmp.lower_bound = rtmp.lower_bound 
        rnewtmp.upper_bound = rtmp.upper_bound
        
    newmod.objective = oldmod.objective
    
    return newmod

 

    