This repository contains scripts needed to assemble microbial communities and analyze their composition, as described in the manuscript "Competition for nutrients increases invasion resistance during assembly of microbial communities" (A. Wagner).

Files include

1. metfuncs_aw_pub.py:  a file with multiple python functions needed for the simulations

2. assembly_cstat_agora_cluscode_pub.py: main simulation routine for agora gut microbial communities. invoke as assembly_cstat_agora_cluscode.py [jobid], where jobid is a positive integer that will determine the random number seed to be used, intended to run on a cluster where multiple jobs can be run simultaneously.

3. asspar_agora_clus_pub.py: parameter file for main routine, also contains paths to data files loaded into the main routine. These comprise the following three files

4. ./metmodels/Agora1_03_without_mucins/agora_modelsample_spec_100_08-09-2021_out.txt: This file contains a list of 100 randomly sampled agora models as represented by their mat file names, together with some basic statistics on these models.

5. ./metmodels/Agora1_03_without_mucins/Magnusdottir_2017_TableS12_diets.csv: This file contains the definition of a Westermn diet, and except for formatting changes, is identical to Table S12 in Magnusdottir et al. (2017).

6. The agora models are loaded as .mat files and should be in the directory ./metmodels/Agora1_03_without_mucins/mat/. They are not provided here, because A. Wagner does not own them. They can be obtained from the Virtual Metabolic Human data base at https://www.vmh.life/#microbes/search as  ‘AGORA 1.03 without mucin’ (see also Noronha et al., Nucleic Acids Research 2019)

7. assembly_cstat_ranvia_cluscode.py: main simulation routine for communities of random viable metabolisms. Call as 'assembly_cstat_ranvia_cluscode_pub.py [jobid]', where jobid is a positive integer that will determine the random number seed to be used, intended to run on a cluster where multiple jobs can be run simultaneously.

8. asspar_ranvia_clus_pub.py: parameter file for main simulation routine for random viable networks, also contains paths to data files to be loaded into the main simulation routine for random viable networks. These comprise the following files

9. metmodels/ranviable/minglc_r871/ranvia_stats_minglc_871_04-11-2021: A list of file names for 100 random viable networks to be loaded into the simulation routine, together with some basic statistics on these models.

10. metmodels/ranviable/minglc_r871/S*.xml: 100 xml files, each containing one of the random viable metabolisms to be loaded into the main simulation routine.

9. 'kegg_compounds_magda.txt':  A map of KEGG comppound identifiers to common names for multiple types of molecules from San Roman et al. PLoS Computational Biology 2018. 



