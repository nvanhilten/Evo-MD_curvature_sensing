import os
import sys
import core.system as System
import getopt
import subprocess
import numpy as np

##############
## Settings ##
##############

# PDB generator #
pdb_gen                             = "PeptideBuilder" # "Modeller"

# Force field #
ff                                  = "martini3001"

# Input #
filename_seq                        = "input.seq"
filename_ssd                        = "input.ssd"

# Output #
name_output_score                   = "fitness"
name_output_peptide                 = "peptide"
name_output_system                  = "system"
name_output_production              = "production"

# Simulation #
timeout_generate_peptide_MDRUN      = 60*2 # In seconds
timeout_Production_relax_MDRUN      = 60*2 # In seconds
timeout_Production_MDRUN            = 60*60*2 # In seconds

# Basefile #
path_basedir_generate_peptide       = System.get_file_dir(__file__) + "/base/generate_peptide/"
path_basedir_insert_peptide         = System.get_file_dir(__file__) + "/base/insert_peptide/"
path_basedir_production             = System.get_file_dir(__file__) + "/base/production/"

# generate peptide
filename_steep_peptide_TOP          = "steep_peptide.top"
if pdb_gen == "Modeller":    
    template_pdb                        = "template.pdb"

if ff == "martini22":
    template_top                        = "template_m22.top"
elif ff == "martini3001":
    template_top                        = "template_m3001.top"
else:
    print("Force field selection not supported")
#filename_NVT_peptide_MDP            = "nvt_peptide.mdp"
filename_MDP_steep_pept             = "steep_pept.mdp"
filename_MDP_eq_pept                = "eq_small_dt_pept.mdp"

# insert peptide
membrane_types                      = ['high_tension', 'low_tension']
filename_membrane_high_tension_GRO  = "membrane_high_tension.gro"
filename_membrane_low_tension_GRO   = "membrane_low_tension.gro"
filename_membrane_TOP               = "membrane.top"
filename_steep_system_MDP           = "steep_system.mdp"

# production
filename_steep_softcore_MDP         = "steep_softcore_system.mdp"
filename_steep_restr_MDP            = "steep_restr_system.mdp"
filename_NPT_equilibrate_system_MDP = "equilibrate_low_timestep.mdp"
filename_NPT_system_MDP             = "npt.mdp"

# tags
tags_steep_peptide_TOP              = ("$PEPTIDE_ITP",)
tags_steep_peptide_2pept_TOP        = ("$PROTEIN_NUM",)
tags_NVT_peptide_MDP                = ("$INCLUDE_CWD",)
tags_steep_system_MDP               = ("$INCLUDE_CWD",)
tags_NPT_equilibrate_system_MDP     = ("$INCLUDE_CWD",)
tags_NPT_system_MDP                 = ("$INCLUDE_CWD",)
tags_steep_softcore_MDP             = ("$INCLUDE_CWD",)
tags_steep_restr_MDP                = ("$INCLUDE_CWD",)
tags_SIM_dt                         = ("$SIM_DT",)
tags_SIM_nsteps                     = ("$SIM_NSTEPS",)


# Miscellaneous: generate_peptide.py #
STEEP_constraint_forceConstant      = "1000000" #"5000"
# append_ITP_noConstraints            = ".noConstraints.itp"
# append_ITP_seq2itp                  = ".seq2itp.itp"
# append_ITP_noConstr_Restr           = ".noConstr_Restr.itp"
# append_TOP_noConstraints            = ".noConstraints.top"
# append_TOP_noConstr_Restr           = ".noConstr_Restr.top"

# Miscellaneous: insert_peptide.py #
resnm_to_indexgroup_dict            = {
                                        "POPC"  : ("POPC", "Lipids"),
                                        "POPE"  : ("POPE", "Lipids"),
                                        "POPS"  : ("POPS", "Lipids"),
                                        "DPSM"  : ("DPSM", "Lipids"),
                                        "CHOL"  : ("CHOL", "Lipids"),
                                        "W"     : "Solvent", "NA"  : "Solvent", "CL"  : "Solvent", "ION"  : "Solvent",
                                        "GLY"   : "Protein", "ALA"  : "Protein", "ASP"  : "Protein", "ASN"  : "Protein", "GLU"  : "Protein",
                                        "GLN"   : "Protein", "VAL"  : "Protein", "LEU"  : "Protein", "ILE"  : "Protein", "MET"  : "Protein",
                                        "THR"   : "Protein", "SER"  : "Protein", "CYS"  : "Protein", "LYS"  : "Protein", "ARG"  : "Protein",
                                        "HIS"   : "Protein", "PHE"  : "Protein", "PRO"  : "Protein", "TRP"  : "Protein", "TYR"  : "Protein"      
                                      }

# Miscellaneous: production.py #
name_positive_ion_GENION            = "NA"
name_negative_ion_GENION            = "CL"
name_solvent_GENION                 = "W"
neutralize                          = True
WRITE_XTC                           = False #True

# Miscellaneous
FILE_WRITE_DEBUGGING                = False

# # Miscellaneous: seq2itp_aminoacids.py #
# generate_bond_length                = 0.350
# generate_sidechain_angle            = 90
# generate_box_factor                 = 50

# # Dependency: Seq2ITP #
# path_Seq2ITP                        = System.get_file_dir(__file__) + "/seq2itp/seq2itp.pl"
# #path_Seq2ITP_2pept                  = System.get_file_dir(__file__) + "/seq2itp/seq2itp_2pept.pl"

# Dependency: Porter5 #
path_Porter5                        = System.get_file_dir(__file__) + "/Porter5/Porter5.py"

# Dependency: modeller #
if pdb_gen == "Modeller":
    path_modeller                        = System.get_file_dir(__file__) + "/modeller/run_modeller.py"

# fitness
SIM_dt                    = 0.03 # ps
SIM_time                  = 400000 # ps
ENERGY_init_time          = 0.95*SIM_time

# TMD fitness correction factor
a = 0.0
b = [1.65646, 1.82942] # thickness values: 3.601, 3.977 (b_value = 0.5*0.92*thickness) [high_tension, low_tension]

# other stuff
selection_sorting                   = True # goes in SelectionMethods, should be True if fitness is maximized, False if minimized (i.e. in case of free energy difference as fitness)             
SSP                                 = False
normalize_by_nr_of_beads            = False # if True, fitness values are devided by the number of beads in the peptide

rep_lengths                         = [24] # possible block lengths
repeats                             = [1] # possible repeat numbers
min_length                          = 24
max_length                          = 24 # maximal sequence length
charge_range                        = [-24,24] # range of allowed charges
repeat_mutations                    = [-1,1] # -1: decrease/increase number of repeats. Set to [0] to switch off.
pept_position                       = [1/3] # position on the membrane (*box_length from the center)
rotate_vector                       = (0,0,90)
equilibration_time                  = 500 # ps
blocktime                           = 500 # ps
bootstrapping_sampling              = 10000
ENERGY_select_groups_Protein1       = ["Coul-SR:Protein1-Lipids", "LJ-SR:Protein1-Lipids", "LJ-SR:Protein1-Solvent"]
ENERGY_select_groups_Protein2       = ["Coul-SR:Protein2-Lipids", "LJ-SR:Protein2-Lipids", "LJ-SR:Protein2-Solvent"]
#Area_free_low_tension               = 42.3853082 # nm2; x*y of /user/base/insert_peptide/membrane_low_tension.gro
Area_free_low_tension               = 6.5104**2 # nm2; x*y of /user/base/insert_peptide/membrane_low_tension.gro
Area_free_high_tension              = 7.02833**2 # nm2; x*y of /user/base/insert_peptide/membrane_high_tension.gro
#SurfTen_free_low_tension            = [2.75821, 1.2] # avg, err; measured on /user/base/production/free_low_tension.edr
SurfTen_free_low_tension            = [1.92718, 0.99] # avg, err; measured on /user/base/production/free_low_tension.edr
SurfTen_free_high_tension           = [314.934, 0.59] # avg, err; measured on /user/base/production/free_high_tension.edr
conversion_factor                   = 6.02214076e-2 # bar*nm to kJ*mol-1*nm-2
plotting                            = False
propensities_initial                = True # false = all AA equal probability (0.05), true = take propensities from ALPS-like peptides (Drin et al 2007)
propensities_mutation               = True # false = all AA equal probability (0.05), true = take propensities from ALPS-like peptides (Drin et al 2007)
propensities                        = {"A": 0.1,
                                       "C": 0.0,
                                       "D": 0.0,
                                       "E": 0.1,
                                       "F": 0.1,
                                       "G": 0.0,
                                       "H": 0.0,
                                       "I": 0.0,
                                       "K": 0.1,
                                       "L": 0.1,
                                       "M": 0.1,
                                       "N": 0.0,
                                       "P": 0.0,
                                       "Q": 0.1,
                                       "R": 0.0,
                                       "S": 0.1,
                                       "T": 0.0,
                                       "V": 0.0,
                                       "W": 0.1,
                                       "Y": 0.1}

# propensities_ALPS                   = {"A": 0.063,
#                                        "C": 0.000,
#                                        "D": 0.025,
#                                        "E": 0.024,
#                                        "F": 0.040,
#                                        "G": 0.132,
#                                        "H": 0.012,
#                                        "I": 0.053,
#                                        "K": 0.037,
#                                        "L": 0.111,
#                                        "M": 0.021,
#                                        "N": 0.036,
#                                        "P": 0.018,
#                                        "Q": 0.033,
#                                        "R": 0.026,
#                                        "S": 0.173,
#                                        "T": 0.101,
#                                        "V": 0.059,
#                                        "W": 0.014,
#                                        "Y": 0.022}

###############
## Functions ##
###############

# # Seq2ITP #
# def command_seq2ITP(path_sequence, path_structure, path_ITP_out):
#     command = "perl "
#     command += path_Seq2ITP
#     command += " --sequence " + path_sequence
#     command += " --structure " + path_structure
#     command += " --itp " + path_ITP_out
#     command += " --noelastic"
#     return command

# Porter5 #
def command_Porter5(path_sequence):
    command = "python3 "
    command += path_Porter5
    command += " -i " + path_sequence
    command += " --cpu 1 "
    command += " --fast"
    return command

# modeller #
def command_modeller(path_sequence, path_structure, path_PDB, name, ali, outdir):
    command = "python "
    command += path_modeller + " "
    command += path_sequence + " "
    command += path_structure + " "
    command += path_PDB + " "
    command += name + " "
    command += ali + " "
    command += outdir
    return command


def command_martinize2(path_PDB_in, structure, path_PDB_out, path_TOP, ff=ff, elastic=False):
    command = "/home/nichilte/software/Python-3.6.9/Python/bin/martinize2"
    command += " -f " + path_PDB_in
    command += " -x " + path_PDB_out
    command += " -o " + path_TOP
    command += " -ss " + structure
    command += " -ff " + ff
    if elastic:
        command += " -elastic"
    return command

def command_whereamI():
    command = "pwd"
    return command

# def command_seq2ITP_2pept(path_sequence, path_structure, path_ITP_out, molname):
#     command = "perl "
#     command += path_Seq2ITP_2pept
#     command += " --sequence " + path_sequence
#     command += " --structure " + path_structure
#     command += " --itp " + path_ITP_out
#     command += " --noelastic"
#     command += " --molname " + molname
#     return command

# Energy Analysis #
def get_ENERGY_select_groups():
    result = []
    result.extend(ENERGY_select_groups)
    return result

def old_correction_factor(x, a, b):
    if x < a or x > 2*(b-a):
            factor = 0.0
    else:
            factor = -0.5 * (np.cos((x-a)*np.pi/(b-a)) - 1)
    return factor

def correction_factor(x, b):
    if x < b/2 or x > 3*b/2:
            factor = 0.0
    else:
            factor = -0.5 * (np.cos((x-b/2)*2*np.pi/b) - 1)
    return factor
