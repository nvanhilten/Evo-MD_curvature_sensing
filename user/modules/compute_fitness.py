"""
Modules that handle computation of fitness. This file can contain multiple modules.
"""

import os
import numpy as np

import core.settings as CoreSettings
import core.system as System
import core.gromacs as Gromacs
import core.gromacs_commands

import user.usersettings as UserSettings


def module_compute_fitness_2mem_dict(simulation):

    class LocalFiles(): # Filenames used within this module
        prod_GRO_dict = {}
        prod_XTC_dict = {}
        prod_TPR_dict = {}
        prod_EDR_dict = {}
        surften_XVG_dict = {}
        distance_XVG_dict = {}
        ENERGY_dict = {}
        DIST_dict = {}
        fitness_dict = {}
        
        for membrane_type in UserSettings.membrane_types:
            prod_GRO_dict[membrane_type] = simulation._share.get_filename("prod-%s-GRO" %membrane_type)
            prod_XTC_dict[membrane_type] = simulation._share.get_filename("prod-%s-XTC" %membrane_type)
            prod_TPR_dict[membrane_type] = simulation._share.get_filename("prod-%s-TPR" %membrane_type)
            prod_EDR_dict[membrane_type] = simulation._share.get_filename("prod-%s-EDR" %membrane_type)
            surften_XVG_dict[membrane_type] = UserSettings.name_output_production + "_" + membrane_type + "_surften.xvg"
            distance_XVG_dict[membrane_type] = UserSettings.name_output_production + "_" + membrane_type + "_distance.xvg"
            ENERGY_dict[membrane_type]   = UserSettings.name_output_production + "_" + membrane_type + ".ENERGY"
            DIST_dict[membrane_type]     = UserSettings.name_output_production + "_" + membrane_type + ".DIST"

        NDX = simulation._share.get_filename("prod-NDX")
        fitness = UserSettings.name_output_score 

        if UserSettings.normalize_by_nr_of_beads:
            sequence = UserSettings.filename_seq
            AA_dict = {"A": 2, "C": 2, "D": 2, "E": 2, "F": 4, "G": 1, "H": 4, "I": 2, "K": 3, "L": 2, "M": 2, "N": 2, "P": 2, "Q": 2, "R": 3, "S": 2, "T": 2, "V": 2, "W": 6, "Y": 5}

        stdout_stderr = None
        if CoreSettings.write_stdout_stderr_to_file:
            stdout_stderr = os.path.join(simulation._output_dir, CoreSettings.filename_output_stdout_stderr)

    def path(filename):
        if filename is None:
            return None
        return os.path.join(simulation._temp_dir, filename)

    def run_ENERGY():
        for membrane_type in UserSettings.membrane_types:
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.ENERGY(
                    path(LocalFiles.prod_EDR_dict[membrane_type]),
                    path(LocalFiles.surften_XVG_dict[membrane_type]),
                    -1,
                ),
                input="#Surf*SurfTen\n".encode(),
                path_stdout=path(LocalFiles.ENERGY_dict[membrane_type]),
                path_stderr=path(LocalFiles.ENERGY_dict[membrane_type])
            )

    def compute_ddF():
        surften_dict = {}
        surften_err_dict = {}
        for membrane_type in UserSettings.membrane_types:
            average, err = Gromacs.GMX_average_ENERGY_output(path(LocalFiles.ENERGY_dict[membrane_type]))
            surften_dict[membrane_type] = average*UserSettings.conversion_factor # convert to kJ mol-1 nm-2
            surften_err_dict[membrane_type] = err*UserSettings.conversion_factor # convert to kJ mol-1 nm-2

        dA = UserSettings.Area_free_high_tension - UserSettings.Area_free_low_tension
        diff_low_tension = UserSettings.conversion_factor*UserSettings.SurfTen_free_low_tension[0] - surften_dict[UserSettings.membrane_types[1]]
        diff_high_tension = UserSettings.conversion_factor*UserSettings.SurfTen_free_high_tension[0] - surften_dict[UserSettings.membrane_types[0]]
        ddF = 0.5 * dA * (diff_low_tension + diff_high_tension)
        ddF_err = np.sqrt(( (0.5*dA)**2 * (UserSettings.SurfTen_free_high_tension[1]*UserSettings.conversion_factor)**2 + surften_err_dict[UserSettings.membrane_types[0]]**2 + (UserSettings.SurfTen_free_low_tension[1]*UserSettings.conversion_factor)**2 + surften_err_dict[UserSettings.membrane_types[1]]**2))
        return ddF, ddF_err

    def run_DISTANCE():
        for membrane_type in UserSettings.membrane_types:
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.DISTANCE(
                    path(LocalFiles.prod_XTC_dict[membrane_type]),
                    path(LocalFiles.prod_TPR_dict[membrane_type]),
                    path(LocalFiles.NDX),
                    path(LocalFiles.distance_XVG_dict[membrane_type]),
                    UserSettings.ENERGY_init_time
                ),
                input="com of group Protein plus com of group Lipids\n".encode(),
                path_stdout=path(LocalFiles.DIST_dict[membrane_type]),
                path_stderr=path(LocalFiles.DIST_dict[membrane_type])
            )

    def compute_correction_factor():
        factor_list = []
        for membrane_type in UserSettings.membrane_types:
            z_list = []
            for line in open(path(LocalFiles.distance_XVG_dict[membrane_type]), "r"):
                if not "@" in line and not "#" in line:
                    z_list.append(float(line.split()[-1]))
            z_av = np.mean(z_list)
            factor = UserSettings.old_correction_factor(z_av, UserSettings.a, UserSettings.b[UserSettings.membrane_types.index(membrane_type)])
            #factor = UserSettings.correction_factor(z_av, UserSettings.b[UserSettings.membrane_types.index(membrane_type)])
            factor_list.append(factor) # apply fuction
            print(z_av, factor)
        return min(factor_list) # return the lowest factor

    def count_beads(seq, AA_dict):
        count = 0
        print(seq)
        for AA in seq:
            count+=AA_dict[AA]
        return count
    
    run_ENERGY()
    ddF, ddF_err = compute_ddF()
    run_DISTANCE()
    if ddF > 0.0:
        factor = compute_correction_factor()
    else:
        factor = 1.0
    fitness = factor*ddF

    if UserSettings.normalize_by_nr_of_beads:
        sequence = System.read_text_file(path(LocalFiles.sequence), strip=True)[1]
        nr_of_beads = count_beads(sequence, LocalFiles.AA_dict)
        fitness = fitness/nr_of_beads


    print(str(ddF)+" +/- "+str(ddF_err)+"\tpenalty="+str(factor)+"\tFitness="+str(fitness)+"\n")

    System.write_text_file(path(LocalFiles.fitness), [str(ddF)+" +/- "+str(ddF_err)+"\tpenalty="+str(factor)+"\tFitness="+str(fitness)])
    
    simulation._output_variables["ddF"] = ddF
    simulation._output_variables["fitness"] = fitness

    simulation._share.add_file("fitness", LocalFiles.fitness, True)
    for membrane_type in UserSettings.membrane_types:
        simulation._share.add_file("%s-ENERGY" %membrane_type, LocalFiles.ENERGY_dict[membrane_type], True)
        simulation._share.add_file("%s-DIST" %membrane_type, LocalFiles.DIST_dict[membrane_type], False)
        simulation._share.add_file("surften-%s-XVG" %membrane_type, LocalFiles.surften_XVG_dict[membrane_type], False)
        simulation._share.add_file("distance-%s-XVG" %membrane_type, LocalFiles.distance_XVG_dict[membrane_type], True)

# KAI END
