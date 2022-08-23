"""
MODULE VERSION:
Combining peptide and membrane structures into a single system.
"""

import os
import numpy as np

import core.settings as CoreSettings
import core.system as System
import core.gromacs as Gromacs
import core.gromacs_commands

import user.usersettings as UserSettings

def module_insert_peptide_2mem(simulation):
    """
    Takes a peptide GRO file and a membrane GRO file, combines the two into a new GRO file and runs a steepest descent energy minimization.
    """
    for membrane_type in UserSettings.membrane_types:
        #print('Inserting peptide into system: %s \n' %membrane_type)
        class LocalFiles(): # Filenames used within this module
            if membrane_type == 'high_tension':
                membrane_GRO                = UserSettings.filename_membrane_high_tension_GRO
            if membrane_type == 'low_tension':
                membrane_GRO                = UserSettings.filename_membrane_low_tension_GRO
            membrane_TOP                = UserSettings.filename_membrane_TOP

            peptide_GRO                 = simulation._share.get_filename("gen-GRO")
            peptide_TOP                 = simulation._share.get_filename("gen-TOP")

            peptide_GRO_tmp             = "peptide_tmp.gro"
            # peptide_TOP_restr           = simulation._share.get_filename("gen-TOP_restr")
            # peptide_TOP_noConstraints   = simulation._share.get_filename("gen-TOP_noConstraints")
            # peptide_TOP_noConstr_Restr  = simulation._share.get_filename("gen-TOP_noConstr_Restr")

            # TOP_restr                   = UserSettings.name_output_system + "_" + membrane_type + ".restr.top"
            # TOP_noConstraints           = UserSettings.name_output_system + "_" + membrane_type + ".noConstraints.top"
            # TOP_noConstr_Restr          = UserSettings.name_output_system + "_" + membrane_type + ".noConstr_Restr.top"

            TOP                         = UserSettings.name_output_system + "_" + membrane_type + ".top"
            TOPOUT                      = "insertion_out_" + membrane_type + ".top"

            GRO                         = UserSettings.name_output_system + "_" + membrane_type + ".gro"
            NDX                         = UserSettings.name_output_system + "_" + membrane_type + ".ndx"
            TPR                         = UserSettings.name_output_system + "_" + membrane_type + ".tpr"
            MDOUT                       = UserSettings.name_output_system + "_" + membrane_type + ".mdout.mdp"
            TRR                         = UserSettings.name_output_system + "_" + membrane_type + ".trr"
            EDR                         = UserSettings.name_output_system + "_" + membrane_type + ".edr"
            LOG                         = UserSettings.name_output_system + "_" + membrane_type + ".log"
            CPT                         = UserSettings.name_output_system + "_" + membrane_type + ".cpt"

            MDP = UserSettings.filename_steep_system_MDP

            stdout_stderr = None
            if CoreSettings.write_stdout_stderr_to_file:
                stdout_stderr = os.path.join(simulation._output_dir, CoreSettings.filename_output_stdout_stderr)

        def path(filename):
            if filename is None:
                return None
            return os.path.join(simulation._temp_dir, filename)

        def center_peptide_to_membrane(num=0):
            box_vector = Gromacs.GRO_read_box(path(LocalFiles.membrane_GRO))
            data_peptide_GRO = System.read_text_file(path(LocalFiles.peptide_GRO))

            for i in reversed(range(0, len(data_peptide_GRO))):
                line = data_peptide_GRO[i].strip()
                if not len(line) == 0:
                    data_peptide_GRO[i] = " ".join(box_vector)
                    break

            #os.remove(path(LocalFiles.peptide_GRO_tmp))
            System.write_text_file(path(LocalFiles.peptide_GRO_tmp), data_peptide_GRO, add_newline=False)

            #if num == 0:
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.EDITCONF_align_X_center(
                    path(LocalFiles.peptide_GRO_tmp),
                    path(LocalFiles.peptide_GRO_tmp)
                ),
                input="0\n".encode(),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )
                # rotate_vector = UserSettings.rotate_vector

                # System.run_command(
                #     core.gromacs_commands.GROMACSRunCommands.EDITCONF_rotate( # ROTATE
                #         path(LocalFiles.peptide_GRO),
                #         path(LocalFiles.peptide_GRO),
                #         rotate_vector
                #     ),
                #     path_stdout=path(LocalFiles.stdout_stderr),
                #     path_stderr=path(LocalFiles.stdout_stderr)
                # )
            ## two peptides, one on each side of the membrane is planned to be added

            # System.run_command(
            #     core.gromacs_commands.GROMACSRunCommands.EDITCONF_center( # CENTER
            #         path(LocalFiles.peptide_GRO_tmp),
            #         path(LocalFiles.peptide_GRO_tmp)
            #     ),
            #     path_stdout=path(LocalFiles.stdout_stderr),
            #     path_stderr=path(LocalFiles.stdout_stderr)
            # )
            
            translate_vector = (0, 0, -1*UserSettings.b[UserSettings.membrane_types.index(membrane_type)]) # thickness = 3.977 (1.82942 = 92% of monolayer); thickness = 3.540 (1.6284 = 92% of monolayer)

            
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.EDITCONF_on_top_of_mem( # TRANSLATE
                    path(LocalFiles.peptide_GRO_tmp),
                    path(LocalFiles.peptide_GRO_tmp),
                    translate_vector
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

        def insert_into_membrane(num=0):
            base_sections = [Gromacs.TOP_HEADER_SECTION, "system", "molecules"]
            append_sections = [Gromacs.TOP_HEADER_SECTION, "molecules"]

            Gromacs.GRO_merge(
                path(LocalFiles.membrane_GRO),
                path(LocalFiles.peptide_GRO_tmp),
                path(LocalFiles.GRO)
            )

            Gromacs.TOP_merge(
                path(LocalFiles.membrane_TOP),
                path(LocalFiles.peptide_TOP),
                path(LocalFiles.TOP),
                base_sections, append_sections
            )

            # Gromacs.TOP_merge(
            #     path(LocalFiles.membrane_TOP),
            #     path(LocalFiles.peptide_TOP_restr),
            #     path(LocalFiles.TOP_restr),
            #     base_sections, append_sections
            # )
            # Gromacs.TOP_merge(
            #     path(LocalFiles.membrane_TOP),
            #     path(LocalFiles.peptide_TOP_noConstraints),
            #     path(LocalFiles.TOP_noConstraints),
            #     base_sections, append_sections
            # )
            # Gromacs.TOP_merge(
            #     path(LocalFiles.membrane_TOP),
            #     path(LocalFiles.peptide_TOP_noConstr_Restr),
            #     path(LocalFiles.TOP_noConstr_Restr),
            #     base_sections, append_sections
            # )
            
            Gromacs.NDX_from_GRO(
                path(LocalFiles.GRO),
                path(LocalFiles.NDX),
                UserSettings.resnm_to_indexgroup_dict
            )

        def steep_system():
            System.replace_in_file(
                path(LocalFiles.MDP),
                UserSettings.tags_steep_system_MDP[0],
                "-I" + os.path.realpath(simulation._temp_dir)
            )

            # Actual simulation #
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX_TOPOUT(
                    path(LocalFiles.MDP),
                    path(LocalFiles.GRO),
                    path(LocalFiles.TOP),
                    path(LocalFiles.TPR),
                    path(LocalFiles.MDOUT),
                    path(LocalFiles.NDX),
                    path(LocalFiles.TOPOUT)
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            System.run_simulation_command(
                core.gromacs_commands.GROMACSRunCommands.MDRUN(
                    path(LocalFiles.TPR),
                    path(LocalFiles.TRR),
                    path(LocalFiles.GRO),
                    path(LocalFiles.EDR),
                    path(LocalFiles.LOG),
                    path(LocalFiles.CPT),
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

        System.copy_dir_files(UserSettings.path_basedir_insert_peptide, path(""))

        center_peptide_to_membrane()
        insert_into_membrane()
        steep_system()

        simulation._share.add_file("ins-%s-peptide_tmp-GRO" %membrane_type,       LocalFiles.peptide_GRO_tmp,     UserSettings.FILE_WRITE_DEBUGGING)
        simulation._share.add_file("ins-%s-GRO" %membrane_type,                   LocalFiles.GRO,                 False)
        simulation._share.add_file("ins-%s-TOP" %membrane_type,                   LocalFiles.TOP,                 True)
        simulation._share.add_file("ins-%s-TOPOUT" %membrane_type,                LocalFiles.TOPOUT,              UserSettings.FILE_WRITE_DEBUGGING)
        # simulation._share.add_file("ins-%s-TOP_restr" %membrane_type,             LocalFiles.TOP_restr,           True)
        # simulation._share.add_file("ins-%s-TOP_noConstraints" %membrane_type,     LocalFiles.TOP_noConstraints,   True)
        # simulation._share.add_file("ins-%s-TOP_noConstr_Restr" %membrane_type,    LocalFiles.TOP_noConstr_Restr,  True)
        simulation._share.add_file("ins-%s-NDX" %membrane_type,                   LocalFiles.NDX,                 True)
