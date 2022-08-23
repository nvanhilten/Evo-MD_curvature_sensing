"""
MODULE VERSION:
Module that handles production simulations and energy calculation LJ+coul CHOL_Prot/(CHOL_Prot+POPC_Prot)
"""

import os

import core.settings as CoreSettings
import core.system      as System
import core.gromacs     as Gromacs
import core.gromacs_commands

import user.usersettings as UserSettings


def module_production_2mem(simulation):
    for membrane_type in UserSettings.membrane_types:
        class LocalFiles(): # Filenames used within this module
            MDP_steep_softcore          = UserSettings.filename_steep_softcore_MDP
            MDP_steep_restr             = UserSettings.filename_steep_restr_MDP
            MDP_NPT_equil               = UserSettings.filename_NPT_equilibrate_system_MDP
            MDP_NPT                     = UserSettings.filename_NPT_system_MDP

            system_GRO                  = simulation._share.get_filename("ins-%s-GRO" %membrane_type)
            system_TOP                  = simulation._share.get_filename("ins-%s-TOP" %membrane_type)

            GRO_neutralized = system_GRO + ".neutralized.gro"

            GRO                         = UserSettings.name_output_production + "_" + membrane_type + ".gro"
            GRO_softcore                = UserSettings.name_output_production + "_" + membrane_type + ".softcore.gro"
            GRO_steep                   = UserSettings.name_output_production + "_" + membrane_type + ".steep.gro"
            GRO_equil                   = UserSettings.name_output_production + "_" + membrane_type + ".equil.gro"
            NDX                         = UserSettings.name_output_production + "_" + membrane_type + ".ndx"
            TPR                         = UserSettings.name_output_production + "_" + membrane_type + ".tpr"
            MDOUT                       = UserSettings.name_output_production + "_" + membrane_type + ".mdout.mdp"
            TRR                         = UserSettings.name_output_production + "_" + membrane_type + ".trr"
            EDR                         = UserSettings.name_output_production + "_" + membrane_type + ".edr"
            LOG                         = UserSettings.name_output_production + "_" + membrane_type + ".log"
            CPT                         = UserSettings.name_output_production + "_" + membrane_type + ".cpt"
            XTC                         = UserSettings.name_output_production + "_" + membrane_type + ".xtc"


            stdout_stderr = None
            if CoreSettings.write_stdout_stderr_to_file:
                stdout_stderr = os.path.join(simulation._output_dir, CoreSettings.filename_output_stdout_stderr)

        def path(filename):
            if filename is None:
                return None
            return os.path.join(simulation._temp_dir, filename)

        def setup():
            nsteps = int(round(UserSettings.SIM_time/UserSettings.SIM_dt))

            # steep softcore
            System.replace_in_file(
                path(LocalFiles.MDP_steep_softcore),
                UserSettings.tags_steep_softcore_MDP[0],
                "-I" + os.path.realpath(simulation._temp_dir)
            )

            # steep restr
            System.replace_in_file(
                path(LocalFiles.MDP_steep_restr),
                UserSettings.tags_steep_restr_MDP[0],
                "-I" + os.path.realpath(simulation._temp_dir)
            )

            # equilibrate
            System.replace_in_file(
                path(LocalFiles.MDP_NPT_equil),
                UserSettings.tags_NPT_equilibrate_system_MDP[0],
                "-I" + os.path.realpath(simulation._temp_dir)
            )

            # npt
            System.replace_in_file(
                path(LocalFiles.MDP_NPT),
                UserSettings.tags_NPT_system_MDP[0],
                "-I" + os.path.realpath(simulation._temp_dir)
            )

            System.replace_in_file(
                path(LocalFiles.MDP_NPT),
                UserSettings.tags_SIM_dt[0],
                str(UserSettings.SIM_dt)
            )

            System.replace_in_file(
                path(LocalFiles.MDP_NPT),
                UserSettings.tags_SIM_nsteps[0],
                str(nsteps)
            )

        def neutralize_system():
            Gromacs.GRO_TOP_neutralize_system(
                path(LocalFiles.system_GRO),
                path(LocalFiles.system_TOP),
                path(LocalFiles.GRO_neutralized),
                path(LocalFiles.system_TOP),
                UserSettings.name_solvent_GENION,
                UserSettings.name_positive_ion_GENION,
                UserSettings.name_negative_ion_GENION,
                [os.path.realpath(simulation._temp_dir)]
            )

            # Gromacs.GRO_TOP_neutralize_system(
            #     path(LocalFiles.system_GRO),
            #     path(LocalFiles.system_TOP_noConstraints),
            #     path(LocalFiles.GRO_neutralized),
            #     path(LocalFiles.system_TOP_noConstraints),
            #     UserSettings.name_solvent_GENION,
            #     UserSettings.name_positive_ion_GENION,
            #     UserSettings.name_negative_ion_GENION,
            #     [os.path.realpath(simulation._temp_dir)]
            # )

            # Gromacs.GRO_TOP_neutralize_system(
            #     path(LocalFiles.system_GRO),
            #     path(LocalFiles.system_TOP_noConstr_Restr),
            #     path(LocalFiles.GRO_neutralized),
            #     path(LocalFiles.system_TOP_noConstr_Restr),
            #     UserSettings.name_solvent_GENION,
            #     UserSettings.name_positive_ion_GENION,
            #     UserSettings.name_negative_ion_GENION,
            #     [os.path.realpath(simulation._temp_dir)]
            # )

            # Gromacs.GRO_TOP_neutralize_system(
            #     path(LocalFiles.system_GRO),
            #     path(LocalFiles.system_TOP_restr),
            #     path(LocalFiles.GRO_neutralized),
            #     path(LocalFiles.system_TOP_restr),
            #     UserSettings.name_solvent_GENION,
            #     UserSettings.name_positive_ion_GENION,
            #     UserSettings.name_negative_ion_GENION,
            #     [os.path.realpath(simulation._temp_dir)]
            # )

            # Gromacs.NDX_from_GRO(
            #     path(LocalFiles.GRO_neutralized),
            #     path(LocalFiles.NDX),
            #     UserSettings.resnm_to_indexgroup_dict
            # )

            # System.run_command(
            #     core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
            #         path(LocalFiles.MDP_steep_genion),
            #         path(LocalFiles.GRO_neutralized),
            #         path(LocalFiles.system_TOP),
            #         path(LocalFiles.TPR),
            #         path(LocalFiles.MDOUT),
            #         path(LocalFiles.NDX),
            #         max_warnings=1
            #     ),
            #     path_stdout=path(LocalFiles.stdout_stderr),
            #     path_stderr=path(LocalFiles.stdout_stderr)
            # )

            # System.run_simulation_command(
            #     core.gromacs_commands.GROMACSRunCommands.MDRUN(
            #         path(LocalFiles.TPR),
            #         path(LocalFiles.TRR),
            #         path(LocalFiles.GRO_neutralized),
            #         path(LocalFiles.EDR),
            #         path(LocalFiles.LOG),
            #         path(LocalFiles.CPT),
            #     ),
            #     path_stdout=path(LocalFiles.stdout_stderr),
            #     path_stderr=path(LocalFiles.stdout_stderr)
            # )

        def relax_softcore():
            if UserSettings.neutralize:
                gropath = path(LocalFiles.GRO_neutralized)
            else:
                gropath = path(LocalFiles.system_GRO)

            Gromacs.NDX_from_GRO(
                gropath,
                path(LocalFiles.NDX),
                UserSettings.resnm_to_indexgroup_dict
            )

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
                    path(LocalFiles.MDP_steep_softcore),
                    gropath,
                    path(LocalFiles.system_TOP),
                    path(LocalFiles.TPR),
                    path(LocalFiles.MDOUT),
                    path(LocalFiles.NDX),
                    max_warnings=1
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            System.run_simulation_command(
                core.gromacs_commands.GROMACSRunCommands.MDRUN(
                    path(LocalFiles.TPR),
                    path(LocalFiles.TRR),
                    path(LocalFiles.GRO_softcore),
                    path(LocalFiles.EDR),
                    path(LocalFiles.LOG),
                    path(LocalFiles.CPT),
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            ) 

        def actual_simulation():
            # Steep with restricted angles, to avoid possible 180 degree angles (cannot be done in previous steep as restricted angles + free energy is not yet available) #
            
            # if UserSettings.neutralize:
            #     System.run_command(
            #         core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
            #             path(LocalFiles.MDP_steep_restr),
            #             path(LocalFiles.GRO_neutralized),
            #             path(LocalFiles.system_TOP),
            #             path(LocalFiles.TPR),
            #             path(LocalFiles.MDOUT),
            #             path(LocalFiles.NDX),
            #             max_warnings=1
            #         ),
            #         path_stdout=path(LocalFiles.stdout_stderr),
            #         path_stderr=path(LocalFiles.stdout_stderr)
            #     )
            # else:
            #     Gromacs.NDX_from_GRO(
            #         path(LocalFiles.system_GRO),
            #         path(LocalFiles.NDX),
            #         UserSettings.resnm_to_indexgroup_dict
            #     )
                
            #     System.run_command(
            #         core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
            #             path(LocalFiles.MDP_steep_restr),
            #             path(LocalFiles.system_GRO),
            #             path(LocalFiles.system_TOP),
            #             path(LocalFiles.TPR),
            #             path(LocalFiles.MDOUT),
            #             path(LocalFiles.NDX),
            #             max_warnings=1
            #         ),
            #         path_stdout=path(LocalFiles.stdout_stderr),
            #         path_stderr=path(LocalFiles.stdout_stderr)
            #     )

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
                    path(LocalFiles.MDP_steep_restr),
                    path(LocalFiles.GRO_softcore),
                    path(LocalFiles.system_TOP),
                    path(LocalFiles.TPR),
                    path(LocalFiles.MDOUT),
                    path(LocalFiles.NDX),
                    max_warnings=1
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            System.run_simulation_command(
                core.gromacs_commands.GROMACSRunCommands.MDRUN(
                    path(LocalFiles.TPR),
                    path(LocalFiles.TRR),
                    path(LocalFiles.GRO_steep),
                    path(LocalFiles.EDR),
                    path(LocalFiles.LOG),
                    path(LocalFiles.CPT),
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            # additional NPT equilibration run, very short run at very short timestep to avoid exploding constraints #
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
                    path(LocalFiles.MDP_NPT_equil),
                    path(LocalFiles.GRO_steep),
                    path(LocalFiles.system_TOP),
                    path(LocalFiles.TPR),
                    path(LocalFiles.MDOUT),
                    path(LocalFiles.NDX),
                    max_warnings=1
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            System.run_simulation_command(
                core.gromacs_commands.GROMACSRunCommands.MDRUN_XTC(
                    path(LocalFiles.TPR),
                    path(LocalFiles.TRR),
                    path(LocalFiles.GRO_equil),
                    path(LocalFiles.EDR),
                    path(LocalFiles.LOG),
                    path(LocalFiles.CPT),
                    path(LocalFiles.XTC),
                ),
                tTimeout=UserSettings.timeout_Production_relax_MDRUN,
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )


            # Actual sim #
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
                    path(LocalFiles.MDP_NPT),
                    path(LocalFiles.GRO_equil),
                    path(LocalFiles.system_TOP),
                    path(LocalFiles.TPR),
                    path(LocalFiles.MDOUT),
                    path(LocalFiles.NDX),
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )
    
            System.run_simulation_command(
                core.gromacs_commands.GROMACSRunCommands.MDRUN_XTC(
                    path(LocalFiles.TPR),
                    path(LocalFiles.TRR),
                    path(LocalFiles.GRO),
                    path(LocalFiles.EDR),
                    path(LocalFiles.LOG),
                    path(LocalFiles.CPT),
                    path(LocalFiles.XTC),
                ),
                tTimeout=UserSettings.timeout_Production_MDRUN,
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

        System.copy_dir_files(UserSettings.path_basedir_production, path(""))
        setup()
        if UserSettings.neutralize:
            neutralize_system()
        relax_softcore()
        actual_simulation()

        simulation._share.add_file("prod-%s-GRO" %membrane_type, LocalFiles.GRO, True)
        simulation._share.add_file("prod-%s-EDR" %membrane_type, LocalFiles.EDR, True)
        simulation._share.add_file("prod-%s-LOG" %membrane_type, LocalFiles.LOG, False)
        simulation._share.add_file("prod-%s-MDOUT" %membrane_type, LocalFiles.MDOUT, True)
        if UserSettings.WRITE_XTC:
            #print("production.py: Writing Trajectories.")
            simulation._share.add_file("prod-%s-XTC" %membrane_type, LocalFiles.XTC, True)
            simulation._share.add_file("prod-%s-TPR" %membrane_type, LocalFiles.TPR, True)
        else:
            simulation._share.add_file("prod-%s-XTC" %membrane_type, LocalFiles.XTC, False)
            simulation._share.add_file("prod-%s-TPR" %membrane_type, LocalFiles.TPR, True)
        try:
            simulation._share.add_file("prod-NDX", LocalFiles.NDX, False)
        except ValueError:
            pass
