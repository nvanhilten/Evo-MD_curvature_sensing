"""
MODULE VERSION:
Combining peptide and membrane structures into a single system.
"""

import os

import core.settings as CoreSettings
import core.system as System
import core.gromacs as Gromacs
import core.gromacs_commands

import user.usersettings as UserSettings

def module_insert_peptide(simulation):
    """
    Takes a peptide GRO file and a membrane GRO file, combines the two into a new GRO file and runs a steepest descent energy minimization.
    """
    class LocalFiles(): # Filenames used within this module
        membrane_GRO                = UserSettings.filename_membrane_GRO
        membrane_TOP                = UserSettings.filename_membrane_TOP

        peptide_GRO                 = simulation._share.get_filename("gen-GRO")
        peptide_TOP_restr           = simulation._share.get_filename("gen-TOP_restr")
        peptide_TOP_noConstraints   = simulation._share.get_filename("gen-TOP_noConstraints")
        peptide_TOP_noConstr_Restr  = simulation._share.get_filename("gen-TOP_noConstr_Restr")

        TOP_restr                   = UserSettings.name_output_system + ".restr.top"
        TOP_noConstraints           = UserSettings.name_output_system + ".noConstraints.top"
        TOP_noConstr_Restr          = UserSettings.name_output_system + ".noConstr_Restr.top"

        GRO                         = UserSettings.name_output_system + ".gro"
        NDX                         = UserSettings.name_output_system + ".ndx"
        TPR                         = UserSettings.name_output_system + ".tpr"
        MDOUT                       = UserSettings.name_output_system + ".mdout.mdp"
        TRR                         = UserSettings.name_output_system + ".trr"
        EDR                         = UserSettings.name_output_system + ".edr"
        LOG                         = UserSettings.name_output_system + ".log"
        CPT                         = UserSettings.name_output_system + ".cpt"

        MDP = UserSettings.filename_steep_system_MDP

        stdout_stderr = None
        if CoreSettings.write_stdout_stderr_to_file:
            stdout_stderr = os.path.join(simulation._output_dir, CoreSettings.filename_output_stdout_stderr)

    def path(filename):
        if filename is None:
            return None
        return os.path.join(simulation._temp_dir, filename)

    def center_peptide_to_membrane():
        box_vector = Gromacs.GRO_read_box(path(LocalFiles.membrane_GRO))
        data_peptide_GRO = System.read_text_file(path(LocalFiles.peptide_GRO))

        for i in reversed(range(0, len(data_peptide_GRO))):
            line = data_peptide_GRO[i].strip()
            if not len(line) == 0:
                data_peptide_GRO[i] = " ".join(box_vector)
                break

        os.remove(path(LocalFiles.peptide_GRO))
        System.write_text_file(path(LocalFiles.peptide_GRO), data_peptide_GRO, add_newline=False)

        System.run_command(
            core.gromacs_commands.GROMACSRunCommands.EDITCONF_center(
                path(LocalFiles.peptide_GRO),
                path(LocalFiles.peptide_GRO)
            ),
            path_stdout=path(LocalFiles.stdout_stderr),
            path_stderr=path(LocalFiles.stdout_stderr)
        )

    def insert_into_membrane():
        base_sections = [Gromacs.TOP_HEADER_SECTION, "system", "molecules"]
        append_sections = [Gromacs.TOP_HEADER_SECTION, "molecules"]

        Gromacs.GRO_merge(
            path(LocalFiles.membrane_GRO),
            path(LocalFiles.peptide_GRO),
            path(LocalFiles.GRO)
        )
        Gromacs.TOP_merge(
            path(LocalFiles.membrane_TOP),
            path(LocalFiles.peptide_TOP_restr),
            path(LocalFiles.TOP_restr),
            base_sections, append_sections
        )
        Gromacs.TOP_merge(
            path(LocalFiles.membrane_TOP),
            path(LocalFiles.peptide_TOP_noConstraints),
            path(LocalFiles.TOP_noConstraints),
            base_sections, append_sections
        )
        Gromacs.TOP_merge(
            path(LocalFiles.membrane_TOP),
            path(LocalFiles.peptide_TOP_noConstr_Restr),
            path(LocalFiles.TOP_noConstr_Restr),
            base_sections, append_sections
        )
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
            core.gromacs_commands.GROMACSRunCommands.GROMPP_NDX(
                path(LocalFiles.MDP),
                path(LocalFiles.GRO),
                path(LocalFiles.TOP_noConstraints),
                path(LocalFiles.TPR),
                path(LocalFiles.MDOUT),
                path(LocalFiles.NDX)
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

    simulation._share.add_file("ins-GRO",                   LocalFiles.GRO,                 False)
    simulation._share.add_file("ins-TOP_restr",             LocalFiles.TOP_restr,           False)
    simulation._share.add_file("ins-TOP_noConstraints",     LocalFiles.TOP_noConstraints,   False)
    simulation._share.add_file("ins-TOP_noConstr_Restr",    LocalFiles.TOP_noConstr_Restr,  False)
