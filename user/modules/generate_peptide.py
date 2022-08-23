"""
MODULE VERSION:
Create peptide GRO and ITP from an amino acid sequence and an ssdump secondary structure.
"""

import subprocess
import shutil
import os

import core.settings as CoreSettings
import core.system as System
import core.gromacs as Gromacs
import core.gromacs_commands

import user.usersettings as UserSettings

def module_generate_peptide(simulation):
    class LocalFiles(): # Filenames used within this module
        _TOP = UserSettings.filename_steep_peptide_TOP

        sequence                = UserSettings.filename_seq
        structure               = UserSettings.filename_ssd

        ITP_noConstraints       = UserSettings.name_output_peptide + ".seq2itp.noconstr.itp"
        ITP_noConstraints_Restr = UserSettings.name_output_peptide + ".seq2itp.noconstr.restr.itp"
        ITP_seq2itp             = UserSettings.name_output_peptide + ".seq2itp.itp"
        ITP_restr               = UserSettings.name_output_peptide + ".restr.itp"

        TOP_restr               = _TOP + ".restr.top"
        TOP_noConstraints       = _TOP + ".noConstraints.top"
        TOP_noConstr_Restr      = _TOP + ".noConstr_Restr.top"

        MDP                     = UserSettings.filename_NVT_peptide_MDP
        TPR                     = UserSettings.name_output_peptide + ".tpr"
        TRR                     = UserSettings.name_output_peptide + ".trr"
        GRO                     = UserSettings.name_output_peptide + ".gro"
        EDR                     = UserSettings.name_output_peptide + ".edr"
        LOG                     = UserSettings.name_output_peptide + ".log"
        CPT                     = UserSettings.name_output_peptide + ".cpt"
        XTC                     = UserSettings.name_output_peptide + ".xtc"
        MDOUT                   = UserSettings.name_output_peptide + ".mdout.mdp"

        stdout_stderr           = None
        if CoreSettings.write_stdout_stderr_to_file:
            stdout_stderr = os.path.join(simulation._output_dir, CoreSettings.filename_output_stdout_stderr)

    def path(filename):
        if filename is None:
            return None
        return os.path.join(simulation._temp_dir, filename)

    def generateITP():
        def run_seq2itp():
            _n_retries = 3
            while True:
                try:
                    System.run_command(
                        UserSettings.command_seq2ITP(
                            path(LocalFiles.sequence),
                            path(LocalFiles.structure),
                            path(LocalFiles.ITP_seq2itp)
                        ),
                        path_stdout=path(LocalFiles.stdout_stderr),
                        path_stderr=path(LocalFiles.stdout_stderr)
                    )
                except OSError as e:
                    _n_retries -= 1
                    print("Seq2itp OSError caught, retry left=" + str(_n_retries))
                    if _n_retries > 0:
                        continue
                    else:
                        raise
                break

        def restrict_backbone_angles_seq2itp():
            data = System.read_text_file(path(LocalFiles.ITP_seq2itp))
            state = False
            for i in range(len(data)):
                line = data[i].strip()
                if state:
                    if line.startswith(";"):
                        state = False
                        continue
                    line_data = line.split()
                    data[i] = "\t".join(line_data[0:3] + ["10"] + line_data[4:] + ["\n"])
                else:
                    if "backbone-backbone-backbone angles" in line:
                        state = True
                        continue
            System.write_text_file(path(LocalFiles.ITP_restr), data)

        run_seq2itp()

        restrict_backbone_angles_seq2itp()

        Gromacs.ITP_constraints_to_bonds(
            path(LocalFiles.ITP_seq2itp),
            path(LocalFiles.ITP_noConstraints),
            UserSettings.STEEP_constraint_forceConstant
        )
        Gromacs.ITP_constraints_to_bonds(
            path(LocalFiles.ITP_restr),
            path(LocalFiles.ITP_noConstraints_Restr),
            UserSettings.STEEP_constraint_forceConstant
        )

    def generate_GRO_from_ITP():
        def relax_structure():
            shutil.copy2(path(LocalFiles._TOP), path(LocalFiles.TOP_restr))
            shutil.copy2(path(LocalFiles._TOP), path(LocalFiles.TOP_noConstraints))
            shutil.copy2(path(LocalFiles._TOP), path(LocalFiles.TOP_noConstr_Restr))

            System.replace_in_file(
                path(LocalFiles.TOP_restr),
                UserSettings.tags_steep_peptide_TOP[0],
                "\"" + LocalFiles.ITP_restr + "\""
            )
            System.replace_in_file(
                path(LocalFiles.TOP_noConstraints),
                UserSettings.tags_steep_peptide_TOP[0],
                "\"" + LocalFiles.ITP_noConstraints + "\""
            )
            System.replace_in_file(
                path(LocalFiles.TOP_noConstr_Restr),
                UserSettings.tags_steep_peptide_TOP[0],
                "\"" + LocalFiles.ITP_noConstraints_Restr + "\""
            )

            System.replace_in_file(
                path(LocalFiles.MDP),
                UserSettings.tags_NVT_peptide_MDP[0],
                "-I" + os.path.realpath(simulation._temp_dir)
            )

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP(
                    path(LocalFiles.MDP),
                    path(LocalFiles.GRO),
                    path(LocalFiles.TOP_noConstraints),
                    path(LocalFiles.TPR),
                    path(LocalFiles.MDOUT)
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
                    path(LocalFiles.XTC)
                ),
                tTimeout=UserSettings.timeout_generate_peptide_MDRUN,
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.TRJCONV_center_PBC(
                    path(LocalFiles.GRO),
                    path(LocalFiles.TPR)
                ),
                input="0\n0\n".encode(),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.EDITCONF_align_Z_center(
                    path(LocalFiles.GRO),
                    path(LocalFiles.GRO)
                ),
                input="0\n".encode(),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

        Gromacs.ITP_convert_to_GRO(
            path(LocalFiles.ITP_restr),
            path(LocalFiles.GRO)
        )

        relax_structure()

    System.copy_dir_files(UserSettings.path_basedir_generate_peptide, path(""))

    generateITP()
    generate_GRO_from_ITP()

    simulation._share.add_file("gen-GRO",                   LocalFiles.GRO,                 False)
    simulation._share.add_file("gen-TOP_restr",             LocalFiles.TOP_restr,           False)
    simulation._share.add_file("gen-TOP_noConstraints",     LocalFiles.TOP_noConstraints,   False)
    simulation._share.add_file("gen-TOP_noConstr_Restr",    LocalFiles.TOP_noConstr_Restr,  False)

