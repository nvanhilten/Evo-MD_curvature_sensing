"""
MODULE VERSION:
Create peptide atomistic PDB with PeptideBuilder (or Modeller in the future) from an amino acid sequence and use martinize2 to generate coarse-grained GRO and ITP using an ssdump secondary structure.
"""

import subprocess
import shutil
import os

import PeptideBuilder

#from Bio import PDB

import core.settings as CoreSettings
import core.system as System
import core.gromacs as Gromacs
import core.gromacs_commands

import user.usersettings as UserSettings
if UserSettings.pdb_gen == "PeptideBuilder":
    try:
        import user.PeptideBuilder.peptide_builder as PB
    except ImportError:
        print("PeptideBuilder was set as PDB generator, but PeptideBuilder import failed.")


def module_generate_peptide(simulation):
    class LocalFiles(): # Filenames used within this module
        if UserSettings.pdb_gen == "Modeller":
            _PDB = UserSettings.template_pdb
            ALI                     = "template.ali"
        _TOP = UserSettings.template_top

        sequence                = UserSettings.filename_seq
        structure               = UserSettings.filename_ssd

        PDB_atomistic           = UserSettings.name_output_peptide + "_atomistic.pdb"
        PDB                     = UserSettings.name_output_peptide + ".pdb"

        TOP_martinize2          = "martinize2_out.top"
        ITP_martinize2          = "molecule_0.itp"

        ITP                     = UserSettings.name_output_peptide + ".itp"
        TOP                     = UserSettings.name_output_peptide + ".top"

        # ITP_noConstraints       = UserSettings.name_output_peptide + ".seq2itp.noconstr.itp"
        # ITP_noConstraints_Restr = UserSettings.name_output_peptide + ".seq2itp.noconstr.restr.itp"
        # ITP_seq2itp             = UserSettings.name_output_peptide + ".seq2itp.itp"
        # ITP_restr               = UserSettings.name_output_peptide + ".restr.itp"

        # TOP_restr               = _TOP + ".restr.top"
        # TOP_noConstraints       = _TOP + ".noConstraints.top"
        # TOP_noConstr_Restr      = _TOP + ".noConstr_Restr.top"

        MDP_steep_pept          = UserSettings.filename_MDP_steep_pept
        MDP_eq_pept             = UserSettings.filename_MDP_eq_pept

        # MDP                     = UserSettings.filename_NVT_peptide_MDP
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

    def generatePDB():
        if UserSettings.pdb_gen == "PeptideBuilder":
            sequence = System.read_text_file(path(LocalFiles.sequence), strip=True)[1]
            #print(sequence)
            PB.generate_peptide_alpha(sequence, path(LocalFiles.PDB_atomistic))
        elif UserSettings.pdb_gen == "Modeller":
            def run_modeller(): # generates a homolgy model (pdb file) from sequence + structure
                System.run_command(
                    UserSettings.command_modeller(
                        path(LocalFiles.sequence),
                        path(LocalFiles.structure),
                        path(LocalFiles.PDB_atomistic),
                        UserSettings.name_output_peptide,
                        path(LocalFiles.ALI),
                        simulation._output_dir
                    ),
                    path_stdout=path(LocalFiles.stdout_stderr),
                    path_stderr=path(LocalFiles.stdout_stderr)
                )
            
            run_modeller()

    def generateCG_PDB_ITP():
        structure = System.read_text_file(path(LocalFiles.structure), strip=False)[1]
        #print(structure)
        try:
            System.run_command(
                UserSettings.command_martinize2(
                    path(LocalFiles.PDB_atomistic),
                    structure,
                    path(LocalFiles.PDB),
                    path(LocalFiles.TOP_martinize2)
                ),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr),
                cwd=simulation._temp_dir,
            )
        except OSError as e:
            print("Martinize2 OSError caught")

        System.replace_in_file(
            path(LocalFiles.ITP_martinize2),
            "molecule_0",
            "Protein"
        )

    def restrict_BBB_angles():
        data = System.read_text_file(path(LocalFiles.ITP_martinize2))
        state = False
        count = 0
        memory = []
        for i in range(len(data)):
            line = data[i].strip()
            if state:
                if line.startswith(";"):
                    state = False
                    new = "#endif\n\n#ifndef RESTRICT_BBB\n"
                    for angle in memory:
                        angle = angle.split()
                        try:
                            if angle[3] == "10":
                                angle[3] = "2"
                                angle[-1] = "1000"
                            new += " ".join(angle) + "\n"
                        except IndexError:
                            pass
                    new += "#endif\n\n"
                    data[i-1] = new
                    continue
                line_data = line.split()
                if count == 0:
                    data[i] = "#ifdef RESTRICT_BBB\n" + "\t".join(line_data[0:3] + ["10"] + line_data[4:] + ["\n"])
                    count+=1
                elif count > 0:
                    data[i] = "\t".join(line_data[0:3] + ["10"] + line_data[4:] + ["\n"])
                memory.append(line)
            else:
                if "BBB angles" in line:
                    state = True
                    continue
        System.write_text_file(path(LocalFiles.ITP), data)

    def generateGRO():
        def reorganize_TOP():
            shutil.copy2(path(LocalFiles._TOP), path(LocalFiles.TOP))

            System.replace_in_file(
                path(LocalFiles.TOP),
                "template.itp",
                LocalFiles.ITP
            )

            System.replace_in_file(
                path(LocalFiles.TOP),
                "molecule_0",
                "Protein"
            )

            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.EDITCONF_pdb2gro(
                    path(LocalFiles.PDB),
                    path(LocalFiles.GRO)
                ),
                input="0\n".encode(),
                path_stdout=path(LocalFiles.stdout_stderr),
                path_stderr=path(LocalFiles.stdout_stderr)
            )

        def relax_peptide():
            # steep pept
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP(
                    path(LocalFiles.MDP_steep_pept),
                    path(LocalFiles.GRO),
                    path(LocalFiles.TOP),
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

            # eq pept
            System.run_command(
                core.gromacs_commands.GROMACSRunCommands.GROMPP(
                    path(LocalFiles.MDP_eq_pept),
                    path(LocalFiles.GRO),
                    path(LocalFiles.TOP),
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
        
        reorganize_TOP()
        relax_peptide()

    # def generateITP():
    #     def run_seq2itp():
    #         _n_retries = 3
    #         while True:
    #             try:
    #                 System.run_command(
    #                     UserSettings.command_seq2ITP(
    #                         path(LocalFiles.sequence),
    #                         path(LocalFiles.structure),
    #                         path(LocalFiles.ITP_seq2itp)
    #                     ),
    #                     path_stdout=path(LocalFiles.stdout_stderr),
    #                     path_stderr=path(LocalFiles.stdout_stderr)
    #                 )
    #             except OSError as e:
    #                 _n_retries -= 1
    #                 print("Seq2itp OSError caught, retry left=" + str(_n_retries))
    #                 if _n_retries > 0:
    #                     continue
    #                 else:
    #                     raise
    #             break

        # def restrict_backbone_angles_seq2itp():
        #     data = System.read_text_file(path(LocalFiles.ITP_seq2itp))
        #     state = False
        #     for i in range(len(data)):
        #         line = data[i].strip()
        #         if state:
        #             if line.startswith(";"):
        #                 state = False
        #                 continue
        #             line_data = line.split()
        #             data[i] = "\t".join(line_data[0:3] + ["10"] + line_data[4:] + ["\n"])
        #         else:
        #             if "backbone-backbone-backbone angles" in line:
        #                 state = True
        #                 continue
        #     System.write_text_file(path(LocalFiles.ITP_restr), data)

        # run_seq2itp()

        # restrict_backbone_angles_seq2itp()

        # Gromacs.ITP_constraints_to_bonds(
        #     path(LocalFiles.ITP_seq2itp),
        #     path(LocalFiles.ITP_noConstraints),
        #     UserSettings.STEEP_constraint_forceConstant
        # )
        # Gromacs.ITP_constraints_to_bonds(
        #     path(LocalFiles.ITP_restr),
        #     path(LocalFiles.ITP_noConstraints_Restr),
        #     UserSettings.STEEP_constraint_forceConstant
        # )

    # def generate_GRO_from_ITP():
    #     def relax_structure():
    #         shutil.copy2(path(LocalFiles._TOP), path(LocalFiles.TOP_restr))
    #         shutil.copy2(path(LocalFiles._TOP), path(LocalFiles.TOP_noConstraints))
    #         shutil.copy2(path(LocalFiles._TOP), path(LocalFiles.TOP_noConstr_Restr))

    #         System.replace_in_file(
    #             path(LocalFiles.TOP_restr),
    #             UserSettings.tags_steep_peptide_TOP[0],
    #             "\"" + LocalFiles.ITP_restr + "\""
    #         )
    #         System.replace_in_file(
    #             path(LocalFiles.TOP_noConstraints),
    #             UserSettings.tags_steep_peptide_TOP[0],
    #             "\"" + LocalFiles.ITP_noConstraints + "\""
    #         )
    #         System.replace_in_file(
    #             path(LocalFiles.TOP_noConstr_Restr),
    #             UserSettings.tags_steep_peptide_TOP[0],
    #             "\"" + LocalFiles.ITP_noConstraints_Restr + "\""
    #         )

    #         System.replace_in_file(
    #             path(LocalFiles.MDP),
    #             UserSettings.tags_NVT_peptide_MDP[0],
    #             "-I" + os.path.realpath(simulation._temp_dir)
    #         )

    #         System.run_command(
    #             core.gromacs_commands.GROMACSRunCommands.GROMPP(
    #                 path(LocalFiles.MDP),
    #                 path(LocalFiles.GRO),
    #                 path(LocalFiles.TOP_noConstraints),
    #                 path(LocalFiles.TPR),
    #                 path(LocalFiles.MDOUT)
    #             ),
    #             path_stdout=path(LocalFiles.stdout_stderr),
    #             path_stderr=path(LocalFiles.stdout_stderr)
    #         )

    #         System.run_simulation_command(
    #             core.gromacs_commands.GROMACSRunCommands.MDRUN_XTC(
    #                 path(LocalFiles.TPR),
    #                 path(LocalFiles.TRR),
    #                 path(LocalFiles.GRO),
    #                 path(LocalFiles.EDR),
    #                 path(LocalFiles.LOG),
    #                 path(LocalFiles.CPT),
    #                 path(LocalFiles.XTC)
    #             ),
    #             tTimeout=UserSettings.timeout_generate_peptide_MDRUN,
    #             path_stdout=path(LocalFiles.stdout_stderr),
    #             path_stderr=path(LocalFiles.stdout_stderr)
    #         )

    #         System.run_command(
    #             core.gromacs_commands.GROMACSRunCommands.TRJCONV_center_PBC(
    #                 path(LocalFiles.GRO),
    #                 path(LocalFiles.TPR)
    #             ),
    #             input="0\n0\n".encode(),
    #             path_stdout=path(LocalFiles.stdout_stderr),
    #             path_stderr=path(LocalFiles.stdout_stderr)
    #         )

    #         System.run_command(
    #             core.gromacs_commands.GROMACSRunCommands.EDITCONF_align_Z_center(
    #                 path(LocalFiles.GRO),
    #                 path(LocalFiles.GRO)
    #             ),
    #             input="0\n".encode(),
    #             path_stdout=path(LocalFiles.stdout_stderr),
    #             path_stderr=path(LocalFiles.stdout_stderr)
    #         )

    #     Gromacs.ITP_convert_to_GRO(
    #         path(LocalFiles.ITP_restr),
    #         path(LocalFiles.GRO)
    #     )

    #     relax_structure()

    System.copy_dir_files(UserSettings.path_basedir_generate_peptide, path(""))

    generatePDB()
    generateCG_PDB_ITP()
    restrict_BBB_angles()
    generateGRO()

    simulation._share.add_file("gen-PDB_atomistic",         LocalFiles.PDB_atomistic,       UserSettings.FILE_WRITE_DEBUGGING)
    simulation._share.add_file("gen-PDB",                   LocalFiles.PDB,                 UserSettings.FILE_WRITE_DEBUGGING)
    simulation._share.add_file("gen-ITP_martinize2",        LocalFiles.ITP_martinize2,      UserSettings.FILE_WRITE_DEBUGGING)
    simulation._share.add_file("gen-ITP",                   LocalFiles.ITP,                 True)
    simulation._share.add_file("gen-GRO",                   LocalFiles.GRO,                 True)
    simulation._share.add_file("gen-TOP",                   LocalFiles.TOP,                 UserSettings.FILE_WRITE_DEBUGGING)
    simulation._share.add_file("gen-TPR",                   LocalFiles.TPR,                 UserSettings.FILE_WRITE_DEBUGGING)
    simulation._share.add_file("gen-XTC",                   LocalFiles.XTC,                 UserSettings.FILE_WRITE_DEBUGGING)


    # simulation._share.add_file("gen-TOP_restr",             LocalFiles.TOP_restr,           False)
    # simulation._share.add_file("gen-TOP_noConstraints",     LocalFiles.TOP_noConstraints,   False)
    # simulation._share.add_file("gen-TOP_noConstr_Restr",    LocalFiles.TOP_noConstr_Restr,  False)

