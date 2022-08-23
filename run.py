"""
Executable wrapper script. 
Initializes genetic algorithm and defines the fitness function.
"""

import sys
import os
import subprocess

if not sys.version_info[0] == 3:
    sys.exit("Script requires Python 3.")

import core.settings            as CoreSettings
import core.simulation          as Sim
import core.system              as System
#import core.gromacs             as Gromacs
import core.genetic_algorithm   as GA

import user.usersettings        as UserSettings

from core.system                import mprint


def _setup_simulation_directory(sequence, path_parent_dir, dirname=None):
    if dirname is None:
        dirname = sequence.upper()
    sequence_dir = os.path.join(path_parent_dir, dirname)
    System.mkdir(sequence_dir)
    System.rm_dir_content(sequence_dir) # Ensure folder is empty in case of a rerun

    System.write_text_file(os.path.join(sequence_dir, UserSettings.filename_seq), [">Test", sequence], add_newline=True)
    if UserSettings.SSP:
        _n_retries = 3
        while True:
            try:
                System.check_exitcode(
                    subprocess.run(
                        UserSettings.command_Porter5(
                            os.path.join(sequence_dir, UserSettings.filename_seq)
                        ).split()
                    ).returncode
                )
            except OSError as e:
                _n_retries -= 1
                print("OSError caught, retry left=" + str(_n_retries))
                if _n_retries > 0:
                    continue
                else:
                    raise
            break
    else:
#        n_tag = len(sequence)-24
#        System.write_text_file(os.path.join(sequence_dir, UserSettings.filename_ssd), [str(len(sequence)), "C"*n_tag+"H"*(len(sequence)-n_tag)], add_newline=True)
        System.write_text_file(os.path.join(sequence_dir, UserSettings.filename_ssd), [str(len(sequence)), "H"*(len(sequence))], add_newline=True)
    return sequence_dir

    # KAI START
def _compute_fitness_simulation_2mem(sequence, dirname=None):
    from user.modules.generate_peptide_PB import module_generate_peptide
    from user.modules.insert_peptide_2mem import module_insert_peptide_2mem 
    from user.modules.production          import module_production_2mem 
    from user.modules.compute_fitness     import module_compute_fitness_2mem_dict
    
    sequence = "".join(sequence)
    sequence_dir = _setup_simulation_directory(sequence, CoreSettings.output_dir, dirname)

    N_retries = 0
    fitness = -1

    while True:
        try:
            simulation = Sim.Simulation(sequence_dir)

            simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_seq))
            simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_ssd))

            simulation.add_module(Sim.Module(module_generate_peptide))
            simulation.add_module(Sim.Module(module_insert_peptide_2mem))
            simulation.add_module(Sim.Module(module_production_2mem))
            simulation.add_module(Sim.Module(module_compute_fitness_2mem_dict))

            simulation.run()

            ddF = simulation.get_output_variable("ddF")
            fitness = simulation.get_output_variable("fitness")

        except System.SimulationError:
            N_retries += 1
            if N_retries > CoreSettings.N_simulation_retry:
                print("Simulation stopped due to an error. Simulation retries exceeded. Exiting.")
                sys.exit(1)
            else:
                print("Simulation stopped due to an error. Retrying simulation. %d retries left." % (CoreSettings.N_simulation_retry - N_retries))
            continue
        break

    return [fitness, ddF]


def _run_specific_sequences_trajectories(sequences):
    if not UserSettings.WRITE_XTC:
        print("Warning. UserSettings.WRITE_XTC set to False.")
#        sys.exit(1)
    MPI = None
    try:
        from mpi4py import MPI
    except ImportError:
        print("ERROR: Missing mpi4py. Exiting.")
        sys.exit(1)

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    buffer = [[] for i in range(size)]
    
    count = 0
    while len(sequences) > 0:
        buffer[count].append(sequences.pop())
        count += 1
        if count >= size:
            count = 0

    queue = buffer[rank]
    count = 0
    for sequence in queue:
        _compute_fitness_simulation_2mem(sequence, "".join(sequence) + str(rank) + "_" + str(count))
        count += 1
    print("Done.")

# def _compute_fitness_simulation_2mem_long(sequence, dirname=None):
#     from user.modules.generate_peptide_PB import module_generate_peptide
#     from user.modules.insert_peptide_2mem import module_insert_peptide_2mem 
#     from user.modules.production_long     import module_production_long_2mem 
#     from user.modules.compute_fitness     import module_compute_fitness_2mem_dict
    
#     sequence = "".join(sequence)
#     sequence_dir = _setup_simulation_directory(sequence, CoreSettings.output_dir, dirname)

#     N_retries = 0
#     fitness = -1

#     while True:
#         try:
#             simulation = Sim.Simulation(sequence_dir)

#             simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_seq))
#             simulation.add_input_file(os.path.join(sequence_dir, UserSettings.filename_ssd))

#             simulation.add_module(Sim.Module(module_generate_peptide))
#             simulation.add_module(Sim.Module(module_insert_peptide_2mem))
#             simulation.add_module(Sim.Module(module_production_long_2mem))
#             simulation.add_module(Sim.Module(module_compute_fitness_2mem_dict))


#             simulation.run()

#             fitness = simulation.get_output_variable("fitness")

#         except System.SimulationError:
#             N_retries += 1
#             if N_retries > CoreSettings.N_simulation_retry:
#                 print("Simulation stopped due to an error. Simulation retries exceeded. Exiting.")
#                 sys.exit(1)
#             else:
#                 print("Simulation stopped due to an error. Retrying simulation. %d retries left." % (CoreSettings.N_simulation_retry - N_retries))
#             continue
#         break

#     return fitness

# def _run_specific_sequences_trajectories_long(sequences):
#     if not UserSettings.WRITE_XTC:
#         print("Warning. UserSettings.WRITE_XTC set to False.")
#         sys.exit(1)
#     MPI = None
#     try:
#         from mpi4py import MPI
#     except ImportError:
#         print("ERROR: Missing mpi4py. Exiting.")
#         sys.exit(1)

#     comm = MPI.COMM_WORLD
#     size = comm.Get_size()
#     rank = comm.Get_rank()

#     buffer = [[] for i in range(size)]
    
#     count = 0
#     while len(sequences) > 0:
#         buffer[count].append(sequences.pop())
#         count += 1
#         if count >= size:
#             count = 0

#     queue = buffer[rank]
#     count = 0
#     for sequence in queue:
#         _compute_fitness_simulation_2mem_long(sequence, "".join(sequence) + str(rank) + "_" + str(count))
#         count += 1
#     print("Done.")

class Constraint:
    """Contains functions that can be passed to GeneticAlgorithm.ptr_constraint_function to constrain the search space."""
    @staticmethod
    def _constrain_peptide_charge(sequence):
        CHARGE_THRESHOLD = 1 # x*2 when using mirror sequences.

        positive = ["R", "K"] # Histidine left out at pH 7, map to one of single-protonated forms in CHARMM36
        negative = ["D", "E"]

        charge = 0
        for x in positive:
            charge += sequence.count(x)

        for x in negative:
            charge -= sequence.count(x)

        return (abs(charge) > CHARGE_THRESHOLD)

    @staticmethod
    def _constrain_hydrophilic_residue_count(sequence):
        MAX_HYDROPHILIC_COUNT = 2 # x*2 when using mirror sequences

        hydrophilic = ["D", "E", "H", "K", "N", "P", "Q", "R", "S", "T", "W", "Y"]

        count = 0
        for x in hydrophilic:
            count += sequence.count(x)

        return count > MAX_HYDROPHILIC_COUNT

    @staticmethod
    def _constrain_hydrophilic_edge(sequence):
        MAX_HYDROPHILIC_COUNT = 2

        hydrophilic = ["D", "E", "H", "K", "N", "P", "Q", "R", "S", "T", "W", "Y"]

        for x in sequence[0:MAX_HYDROPHILIC_COUNT]:
            if x not in hydrophilic:
                return True

        for x in sequence[MAX_HYDROPHILIC_COUNT:]:
            if x in hydrophilic:
                return True

        return False

    @staticmethod
    def _constrain_length(sequence):
        if len(sequence) < UserSettings.min_length:
            return True
        elif len(sequence) > UserSettings.max_length:
            return True

        return False

###########################
## GA setup and settings ##
###########################

def main():
    g = GA.GeneticAlgorithm()
    g.description = """Default GA run."""

    #g.sequence_length = 10 # Resulting peptides will be length=20, implemented in _compute_fitness_simulation_mirror

    g.ptr_scoring_function      = _compute_fitness_simulation_2mem
    g.population_size           = 144 # 128 #128
    g.iterations                = 2 # 100 #100
    g.N_parents                 = 36 # 16  #16
    g.ptr_selection_function    = GA.SelectionMethods.best_N
    g.selfun_exponential_rank_C = 0.5
    g.N_elites                  = 1   #2
    g.N_rerun_elites            = 1   #2
    g.N_rerun_threshold         = 4
    g.filepath_progress_output  = os.path.join(CoreSettings.output_dir, CoreSettings.name_output_progress)
    g.filepath_checkpoint       = os.path.join(CoreSettings.output_dir, CoreSettings.name_checkpoint_file)

    g.initialize(path_restart_file=CoreSettings.filepath_checkpoint)
    g.run()
    result = g.get_result()
    print(result[0][0] + "\t" + str(result[0][1]))

##########
## MAIN ##
##########

if __name__ == "__main__":
    CoreSettings.parse_user_input()

    STATE = 0 # Allow switch between actual genetic algorithm and 'test' states
    if not STATE == 0:
        mprint("WARNING:\tRunning non-default GA script, see run.py for details.")

    if      STATE == 0: # Genetic Algorithm
        if UserSettings.WRITE_XTC:
            mprint("Warning. UserSettings.WRITE_XTC set to True. This will result in a very large amount of data on disk.")
            sys.exit(1)
        main()
    elif    STATE == 2:
        mprint("Running specific sequences buffer.")
        if not UserSettings.WRITE_XTC:
            mprint("Warning. UserSettings.WRITE_XTC set to False.")
            #sys.exit(1)

        REPEAT = 2 # 100

        seq_file = open('seq.txt')
        seq_lines = seq_file.readlines()
        SEQUENCES = []
        for line in seq_lines:
            SEQUENCES.append(line[:-1])


        _run_specific_sequences_trajectories(SEQUENCES*REPEAT)
    elif    STATE == 3: # run list of sequences using production_long module
        mprint("Running LONG specific sequences buffer.")
        if not UserSettings.WRITE_XTC:
            mprint("Warning. UserSettings.WRITE_XTC set to False.")
            sys.exit(1)

        # change mdp files
        UserSettings.filename_NPT_system_MDP = UserSettings.filename_NPT_long_system_MDP
        if UserSettings.useGPU:
            UserSettings.filename_NPT_gpu_system_MDP = "npt_long_gpu.mdp"

        seq_file = open('seq.txt')
        seq_lines = seq_file.readlines()
        SEQUENCES = []
        for line in seq_lines:
            SEQUENCES.append(line[:-1])

        _run_specific_sequences_trajectories(SEQUENCES)
