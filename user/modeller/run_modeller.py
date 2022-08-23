import sys
import os
import shutil

from modeller import *
from modeller.automodel import *

seq = sys.argv[1]
ssd = sys.argv[2]
pdb = sys.argv[3]
name = sys.argv[4]
ali = sys.argv[5]
outdir = sys.argv[6]

# set working directory
dir = os.path.join(outdir, name + "_modeller")
if not os.path.exists(dir):
    os.makedirs(dir)
os.chdir(dir)

# get sequence
infile = open(seq, "r")
for line in infile:
    line = line.strip()
    if not ">" in line:
        sequence = line
infile.close()

# get structure
infile = open(ssd, "r")
for line in infile:
    line = line.strip()
    if not ">" in line:
        structure = line
infile.close()

env = Environ()
log.level(output=0, notes=0, warnings=0, errors=0, memory=0)
env.io.atom_files_directory = [os.path.dirname(pdb)]

# write ali file
length = len(sequence)
temp_seq = "IGDKVNLRQKLLNM"
temp_length = len(temp_seq)
outfile = open(ali, "w")
outfile.write(">P1;templateA\nstructureX:template:81:A:94:A::::\n" + temp_seq + "-"*length + "*\n>P1;" + name + "\nsequence:" + name + "::::::::\n" + "-"*temp_length + sequence + "*")
outfile.close()



class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms

        def handle_SS(structure):
            prev = "A"
            count = 1
            helix_list = []
            strand_list = []
            h = -1
            e = -1
            for i in structure:
                if i != prev: # when entering a new SS region
                    if i == "C":
                        if prev == "H": # finish the previous helix
                            helix_list[h].append(str(count-1) + ":A")
                        if prev == "E": # finish the previous strand
                            strand_list[e].append(str(count-1) + ":A")
                    if i == "H":
                        helix_list.append([str(count) + ":A"]) # begin the new helix
                        h+=1
                        if prev == "E": # finish the previous strand
                            strand_list[e].append(str(count-1) + ":A")
                    if i == "E":
                        strand_list.append([str(count) + ":A"]) # begin the new strand
                        e+=1
                        if prev == "H": # finish the previous helix
                            helix_list[h].append(str(count-1) + ":A")
                if count == len(structure): # end of sequence
                    if i == "H":
                        helix_list[h].append(str(count)+":A")
                    if i == "E":
                        strand_list[e].append(str(count)+":A")
                prev = i
                count+=1
            return(helix_list, strand_list)

        helix_list, strand_list = handle_SS(structure)

        # Helices
        if len(helix_list) > 0:
            for helix in helix_list:
                rsr.add(secondary_structure.alpha(self.residue_range(helix[0], helix[1])))

        # Strands
        if len(strand_list) > 0:
            for strand in strand_list:
                rsr.add(secondary_structure.strand(self.residue_range(strand[0], strand[1])))

        # Anti-parallel sheet (-E3_C2-5_E3-)
        i_list = []
        m_list = []
        for c in [2,3,4,5]:
            motif = "E"*3 + "C"*c + "E"*3
            if motif in structure:
                for i in range(len(structure)):
                    if structure[i:i+len(motif)] == motif:
                        i_list.append(i)
                        m_list.append(len(motif))
        
        sheet_list = []
        start_switch = False
        l = 0
        for i in i_list:
            c = 0
            for strand in strand_list:
                start = int(strand[0].split(":")[0]) - 1
                try:
                    end = int(strand_list[c+1][1].split(":")[0]) -1
                    if start <= i and end >= i + m_list[l] - 1:
                        sheet_list.append(["N:"+str(start+1)+":A",  "O:"+str(end+1)+":A"])
                    c+=1
                except IndexError:
                    pass
            l+=1
        
        if len(sheet_list) > 0:
            for sheet in sheet_list:
                rsr.add(secondary_structure.sheet(at[sheet[0]], at[sheet[1]], sheet_h_bonds=-5))
        

a=MyModel(env, alnfile=ali,
                knowns='templateA', sequence=name,
                assess_methods=(assess.DOPE, assess.GA341))
a.starting_model=1
a.ending_model=50
a.write_intermediates = False
a.make()

#   Get a list of all successfully built models from a.outputs
ok_models = [x for x in a.outputs if x['failure'] is None]

#   Rank the models by DOPE score
key = 'DOPE score'

#   ok_models.sort(lambda a,b: cmp(a[key], b[key]))
ok_models.sort(key=lambda a: a[key])

#   Get top model
m = ok_models[0]
print ("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))

shutil.copy(os.path.join(dir, m['name']), pdb)
