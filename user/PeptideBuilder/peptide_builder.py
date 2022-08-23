from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB



def generate_peptide_alpha(sequence, filepath):
    alpha_phi = -60
    alpha_psi_im1 = -40

    # create a geometry for initialization of structure
    geo = Geometry.geometry(sequence[0])
    geo.phi = alpha_phi
    geo.psi_im1 = alpha_psi_im1
    structure = PeptideBuilder.initialize_res(geo)
    # add remaining AAs
    for aa in sequence[1:]:
        PeptideBuilder.add_residue(structure, aa, alpha_phi, alpha_psi_im1)
    # add terminal oxygen (OXT) to the final aa
    PeptideBuilder.add_terminal_OXT(structure)

    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(filepath) 