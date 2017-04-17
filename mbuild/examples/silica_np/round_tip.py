from mbuild.compound import Compound
from mbuild.tools.tiled_compound import TiledCompound
from mbuild.tools.mask import hemisphere_mask
from mbuild.tools.mask import apply_mask
import numpy as np

from mbuild.examples.alkane_monolayer.alkylsilane import AlkylSilane
from tip_2nm import Tip2nm
from mbuild.components.atoms.H import H
from mbuild.coordinate_transform import _extract_atom_positions, _write_back_atom_positions


class RoundTip(Compound):
    """ """

    def __init__(self, tile_x=1, tile_y=1, chain_length=10, n_chains=38):
        """

        Args:
            chain_length (int): Number of carbon atoms per chain.
            mask (np.ndarray): A 2D array of binding locations.
        """
        super(RoundTip, self).__init__()

        surface = Tip2nm()
        # Replicate the surface.
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")
        self.add(tc, 'tiled_surface')
        for atom in self.yield_atoms():
            if atom.kind == 'O' or atom.kind == 'OB':
                atom.kind = 'OT'
            if atom.kind == 'Si':
                atom.kind == 'SiT'
        alkylsilane = AlkylSilane(chain_length)
        hydrogen = H()

        # Attach chains to specified binding sites. Other sites get a hydrogen.
        mask = hemisphere_mask(n_chains*2)
        mask *= 2.0
        atom_positions = _extract_atom_positions(self)
        x = atom_positions[:,0]
        y = atom_positions[:,1]
        z = atom_positions[:,2]
        check = open('test.xyz','w')
        check.write(str(len(mask)+len(z))+'\n\n')
        for a in range(len(mask)):
            check.write('C\t')
            for b in range(0,3):
                check.write(str(mask[a][b])+'\t')
            check.write('\n')
        for a in range(len(atom_positions)):
            check.write('Si\t')
            for b in range(0,3):
                check.write(str(atom_positions[a][b])+'\t')
            check.write('\n')
        apply_mask(self.tiled_surface, alkylsilane, mask, backfill=hydrogen)

def main():
    tip = RoundTip(chain_length=18,n_chains=38)

    tip = tip.to_trajectory()  # Convert from mbuild to mdtraj
    tip.topology.find_forcefield_terms()  # Find angles/dihedrals from bonds.

    tip.save(filename='data.round_tip_c18')  # Print a LAMMPS data file.
    tip.save(filename='round_tip_c18.pdb')  # Print a .pdb file.

if __name__ == "__main__":
    main()
