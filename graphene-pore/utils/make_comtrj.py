import mdtraj as md
import numpy as np
from mdtraj.core.element import virtual_site


def make_comtrj(trj):
    """
    *** THIS FUNCTION IS FROM MTOOLS: https://github.com/mattwthompson/mtools ***
    Takes a trj and returns a trj with COM positions as atoms
    """

    comtop = md.Topology()
    coords = np.ndarray(shape=(trj.n_frames, trj.n_residues, 3))

    for j, res in enumerate(trj.topology.residues):
        comtop.add_atom(res.name, virtual_site, comtop.add_residue(res.name,     comtop.add_chain()))
        res_frame = trj.atom_slice([at.index for at in res.atoms])
        coords[:, j, :] = md.compute_center_of_mass(res_frame)

    comtrj = md.Trajectory(xyz=coords,
                           topology=comtop,
                           time=trj.time,
                           unitcell_angles=trj.unitcell_angles,
                           unitcell_lengths=trj.unitcell_lengths)

    return comtrj
