import numpy as np
import mdtraj as md
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from utils.make_comtrj import make_comtrj


def calc_number_density(gro_file, trj_file, top_file, area,
        dim, box_range, n_bins, frame_range=None,maxs=None, mins=None):
    """
    Calculate a 1-dimensional number density profile for each residue

    Parameters
    ----------
    gro_file: str
        GROMACS '.gro' file to load 
    trj_file: str
        Trajectory to load
    top_file: str
        GROMACS '.top' file to load
    bin_width: int
        Width (nm) of numpy histogram bins
    dim: int
        Dimension to calculate number density profile (0,1 or 2)
    box_range: array
        Range of coordinates in 'dim' to evaluate
    frame_range: Python range() (optional)
        Range of frames to calculate number density function over
    maxs: array (optional)
        Maximum coordinate to evaluate 
    mins: array (optional)
        Minimum coordinate to evalute
    
    Attributes
    ----------
    """
    trj = md.load(trj_file, top=gro_file)
    com_trj = make_comtrj(trj)
    resnames = np.unique([x.name for x in
               com_trj.topology.residues])
    rho_list = list()
    
    for resname in resnames:
        sliced = com_trj.topology.select('resname {}'.format(resname))
        trj_slice = com_trj.atom_slice(sliced)
        if frame_range:
            trj_slice = trj_slice[frame_range]
        for i,frame in enumerate(trj_slice):
            if maxs is None:
                indices = [[atom.index for atom in compound.atoms]
                          for compound in
                          list(frame.topology.residues)]
            else:
                indices = np.intersect1d(
                          np.intersect1d(np.where(frame.xyz[-1, :, 0]
                              > mins[0]),
                          np.where(frame.xyz[-1, :, 0] < maxs[0])),
                          np.intersect1d(np.where(frame.xyz[-1, :, 1] 
                              > box_range[0]),
                          np.where(frame.xyz[-1, :, 1] < box_range[1])))

            if frame_range:
                if i == 0:
                    x = np.histogram(frame.xyz[0,indices,dim].flatten(), 
                        bins=n_bins, range=(box_range[0], box_range[1]))
                    rho = x[0]
                    bins = x[1]
                else:
                    rho += np.histogram(frame.xyz[0, indices, dim].
                            flatten(),bins=n_bins, range=(box_range[0],
                                box_range[1]))[0]
            else:
                if i == 0:
                    x = np.histogram(frame.xyz[0,indices,dim].flatten(), 
                        bins=n_bins, range=(box_range[0], box_range[1]))
                    rho = x[0]
                    bins = x[1]
                else:
                    rho += np.histogram(frame.xyz[0, indices, dim].
                            flatten(),bins=n_bins, range=(box_range[0],
                                box_range[1]))[0]
        rho = np.divide(rho, trj_slice.n_frames * area *
                2 / n_bins)
        rho_list.append(rho)

    bin_list = bins[:-1]
    
    return(rho_list, bin_list)
