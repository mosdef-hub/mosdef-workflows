#!/usr/bin/env python3
import numpy as np
import mdtraj as md
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import InsetPosition
from project import Project


matplotlib.use('Agg')
START_FRAME = 100
LONG_STRIDE = 10


def potential_energy(jobs):
    fig, ax = plt.subplots()

    for job in jobs:
        data = np.loadtxt(job.fn('energy.xvg'), comments=['#', '@'])

        ax.plot(data[:, 0], data[:, 1], label=job.sp.model_name)
        axins = plt.axes([0, 0, 1, 1])

        ip = InsetPosition(ax, [0.2, 0.25, 0.7, 0.5])
        axins.set_axes_locator(ip)

        axins.plot(data[START_FRAME:, 0], data[START_FRAME:, 1])

    ax.legend(loc=0)
    ax.set_xlabel('Time, ps')
    ax.set_ylabel('Potential energy, kJ/mol')

    fig.savefig('potential_energy.png', bbox_inches='tight')


def density(jobs):
    fig, ax = plt.subplots()

    for job in jobs:
        trj = md.load(
            job.fn('water.xtc'),
            top=job.fn('water.gro'),
        )

        rho = md.density(trj)

        ax.plot(trj.time, rho, label=job.sp.model_name)
        axins = plt.axes([0, 0, 1, 1])

        ip = InsetPosition(ax, [0.2, 0.3, 0.7, 0.4])
        axins.set_axes_locator(ip)

        axins.plot(trj.time[START_FRAME:], rho[START_FRAME:])

    ax.legend(loc='lower right')
    ax.set_xlabel('Time, ps')
    ax.set_ylabel('Density, $kg/m^3$')
    fig.savefig('density.png')


def o_o_rdf(jobs):
    fig, ax = plt.subplots()

    for job in jobs:
        trj = md.load(
            job.fn('water.xtc'),
            top=job.fn('water.gro'),
            stride=LONG_STRIDE,
        )

        pairs = trj.top.select_pairs('name O', 'name O')

        r, g_r = md.compute_rdf(trj[START_FRAME:], pairs)

        plt.plot(r, g_r, label=job.sp.model_name)

    ax.legend()
    ax.set_title('O-O RDF')
    ax.set_xlabel('Interatomic separation, $r$')
    ax.set_ylabel('Pair distribution function, unitless')
    fig.savefig('o_o_rdf.png')


def h_h_rdf(jobs):
    fig, ax = plt.subplots()

    for job in jobs:
        trj = md.load(
            job.fn('water.xtc'),
            top=job.fn('water.gro'),
            stride=LONG_STRIDE,
        )

        pairs = trj.top.select_pairs('name H', 'name H')

        r, g_r = md.compute_rdf(trj[START_FRAME:], pairs,
                                r_range=(0.145, 0.17), n_bins=200)

        plt.plot(r, g_r, label=job.sp.model_name)

    ax.plot([0.1633, 0.1633], [0, 60], 'k-')
    ax.plot([0.1513, 0.1513], [0, 60], 'k-')

    ax.legend()
    ax.set_title('Intramolecular H-H RDF')
    ax.set_xlim((0.145, 0.17))
    ax.set_xlabel('Interatomic separation, $r$')
    ax.set_ylabel('Pair distribution function, unitless')
    fig.savefig('h_h_rdf.png')


def o_h_rdf(jobs):
    fig, ax = plt.subplots()

    for job in jobs:
        trj = md.load(
            job.fn('water.xtc'),
            top=job.fn('water.gro'),
            stride=LONG_STRIDE,
        )

        pairs = trj.top.select_pairs('name O', 'name H')

        r, g_r = md.compute_rdf(trj[START_FRAME:], pairs,
                                r_range=(0.085, 0.105), n_bins=200)

        plt.plot(r, g_r, label=job.sp.model_name)

    ax.plot([0.09572, 0.09572], [0, 250], 'k-')
    ax.plot([0.1, 0.1], [0, 250], 'k-')

    ax.legend()
    ax.set_title('Intramolecular O-H RDF')
    ax.set_xlim((0.085, 0.105))
    ax.set_xlabel('Interatomic separation, $r$')
    ax.set_ylabel('Pair distribution function, unitless')
    fig.savefig('o_h_rdf.png')


if __name__ == '__main__':
    project = Project()
    potential_energy(project)
    density(project)
    o_o_rdf(project)
    h_h_rdf(project)
    o_h_rdf(project)
