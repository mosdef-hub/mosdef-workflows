#!/usr/bin/env python3
from flow import FlowProject, cmd, with_job
import signac

# Water models
# Here we will demonstrate how MoSDeF can be used to compare two different
# force fields; in this case, two atomistic water models. Here we are using two
# rigid, classical 3-site models,
# [SPC/E](https://dx.doi.org/doi/10.1021/j100308a038) and
# [tip3p](https://dx.doi.org/10.1063/1.445869).

root_directory = signac.get_project().root_directory()


class Project(FlowProject):
    pass


@Project.operation
@Project.post.isfile('water.gro')
@Project.post.isfile('water.top')
def initialize(job):
    import mbuild as mb
    from foyer import Forcefield
    # First, we will use mBuild and foyer to generate initial configurations

    # Generate a water compound from a SMILES string
    water = mb.load('O', smiles=True)
    water.name = 'water'

    # Fill a simulation box with water molecules
    system = mb.fill_box(compound=water, n_compounds=1000, density=993)

    # Load and apply force field
    ff = Forcefield(job.fn('ff.xml'))
    system_ff = ff.apply(system.to_parmed(residues=['water']))

    # Save GROMACS files for system
    system_ff.save(job.fn('water.gro'))
    system_ff.save(job.fn('water.top'))


@Project.operation
@Project.pre.after(initialize)
@Project.post.isfile('water.edr')
@with_job
@cmd
def simulate(job):
    # Now that we have the GRO and TOP files, we can prepare and run GROMACS
    # simulations. These should take approximately an hour each on a typical
    # quad-core laptop.
    return ('gmx grompp -c water.gro -p water.top -f {root}/files/npt.mdp -o '
            'water.tpr && gmx mdrun -v -deffnm water').format(
                root=root_directory)


@Project.operation
@Project.pre.after(simulate)
@Project.post.isfile('energy.xvg')
@with_job
@cmd
def energy(job):
    # Now that the simulations are complete, we begin our analysis. We will
    # first do some quick time series analyses and then compute a few RDFs. The
    # simulations equilibrate quickly, so we will only throw out the first 100
    # ps of the trajectories and run analyses on the last 900 ps.
    return 'echo 4 | gmx energy -f water.edr -o energy.xvg'


@Project.operation
@Project.pre.after(energy)
@Project.post.isfile('potential_energy.png')
@Project.post.isfile('density.png')
@Project.post.isfile('o_o_rdf.png')
@Project.post.isfile('h_h_rdf.png')
@Project.post.isfile('o_h_rdf.png')
@with_job
def analyze(job):
    from plots import potential_energy, density, o_o_rdf, h_h_rdf, o_h_rdf
    potential_energy([job])
    density([job])
    o_o_rdf([job])
    h_h_rdf([job])
    o_h_rdf([job])


if __name__ == '__main__':
    Project().main()
