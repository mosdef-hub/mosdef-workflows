# **signac** demo of water models

[![signac](http://docs.signac.io/en/latest/_images/signac.io-blue.svg)](https://signac.io)

This example demonstrates a workflow using the **signac** data management framework.
This **signac** project manages a data space simulating SPC/E and tip3p water models.

## Setup

This demo requires all of the dependencies of the `water-models` demonstration (including `mbuild`, `foyer`, `gromacs`, and others).

Additionally, this demo requires the **signac** and **signac-flow** packages, which can be installed via `pip install signac signac-flow` or `conda install -c conda-forge signac signac-flow`.

## Initialization

Running `python init.py` initializes the **signac** workspace.

It creates the `signac.rc` file, which shows this directory contains a **signac** project, and initializes jobs.

The *state point* is a key-value mapping defining the essential and uniquely identifying parameters of a job. In some file-based workflows, it is easier to use the state point as a simple "unique identifier" and rely on files for the actual parameter definitions.

In this demo, we use the "model name" (for example, `spce`) as the only state point parameter, because it uniquely identifies each job.
We store the force field file `ff.xml` in the job directory, which actually defines the force field parameters.
*Note: An alternative approach would be to use state points defining each part of the water model's potential, like bond strengths and atom charges.*

## Running the Project

The `project.py` script uses **signac-flow** to prepare and execute the workflow.

Running `python project.py status --detailed` will show that there are two operations ready to run:

```
job_id                            operation                     labels
--------------------------------  ----------------------------  --------
70113b8988d2521c6738a034646e8abf  initialize_configuration [U]
a7593ab3ddff7f6e49039a16f68f2dbd  initialize_configuration [U]

[U]:unknown [R]:registered [Q]:queued [A]:active [I]:inactive [!]:requires_attention
```

Looking at `project.py`, we see a set of operations, which define the *workflow*.
Each operation has a set of *preconditions* and *postconditions*. These determine when each operation can be run.
For example, the *postconditions* for the `initialize` operation show that it is "complete" if the files `water.gro` and `water.top` exist.

The precondition  `@Project.pre.after(initialize)` on the `simulate` operation means that it will be eligible only if the `initialize` operation's postconditions are met.
This makes sense because we need the files `water.gro` and `water.top` to run the `simulate` operation.

Looking at the conditions in `project.py`, the operations are constructed in a linear fashion.
The operations will execute as follows:
```
initialize --> simulate --> energy --> analyze
```

We can execute these operations locally with `python project.py run`.
This will execute all eligible operations, and loop until no new operations are eligible.

If we are on a supercomputer, e.g. a university cluster with a SLURM scheduler, we can use `python project.py submit` to submit these operations as jobs to the cluster.

For further details about all features of the **signac** framework, please refer to the [signac documentation](https://docs.signac.io).
