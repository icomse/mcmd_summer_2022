# Full documentation at http://physical-validation.readthedocs.io/
import physical_validation as pv

# Create a GROMACS parser, needs the location of the GROMACS executable
# and the location of the topology include folder.
# Here, we assume that `gmx` is in the PATH, and that the topology
# folder is in its standard location.
gmx = pv.data.GromacsParser(exe='gmx',
                            includepath='/usr/local/share/gromacs/top/')

# The directories containing the simulation results.
# 'vr' stands for velocity-rescale, 'be' for Berendsen thermostat.
# '_1' denotes a simulation using a single thermostat for the entire system.
# '_2' denotes a simulation using two separate thermostats, one for the
# protein and one for the solvent.
dirs = ['vr_1', 'vr_2', 'be_1', 'be_2']

for d in dirs:
    # Read in the simulation results using the GROMACS parser.
    # This uses the `mdp` parameter file and the `top` topology file
    # to gather information about the system and the simulation settings,
    # and read the results from the `edr` file (trajectory of energy /
    # volume / pressure / ...) and the `trr` file (position and velocity
    # trajectory)
    # Note that the topology and the trajectory have been modified to
    # only contain the protein, as we are only interested in the equiparition
    # of the solute, not the solvent. This reduces file size and execution
    # time considerably.
    results_protein = gmx.get_simulation_data(mdp=str(d) + '/protein.mdp',
                                              top=str(d) + '/trp-cage.top',
                                              edr=str(d) + '/run.edr',
                                              trr=str(d) + '/protein.trr')
    # In the simulations, the center of mass was artificially kept immobile
    # to avoid the build up of numerical errors, effectively reducing the
    # number of translational degrees of freedom of the system by 3. As we
    # are only looking at the protein here, we're using a mass-dependent
    # fraction of these three degrees of freedom (weight of protein:
    # 2170.4375 amu; total weight of the system: 96390.4017 amu).
    results_protein.system.ndof_reduction_tra *= 2170.4375 / 96390.4017
    print('===> ' + d)

    # Here, we run the equipartition test.
    # The first input is the simulation results read in earlier,
    # the second denotes the filename for the plotting,
    # and the third sets the level of detail of the output (with
    # verbosity=0 being quiet and verbosity=3 being the most chatty).
    pv.kinetic_energy.equipartition(results_protein,strict=False,
                                    filename='equipartition_'+str(d),
                                    verbosity=3)
