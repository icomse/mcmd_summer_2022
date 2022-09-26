# Full documentation at http://physical-validation.readthedocs.io/
import physical_validation as pv

# Create a GROMACS parser, needs the location of the GROMACS executable
# and the location of the topology include folder.
# Here, we assume that `gmx` is in the PATH, and that the topology
# folder is in its standard location.
gmx = pv.data.GromacsParser(exe='gmx',
                            includepath='/usr/local/share/gromacs/top/')

# We'll test simulations ran with two thermostats:
# 'vr' stands for velocity-rescale, 'be' for Berendsen thermostat.
algos = ['vr', 'be']
# We'll test simulations performed in two ensembles: NVT and NPT
# Note that in NPT, the 'vr' thermostat was complemented with a
# Parinello-Rahman barostat, while the 'be' thermostat was complemented
# with a Berendsen barostat.
ensembles = ['NVT', 'NPT']
# The number of simulations for each ensemble
sims = {
    'NPT': 4, # We use 4 simulations for NPT: (T,P), (T+dT,P), (T,P+dP), (T+dT,P+dP)
    'NVT': 2  # We need 2 simulations for NVT: T and T+dT
}

# Dictionary we will store the parsed data in
simulations = {}

# Loop over the different thermostats and ensembles
for a in algos:
    simulations[a] = {}
    for e in ensembles:
        simulations[a][e] = []
        # Parse 4 simulations for NPT, 2 simulations for NVT
        for n in range(1, sims[e]+1):
            # Set directory
            d = 'md_' + e + '_' + a + '_' + str(n) + '/'
            # Read in the simulation results using the GROMACS parser.
            # This uses the `mdp` parameter file and the `top` topology file
            # to gather information about the system and the simulation settings,
            # and read the results from the `edr` file (trajectory of energy /
            # volume / pressure / ...) and the `gro` file (position and velocity
            # snapshot - used to read the box volume in NVT)
            simulations[a][e].append(
                gmx.get_simulation_data(
                    mdp=d + 'mdout.mdp',
                    top='top/system.top',
                    edr=d + 'system.edr',
                    gro=d + 'system.gro'
                )
            )
            
            # Test the kinetic energy distribution of the simulation result
            # read in last.
            # The first input is the simulation results read in,
            # `strict` determines whether we test the full distribution (True)
            # or only determine the mean and the variance of the distribution (False),
            # `verbosity` sets the level of detail of the output  (with verbosity=0
            # being quiet and verbosity=3 being the most chatty),
            # and the filename is being used to plot the resulting distribution for
            # visual inspection.
            print('==> Kinetic energy test of simulation' + e + '_' + a + '_' + str(n))
            pv.kinetic_energy.distribution(simulations[a][e][-1], strict=True, verbosity=2)
            pv.kinetic_energy.distribution(simulations[a][e][-1], strict=False, verbosity=2,
                                           filename='_'.join(['ke', a, e, str(n)]))

        # Now that we have all simulations of the current thermostat and ensemble
        # read in, we can test the distribution of the potential energy and (for NPT)
        # the volume. While the first two inputs to the tests are the parsed simulation
        # results, `verbosity` sets the level of detail of the output  (with verbosity=0
        # being quiet and verbosity=3 being the most chatty), and the filename is being
        # used to plot the resulting distribution for visual inspection.
        if e == 'NVT':
            pv.ensemble.check(simulations[a][e][0], simulations[a][e][1],
                              verbosity=2, filename='_'.join(['pe', a, e]))
        else:
            # There are three checks we can do: P(E_1)/P(E_2), P(V_1)/P(V_2) and P(E_1,V_1)/P(E_2,V_2)
            pv.ensemble.check(simulations[a][e][0], simulations[a][e][1],
                              verbosity=2, filename='_'.join(['pe', a, e, 'dT']))
            pv.ensemble.check(simulations[a][e][0], simulations[a][e][2],
                              verbosity=2, filename='_'.join(['pe', a, e, 'dP']))
            pv.ensemble.check(simulations[a][e][0], simulations[a][e][3],
                              verbosity=2) # Plotting for the 2D case is not supported
